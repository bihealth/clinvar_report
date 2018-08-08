# -*- coding: utf-8 -*-

import argparse
import collections
import contextlib
import gzip
import logging
import re
import sys
import typing
from urllib.parse import unquote

import attr
import intervaltree
import vcfpy
from clinvar_report import __version__


Individual = collections.namedtuple(
    "Individual", ["family", "name", "father", "mother", "gender", "affected"]
)

SEX = {"0": "unknown", "1": "male", "2": "female"}

AFFECTED = {"0": "unknown", "1": "unaffected", "2": "affected"}


class Pedigree:
    """Represents a pedigree"""

    @classmethod
    def load(klass, file):
        entries = []
        for line in file:
            arr = line.strip().split("\t")
            entries.append(
                Individual(
                    arr[0],
                    "-".join(arr[1].split("-")[:-3]),
                    "-".join(arr[2].split("-")[:-3]),
                    "-".join(arr[3].split("-")[:-3]),
                    SEX[arr[4]],
                    AFFECTED[arr[5]],
                )
            )
        return Pedigree(entries)

    def __init__(self, entries):
        self.entries = entries
        self._father_of = {e.name: e.father for e in self.entries}
        self._mother_of = {e.name: e.mother for e in self.entries}
        self._by_name = {e.name: e for e in self.entries}
        self.by_family = {}
        for entry in self.entries:
            self.by_family.setdefault(entry.family, []).append(entry)

    def get_individual(self, name):
        return self._by_name.get(name)

    def get_father(self, name):
        """Return id of father, if any, otherwise None"""
        result = self._father_of.get(name)
        if result == "0":
            return None
        return result

    def get_mother(self, name):
        """Return id of mother, if any, otherwise None"""
        result = self._mother_of.get(name)
        if result == "0":
            return None
        return result

    def print_ped(self, file):
        for e in self.entries:
            print(
                "\t".join(
                    map(
                        str,
                        [e.family, e.name, e.father, e.mother, e.gender, e.affected],
                    )
                ),
                file=file,
            )


class Annotation(typing.NamedTuple):
    """Jannovar ANN"""

    allele: str
    impact: str
    effect: str
    gene_name: str
    gene_id: str
    transcript: str
    hgvs_c: str
    hgvs_p: str


class BasicInfo(typing.NamedTuple):
    """Represents a ``BASIC_INFO`` entry."""

    allele: str
    hgvs_string: str
    origin: str


class VarInfo(typing.NamedTuple):
    """Represents a ``VAR_INFO`` entry."""

    allele: str
    database_name: str
    id_in_db: str
    origins: typing.Tuple[str]


class GenotypeInfo(typing.NamedTuple):
    """Represent genotype information."""

    sample: str
    gt: str
    ad: typing.List[int]

    @property
    def aaf(self):
        if sum(self.ad) == 0:
            return 0.0
        else:
            return 100 * self.ad[1] // sum(self.ad)

@attr.s
class Variant:
    """Represents an annotated variant."""

    chrom = attr.ib()
    pos = attr.ib()
    ref = attr.ib()
    alt = attr.ib()
    affected_start = attr.ib()
    affected_end = attr.ib()
    exac_freq = attr.ib()
    exac_hom = attr.ib()
    annotations = attr.ib()
    gt_infos = attr.ib()
    inheritance_modes = attr.ib()
    common = attr.ib()
    candidate = attr.ib()
    clinvars = attr.ib(default=[])
    clinvars_ovl = attr.ib(default=[])

    gene_name = attr.ib(init=False)
    transcript = attr.ib(init=False)
    hgvs_c = attr.ib(init=False)
    hgvs_p = attr.ib(init=False)
    effects = attr.ib(init=False)

    def __attrs_post_init__(self):
        allele = self.alt[0]
        if allele in self.annotations:
            self.gene_name = self.annotations[allele][0].gene_name
        else:
            self.gene_name = "-"

        if allele in self.annotations:
            self.transcript = self.annotations[allele][0].transcript
        else:
            self.transcript = "-"
        if allele in self.annotations:
            self.hgvs_c = self.annotations[allele][0].hgvs_c
        else:
            self.hgvs_c = None
        if allele in self.annotations:
            self.hgvs_p = self.annotations[allele][0].hgvs_p
        else:
            self.hgvs_p = None
        if allele in self.annotations:
            self.effects = self.annotations[allele][0].effect
        else:
            self.effects = []

    @property
    def benign(self):
        return all(
            "benign" in clinvar.clinical_significance.lower()
            or "drug response" in clinvar.clinical_significance.lower()
            or "protective" in clinvar.clinical_significance.lower()
            for clinvar in self.clinvars
        )

    @property
    def uncertain_conflict(self):
        return all(
            "conflict" in clinvar.clinical_significance.lower()
            or "uncertain" in clinvar.clinical_significance.lower()
            or "not provided" in clinvar.clinical_significance.lower()
            for clinvar in self.clinvars
        )

    @property
    def variation_ids(self):
        result = []
        for clinvar in self.clinvars:
            result += clinvar.variation_id
        return result

    @property
    def significances(self):
        result = []
        for clinvar in self.clinvars:
            result.append(clinvar.clinical_significance)
        return result

    @property
    def conflict(self):
        return any(var.conflicted for var in self.clinvars)

    @property
    def diseases(self):
        result = []
        result_lower = ["not specified"]
        for var in self.clinvars:
            for trait in var.all_traits:
                if not trait.lower() in result_lower:
                    result.append(trait)
                    result_lower.append(trait.lower())
        return result

    @property
    def locus(self):
        return "{}:{}-{}".format(self.chrom, self.pos, self.pos)

    @property
    def pretty(self):
        import simplejson as json
        import pprint

        return pprint.pformat(json.loads(json.dumps(attr.asdict(self))))

    @property
    def gold_stars(self):
        """Return maximal number of gold stars for all ClinVar entries, including overlaps."""
        if not self.clinvars and not self.clinvars_ovl:
            return 0
        else:
            return max(
                self.clinvars + self.clinvars_ovl, key=lambda x: x.gold_stars
            ).gold_stars

    @property
    def gold_stars_title(self):
        """Return review status title to display"""
        if not self.clinvars and not self.clinvars_ovl:
            return 0
        else:
            return max(
                self.clinvars + self.clinvars_ovl, key=lambda x: x.gold_stars
            ).review_status


def variant_from_vcf_record(var, ped, args):
    # Parse ANN
    annotations = collections.OrderedDict()
    if "ANN" in var.INFO:
        for strval in var.INFO["ANN"].split(","):
            values = strval.split("|")
            allele = values[0]
            effect = values[1].split("&")
            impact = values[2]
            gene_name = values[3]
            gene_id = values[4]
            transcript = values[6]
            hgvs_c = values[9]
            hgvs_p = values[10]
            annotations.setdefault(allele, []).append(
                Annotation(
                    allele,
                    impact,
                    effect,
                    gene_name,
                    gene_id,
                    transcript,
                    hgvs_c,
                    hgvs_p.replace("(", "").replace(")", ""),
                )
            )

    # Get out GT information
    gt_infos = collections.OrderedDict()
    for call in var.calls:
        gt_infos[call.sample] = GenotypeInfo(
            call.sample, call.data["GT"], call.data.get("AD")
        )

    # Get compatible modes of inheritance
    inheritance_modes = var.INFO.get("INHERITANCE", [])

    # Any affected individual has a het. or hom. alt. call
    affecteds = [i.name for i in ped.entries if i.affected == AFFECTED["2"]]
    candidate = False
    for call in var.calls:
        sample_short = "-".join(call.sample.split("-")[:-3])
        if sample_short in affecteds and call.called and call.gt_type != vcfpy.HOM_REF:
            candidate = True

    exac_freq = var.INFO.get('EXAC_AF_POPMAX', [0.0])
    exac_hom = var.INFO.get('EXAC_HOM_ALL', [0])

    alts = [a.serialize() for a in var.ALT]
    return Variant(
        var.CHROM,
        var.POS,
        var.REF,
        alts,
        var.affected_start,
        var.affected_end,
        exac_freq,
        exac_hom,
        annotations,
        gt_infos,
        inheritance_modes,
        var.INFO.get("DBSNP_COMMON", False),
        candidate,
    )


def setup_logging(args):
    """Setup logger."""
    logging.basicConfig(
        format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
        datefmt="%m-%d %H:%M",
    )
    logger = logging.getLogger("")
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)


def var_group(variant):
    """Returns group for the variant"""
    if not variant.candidate:
        return "bad_affecteds"
    elif variant.common:
        return "bad_common"
    elif variant.benign:
        return "bad_benign"
    elif variant.uncertain_conflict:
        return "bad_uncertain_conflict"
    else:
        return "good"


class ClinVarRecord(typing.NamedTuple):
    """One record from MacArthur ClinVar TSV file."""

    chrom: str
    pos: int
    ref: str
    alt: str
    start: int
    stop: int
    strand: str
    variation_type: str
    variation_id: typing.List[int]
    rcv: typing.List[str]
    scv: typing.List[str]
    allele_id: str
    symbol: str
    hgvs_c: str
    hgvs_p: str
    molecular_consequence: str
    clinical_significance: str
    clinical_significance_ordered: typing.List[str]
    pathogenic: int
    likely_pathogenic: int
    uncertain_significance: int
    likely_benign: int
    benign: int
    review_status: str
    review_status_ordered: typing.List[str]
    last_evaluated: str
    all_submitters: typing.List[str]
    submitters_ordered: str
    all_traits: typing.List[str]
    all_pmids: typing.List[str]
    inheritance_modes: str
    age_of_onset: str
    prevalence: str
    disease_mechanism: str
    origin: str
    xrefs: typing.List[str]
    dates_ordered: typing.List[str]
    gold_stars: int
    conflicted: int


class ClinVarContig:
    """Database of ClinVar variants on one contig."""

    @staticmethod
    def _coerce(val_type):
        val, t = val_type
        if getattr(t, "__origin__", None) == typing.List:
            f = t.__args__[0]
            return list(map(f, val.split(";")))
        else:
            if val == "" or (t is int and val == "-"):
                return None
            else:
                return t(val)

    @classmethod
    def load(klass, tsv_readers, contig):
        """Load ClinVar variants from the given contig.

        NB: Because pysam.TabixFile cannot parse the MacArthur ClinVar TSV file, we are simply
        reading the whole file on each iteration which has OK performance...
        """
        records = []
        for reader in tsv_readers:
            reader.seek(0)
            for line in reader:
                if line.startswith(contig + "\t"):
                    vals = list(
                        map(
                            klass._coerce,
                            zip(
                                line.split("\t"), ClinVarRecord.__annotations__.values()
                            ),
                        )
                    )
                    records.append(ClinVarRecord(*vals))
        logging.info("Found %d ClinVar records on contig %s", len(records), contig)
        return ClinVarContig(contig, records)

    def __init__(self, contig, records):
        #: Name of the contig to query.
        self.contig = contig
        #: Flat list of all records.
        self.records = records
        #: Interval tree of all records.
        self.itree = intervaltree.IntervalTree.from_tuples(
            (record.start - 1, record.stop, record) for record in self.records
        )

    def query(self, variant):
        """Query for matches and overlaps with ``variant``.

        Args:

        ``variant`` - The ``Variant`` to query.
        """
        matches, overlaps = [], []
        for itv in self.itree.search(variant.affected_start, variant.affected_end):
            record = itv.data
            if (
                record.start - 1 == variant.affected_start
                and record.stop == variant.affected_end
                and record.alt in variant.alt
            ):
                matches.append(record)
            else:
                overlaps.append(record)
        return matches, overlaps


def process_contig(contig, length, vcf_reader, clinvar_readers, ped, args):
    """Process the given contig and yield ``Variant`` objects."""
    logging.info("Processing contig %s", contig)
    yielded = 0
    clinvar_contig = ClinVarContig.load(clinvar_readers, contig)
    try:
        vcf_it = vcf_reader.fetch(contig, 0, int(length))
    except ValueError as e:
        logging.warn("Problem fetching VCF records: %s", e)
        vcf_it = []
    for record in vcf_it:
        variant = variant_from_vcf_record(record, ped, args)
        variant.clinvars, variant.clinvars_ovl = clinvar_contig.query(variant)
        variant.clinvars_ovl = []  # TODO: ignore for now --> activate?
        if variant.clinvars or variant.clinvars_ovl:
            yielded += 1
            yield variant
    logging.info("Annotated %d variants on contig %s", yielded, contig)


def get_contigs(vcf_reader):
    """Get list of ``(name, length)`` of contigs from ``vcf_reader``."""
    contigs = [(line.id, line.length) for line in vcf_reader.header.get_lines("contig")]
    return contigs


def run(args):
    """Main program entry point after parsing command line arguments."""
    setup_logging(args)

    with open(args.input_ped, "rt") as f:
        ped = Pedigree.load(f)

    logging.info("ClinVar-annotating input VCF file %s", args.input_vcf)
    variants = []
    with vcfpy.Reader.from_path(args.input_vcf) as vcf_reader:
        with contextlib.ExitStack() as stack:
            # Open all ClinVar tsv readers.
            tsv_readers = [
                stack.enter_context(gzip.open(clinvar_tsv, "rt"))
                for clinvar_tsv in args.clinvar_tsvs
            ]
            # Get samples, and short sample names.
            samples = vcf_reader.header.samples.names
            short_samples = ["-".join(s.split("-")[:-3]) for s in samples]
            # Process all contigs.
            for contig, length in get_contigs(vcf_reader):
                if re.match(args.contig_regex, contig):
                    variants += list(
                        process_contig(
                            contig, length, vcf_reader, tsv_readers, ped, args
                        )
                    )
    logging.info(" => done, have %d records", len(variants))

    logging.info("Grouping variants...")
    grouped_vars = {}
    for var in variants:
        grouped_vars.setdefault(var_group(var), []).append(var)
    for key in grouped_vars:
        grouped_vars[key].sort(key=lambda v: v.gold_stars, reverse=True)
    logging.info(" => done.")

    logging.info("Grouped variants")
    for key, lst in sorted(grouped_vars.items()):
        logging.info("%s: %d", key, len(lst))

    values = {
        "variants": grouped_vars,
        "samples": samples,
        "short_samples": short_samples,
        "pedigree": ped,
        "version": __version__,
    }

    logging.info("Writing report file %s", args.output_html)
    from jinja2 import Environment, PackageLoader, select_autoescape

    env = Environment(
        loader=PackageLoader("clinvar_report", "templates"),
        autoescape=select_autoescape(["html", "xml"]),
    )
    template = env.get_template("report.html")
    with open(args.output_html, "wt") as f:
        print(template.render(**values), file=f)
    logging.info(" => done writing ClinVar report")


def main(argv=None):
    """Main program entry point, starts parsing command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--contig-regex",
        default="(chr)?([0-9]+|X|Y|M|MT)",
        help="Regular expression for selecting contigs",
    )
    parser.add_argument(
        "--verbose", action="store_true", default=False, help="Enable verbose mode"
    )
    parser.add_argument(
        "--input-ped",
        required=True,
        help="Path to pedigree to use for affected filtering",
    )
    parser.add_argument("--input-vcf", required=True, help="VCF file to annotate")
    parser.add_argument(
        "--clinvar-tsv",
        dest="clinvar_tsvs",
        type=str,
        required=True,
        default=[],
        action="append",
        help="Path to tabix-indexed clinvar TSV",
    )
    parser.add_argument(
        "--output-html", type=str, required=True, help="Path to output HTML file"
    )

    args = parser.parse_args(argv)
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
