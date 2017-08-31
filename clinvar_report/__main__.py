# -*- coding: utf-8 -*-

import argparse
import collections
import logging
import re
import sys
from urllib.parse import unquote

import vcfpy

#: Constant with INFO key for ``BASIC_INFO``.
BASIC_INFO = 'BASIC_INFO'
#: Constant with INFO key for ``DISEASE_INFO``.
DISEASE_INFO = 'DISEASE_INFO'
#: Constant with INFO key for ``VAR_INFO``.
VAR_INFO = 'VAR_INFO'
#: DB base URLs
DB_BASE_URLS = {
    'MedGen': 'https://www.ncbi.nlm.nih.gov/gtr/conditions/',
    'OMIM': 'https://www.omim.org/entry/',
    'SNOMED_CT': 'https://phinvads.cdc.gov/vads/ViewValuhttps:/phinvads.cdc.gov/vads/ViewCodeSystemConcept.action?oid=2.16.840.1.113883.6.96&code='
}
#: Pathogenicity to level
SIG_LEVELS = {
    'histocompatibility': 7,
    'drug_response': 6,
    'pathogenic': 5,
    'likely_pathogenic': 4,
    'likely_benign': 3,
    'benign': 2,
    'uncertain': 1,
    'not_provided': 0,
    'other': 255,
}

def shorten_aa(x):
    d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
         'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
         'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    prefix = x.startswith('p.')
    x = x.replace('p.(', '').replace(')', '').upper()
    for k, v in d.items():
        x = x.replace(k, v)
    if prefix:
        return 'p.' + x
    else:
        return x

def fixhex(text):
    """Fix hexadecimal encoding of utf 8 characters."""
    return re.sub(r'\\x[0-9a-f][0-9a-f]', lambda m: chr(int(m.group(0)[2:], 16)), text)


class Annotation:
    """Jannovar ANN"""

    # Allele|Annotation|Annotation_Impact|Gene_Name|Gene_ID|Feature_Type|Feature_ID|Transcript_BioType|Rank|HGVS.c|HGVS.p|cDNA.pos / cDNA.length|CDS.pos / CDS.length|AA.pos / AA.length|Distance|ERRORS / WARNINGS / INFO
    def __init__(
            self, allele, impact, effect, gene_name, gene_id, transcript,
            hgvs_c, hgvs_p):
        self.allele = allele
        self.impact = impact
        self.effect = effect
        self.gene_name = gene_name
        self.gene_id = gene_id
        self.transcript = transcript
        self.hgvs_c = hgvs_c
        self.hgvs_p = shorten_aa(hgvs_p)


class BasicInfo:
    """Represents a ``BASIC_INFO`` entry."""

    def __init__(self, allele, hgvs_string, origin):
        #: The annotated allele
        self.allele = allele
        #: The annotating HGVS string
        self.hgvs_string = hgvs_string
        #: The variant origin
        self.origin = origin

    def __repr__(self):
        return 'BasicInfo({})'.format(', '.join(map(
            str, (self.allele, self.hgvs_string, self.origin))))


# "Annotation of disease information of the form 'allele | significance | disease db | id in disease db | name in disease db | revision status | clinical accession'

class DiseaseInfo:
    """Represents a ``DISEASE_INFO`` entry."""

    def __init__(
            self, allele, significance, disease_db, id_in_db, name_in_db,
            revision_status, clinical_accession):
        self.allele = allele
        self.significance = significance
        self.db_infos = []
        ids = list(map(lambda x: x.replace('HP__', 'HP:'),
                       id_in_db.replace('HP:', 'HP__').split(':')))
        for db, id_ in zip(disease_db.split(':'), ids):
            self.db_infos.append({
                'db': db,
                'id': id_,
                'name': name_in_db,
                'url': DB_BASE_URLS.get(db, '') + id_,
            })
        self.name_in_db = name_in_db
        self.revision_status = revision_status
        self.clinical_accession = clinical_accession

    def __repr__(self):
        return 'DiseaseInfo({})'.format(', '.join(map(
            str, (self.allele, self.significance, self.disease_db,
                  self.id_in_db, self.name_in_db, self.revision_status,
                  self.clinical_accession))))


class VarInfo:
    """Represents a ``VAR_INFO`` entry."""

    def __init__(self, allele, database_name, id_in_db, origins):
        #: The annotated allele
        self.allele = allele
        #: The database name
        self.database_name = database_name
        #: The variant identifier in database
        self.id_in_db = id_in_db
        #: The variant origins
        self.origins = tuple(origins)

    def __repr__(self):
        return 'VarInfo({})'.format(', '.join(map(
            str, (self.allele, self.database_name, self.id_in_db,
                  self.origins))))


class GenotypeInfo:
    """Represent genotype information."""

    def __init__(self, sample, gt, ad):
        self.sample = sample
        self.gt = gt
        self.ad = ad

    def __repr__(self):
        return 'GenotypeInfo({})'.format(', '.join(map(
            str, (self.sample, self.gt, self.ad))))


class Variant:
    """Represents an annotated variant."""

    def __init__(
            self, chrom, pos, ref, alt, annotations, basic_infos,
            disease_infos, variant_infos, gt_infos):
        #: Chromosome with variant
        self.chrom = chrom
        #: Chromosomal variant position
        self.pos = pos
        #: Reference allele
        self.ref = ref
        #: Alternative alleles
        self.alt = alt
        #: Annotations
        self.annotations = annotations
        #: Basic infos for each allele
        self.basic_infos = basic_infos
        #: Disease-related information
        self.disease_infos = disease_infos
        #: Variant-specific information
        self.variant_infos = variant_infos
        #: Genotype information
        self.gt_infos = gt_infos
        #: Gene name
        allele = self.alt[0]
        if allele in self.annotations:
            self.gene_name = self.annotations[allele][0].gene_name
        else:
            self.gene_name = '-'
        #: Transcript
        if allele in self.annotations:
            self.transcript = self.annotations[allele][0].transcript
        else:
            self.transcript = '-'
        #: HGVS.c
        if allele in self.annotations:
            self.hgvs_c = self.annotations[allele][0].hgvs_c
        else:
            self.hgvs_c = None
        #: HGVS.p
        if allele in self.annotations:
            self.hgvs_p = self.annotations[allele][0].hgvs_p
        else:
            self.hgvs_p = None
        #: effect
        if allele in self.annotations:
            self.effects = self.annotations[allele][0].effect
        else:
            self.effects = []
        #: ClinVar accessions
        self.clinical_accessions = [
            i.clinical_accession for i in self.disease_infos[allele]]
        #: Significance
        self.significances = {}
        for i in self.disease_infos[allele]:
            if i.significance != 'other':
                self.significances.setdefault(SIG_LEVELS[i.significance], 0)
                self.significances[SIG_LEVELS[i.significance]] += 1
        #: Diseases
        self.diseases = {}
        for i in self.disease_infos[allele]:
            for db_info in i.db_infos:
                if db_info['name'] not in ('not specified', 'not provided'):
                    self.diseases[db_info['id']] = db_info['name']
        self.diseases = list(sorted(set(self.diseases.values())))

    def __repr__(self):
        return 'VarInfo({})'.format(', '.join(map(
            str, (self.chrom, self.pos, self.ref, self.alt, self.annotations,
                  self.basic_infos, self.disease_infos, self.variant_infos,
                  self.gt_infos))))


def parse_record(var, args):
    # Parse ANN
    annotations = collections.OrderedDict()
    if 'ANN' in var.INFO:
        for strval in var.INFO['ANN'].split(','):
            values = strval.split('|')
            allele = values[0]
            effect = values[1].split('&')
            impact = values[2]
            gene_name = values[3]
            gene_id = values[4]
            transcript = values[6]
            hgvs_c = values[9]
            hgvs_p = values[10]
            annotations.setdefault(allele, []).append(
                Annotation(
                    allele, impact, effect, gene_name, gene_id, transcript,
                    hgvs_c, hgvs_p))

    # Parse BASIC_INFO
    basic_infos = collections.OrderedDict()
    for basic_info in var.INFO.get(args.clinvar_prefix + BASIC_INFO, []):
        allele, hgvs_string, origin = basic_info.split('|')
        hgvs_string = unquote(hgvs_string)
        basic_infos.setdefault(allele, []).append(
            BasicInfo(allele, hgvs_string, origin))

    # Parse DISEASE_INFO
    disease_infos = collections.OrderedDict()
    for disease_info in var.INFO.get(args.clinvar_prefix + DISEASE_INFO, []):
        (allele, significance, disease_db, id_in_db, name_in_db,
         revision_status, clinical_accession) = disease_info.split('|')
        hgvs_string = unquote(hgvs_string)
        name_in_db = fixhex(unquote(name_in_db)).replace('_', ' ')
        disease_infos.setdefault(allele, []).append(
            DiseaseInfo(allele, significance, disease_db, id_in_db, name_in_db,
                        revision_status, clinical_accession))

    # Parse VAR_INFO
    var_infos = collections.OrderedDict()
    for var_info in var.INFO.get(args.clinvar_prefix + VAR_INFO, []):
        # allele | db name | id in db
        allele, db_name, id_in_db, origins = var_info.split('|')
        db_name = unquote(db_name)
        id_in_db = unquote(id_in_db)
        origins = origins.split('&')
        var_infos.setdefault(allele, []).append(
            VarInfo(allele, db_name, id_in_db, origins))

    # Get out GT information
    gt_infos = collections.OrderedDict()
    for call in var.calls:
        gt_infos[call.sample] = GenotypeInfo(
            call.sample, call.data['GT'], call.data['AD'])

    alts = [a.serialize() for a in var.ALT]
    alts = [a for a in alts
            if a in basic_infos or a in disease_infos or a in var_infos]
    return Variant(
        var.CHROM, var.POS, var.REF, alts, annotations, basic_infos,
        disease_infos, var_infos, gt_infos)


def setup_logging(args):
    """Setup logger."""
    logging.basicConfig(
        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
        datefmt='%m-%d %H:%M')
    logger = logging.getLogger('')
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)


def run(args):
    """Main program entry point after parsing command line arguments."""
    setup_logging(args)
    logging.info('Reading input VCF file %s', args.input_vcf)
    with vcfpy.Reader.from_path(args.input_vcf) as reader:
        samples = reader.header.samples.names
        short_samples = ['-'.join(s.split('-')[:-3]) for s in samples]
        records = [r for r in reader]
    logging.info(' => done reading %d records', len(records))

    logging.info('Parsing ClinVar VCF')
    values = {
        'variants': [parse_record(r, args) for r in records],
        'samples': samples,
        'short_samples': short_samples,
    }
    logging.info(' => done parsing ClinVar VCF')

    logging.info('Writing report file %s', args.output_html)
    from jinja2 import Environment, PackageLoader, select_autoescape
    env = Environment(
        loader=PackageLoader('clinvar_report', 'templates'),
        autoescape=select_autoescape(['html', 'xml']))
    template = env.get_template('report.html')
    with open(args.output_html, 'wt') as f:
        print(template.render(**values), file=f)
    logging.info(' => done writing ClinVar report')


def main(argv=None):
    """Main program entry point, starts parsing command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--verbose', action='store_true', default=False,
        help='Enable verbose mode')
    parser.add_argument(
        '--input-vcf', type=str, required=True,
        help='Path to input VCF file')
    parser.add_argument(
        '--output-html', type=str, required=True,
        help='Path to output HTML file')
    parser.add_argument(
        '--clinvar-prefix', default='CLINVAR_',
        help='Prefix of clinvar variants in VCF file')

    args = parser.parse_args(argv)
    return run(args)


if __name__ == '__main__':
    sys.exit(main())