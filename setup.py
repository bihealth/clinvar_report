#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

from setuptools import setup, find_packages
import pip


def parse_requirements(path):
    """Parse ``requirements.txt`` at ``path``."""
    requirements = []
    with open(path, "rt") as reqs_f:
        for line in reqs_f:
            line = line.strip()
            if line.startswith("-r"):
                fname = line.split()[1]
                inner_path = os.path.join(os.path.dirname(path), fname)
                requirements += parse_requirements(inner_path)
            elif line != "" and not line.startswith("#"):
                requirements.append(line)
    return requirements


setup(
    name="clinvar-report",
    description=(
        "Utility for building Clinvar reports from " "Jannovar-annotated VCF files"
    ),
    packages=find_packages(),
    package_dir={"clinvar_report": "clinvar_report"},
    entry_points={"console_scripts": ["clinvar-report = clinvar_report.__main__:main"]},
    include_package_data=True,
    license="MIT license",
    zip_safe=False,
)
