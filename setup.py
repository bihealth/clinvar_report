#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

from setuptools import setup, find_packages
import pip
from pip.req import parse_requirements

setup(
    name='clinvar-report',
    description=(
        'Utility for building Clinvar reports from '
        'Jannovar-annotated VCF files'),
    packages=find_packages(),
    package_dir={
        'clinvar_report': 'clinvar_report',
    },
    entry_points={
        'console_scripts': [
            'clinvar-report = clinvar_report.__main__:main',
        ]
    },
    include_package_data=True,
    license='MIT license',
    zip_safe=False,
)
