#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
psi4_step
A step in SEAMM to setup, run and analyze quantum chemistry calculations
using Psi4
"""
import sys
from setuptools import setup, find_packages
import versioneer

short_description = __doc__.split("\n")

# from https://github.com/pytest-dev/pytest-runner#conditional-requirement
needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)
pytest_runner = ['pytest-runner'] if needs_pytest else []

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

with open('requirements_install.txt') as fd:
    requirements = fd.read()

setup(
    name='psi4_step',
    author="Paul Saxe",
    author_email='psaxe@molssi.org',
    description=short_description[1],
    long_description=readme + '\n\n' + history,
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license="BSD-3-Clause",
    url='https://github.com/molssi-seamm/psi4_step',

    # Which Python importable modules should be included when your package is
    # installed, handled automatically by setuptools. Use 'exclude' to prevent
    # some specific subpackage(s) from being added, if needed
    packages=find_packages(include=['psi4_step']),

    # Optional include package data to ship with your package. Customize
    # MANIFEST.in if the general case does not suit your needs. Comment out
    # this line to prevent the files from being packaged with your software
    include_package_data=True,

    # Allows `setup.py test` to work correctly with pytest
    setup_requires=[] + pytest_runner,

    # Required packages, pulls from pip if needed; do not use for Conda
    # deployment
    install_requires=requirements,

    test_suite='tests',

    # Valid platforms your code works on, adjust to your flavor
    platforms=['Linux',
               'Mac OS-X',
               'Unix',
               'Windows'],

    # Manual control if final package is compressible or not, set False to
    # prevent the .egg from being made
    zip_safe=True,

    keywords='psi4_step',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    entry_points={
        'org.molssi.seamm': [
            'Psi4 = psi4_step:Psi4Step',
        ],
        'org.molssi.seamm.tk': [
            'Psi4 = psi4_step:Psi4Step',
        ],
        'org.molssi.seamm.psi4': [
            'Energy = psi4_step:EnergyStep',
            'Initialization = psi4_step:InitializationStep',
        ],
        'org.molssi.seamm.psi4.tk': [
            'Energy = psi4_step:EnergyStep',
            'Initialization = psi4_step:InitializationStep',
        ],
    }
)
