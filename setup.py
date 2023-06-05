#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""psi4_step
A SEAMM plug-in to setup, run and analyze quantum chemistry calculations using Psi4
"""
import sys
from setuptools import setup, find_packages
import versioneer

short_description = __doc__.splitlines()[1]

# from https://github.com/pytest-dev/pytest-runner#conditional-requirement
needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)
pytest_runner = ['pytest-runner'] if needs_pytest else []

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

with open('requirements.txt') as fd:
    requirements = fd.read()

setup(
    name='psi4_step',
    author="Paul Saxe",
    author_email='psaxe@molssi.org',
    description=short_description,
    long_description=readme + '\n\n' + history,
    long_description_content_type='text/x-rst',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license="BSD-3-Clause",
    url='https://github.com/molssi-seamm/psi4_step',
    packages=find_packages(include=['psi4_step']),
    include_package_data=True,
    setup_requires=[] + pytest_runner,
    install_requires=requirements,
    test_suite='tests',
    platforms=['Linux',
               'Mac OS-X',
               'Unix',
               'Windows'],
    zip_safe=True,

    keywords=['SEAMM', 'plug-in', 'flowchart', 'Psi4', 'DFT',
              'quantum chemistry', 'CCSD'],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    entry_points={
        'console_scripts': [
            'psi4-step-installer=psi4_step.__main__:run',
        ],
        'org.molssi.seamm': [
            'Psi4 = psi4_step:Psi4Step',
        ],
        'org.molssi.seamm.tk': [
            'Psi4 = psi4_step:Psi4Step',
        ],
        'org.molssi.seamm.psi4': [
            'Thermochemistry = psi4_step:ThermochemistryStep',
            'Energy = psi4_step:EnergyStep',
            'Initialization = psi4_step:InitializationStep',
            'Optimization = psi4_step:OptimizationStep',
        ],
        'org.molssi.seamm.psi4.tk': [
            'Thermochemistry = psi4_step:ThermochemistryStep',
            'Energy = psi4_step:EnergyStep',
            'Initialization = psi4_step:InitializationStep',
            'Optimization = psi4_step:OptimizationStep',
        ],
    }
)