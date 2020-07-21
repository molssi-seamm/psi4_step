# -*- coding: utf-8 -*-
"""Global control parameters for Psi4
"""

import logging

from psi4_step import methods, dft_functionals
import seamm

logger = logging.getLogger(__name__)


class EnergyParameters(seamm.Parameters):
    """The control parameters for the energy."""

    parameters = {
        "level": {
            "default": "recommended",
            "kind": "string",
            "format_string": "s",
            "enumeration": (
                "recommended",
                "advanced"
            ),
            "description": "The level of disclosure in the interface",
            "help_text": (
                "How much detail to show in the GUI. Currently 'recommended' "
                "or 'advanced', which shows everything."
            )
        },
        "method": {
            "default": "Kohn-Sham (KS) density functional theory (DFT)",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": [
                x for x in methods if methods[x]['level'] == 'normal'
            ],
            "format_string": "s",
            "description": "Method:",
            "help_text": ("The computational method to use.")
        },
        "advanced_method": {
            "default": "Kohn-Sham (KS) density functional theory (DFT)",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": [x for x in methods],
            "format_string": "s",
            "description": "Method:",
            "help_text": ("The computational method to use.")
        },
        "functional": {
            "default": "B3LYP Hyb-GGA Exchange-Correlation Functional",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": [
                x for x in dft_functionals
                if dft_functionals[x]['level'] == 'normal'
            ],
            "format_string": "s",
            "description": "DFT Functional:",
            "help_text": ("The exchange-correlation functional to use.")
        },
        "advanced_functional": {
            "default": "B3LYP Hyb-GGA Exchange-Correlation Functional",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": [x for x in dft_functionals],
            "format_string": "s",
            "description": "DFT Functional:",
            "help_text": ("The exchange-correlation functional to use.")
        },
        "dispersion": {
            "default": "d3bj",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": ['none', 'd3bj', 'd3mbj', 'nl'],
            "format_string": "s",
            "description": "Dispersion correction:",
            "help_text": ("The dispersion correction to use.")
        },
        "spin-restricted": {
            "default": "no",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": ('yes', 'no'),
            "format_string": "s",
            "description": "Spin-restricted:",
            "help_text": (
                "Whether to restrict the spin (RHF, ROHF, RKS) or not "
                "(UHF, UKS)"
            )
        },
        "results": {
            "default": {},
            "kind": "dictionary",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "results",
            "help_text": ("The results to save to variables or in "
                          "tables. ")
        },
        "create tables": {
            "default": "yes",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ('yes', 'no'),
            "format_string": "",
            "description": "Create tables as needed:",
            "help_text": ("Whether to create tables as needed for "
                          "results being saved into tables.")
        },
    }

    def __init__(self, defaults={}, data=None):
        """Initialize the instance, by default from the default
        parameters given in the class"""

        super().__init__(
            defaults={**EnergyParameters.parameters, **defaults},
            data=data
        )
