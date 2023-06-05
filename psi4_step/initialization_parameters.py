# -*- coding: utf-8 -*-
"""Global control parameters for Psi4
"""

import logging
import seamm

logger = logging.getLogger(__name__)


class InitializationParameters(seamm.Parameters):
    """The control parameters for initializing Psi4"""

    parameters = {
        "basis": {
            "default": "6-31G**",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": (
                "6-31G",
                "6-31G(d)",
                "6-31G(d,p)",
                "6-31+G",
                "6-31+G(d)",
                "6-31+G(d,p)",
                "6-311G",
                "6-311G(d)",
                "6-311G(d,p)",
                "6-311+G",
                "6-311+G(d)",
                "6-311+G(d,p)",
                "cc-pVDZ",
                "cc-pVTZ",
                "cc-pVQZ",
                "def2-SV(P)",
                "def2-SVP",
                "def2-TZVP",
                "def2-TZVPP",
                "def2-QZVP",
                "def2-QZVPP",
            ),
            "format_string": "s",
            "description": "Basis:",
            "help_text": ("The basis set to use."),
        },
        "symmetry_tolerance": {
            "default": "0.00001",
            "kind": "float",
            "default_units": None,
            "enumeration": tuple(),
            "format_string": ".1e",
            "description": "Symmetry tolerance:",
            "help_text": "The tolerance used when determining the symmetry.",
        },
    }

    def __init__(self, defaults={}, data=None):
        """Initialize the instance, by default from the default
        parameters given in the class"""

        super().__init__(
            defaults={**InitializationParameters.parameters, **defaults}, data=data
        )
