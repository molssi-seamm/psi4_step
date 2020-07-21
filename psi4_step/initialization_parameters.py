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
                '6-31G',
                '6-31G*',
                '6-31G**',
                'cc-pVDZ',
                'cc-pVTZ',
                'cc-pVQZ',
            ),
            "format_string": "s",
            "description": "Basis:",
            "help_text": ("The basis set to use.")
        },
        "symmetry_tolerance": {
            "default": "0.05",
            "kind": "float",
            "default_units": None,
            "enumeration": tuple(),
            "format_string": ".1e",
            "description": "Symmetry tolerance:",
            "help_text": "The tolerance used when determining the symmetry."
        },
    }

    def __init__(self, defaults={}, data=None):
        """Initialize the instance, by default from the default
        parameters given in the class"""

        super().__init__(
            defaults={**InitializationParameters.parameters, **defaults},
            data=data
        )
