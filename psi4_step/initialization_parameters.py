# -*- coding: utf-8 -*-
"""Global control parameters for Psi4"""

import logging
import seamm

logger = logging.getLogger(__name__)


class InitializationParameters(seamm.Parameters):
    """The control parameters for initializing Psi4"""

    parameters = {
        "basis source": {
            "default": "library",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": ("library", "basis file"),
            "format_string": "s",
            "description": "Basis set from:",
            "help_text": (
                "Where the basis set comes from: Psi4's library, by name, or a "
                "basis file of Psi4 basis blocks, which can also give the "
                "density-fitting basis sets."
            ),
        },
        "basis file": {
            "applies_when": {"basis source": "basis file"},
            "default": "cc-pwCVTZ-ORCA",
            "kind": "string",
            "default_units": "",
            "enumeration": ("cc-pwCVTZ-ORCA",),
            "format_string": "s",
            "description": "Basis file:",
            "help_text": (
                "A basis file of Psi4 input blocks (basis {...}, and optionally "
                "df_basis_cc {...} and df_basis_scf {...}), copied into the input. "
                "Either the name of a file shipped with the step (e.g. "
                "cc-pwCVTZ-ORCA, ORCA 6.1.1's cc-pwCVTZ with explicit fitting "
                "sets) or the path of a file on the machine running the flowchart. "
                "The name of the file, without its extension, is the basis set's "
                "name in the model."
            ),
        },
        "basis": {
            "applies_when": {"basis source": "library"},
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
