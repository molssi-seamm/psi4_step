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
            "enumeration": ("recommended", "advanced"),
            "description": "The level of disclosure in the interface",
            "help_text": (
                "How much detail to show in the GUI. Currently 'recommended' "
                "or 'advanced', which shows everything."
            ),
        },
        "method": {
            "default": "Kohn-Sham (KS) density functional theory (DFT)",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": [x for x in methods if methods[x]["level"] == "normal"],
            "format_string": "s",
            "description": "Method:",
            "help_text": ("The computational method to use."),
        },
        "advanced_method": {
            "default": "Kohn-Sham (KS) density functional theory (DFT)",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": [x for x in methods],
            "format_string": "s",
            "description": "Method:",
            "help_text": ("The computational method to use."),
        },
        "functional": {
            "default": "B3LYP Hyb-GGA Exchange-Correlation Functional",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": [
                x for x in dft_functionals if dft_functionals[x]["level"] == "normal"
            ],
            "format_string": "s",
            "description": "DFT Functional:",
            "help_text": ("The exchange-correlation functional to use."),
        },
        "advanced_functional": {
            "default": "B3LYP Hyb-GGA Exchange-Correlation Functional",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": [x for x in dft_functionals],
            "format_string": "s",
            "description": "DFT Functional:",
            "help_text": ("The exchange-correlation functional to use."),
        },
        "dispersion": {
            "default": "d3bj",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": ["none", "d3bj", "d3mbj", "nl"],
            "format_string": "s",
            "description": "Dispersion correction:",
            "help_text": ("The dispersion correction to use."),
        },
        "spin-restricted": {
            "default": "default",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": ("default", "yes", "no"),
            "format_string": "s",
            "description": "Spin-restricted:",
            "help_text": (
                "Whether to restrict the spin (RHF, ROHF, RKS) or not "
                "(UHF, UKS)."
                " Default is restricted for singlets, unrestricted otherwise."
            ),
        },
        "use damping": {
            "default": "no",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "s",
            "description": "Damp the iterations:",
            "help_text": (
                "Whether to damp the iterations using a fraction of the previous "
                "density to damp oscillations."
            ),
        },
        "damping percentage": {
            "default": 20.0,
            "kind": "float",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "s",
            "description": "Percent damping:",
            "help_text": "Percent of previous density to use to damp oscillations.",
        },
        "damping convergence": {
            "default": 0.0,
            "kind": "float",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "s",
            "description": "Damp to convergence:",
            "help_text": "Convergence level to stop damping. 0 = always damp.",
        },
        "use level shift": {
            "default": "no",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "s",
            "description": "Use a level shift:",
            "help_text": ("Whether to use a level shift to help convergence."),
        },
        "level shift": {
            "default": 5.0,
            "kind": "float",
            "default_units": "E_h",
            "enumeration": tuple(),
            "format_string": "s",
            "description": "Level shift:",
            "help_text": "The amount to shift the occupied orbitals down.",
        },
        "level shift convergence": {
            "default": 0.01,
            "kind": "float",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "s",
            "description": "Stop when converged to:",
            "help_text": "Convergence level to stop level shifting.",
        },
        "use soscf": {
            "default": "no",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "s",
            "description": "Use second order SCF:",
            "help_text": ("Whether to use the second order SCF to help convergence."),
        },
        "soscf starting convergence": {
            "default": 0.01,
            "kind": "float",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "s",
            "description": "Start when converged to:",
            "help_text": "Convergence level to start second order SCF.",
        },
        "soscf convergence": {
            "default": 0.001,
            "kind": "float",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "s",
            "description": "Microiteration convergence:",
            "help_text": "Convergence level for SOSCF microiterations.",
        },
        "soscf max iterations": {
            "default": 5,
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "s",
            "description": "Maximum microiterations:",
            "help_text": "Maximum number of SOSCF microiterations.",
        },
        "soscf print iterations": {
            "default": "no",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "s",
            "description": "Print microiterations:",
            "help_text": "Print information about the SOSCF microiterations.",
        },
        "density convergence": {
            "default": "default",
            "kind": "float",
            "default_units": "",
            "enumeration": ("default",),
            "format_string": "s",
            "description": "Density convergence criterion:",
            "help_text": (
                "Criterion for convergence of the density, default 10^-6 for "
                "SCF, 10^-8 optimization."
            ),
        },
        "energy convergence": {
            "default": "default",
            "kind": "float",
            "default_units": "",
            "enumeration": ("default",),
            "format_string": "s",
            "description": "Energy convergence criterion:",
            "help_text": (
                "Criterion for convergence of the energy, default 10^-6 for "
                "SCF, 10^-8 optimization."
            ),
        },
        "convergence error": {
            "default": "yes",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "s",
            "description": "Error if not converged:",
            "help_text": "Whether to throw an error if not converged.",
        },
        "maximum iterations": {
            "default": 100,
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "s",
            "description": "Maximum iterations:",
            "help_text": "Maximum number of SCF iterations.",
        },
        "freeze-cores": {
            "default": "yes",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "s",
            "description": "Freeze core orbitals:",
            "help_text": (
                "Whether to freeze the core orbitals in correlated " "methods"
            ),
        },
        "stability analysis": {
            "default": "no",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "s",
            "description": "Do stability analysis:",
            "help_text": ("Analyze the stability of the SCF/DFT wavefunction."),
        },
        "results": {
            "default": {},
            "kind": "dictionary",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "results",
            "help_text": ("The results to save to variables or in " "tables. "),
        },
        "create tables": {
            "default": "yes",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "",
            "description": "Create tables as needed:",
            "help_text": (
                "Whether to create tables as needed for "
                "results being saved into tables."
            ),
        },
    }

    output = {
        "density": {
            "default": "no",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "",
            "description": "Plot total density:",
            "help_text": "Whether to plot the total charge density.",
        },
        "orbitals": {
            "default": "no",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "",
            "description": "Plot orbitals:",
            "help_text": "Whether to plot orbitals.",
        },
        "selected orbitals": {
            "default": "-1, HOMO, LUMO, +1",
            "kind": "string",
            "default_units": "",
            "enumeration": ("all", "-1, HOMO, LUMO, +1"),
            "format_string": "",
            "description": "Selected orbitals:",
            "help_text": "Which orbitals to plot.",
        },
    }

    def __init__(self, defaults={}, data=None):
        """Initialize the instance, by default from the default
        parameters given in the class"""

        super().__init__(
            defaults={
                **EnergyParameters.parameters,
                **EnergyParameters.output,
                **defaults,
            },
            data=data,
        )
