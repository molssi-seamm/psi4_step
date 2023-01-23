# -*- coding: utf-8 -*-
"""Global control parameters for Psi4
"""

import logging

import psi4_step
from seamm.standard_parameters import structure_handling_parameters

# import seamm

logger = logging.getLogger(__name__)


class AcceleratedOptimizationParameters(psi4_step.EnergyParameters):
    """The control parameters for the MOPAC accelerated optimization."""

    parameters = {
        "optimization method": {
            "default": "RFO",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": (
                "RFO",
                "P_RFO",
                "NR",
                "SD",
                "LINESEARCH",
            ),
            "format_string": "s",
            "description": "Method:",
            "help_text": "The optimization method to use.",
        },
        "max geometry steps": {
            "default": "default",
            "kind": "integer",
            "default_units": "",
            "enumeration": ("default",),
            "format_string": "",
            "description": "Maximum steps:",
            "help_text": (
                "The maximum number of steps to take in the optimization. "
                "'default' is based on the system size, giving a reasonable "
                "limit in most cases."
            ),
        },
        "geometry convergence": {
            "default": "QCHEM",
            "kind": "float",
            "default_units": "",
            "enumeration": (
                "QCHEM",
                "MOLPRO",
                "GAU",
                "GAU_LOOSE",
                "GAU_TIGHT",
                "GAU_VERYTIGHT",
                "TURBOMOLE",
                "CFOUR",
                "NWCHEM_LOOSE",
            ),
            "format_string": "",
            "description": "Convergence criteria:",
            "help_text": "The criteria to use for convergence.",
        },
        "recalc hessian": {
            "default": "every step",
            "kind": "integer",
            "default_units": "",
            "enumeration": ("every step", "at beginning", "never"),
            "format_string": "",
            "description": "Recalculate Hessian:",
            "help_text": (
                "How often to recalculate the Hessian (in steps). Smaller "
                "values help convergence but are expensive."
            ),
        },
        "hessian update": {
            "default": "bfgs",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": (
                "bfgs",
                "ms",
                "powell",
                "none",
            ),
            "format_string": "s",
            "description": "Hessian update:",
            "help_text": "The algorithm for updating the Hessian.",
        },
        "coordinates": {
            "default": "Internal",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": (
                "Internal",
                "Delocalized",
                "Natural",
                "Cartesian",
                "Both",
            ),
            "format_string": "s",
            "description": "Type of coordinates:",
            "help_text": "The type of coordinates to use in the minimization.",
        },
    }

    def __init__(self, defaults={}, data=None):
        """Initialize the instance, by default from the default
        parameters given in the class"""

        super().__init__(
            defaults={
                **AcceleratedOptimizationParameters.parameters,
                **structure_handling_parameters,
                **defaults,
            },
            data=data,
        )
