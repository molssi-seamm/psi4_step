# -*- coding: utf-8 -*-
"""Global control parameters for Psi4
"""

import logging

import psi4_step

# import seamm

logger = logging.getLogger(__name__)


class OptimizationParameters(psi4_step.EnergyParameters):
    """The control parameters for the energy."""

    parameters = {
        "optimization method": {
            "default": "Rational Function Optimization",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": [x for x in psi4_step.optimization_methods],
            "format_string": "s",
            "description": "Method:",
            "help_text": ("The optimization method to use."),
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
            "default": "QChem",
            "kind": "string",
            "default_units": "",
            "enumeration": [x for x in psi4_step.optimization_convergence],
            "format_string": "",
            "description": "Convergence criteria:",
            "help_text": "The criteria to use for convergence.",
        },
        "recalc hessian": {
            "default": "never",
            "kind": "integer",
            "default_units": "",
            "enumeration": ("never",),
            "format_string": "",
            "description": "Recalculate Hessian:",
            "help_text": (
                "How often to recalculate the Hessian (in steps). Smaller "
                "values help convergence but are expensive."
            ),
        },
    }

    def __init__(self, defaults={}, data=None):
        """Initialize the instance, by default from the default
        parameters given in the class"""

        super().__init__(
            defaults={**OptimizationParameters.parameters, **defaults}, data=data
        )
