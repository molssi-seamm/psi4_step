# -*- coding: utf-8 -*-
"""
Control parameters for the Thermochemistry step in a Psi4 flowchart
"""

import logging

from .energy_parameters import EnergyParameters

logger = logging.getLogger(__name__)


class ThermochemistryParameters(EnergyParameters):
    """
    The control parameters for Thermochemistry.

    You need to replace the "time" entry in dictionary below these comments with the
    definitions of parameters to control this step. The keys are parameters for the
    current plugin,the values are dictionaries as outlined below.

    Examples
    --------
    ::

        parameters = {
            "time": {
                "default": 100.0,
                "kind": "float",
                "default_units": "ps",
                "enumeration": tuple(),
                "format_string": ".1f",
                "description": "Simulation time:",
                "help_text": ("The time to simulate in the dynamics run.")
            },
        }

    parameters : {str: {str: str}}
        A dictionary containing the parameters for the current step.
        Each key of the dictionary is a dictionary that contains the
        the following keys:

    parameters["default"] :
        The default value of the parameter, used to reset it.

    parameters["kind"] : enum()
        Specifies the kind of a variable. One of  "integer", "float", "string",
        "boolean", or "enum"

        While the "kind" of a variable might be a numeric value, it may still have
        enumerated custom values meaningful to the user. For instance, if the parameter
        is a convergence criterion for an optimizer, custom values like "normal",
        "precise", etc, might be adequate. In addition, any parameter can be set to a
        variable of expression, indicated by having "$" as the first character in the
        field. For example, $OPTIMIZER_CONV.

    parameters["default_units"] : str
        The default units, used for resetting the value.

    parameters["enumeration"]: tuple
        A tuple of enumerated values.

    parameters["format_string"]: str
        A format string for "pretty" output.

    parameters["description"]: str
        A short string used as a prompt in the GUI.

    parameters["help_text"]: str
        A longer string to display as help for the user.

    See Also
    --------
    Thermochemistry, TkThermochemistry, Thermochemistry
    ThermochemistryParameters, ThermochemistryStep
    """

    parameters = {
        "use existing parameters": {
            "default": "yes",
            "kind": "boolean",
            "default_units": None,
            "enumeration": ("yes", "no"),
            "format_string": "",
            "description": "Use existing parameters",
            "help_text": "Whether to use existing parameters from the previous step",
        },
        "T": {
            "default": 298.15,
            "kind": "float",
            "default_units": "K",
            "enumeration": tuple(),
            "format_string": ".2f",
            "description": "Temperature:",
            "help_text": "The temperature for the thermochemistry.",
        },
        "P": {
            "default": 1.0,
            "kind": "float",
            "default_units": "bar",
            "enumeration": tuple(),
            "format_string": ".2f",
            "description": "Pressure:",
            "help_text": "The pressure for the thermochemistry.",
        },
        # # Results handling ... uncomment if needed
        # "results": {
        #     "default": {},
        #     "kind": "dictionary",
        #     "default_units": "",
        #     "enumeration": tuple(),
        #     "format_string": "",
        #     "description": "results",
        #     "help_text": "The results to save to variables or in tables.",
        # },
    }

    def __init__(self, defaults={}, data=None):
        """
        Initialize the parameters, by default with the parameters defined above

        Parameters
        ----------
        defaults: dict
            A dictionary of parameters to initialize. The parameters
            above are used first and any given will override/add to them.
        data: dict
            A dictionary of keys and a subdictionary with value and units
            for updating the current, default values.

        Returns
        -------
        None
        """

        logger.debug("ThermochemistryParameters.__init__")

        super().__init__(
            defaults={
                **self.energy_settings(),
                **ThermochemistryParameters.parameters,
                **defaults,
            },
            data=data,
        )

    @staticmethod
    def energy_settings():
        """The Energy step's calculation settings, which apply only when not using
        the parameters of the previous step (a rule shared by the dialog and the
        flowchart builder, see seamm.Parameters)."""
        result = {}
        for key, definition in EnergyParameters.parameters.items():
            if key in ("results", "create tables"):
                continue
            conditions = {
                "use existing parameters": "no",
                **definition.get("applies_when", {}),
            }
            result[key] = {**definition, "applies_when": conditions}
        return result

    def not_applicable_reason(self, key, values=None):
        """Why a parameter does not apply: first, because the previous step's
        parameters are used."""
        if values is None:
            values = self.current_values()
        if (
            key in self.energy_settings()
            and values.get("use existing parameters") == "yes"
        ):
            return (
                "the parameters of the previous step are used ('use existing "
                "parameters' is 'yes')"
            )
        return super().not_applicable_reason(key, values)
