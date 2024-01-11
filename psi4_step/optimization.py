# -*- coding: utf-8 -*-

"""Setup and run Psi4"""

import logging
from pathlib import Path

import psi4_step
import seamm
import seamm.data
from seamm_util import units_class
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter("psi4")


class Optimization(psi4_step.Energy):
    def __init__(
        self,
        flowchart=None,
        title="Optimization",
        extension=None,
        logger=logger,
    ):
        """Initialize the node"""

        logger.debug("Creating Optimization {}".format(self))

        super().__init__(flowchart=flowchart, title=title, extension=extension)

        self._calculation = "optimization"
        self._model = None
        self._metadata = psi4_step.metadata
        self.parameters = psi4_step.OptimizationParameters()

        self.description = "A geometry optimization"

    def description_text(
        self,
        P=None,
        calculation_type="Geometry optimization",
        configuration=None,
    ):
        """Prepare information about what this node will do"""

        if not P:
            P = self.parameters.values_to_dict()

        text = super().description_text(
            P=P, calculation_type=calculation_type, configuration=configuration
        )

        added = "\nThe geometry optimization will use the {optimization method} "
        if P["max geometry steps"] == "default":
            added += "method, using the default maximum number of steps, which"
            added += " is based on the system size."
        else:
            added += "method, with no more than {max geometry steps} steps."

        if P["geometry convergence"] == "Custom":
            added += " The convergence criterion is"
        else:
            added += " The convergence criterion is '{geometry convergence}'."

        if P["recalc hessian"] != "never":
            added += " The Hessian will be recalculated every {recalc hessian}"
            added += " steps. Note that calculating the second derivatives is "
            added += "quite expensive!"

        return text + "\n" + __(added, **P, indent=4 * " ").__str__()

    def get_input(self, calculation_type="optimize"):
        """Get the input for an optimization calculation for Psi4"""
        _, configuration = self.get_system_configuration()

        # Create the directory
        directory = Path(self.directory)
        directory.mkdir(parents=True, exist_ok=True)

        # Check for atoms and bypass optimization
        if configuration.n_atoms == 1:
            return super().get_input()

        # references = self.parent.references

        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )
        # Have to fix formatting for printing...
        PP = dict(P)
        for key in PP:
            if isinstance(PP[key], units_class):
                PP[key] = "{:~P}".format(PP[key])

        self.description = []
        self.description.append(
            __(
                self.description_text(PP, configuration=configuration),
                **PP,
                indent=self.indent,
            )
        )

        lines = []
        lines.append("")
        lines.append("#" * 80)
        lines.append(f"# {self.header}")
        lines.append("#" * 80)
        lines.append("")
        lines.append("set opt_type min")
        if P["optimization method"] in psi4_step.optimization_methods:
            opt_method = psi4_step.optimization_methods[P["optimization method"]]
        else:
            opt_method = P["optimization method"]
        lines.append(f"set step_type {opt_method}")
        lines.append(f"set opt_coordinates {P['coordinates']}")
        max_steps = P["max geometry steps"]
        if max_steps == "default":
            max_steps = 6 * configuration.n_atoms
        if "nAtoms" in max_steps:
            n_atoms = configuration.n_atoms
            max_steps = max_steps.replace("nAtoms", str(n_atoms))
            max_steps = eval(max_steps)

        lines.append(f"set geom_maxiter {max_steps}")
        if P["geometry convergence"] == "Custom":
            pass
        else:
            if P["geometry convergence"] in psi4_step.optimization_convergence:
                convergence = psi4_step.optimization_convergence[
                    P["geometry convergence"]
                ]
            else:
                convergence = P["geometry convergence"]
            lines.append(f"set g_convergence {convergence}")

        if P["recalc hessian"] == "every step":
            lines.append("set full_hess_every 1")
        elif P["recalc hessian"] == "at beginning":
            lines.append("set full_hess_every 0")
        elif P["recalc hessian"] == "never":
            lines.append("set full_hess_every -1")
        else:
            lines.append(f"set full_hess_every {P['recalc hessian']}")
        lines.append(f"set hess_update {P['hessian update']}")
        lines.append("")

        # Add in the input from the energy part of things
        lines.append(super().get_input(calculation_type=calculation_type))

        return "\n".join(lines)
