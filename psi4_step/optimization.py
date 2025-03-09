# -*- coding: utf-8 -*-

"""Setup and run Psi4"""

import json
import logging
from pathlib import Path

import cclib

from molsystem import RMSD
import psi4_step
import seamm
import seamm.data
from seamm_util import units_class
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __


try:
    from itertools import batched
except ImportError:
    from itertools import islice

    def batched(iterable, n):
        "Batch data into tuples of length n. The last batch may be shorter."
        # batched('ABCDEFG', 3) --> ABC DEF G
        if n < 1:
            raise ValueError("n must be at least one")
        it = iter(iterable)
        while batch := tuple(islice(it, n)):
            yield batch


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

    def analyze(self, indent="", data=None, out=[], table=None):
        """Parse the output and generating the text output and store the
        data in variables for other stages to access
        """
        directory = Path(self.directory)
        text = ""

        # Get the results
        if data is None:
            # Use cclib to get results. Note the file is one level up.
            path = directory.parent / "output.dat"
            if path.exists():
                data = vars(cclib.io.ccread(path))
                data = self.process_data(data)
            else:
                data = {}

            # Read in the results from json
            directory = Path(self.directory)
            json_file = directory / "properties.json"
            if not json_file.exists():
                data = {}
                tmp = str(json_file)
                text = (
                    "\nThere are no results from Psi4. Perhaps it "
                    f"failed? Looking for {tmp}."
                )
                printer.normal(__(text, **data, indent=self.indent + 4 * " "))
                raise RuntimeError(text)

            with json_file.open() as fd:
                tmp = json.load(fd)
            data.update(**tmp)

        # Get the structure, if changed
        method, functional, extended_functional, method_string = self.get_method()

        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )
        _, initial_configuration = self.get_system_configuration()
        system, configuration = self.get_system_configuration(P)

        update_structure = P["structure handling"] != "Discard the structure"

        initial = initial_configuration.to_RDKMol()
        if update_structure:
            final = configuration.to_RDKMol()
        else:
            final = initial_configuration.to_RDKMol()

        # Get the updated coordinates
        structure_file = directory / "final_structure.json"
        if structure_file.exists():
            with structure_file.open(mode="r") as fd:
                structure = json.load(fd)
            if "geom" in structure:
                XYZs = [batched(structure["geom"], 3)]
                final.GetConformer(0).SetPositions(XYZs)

        result = RMSD(final, initial, symmetry=True, include_h=True, align=True)
        data["RMSD with H"] = result["RMSD"]
        data["displaced atom with H"] = result["displaced atom"]
        data["maximum displacement with H"] = result["maximum displacement"]

        # Align the structure
        if update_structure:
            configuration.from_RDKMol(final)

            # And the name of the configuration.
            text = seamm.standard_parameters.set_names(
                system,
                configuration,
                P,
                _first=True,
                model=self.parent.model,
            )
            printer.normal(__(text, **data, indent=self.indent + 4 * " "))
            printer.normal("")

        result = RMSD(final, initial, symmetry=True)
        data["RMSD"] = result["RMSD"]
        data["displaced atom"] = result["displaced atom"]
        data["maximum displacement"] = result["maximum displacement"]

        # Create the optimization part of the output
        if table is None:
            table = {
                "Property": [],
                "Value": [],
                "Units": [],
            }

        if "optdone" in data:
            table["Property"].append("Converged")
            table["Value"].append(str(data["optdone"]))
            table["Units"].append("")
        else:
            table["Property"].append("Converged")
            table["Value"].append("unknown")
            table["Units"].append("")

        if "optstatus" in data:
            table["Property"].append("# Steps")
            table["Value"].append(str(len(data["optstatus"])))
            table["Units"].append("")
        else:
            table["Property"].append("# Steps")
            table["Value"].append("unknown")
            table["Units"].append("")

        if "RMSD" in data:
            tmp = data["RMSD"]
            table["Property"].append("RMSD in Geometry")
            table["Value"].append(f"{tmp:.2f}")
            table["Units"].append("Å")

        if "maximum displacement" in data:
            tmp = data["maximum displacement"]
            table["Property"].append("Largest Displacement")
            table["Value"].append(f"{tmp:.2f}")
            table["Units"].append("Å")

        if "displaced atom" in data:
            tmp = data["displaced atom"]
            table["Property"].append("Displaced Atom")
            table["Value"].append(f"{tmp + 1}")
            table["Units"].append("")

        super().analyze(indent=indent, data=data, out=out, table=table)

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

        if self.parent.model is None:
            kwargs = {"model": "<model>"}
        else:
            kwargs = {"model": self.model}
        added += seamm.standard_parameters.structure_handling_description(P, **kwargs)

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
