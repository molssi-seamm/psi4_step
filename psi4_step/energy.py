# -*- coding: utf-8 -*-

"""Setup and run Psi4"""

import json
import logging
from pathlib import Path

from openbabel import openbabel

import psi4_step
import seamm
import seamm.data
from seamm_util import units_class
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter("psi4")


class Energy(seamm.Node):
    def __init__(self, flowchart=None, title="Energy", extension=None, logger=logger):
        """Initialize the node"""

        logger.debug("Creating Energy {}".format(self))

        super().__init__(
            flowchart=flowchart, title=title, extension=extension, logger=logger
        )

        self._calculation = "energy"
        self._model = None
        self._metadata = psi4_step.metadata
        self.parameters = psi4_step.EnergyParameters()

        self.description = "A single point energy calculation"

    @property
    def header(self):
        """A printable header for this section of output"""
        return "Step {}: {}".format(".".join(str(e) for e in self._id), self.title)

    @property
    def version(self):
        """The semantic version of this module."""
        return psi4_step.__version__

    @property
    def git_revision(self):
        """The git version of this module."""
        return psi4_step.__git_revision__

    def description_text(
        self,
        P=None,
        calculation_type="Single-point energy",
        configuration=None,
    ):
        """Prepare information about what this node will do"""

        if P is not None:
            P = self.parameters.values_to_dict()

        if P["level"] == "recommended":
            method = P["method"]
        else:
            method = P["advanced_method"]

        if self.is_expr(method):
            text = f"{calculation_type} using method given by {method}."
        elif (
            method in psi4_step.methods and psi4_step.methods[method]["method"] == "dft"
        ):
            if P["level"] == "recommended":
                functional = P["functional"]
            else:
                functional = P["advanced_functional"]
            text = f"{calculation_type} using {method} with an "
            text += f"exchange-correlation potential of {functional}"
            if (
                functional in psi4_step.dft_functionals
                and len(psi4_step.dft_functionals[functional]["dispersion"]) > 1
                and P["dispersion"] != "none"
            ):
                text += f" with the {P['dispersion']} dispersion correction."
            else:
                text += " with no dispersion correction."
        else:
            text = f"{calculation_type} using {method}."

        # Spin
        if P["spin-restricted"] == "yes":
            text += " The spin will be restricted to a pure eigenstate."
        elif P["spin-restricted"] == "default":
            if configuration is not None:
                if configuration.spin_multiplicity == 1:
                    text += " The spin will be restricted to a pure eigenstate."
                else:
                    text += " The spin will not be restricted and may not be a "
                    text += "proper eigenstate."
            else:
                text += " The spin will be restricted to a pure "
                text += "eigenstate for singlet states. Otherwise it will not "
                text += "be restricted and may not be a proper eigenstate."
        elif self.is_expr(P["spin-restricted"]):
            text += " Whether the spin will be restricted to a pure "
            text += "eigenstate will be determined by {P['spin-restricted']}"
        else:
            text += " The spin will not be restricted and may not be a "
            text += "proper eigenstate."

        # Plotting
        if isinstance(P["density"], str):
            density = P["density"] != "no"
        else:
            density = P["density"]

        if isinstance(P["orbitals"], str):
            orbitals = P["orbitals"] != "no"
        else:
            orbitals = P["orbitals"]

        if density:
            if orbitals:
                text += "\nThe alpha and beta electron, total, and spin densities, "
                text += f"and orbitals {P['selected orbitals']} will be plotted."
        elif orbitals:
            text += f"\nThe orbitals {P['selected orbitals']} will be plotted."

        return self.header + "\n" + __(text, **P, indent=4 * " ").__str__()

    def get_input(self, calculation_type="energy", restart=None):
        """Get the input for an energy calculation for Psi4"""
        # Create the directory
        directory = Path(self.directory)
        directory.mkdir(parents=True, exist_ok=True)

        _, configuration = self.get_system_configuration(None)

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
        if calculation_type == "energy":
            lines.append("")
            lines.append("#" * 80)
            lines.append(f"# {self.header}")
            lines.append("#" * 80)

        # Figure out what we are doing!
        if P["level"] == "recommended":
            method_string = P["method"]
        else:
            method_string = P["advanced_method"]

        # Allow the full name, or the short name, or just pray.
        if method_string in psi4_step.methods:
            method = psi4_step.methods[method_string]["method"]
        else:
            method = method_string.lower()
            for key in psi4_step.methods:
                if psi4_step.methods[key]["method"] == method:
                    break
            else:
                method = method_string

        # lines.append('set scf_type df')
        # lines.append('set guess sad')

        multiplicity = configuration.spin_multiplicity
        spin_restricted = P["spin-restricted"]
        if spin_restricted == "default":
            if multiplicity == 1:
                restricted = True
            else:
                restricted = False
        elif spin_restricted == "yes":
            restricted = True
        else:
            restricted = False

        if method == "dft":
            if restricted:
                lines.append("set reference rks")
            else:
                lines.append("set reference uks")
        else:
            if restricted:
                if multiplicity == 1:
                    lines.append("set reference rhf")
                else:
                    lines.append("set reference rohf")
            else:
                lines.append("set reference uhf")

        if "freeze-cores" in P and P["freeze-cores"] == "yes":
            lines.append("set freeze_core True")
        else:
            lines.append("set freeze_core False")

        if P["stability analysis"]:
            lines.append("set stability_analysis True")

        # Convergence parameters and methods
        lines.append("set fail_on_maxiter False")  # Psi4 hangs! Cope in the plug-in

        lines.append(f"set maxiter {P['maximum iterations']}")

        if P["density convergence"] != "default":
            lines.append(f"set d_convergence {P['density convergence']}")
        if P["energy convergence"] != "default":
            lines.append(f"set e_convergence {P['energy convergence']}")

        if P["use damping"]:
            lines.append(f"set damping_percentage {P['damping percentage']}")
            lines.append(f"set damping_convergence {P['damping convergence']}")

        if P["use level shift"]:
            lines.append(f"set level_shift {P['level shift']}")
            lines.append(f"set level_shift_cutoff {P['level shift convergence']}")

        if P["use soscf"]:
            lines.append("set soscf True")
            lines.append(
                f"set soscf_start_convergence {P['soscf starting convergence']}"
            )
            lines.append(f"set soscf_conv {P['soscf convergence']}")
            lines.append(f"set soscf_max_iter {P['soscf max iterations']}")
            if P["soscf print iterations"]:
                lines.append("set soscf_print True")

        lines.append("")
        if method == "dft":
            if P["level"] == "recommended":
                functional_string = P["functional"]
            else:
                functional_string = P["advanced_functional"]

            # Allow the full name, or the short name, or just pray.
            if functional_string in psi4_step.dft_functionals:
                functional = psi4_step.dft_functionals[functional_string]["name"]
            else:
                functional = functional_string.lower()
                for key in psi4_step.dft_functionals:
                    if psi4_step.dft_functionals[key]["name"] == functional:
                        break
                else:
                    functional = functional_string

            if (
                P["dispersion"] != "none"
                and len(psi4_step.dft_functionals[functional_string]["dispersion"]) > 1
            ):
                functional = functional + "-" + P["dispersion"]
            if restart is None:
                lines.append(
                    f"Eelec, wfn = {calculation_type}('{functional}', return_wfn=True)"
                )
            else:
                if calculation_type == "gradient":
                    lines.append(
                        f"Eelec, wfn = energy('{functional}', return_wfn=True,"
                        f" restart_file='{restart}')"
                    )
                    lines.append(f"G = gradient('{functional}', ref_wfn=wfn)")
                else:
                    lines.append(
                        f"Eelec, wfn = {calculation_type}('{functional}', "
                        f"return_wfn=True, restart_file='{restart}')"
                    )
        else:
            if restart is None:
                lines.append(
                    f"Eelec, wfn = {calculation_type}('{method}', return_wfn=True)"
                )
            else:
                if calculation_type == "gradient":
                    lines.append(
                        f"Eelec, wfn = energy('{method}', return_wfn=True, "
                        f" restart_file='{restart}')"
                    )
                    lines.append(f"G = gradient('{method}', ref_wfn=wfn)")
                else:
                    lines.append(
                        f"Eelec, wfn = {calculation_type}('{method}', return_wfn=True, "
                        f" restart_file='{restart}')"
                    )

        if calculation_type != "gradient":
            # Dump the properties to a json file
            filename = f"@{self._id[-1]}+properties.json"
            lines.append(
                f"""
try:
    oeprop(
        wfn,
        "MULTIPOLE(5)",
        "ESP_AT_NUCLEI",
        "LOWDIN_CHARGES",
        "MULLIKEN_CHARGES",
        "WIBERG_LOWDIN_INDICES",
        "MAYER_INDICES",
        "NO_OCCUPATIONS",
        title="PROP"
    )
except Exception:
    print("Failed to calculate properties")

variables = scalar_variables()
variables.update(wfn.scalar_variables())
arrays = array_variables()
for item in arrays:
    variables[item] = array_variable(item).np.tolist()
arrays = wfn.array_variables()
for item in arrays:
    variables[item] = wfn.array_variable(item).np.tolist()
variables["Eelec"] = Eelec
variables["_method"] = "{method}"
variables["_method_string"] = "{method_string}"


with open("{filename}", "w") as fd:
    json.dump(fix_multipoles(variables), fd, sort_keys=True, indent=3)
"""
            )

        # Orbital plots
        lines.append(self.plot_input())

        return "\n".join(lines)

    def analyze(self, indent="", data={}, out=[]):
        """Parse the output and generating the text output and store the
        data in variables for other stages to access
        """
        # Read in the results from json
        directory = Path(self.directory)
        json_file = directory / "properties.json"
        if json_file.exists():
            with json_file.open() as fd:
                data = json.load(fd)

            # Put any requested results into variables or tables
            self.store_results(
                data=data,
                create_tables=self.parameters["create tables"].get(),
            )

            text = "The calculated energy is {Eelec:.6f} E_h."
        else:
            data = {}
            tmp = str(json_file)
            text = (
                "\nThere are no results from Psi4. Perhaps it "
                f"failed? Looking for {tmp}."
            )
            printer.normal(__(text, **data, indent=self.indent + 4 * " "))
            raise RuntimeError(text)

        # Write the structure locally for use in density and orbital plots
        system, configuration = self.get_system_configuration()
        obConversion = openbabel.OBConversion()
        obConversion.SetOutFormat("sdf")
        obMol = configuration.to_OBMol(properties="all")
        title = f"SEAMM={system.name}/{configuration.name}"
        obMol.SetTitle(title)
        sdf = obConversion.WriteString(obMol)
        path = directory / "structure.sdf"
        path.write_text(sdf)

        printer.normal(__(text, **data, indent=self.indent + 4 * " "))

    def plot_input(self):
        """Generate the input for plotting to cube files."""
        _, configuration = self.get_system_configuration(None)

        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )

        tasks = []
        if P["density"]:
            tasks.append("density")

        orbitals = []
        if P["orbitals"]:
            tasks.append("orbitals")
            # and work out the orbitals
            txt = P["selected orbitals"]
            if txt == "all":
                pass
            else:
                # Which is the HOMO orbital?
                # This will not work with ECPs.
                n_electrons = (
                    sum(configuration.atoms.atomic_numbers) - configuration.charge
                )
                multiplicity = configuration.spin_multiplicity
                homo = (n_electrons - (multiplicity - 1)) // 2 + (multiplicity - 1)

                for chunk in txt.split(","):
                    chunk = chunk.strip()
                    if ":" in chunk or ".." in chunk:
                        if ":" in chunk:
                            first, last = chunk.split(":")
                        elif ".." in chunk:
                            first, last = chunk.split("..")
                        first = first.strip().upper()
                        last = last.strip().upper()

                        if first == "HOMO":
                            first = homo
                        elif first == "LUMO":
                            first = homo + 1
                        else:
                            first = int(first.removeprefix("HOMO").removeprefix("LUMO"))
                            if first < 0:
                                first = homo + first
                            else:
                                first = homo + 1 + first

                        if last == "HOMO":
                            last = homo
                        elif last == "LUMO":
                            last = homo + 1
                        else:
                            last = int(last.removeprefix("HOMO").removeprefix("LUMO"))
                            if last < 0:
                                last = homo + last
                            else:
                                last = homo + 1 + last

                        orbitals.extend(range(first, last + 1))
                    else:
                        first = chunk.strip().upper()

                        if first == "HOMO":
                            first = homo
                        elif first == "LUMO":
                            first = homo + 1
                        else:
                            first = int(first.removeprefix("HOMO").removeprefix("LUMO"))
                            if first < 0:
                                first = homo + first
                            else:
                                first = homo + 1 + first
                        orbitals.append(first)

        if len(tasks) == 0:
            return ""

        lines = []
        txt = "', '".join(tasks)
        lines.append("")
        lines.append("# Cube files for density and orbitals")
        lines.append(f"set cubeprop_tasks ['{txt}']")
        if len(orbitals) > 0:
            txt = ", ".join([f"{i}, {-i}" for i in orbitals])
            lines.append(f"set cubeprop_orbitals [{txt}]")

        lines.append("")
        lines.append("cubeprop(wfn)")
        lines.append(
            f"""
# Prefix the files with the substep number
paths = Path.cwd().glob('*.cube')
for path in paths:
    name = path.name
    newpath = path.with_name('@{self._id[-1]}+' + name)
    path.rename(newpath)
"""
        )

        return "\n".join(lines)
