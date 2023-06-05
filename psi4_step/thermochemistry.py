# -*- coding: utf-8 -*-

"""Non-graphical part of the Thermochemistry step in a Psi4 flowchart
"""

import json
import logging
from pathlib import Path
import pkg_resources

from tabulate import tabulate

import psi4_step
import molsystem
import seamm
from seamm_util import Q_, units_class
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

# In addition to the normal logger, two logger-like printing facilities are
# defined: "job" and "printer". "job" send output to the main job.out file for
# the job, and should be used very sparingly, typically to echo what this step
# will do in the initial summary of the job.
#
# "printer" sends output to the file "step.out" in this steps working
# directory, and is used for all normal output from this step.

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter("Psi4")

# Add this module's properties to the standard properties
path = Path(pkg_resources.resource_filename(__name__, "data/"))
csv_file = path / "properties.csv"
if path.exists():
    molsystem.add_properties_from_file(csv_file)


class Thermochemistry(psi4_step.Energy):
    """
    The non-graphical part of a Thermochemistry step in a flowchart.

    Attributes
    ----------
    parser : configargparse.ArgParser
        The parser object.

    options : tuple
        It contains a two item tuple containing the populated namespace and the
        list of remaining argument strings.

    subflowchart : seamm.Flowchart
        A SEAMM Flowchart object that represents a subflowchart, if needed.

    parameters : ThermochemistryParameters
        The control parameters for Thermochemistry.

    See Also
    --------
    TkThermochemistry,
    Thermochemistry, ThermochemistryParameters
    """

    def __init__(
        self,
        flowchart=None,
        title="Thermochemistry",
        extension=None,
        logger=logger,
    ):
        """A substep for Thermochemistry in a subflowchart for Psi4.

        You may wish to change the title above, which is the string displayed
        in the box representing the step in the flowchart.

        Parameters
        ----------
        flowchart: seamm.Flowchart
            The non-graphical flowchart that contains this step.

        title: str
            The name displayed in the flowchart.
        extension: None
            Not yet implemented
        logger : Logger = logger
            The logger to use and pass to parent classes

        Returns
        -------
        None
        """
        logger.debug(f"Creating Thermochemistry {self}")

        super().__init__(
            flowchart=flowchart,
            title=title,
            extension=extension,
            logger=logger,
        )

        self._calculation = "thermochemistry"
        self._model = None
        self._metadata = psi4_step.metadata
        self.parameters = psi4_step.ThermochemistryParameters()

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

    def description_text(self, P=None):
        """Create the text description of what this step will do.
        The dictionary of control values is passed in as P so that
        the code can test values, etc.

        Parameters
        ----------
        P: dict
            An optional dictionary of the current values of the control
            parameters.
        Returns
        -------
        str
            A description of the current step.
        """
        if not P:
            P = self.parameters.values_to_dict()

        text = super().description_text(P=P, calculation_type="Thermochemistry")

        added = (
            "\nThe thermodynamic functions will be calculated at temperature {T} and "
            "pressure {P}."
        )

        return text + "\n" + __(added, **P, indent=4 * " ").__str__()

    def get_input(self, calculation_type="frequency"):
        """Get the input for an optimization calculation for Psi4"""
        _, configuration = self.get_system_configuration()

        # Create the directory
        directory = Path(self.directory)
        directory.mkdir(parents=True, exist_ok=True)

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
        self.description.append(__(self.description_text(PP), **PP, indent=self.indent))
        # Figure out what we are doing! The method is HF, B3LYP, CCSD, etc.
        if P["level"] == "recommended":
            method_string = P["method"]
        else:
            method_string = P["advanced_method"]

        if method_string in psi4_step.methods:
            method = psi4_step.methods[method_string]["method"]
        else:
            method = method_string

        if method == "dft":
            if P["level"] == "recommended":
                functional_string = P["functional"]
            else:
                functional_string = P["advanced_functional"]
            method = psi4_step.dft_functionals[functional_string]["name"]
            if (
                P["dispersion"] != "none"
                and len(psi4_step.dft_functionals[functional_string]["dispersion"]) > 1
            ):
                method = method + "-" + P["dispersion"]

        lines = []
        lines.append("")
        lines.append("#" * 80)
        lines.append(f"# {self.header}")
        lines.append("#" * 80)
        lines.append("")
        # lines.append("initial.find_point_group(tolerance=1.0e-5)")
        lines.append("initial.symmetrize(1.0e-5)")
        lines.append("point_group = initial.point_group().symbol()")
        lines.append("")

        if not P["use existing parameters"]:
            # Add in the input from the energy part of things
            lines.append(super().get_input(calculation_type=calculation_type))
        else:
            lines.append(f"Eelec, wfn = frequency('{method}', return_wfn=True)")

            # Orbital plots
            lines.append(self.plot_input())

        lines.append(
            f"""
set writer_file_label thermo
set hessian_write on

thermo_vibinfo = vibanal_wfn(wfn)
path = Path(core.get_writer_file_prefix('thermo') + '.vibrec')
with path.open() as fd:
    tmp = json.load(fd)
tmp2 = dict()
for key, value in tmp.items():
    if type(value) == str:
        tmp2[key] = json.loads(value)
    else:
        tmp2[key] = value

with path.with_name('@{self._id[-1]}+thermochemistry.json').open('w') as fd:
    json.dump(tmp2, fd, sort_keys=True, indent=3)
"""
        )

        return "\n".join(lines)

    def analyze(self, indent="", **kwargs):
        """Do any analysis of the output from this step.

        Also print important results to the local step.out file using
        "printer".

        Parameters
        ----------
        indent: str
            An extra indentation for the output
        """
        # Read in the results from json
        directory = Path(self.directory)
        json_file = directory / "thermochemistry.json"
        if json_file.exists():
            with json_file.open() as fd:
                tmp = json.load(fd)

            # Process the data, changing units
            d = {}
            metadata = self.metadata["results"]
            for key, value in tmp.items():
                if key in metadata:
                    from_units = value["units"].replace("Eh", "E_h")
                    to_units = metadata[key]["units"]
                    if from_units == to_units:
                        d[key] = value["data"]
                    else:
                        d[key] = Q_(value["data"], from_units).m_as(to_units)

            # Put any requested results into variables or tables
            self.store_results(data=d, create_tables=True)

            # And the output
            table = {
                "": [
                    "Units",
                    "Electronic",
                    "Translational",
                    "Rotational",
                    "Vibrational",
                    "Total Correction",
                    "Total (E_h)",
                ],
                "S": [
                    "J/mol/K",
                    round(d["S_elec"], 2),
                    round(d["S_trans"], 2),
                    round(d["S_rot"], 2),
                    round(d["S_vib"], 2),
                    round(d["S_tot"], 2),
                ],
                "Cv": [
                    "J/mol/K",
                    round(d["Cv_elec"], 2),
                    round(d["Cv_trans"], 2),
                    round(d["Cv_rot"], 2),
                    round(d["Cv_vib"], 2),
                    round(d["Cv_tot"], 2),
                ],
                "Cp": [
                    "J/mol/K",
                    round(d["Cp_elec"], 2),
                    round(d["Cp_trans"], 2),
                    round(d["Cp_rot"], 2),
                    round(d["Cp_vib"], 2),
                    round(d["Cp_tot"], 2),
                ],
                "E": [
                    "kJ/mol",
                    round(d["E_elec"], 2),
                    round(d["E_trans"], 2),
                    round(d["E_rot"], 2),
                    round(d["E_vib"], 2),
                    round(d["E_corr"], 2),
                    round(d["E_tot"], 6),
                ],
                "H": [
                    "kJ/mol",
                    round(d["H_elec"], 2),
                    round(d["H_trans"], 2),
                    round(d["H_rot"], 2),
                    round(d["H_vib"], 2),
                    round(d["H_corr"], 2),
                    round(d["H_tot"], 6),
                ],
                "G": [
                    "kJ/mol",
                    round(d["G_elec"], 2),
                    round(d["G_trans"], 2),
                    round(d["G_rot"], 2),
                    round(d["G_vib"], 2),
                    round(d["G_corr"], 2),
                    round(d["G_tot"], 6),
                ],
            }

            text = ""
            tmp = tabulate(
                table,
                headers="keys",
                tablefmt="simple_outline",
                disable_numparse=True,
                numalign="decimal",
                stralign="center",
                colalign=(
                    "center",
                    "decimal",
                    "decimal",
                    "decimal",
                    "decimal",
                    "decimal",
                    "decimal",
                ),
            )
            length = len(tmp.splitlines()[0])
            text += "\n"
            parameters = self.parameters.current_values_to_dict(
                context=seamm.flowchart_variables._data
            )
            T = parameters["T"]
            P = parameters["P"]
            text += f"Thermodynamic Functions at {T:.2f~P} and {P:.2f~P}".center(length)
            text += "\n"
            text += tmp
            text += "\n"
            printer.normal(__(text, indent=8 * " ", wrap=False, dedent=False))
        else:
            text = (
                "\nThere are no thermochemistry results from Psi4. Perhaps it "
                f"failed? Looking for {str(json_file)}."
            )
            printer.normal(__(text, indent=self.indent + 4 * " "))
            raise RuntimeError(text)
