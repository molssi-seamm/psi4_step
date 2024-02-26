# -*- coding: utf-8 -*-

"""Non-graphical part of the Psi4 step in a SEAMM flowchart
"""

import configparser
import importlib
import json
import logging
from pathlib import Path
import pprint
import shutil

import psi4_step
import seamm
import seamm_exec
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

# In addition to the normal logger, two logger-like printing facilities are
# defined: 'job' and 'printer'. 'job' send output to the main job.out file for
# the job, and should be used very sparingly, typically to echo what this step
# will do in the initial summary of the job.
#
# 'printer' sends output to the file 'step.out' in this steps working
# directory, and is used for all normal output from this step.

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter("Psi4")

pre_code = """\
def fix_multipoles(data):
    result = {}
    esp = []
    it = iter(data.items())
    for key, value in it:
        if 'PROP ' == key[0:5]:
            result[key[6:]] = value[0]
        elif 'SCF DIPOLE' == key:
            result[key] = value[0]
        elif 'CURRENT DIPOLE' == key:
            result[key] = value[0]
        elif '32-POLE' in key:
            tmp = []
            while True:
                value = 0.0 if abs(value) < 1.0e-10 else value
                tmp.append(value)
                if 'ZZZZZ' in key:
                    break
                key, value = next(it)
            result['32-POLE'] = tmp
        elif 'HEXADECAPOLE' in key:
            tmp = []
            while True:
                value = 0.0 if abs(value) < 1.0e-10 else value
                tmp.append(value)
                if 'ZZZZ' in key:
                    break
                key, value = next(it)
            result['HEXADECAPOLE'] = tmp
        elif 'OCTUPOLE' in key:
            tmp = []
            while True:
                value = 0.0 if abs(value) < 1.0e-10 else value
                tmp.append(value)
                if 'ZZZ' in key:
                    break
                key, value = next(it)
            result['OCTUPOLE'] = tmp
        elif 'QUADRUPOLE' in key:
            tmp = []
            while True:
                value = 0.0 if abs(value) < 1.0e-10 else value
                tmp.append(value)
                if 'ZZ' in key:
                    break
                key, value = next(it)
            result['QUADRUPOLE'] = tmp
        elif 'DIPOLE' in key:
            tmp = []
            while True:
                value = 0.0 if abs(value) < 1.0e-10 else value
                tmp.append(value)
                result[key] = value
                if 'Z' in key:
                    break
                key, value = next(it)
            result[key[0:-2]] = tmp
        elif 'ESP AT CENTER' in key:
            esp.append(value)
        else:
            result[key] = value

    if len(esp) > 0:
        result['ELECTROSTATIC POTENTIAL'] = esp
    return result
"""


def humanize(memory, suffix="B", kilo=1024):
    """
    Scale memory to its proper format e.g:

        1253656 => '1.20 MiB'
        1253656678 => '1.17 GiB'
    """
    if kilo == 1000:
        units = ["", "k", "M", "G", "T", "P"]
    elif kilo == 1024:
        units = ["", "Ki", "Mi", "Gi", "Ti", "Pi"]
    else:
        raise ValueError("kilo must be 1000 or 1024!")

    for unit in units:
        if memory < kilo:
            return f"{memory:.2f} {unit}{suffix}"
        memory /= kilo


def dehumanize(memory, suffix="B"):
    """
    Unscale memory from its human readable form e.g:

        '1.20 MB' => 1200000
        '1.17 GB' => 1170000000
    """
    units = {
        "": 1,
        "k": 1000,
        "M": 1000**2,
        "G": 1000**3,
        "P": 1000**4,
        "Ki": 1024,
        "Mi": 1024**2,
        "Gi": 1024**3,
        "Pi": 1024**4,
    }

    tmp = memory.split()
    if len(tmp) == 1:
        return memory
    elif len(tmp) > 2:
        raise ValueError("Memory must be <number> <units>, e.g. 1.23 GB")

    amount, unit = tmp
    amount = float(amount)

    for prefix in units:
        if prefix + suffix == unit:
            return int(amount * units[prefix])

    raise ValueError(f"Don't recognize the units on '{memory}'")


class Psi4(seamm.Node):
    """
    The non-graphical part of a Psi4 step in a flowchart.

    Attributes
    ----------
    options : tuple
        It contains a two item tuple containing the populated namespace and the
        list of remaining argument strings.

    subflowchart : seamm.Flowchart
        A SEAMM Flowchart object that represents a subflowchart, if needed.

    parameters : Psi4Parameters
        The control parameters for Psi4.

    See Also
    --------
    TkPsi4,
    Psi4, Psi4Parameters
    """

    def __init__(
        self,
        flowchart=None,
        title="Psi4",
        namespace="org.molssi.seamm.psi4",
        extension=None,
    ):
        """A step for Psi4 in a SEAMM flowchart.

        You may wish to change the title above, which is the string displayed
        in the box representing the step in the flowchart.

        Parameters
        ----------
            flowchart: seamm.Flowchart
                The non-graphical flowchart that contains this step.

            title: str
                The name displayed in the flowchart.
            namespace: The namespace for the plugins of the subflowchart
            extension: None
                Not yet implemented
        """
        logger.debug("Creating Psi4 {}".format(self))

        self.subflowchart = seamm.Flowchart(
            parent=self, name="Psi4", namespace=namespace
        )

        super().__init__(
            flowchart=flowchart, title="Psi4", extension=extension, logger=logger
        )

        self.parameters = psi4_step.Psi4Parameters()

    @property
    def version(self):
        """The semantic version of this module."""
        return psi4_step.__version__

    @property
    def git_revision(self):
        """The git version of this module."""
        return psi4_step.__git_revision__

    def create_parser(self):
        """Setup the command-line / config file parser"""
        # parser_name = 'psi4-step'
        parser_name = self.step_type
        parser = self.flowchart.parser

        # Remember if the parser exists ... this type of step may have been
        # found before
        parser_exists = parser.exists(parser_name)

        # Create the standard options, e.g. log-level
        result = super().create_parser(name=parser_name)

        if parser_exists:
            return result

        # Options for Psi4
        parser.add_argument(
            parser_name,
            "--ncores",
            default="available",
            help="How many threads to use in Psi4",
        )

        parser.add_argument(
            parser_name,
            "--memory",
            default="available",
            help=(
                "The maximum amount of memory to use for Psi4, which can be "
                "'all' or 'available', or a number, which may use k, Ki, "
                "M, Mi, etc. suffixes. Default: available."
            ),
        )

        parser.add_argument(
            parser_name,
            "--memory-factor",
            default="90%",
            help=(
                "The amount of possible memory to use, typically about 90%% to"
                " allow for Psi4 not keeping track of all the memory used."
            ),
        )

        return result

    def set_id(self, node_id):
        """Set the id for node to a given tuple"""
        self._id = node_id

        # and set our subnodes
        self.subflowchart.set_ids(self._id)

        return self.next()

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
            description : str
                A description of the current step.
        """

        self.subflowchart.root_directory = self.flowchart.root_directory

        # Get the first real node
        node = self.subflowchart.get_node("1").next()

        text = self.header + "\n\n"
        while node is not None:
            text += __(node.description_text(), indent=4 * " ").__str__()
            text += "\n"
            node = node.next()

        return text

    def run(self):
        """Run a Psi4 step.

        Returns
        -------

        next_node : seamm.Node
            The next node object in the flowchart.

        """
        printer.important(self.header)
        printer.important("")

        # Add the main citation for Psi4
        self.references.cite(
            raw=self._bibliography["doi:10.1063/5.0006002"],
            alias="psi4",
            module="psi4 step",
            level=1,
            note="The principle Psi4 citation.",
        )

        system_db = self.get_variable("_system_db")
        configuration = system_db.system.configuration
        n_atoms = configuration.n_atoms
        if n_atoms == 0:
            self.logger.error("Psi4 run(): there is no structure!")
            raise RuntimeError("Psi4 run(): there is no structure!")

        # Access the options
        options = self.options
        seamm_options = self.global_options

        # Get the computational environment and set limits
        ce = seamm_exec.computational_environment()

        # Maximum number of threads
        n_threads = ce["NTASKS"]
        if options["ncores"] == "available":
            pass
        else:
            n_threads = int(options["ncores"])
            if n_threads < 1:
                n_threads = 1
        if seamm_options["ncores"] != "available":
            tmp = int(seamm_options["ncores"])
            if tmp < n_threads:
                n_threads = tmp
        ce["NTASKS"] = n_threads

        # And memory
        if seamm_options["memory"] == "all":
            mem_limit = ce["MEM_PER_NODE"]
        elif seamm_options["memory"] == "available":
            # For the default, 'available', use in proportion to number of
            # cores used
            mem_limit = n_threads * ce["MEM_PER_CPU"]
        else:
            mem_limit = dehumanize(seamm_options["memory"])

        if options["memory"] == "all":
            memory = ce["MEM_PER_NODE"]
        elif options["memory"] == "available":
            # For the default, 'available', use in proportion to number of
            # cores used
            memory = n_threads * ce["MEM_PER_CPU"]
        else:
            memory = dehumanize(options["memory"])

        memory = min(memory, mem_limit)

        if "%" in options["memory_factor"]:
            factor = float(options["memory_factor"].rstrip("%")) / 100.0
        else:
            factor = float(options["memory_factor"])

        memory *= factor

        # Psi applies a minimum of 250 MiB
        min_memory = dehumanize("250 MiB")
        if min_memory > memory:
            memory = min_memory
        ce["MEM_PER_NODE"] = memory
        memory = humanize(memory, kilo=1024)

        # Work through the subflowchart to find out what to do.
        self.subflowchart.root_directory = self.flowchart.root_directory

        next_node = super().run(printer)

        # Get the first real node
        node0 = self.subflowchart.get_node("1").next()

        # Start the input data
        input_data = []
        input_data.append("import json")
        input_data.append("import numpy as np")
        input_data.append("from pathlib import Path")
        input_data.append("import pprint")
        input_data.append("")
        input_data.append(f"memory {memory}")
        input_data.append("")
        input_data.append(pre_code)
        input_data.append("")

        # Put the structure into the input
        input_data.append(self._convert_structure(name="initial"))

        node = node0
        while node is not None:
            text = node.get_input()
            input_data.append(text)

            input_data.append("clean()")
            input_data.append("clean_variables()")
            # input_data.append('clean_options()')

            node = node.next()

        # Write out the final structure
        input_data.append("")
        input_data.append("# Write the final structure to disk")
        input_data.append("molecule = get_active_molecule()")
        input_data.append("tmp = molecule.to_dict()")
        input_data.append("for item, value in tmp.items():")
        input_data.append("    if isinstance(value, np.ndarray):")
        input_data.append("        tmp[item] = value.tolist()")
        input_data.append("")
        input_data.append("with open('final_structure.json', 'w') as fd:")
        input_data.append("    json.dump(tmp, fd, sort_keys=True, indent=3)")

        files = {"input.dat": "\n".join(input_data)}
        self.logger.info("input.dat:\n" + files["input.dat"])

        directory = Path(self.directory)
        directory.mkdir(parents=True, exist_ok=True)

        return_files = ["output.dat", "*.json", "*.cube"]

        # Check for already having run
        path = Path(self.directory) / "success.dat"
        if path.exists():
            result = {}
            path = Path(self.directory) / "stdout.txt"
            if path.exists():
                result["stdout"] = path.read_text()
            result["stderr"] = ""
            failed = False
        else:
            executor = self.flowchart.executor

            # Read configuration file for Psi4 if it exists
            executor_type = executor.name
            full_config = configparser.ConfigParser()
            ini_dir = Path(seamm_options["root"]).expanduser()
            path = ini_dir / "psi4.ini"

            if path.exists():
                full_config.read(ini_dir / "psi4.ini")

            # If the section we need doesn't exists, get the default
            if not path.exists() or executor_type not in full_config:
                resources = importlib.resources.files("psi4_step") / "data"
                ini_text = (resources / "psi4.ini").read_text()
                full_config.read_string(ini_text)

            # Getting desperate! Look for an executable in the path
            if executor_type not in full_config:
                path = shutil.which("psi4")
                if path is None:
                    raise RuntimeError(
                        f"No section for '{executor_type}' in Psi4 ini file "
                        f"({ini_dir / 'psi4.ini'}), nor in the defaults, nor "
                        "in the path!"
                    )
                else:
                    full_config[executor_type] = {
                        "installation": "local",
                        "code": f"{path} -n {{NTASKS}}",
                    }

            config = dict(full_config.items(executor_type))
            # Use the matching version of the seamm-psi4 image by default.
            config["version"] = self.version

            printer.important(
                self.indent
                + f"    Psi4 will use {ce['NTASKS']} threads and {memory} memory\n"
            )

            result = executor.run(
                cmd=["{code}"],
                config=config,
                directory=self.directory,
                files=files,
                return_files=return_files,
                in_situ=True,
                shell=True,
                ce=ce,
            )

            if not result:
                self.logger.error("There was an error running Psi4")
                raise RuntimeError("There was an error running Psi4")

            self.logger.debug("\n" + pprint.pformat(result))

            failed = False
            if "output.dat" in result["files"]:
                if result["output.dat"]["data"] is not None:
                    if (
                        "*** Psi4 exiting successfully."
                        not in result["output.dat"]["data"]
                    ):
                        self.logger.warning("Psi4 did not complete successfully.")
                        failed = True
            if not failed:
                # Write a small file to say that LAMMPS ran successfully, so cancel
                # skip if rerunning.
                path = Path(self.directory) / "success.dat"
                path.write_text("success")

        # Analyze the results
        self.analyze()

        if failed:
            raise RuntimeError("Psi4 did not complete successfully.")

        return next_node

    def analyze(self, indent="", **kwargs):
        """Do any analysis of the output from this step.

        Also print important results to the local step.out file using
        'printer'.

        Parameters
        ----------
            indent: str
                An extra indentation for the output
        """
        # Update the structure
        directory = Path(self.directory)
        structure_file = directory / "final_structure.json"
        if structure_file.exists():
            with structure_file.open(mode="r") as fd:
                structure = json.load(fd)
            if "geom" in structure:
                system_db = self.get_variable("_system_db")
                configuration = system_db.system.configuration
                xs = []
                ys = []
                zs = []
                it = iter(structure["geom"])
                for x in it:
                    xs.append(x)
                    ys.append(next(it))
                    zs.append(next(it))
                configuration.atoms["x"][0:] = xs
                configuration.atoms["y"][0:] = ys
                configuration.atoms["z"][0:] = zs
                printer.important(
                    self.indent + "    Updated the system with the structure from Psi4",
                )
                printer.important("")

        # Get the first real node
        node = self.subflowchart.get_node("1").next()

        # Loop over the subnodes, asking them to do their analysis
        while node is not None:
            for value in node.description:
                printer.important(value)

            node.analyze()

            printer.normal("")

            node = node.next()

    def _convert_structure(self, name=None, no_com=True, no_reorient=True):
        """Convert the structure to the input for Psi4."""

        system_db = self.get_variable("_system_db")
        configuration = system_db.system.configuration

        structure = []
        if name is None:
            structure.append("molecule {")
        else:
            structure.append("molecule " + name + " {")

        # Charge and multiplicity
        structure.append(f"  {configuration.charge} {configuration.spin_multiplicity}")

        elements = configuration.atoms.symbols
        coordinates = configuration.atoms.coordinates

        if "freeze" in configuration.atoms:
            freeze = configuration.atoms["freeze"]
        else:
            freeze = [""] * len(elements)

        for element, xyz, frz in zip(elements, coordinates, freeze):
            x, y, z = xyz
            structure.append(
                f"    {element:2s} {float(x): 12.8f} {float(y): 12.8f} "
                f"{float(z): 12.8f}"
            )
        structure.append("")
        if no_com:
            structure.append("no_com")
        if no_reorient:
            structure.append("no_reorient")
        structure.append("}")

        return "\n".join(structure) + "\n"
