# -*- coding: utf-8 -*-

"""Acclerated optimization using MOPAC for Hessians."""

import json
import logging
from pathlib import Path
import pkg_resources
import sys
import textwrap
import traceback

import numpy as np
import optking
from scipy.optimize import minimize, OptimizeResult
from scipy.optimize._optimize import (
    _prepare_scalar_function,
    _status_message,
    _epsilon,
    vecnorm,
    _line_search_wolfe12,
    _LineSearchError,
)

from .optimizer import Optimizer
import psi4_step
import seamm
import seamm.data
from seamm_util import units_class, Q_
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter("psi4")

optking_logger = logging.getLogger("optking")
optking_logger.setLevel(logging.WARNING)


def norm(V):
    return np.linalg.norm(V)


def abs_max(V):
    return max(abs(elem) for elem in V)


def abs_min(V):
    return min(abs(elem) for elem in V)


def rms(V):
    return np.sqrt(np.mean(V**2))


class AcceleratedOptimization(psi4_step.Energy):
    def __init__(
        self,
        flowchart=None,
        title="Accelerated Optimization",
        extension=None,
    ):
        """Initialize the node"""

        logger.debug("Creating AcceleratedOptimization {}".format(self))

        super().__init__(flowchart=flowchart, title=title, extension=extension)

        self._calculation = "optimization"
        self._model = None
        self._metadata = psi4_step.metadata
        self.parameters = psi4_step.AcceleratedOptimizationParameters()

        self.description = "A geometry optimization"
        self._E = None
        self._derivatives = None
        self._triangle = None  # Triangule Hessian matrix
        self.prolog = None
        self.epilog = None
        self.iteration = 0
        self._restart_wfn = None

        self._counters = {}

    def description_text(self, P=None):
        """Prepare information about what this node will do"""

        if not P:
            P = self.parameters.values_to_dict()

        text = super().description_text(P=P, calculation_type="Geometry optimization")

        added = "The geometry optimization will use the {optimization method} "
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
            added += (
                " The Hessian will be recalculated with MOPAC every {recalc hessian}"
                " steps. Note that calculating the second derivatives with MOPAC is "
                "not expensive!"
            )

        return text + "\n" + __(added, **P, indent=4 * " ").__str__()

    def input_text(self, node0, memory, restart=None):
        """Create the first and last parts of the Psi4 input.

        Parameters
        ----------
        node0 : seamm.Node
            The first node in the Psi4 subflowchart

        memory : str
            Amount of memory to allocate Psi4

        n_threads : int
            Number of threads to use in Psi4

        restart : str = None
            A wavefunction file to restart from
        """
        prolog = textwrap.dedent(
            f"""
            import json
            import numpy as np
            import pprint

            memory {memory}

            """
        )

        node = node0
        tmp = []
        while node is not None:
            if node == self:
                break
            text = node.get_input()
            tmp.append(text)

            tmp.append("clean()")
            tmp.append("clean_variables()")
            # tmp.append('clean_options()')

            node = node.next()

        # Add our input
        tmp.append("")
        tmp.append("#" * 80)
        tmp.append(f"# {self.header}")
        tmp.append("#" * 80)
        tmp.append("")

        # Add in the input from the energy part of things
        tmp.append(super().get_input(calculation_type="gradient", restart=restart))

        # Write out the wavefunction
        tmp.append("")
        tmp.append("wfn.to_file('wfn')")

        # Write out the final structure
        tmp.append("")
        tmp.append("# Write the final structure to disk")
        tmp.append("molecule = get_active_molecule()")
        tmp.append("tmp = molecule.to_dict()")
        tmp.append("for item, value in tmp.items():")
        tmp.append("    if isinstance(value, np.ndarray):")
        tmp.append("        tmp[item] = value.tolist()")
        tmp.append("")
        tmp.append("# Add the energy and gradients")
        tmp.append("tmp['energy'] = wfn.energy()")
        tmp.append("tmp['gradient'] = wfn.gradient().to_array().tolist()")
        tmp.append("")
        tmp.append("with open('final_structure.json', 'w') as fd:")
        tmp.append("    json.dump(tmp, fd, sort_keys=True, indent=3)")

        epilog = "\n".join(tmp)

        return prolog, epilog

    def run(self, node0, memory, n_threads):
        """Run the accelerated optimization.

        Parameters
        ----------
        node0 : seamm.Node
            The first node in the Psi4 subflowchart

        memory : str
            Amount of memory to allocate Psi4

        n_threads : int
            Number of threads to use in Psi4
        """
        # Create the directory
        directory = Path(self.directory)
        directory.mkdir(parents=True, exist_ok=True)

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

        # Handle the system
        _, configuration = self.get_system_configuration(P)

        # Get the input other than the structure
        self.n_threads = n_threads
        self._node0 = node0
        self._memory = memory
        if self._restart_wfn is not None:
            self.prolog, self.epilog = self.input_text(
                node0, memory, restart="old_wfn.npy"
            )
        else:
            self.prolog, self.epilog = self.input_text(node0, memory)

        max_iterations = P["max geometry steps"]
        if max_iterations == "default":
            max_iterations = 2 * 3 * configuration.n_atoms

        self._counters["surrogate"] = 0
        self._counters["hessian"] = 0
        self._counters["main"] = 0

        # self._test(directory)
        # return

        use_surrogate = False
        if use_surrogate:
            if True:
                optimizer = Optimizer(
                    energy_callback=self.runPsi4, surrogate_callback=self.runMOPAC
                )
                x0 = configuration.atoms.get_coordinates()
                optimizer.optimize(x0)
            else:
                for self.iteration in range(max_iterations + 1):
                    converged = self._surrogate_optimization(
                        directory / f"step_{self.iteration:03d}"
                    )
                    if converged:
                        printer.normal("Converged!!!!!!")
                        break
                else:
                    printer.normal("NOT Converged!!!!!!")
        else:
            # Iterate
            molecule = configuration.to_qcschema_dict()
            optking_options = {
                "step_type": P["optimization method"],
                "g_convergence": P["geometry convergence"],
                "opt_coordinates": P["coordinates"],
                "hess_update": P["hessian update"],
            }
            if P["recalc hessian"] == "every step":
                optking_options["full_hess_every"] = 1
            elif P["recalc hessian"] == "at beginning":
                optking_options["full_hess_every"] = 0
            elif P["recalc hessian"] == "never":
                optking_options["full_hess_every"] = -1
            else:
                optking_options["full_hess_every"] = P["recalc hessian"]

            # "max_force_g_convergence": 4.5e-04,
            # "rms_force_g_convergence": 3.0e-04,
            # "intrafrag_hess": "SIMPLE",

            optimizer = optking.CustomHelper(molecule, optking_options)

            # The initial coordinates as NUMPY array. Psi4 may reorient, so need
            # to run and cache the first results
            for self.iteration in range(max_iterations + 1):
                subdirectory = directory / f"step_{self.iteration:03d}"
                subdirectory.mkdir(parents=True, exist_ok=True)

                needed = optimizer.calculations_needed()
                factor = Q_(1.0, "a_0").m_as("Å")
                xyz = []
                for row in optimizer.geom:
                    xyz.append([val * factor for val in row])

                configuration.atoms.set_coordinates(xyz)

                if "hessian" in needed:
                    # Run MOPAC
                    self._triangle = self._run_MOPAC_Hessian(subdirectory / "MOPAC")
                    H = np.array(self._hessian(None))
                    optimizer.HX = H

                E, dE = self._run_psi4(subdirectory)
                optimizer.E = E
                optimizer.gX = dE

                optimizer.compute()
                optimizer.take_step()

                lines = optimizer.test_convergence(str_mode="table").splitlines()
                if self.iteration == 0:
                    printer.normal("\n".join(lines[0:-2]))
                else:
                    printer.normal(lines[-3])
                converged = optimizer.test_convergence()
                if converged is True:
                    printer.normal("Optimization SUCCESS:")
                    break
            else:
                printer.normal("Optimization FAILURE:\n")

            result = optimizer.close()
            if False:
                printer.normal(result)

        printer.normal(f"   Calls to Psi4: {self._counters['main']:4d}")
        printer.normal(f"      to Hessian: {self._counters['hessian']:4d}")
        printer.normal(f"    to surrogate: {self._counters['surrogate']:4d}")

    def _test(self, directory):
        """Calculate derivatives as a check"""

        # Handle the system
        _, configuration = self.get_system_configuration()
        directory.mkdir(parents=True, exist_ok=True)

        xyz0 = configuration.atoms.get_coordinates(as_array=True)

        # Run MOPAC and Psi4
        Emopac, gradients = self._run_MOPAC_SPE(directory / "MOPAC")
        E0, dE0 = self._run_psi4(directory / "Psi4")

        print(f"{E0=:9.6f}")

        # And the correction to MOPAC dE
        delta_E = E0 - Emopac
        delta_dE = dE0 - gradients

        angstrom2bohr = Q_(1.0, "Å").m_as("a_0")
        direction = np.array(dE0) / norm(dE0)
        direction.reshape(configuration.n_atoms, 3)
        for factor in range(-3, 4, 1):
            dxyz = direction * (factor / 10)
            xyz = xyz0 + dxyz.reshape(configuration.n_atoms, 3) / angstrom2bohr
            configuration.atoms.set_coordinates(xyz)
            E, dE = self._run_MOPAC_SPE(directory / f"step{factor}")
            Ep, dEp = self._run_psi4(directory / f"psi4_{factor}")

            delta_E1 = np.dot(dxyz, delta_dE)
            E = E + delta_E + delta_E1
            dE += delta_dE

            print(
                f"{factor:2d} {E:9.6f} {np.dot(dE, direction.flatten()):9.6f} "
                f"{Ep:9.6f} {np.dot(dEp, direction.flatten()):9.6f}"
            )

    def _surrogate_optimization(self, directory):
        """Optimize using MOPAC gradients shifted by real ones."""
        debug = False

        # Handle the system
        _, configuration = self.get_system_configuration()

        molecule = configuration.to_qcschema_dict()
        # "g_convergence": "gau",
        optking_options = {
            "opt_coordinates": "cartesian",
            "g_convergence": "cfour",
            "max_force_g_convergence": 4.5e-04,
            "rms_force_g_convergence": 3.0e-04,
            "FULL_HESS_EVERY": 1,
            "intrafrag_hess": "SIMPLE",
        }
        # "hess_update": "none",
        optimizer = optking.CustomHelper(molecule, optking_options)

        # First calculate the derivatives with both Psi4 and MOPAC
        xyz0 = np.array(optimizer.geom).flatten()
        factor = Q_(1.0, "a_0").m_as("Å")

        xyz = xyz0 * factor
        configuration.atoms.set_coordinates(xyz.reshape(configuration.n_atoms, 3))

        directory.mkdir(parents=True, exist_ok=True)

        # Run MOPAC and Psi4
        Emopac, gradients = self._run_MOPAC_SPE(directory / "MOPAC")
        E, dE = self._run_psi4(directory / "Psi4")

        tmp = gradients.reshape(configuration.n_atoms, 3)
        tmp2 = np.sum(tmp, axis=0)
        tmp3 = " ".join([f"{val:12.9f}" for val in tmp2])
        print(f"sum dEm {tmp3}")

        tmp = dE.reshape(configuration.n_atoms, 3)
        tmp2 = np.sum(tmp, axis=0)
        tmp3 = " ".join([f"{val:12.9f}" for val in tmp2])
        print(f"sum dEp {tmp3}")

        tmp -= tmp2 / configuration.n_atoms
        dE = tmp.flatten()

        tmp = dE.reshape(configuration.n_atoms, 3)
        tmp2 = np.sum(tmp, axis=0)
        tmp3 = " ".join([f"{val:12.9f}" for val in tmp2])
        print(f"sum dEp {tmp3}")

        if debug:
            print("MOPAC gradients")
            print(gradients)
            print("Psi4 derivatives")
            print(dE)

        # And the correction to MOPAC dE
        delta_E = E - Emopac
        delta_dE = dE - gradients

        if debug:
            m0 = abs_max(delta_dE)
            m1 = abs_max(gradients)
            m2 = abs_max(dE)
            rms2 = rms(dE)
            print(f"maxes: {m0:9.6f} {m1:9.6f} {m2:9.6f} {rms2:9.6f} {E:12.7f}")

        if debug:
            print("Correction")
            print(delta_dE)

        P = self.parameters.values_to_dict()
        max_iterations = P["max geometry steps"]
        if max_iterations == "default":
            max_iterations = 2 * 3 * configuration.n_atoms

        max_iterations = 2
        converged = False
        for self.iteration in range(max_iterations + 1):
            subdirectory = directory / f"step_{self.iteration:03d}"
            subdirectory.mkdir(parents=True, exist_ok=True)

            needed = optimizer.calculations_needed()
            xyz1 = np.array(optimizer.geom).flatten()
            factor = Q_(1.0, "a_0").m_as("Å")
            xyz = xyz1 * factor
            configuration.atoms.set_coordinates(xyz.reshape(configuration.n_atoms, 3))

            if "hessian" in needed:
                # Run MOPAC
                self._triangle = self._run_MOPAC_Hessian(subdirectory / "MOPAC")
                H = np.array(self._hessian(None))
                Unit = np.identity(configuration.n_atoms * 3) * 0.0
                optimizer.HX = H + Unit

                eigvals = np.linalg.eigvalsh(H + Unit).tolist()
                txt = " ".join([f"{val:6.3f}" for val in eigvals])
                print(f"Eigenvalues {txt}")

            E, dE = self._run_MOPAC_SPE(subdirectory)

            Esv = E
            dEsv = np.array(dE)

            if debug:
                # And Psi4 for testing
                Ex, dEx = self._run_psi4(subdirectory / "Psi4")

            dxyz = xyz1 - xyz0
            delta_E1 = np.dot(dxyz, delta_dE)

            if debug and optimizer.HX is not None:
                tmp = np.matmul(optimizer.HX, dxyz)
                m1 = tmp.max()
                m0 = delta_dE.max()
                e1 = np.dot(dxyz, tmp)
                e0 = np.dot(dxyz, delta_dE)
                print(f"max = {m1:9.6f} {m0:9.6f} delta E {e1:9.6f} {e0:9.6f}")

            E += delta_E + delta_E1
            dE += delta_dE

            if debug:
                m0 = abs_max(delta_dE)
                m1 = abs_max(dEsv)
                m2 = abs_max(dE)
                rms2 = rms(dE)
                print(f"-----> {m0:9.6f} {m1:9.6f} {m2:9.6f} {rms2:9.6f} {E:12.7f}")

            if debug:
                print(f"{Esv=} {delta_E=} {delta_E1=} --> {E} vs {Ex}")
                tmp = [f"{val:10.6f}" for val in dE.tolist()]
                print(" ".join(tmp))
                tmp = [f"{val:10.6f}" for val in dEx.tolist()]
                print(" ".join(tmp))

            optimizer.E = E
            optimizer.gX = dE

            optimizer.compute()
            optimizer.take_step()

            lines = optimizer.test_convergence(str_mode="table").splitlines()
            if self.iteration == 0:
                print("\n".join(lines[0:-2]))
            else:
                print(lines[-3])
            if optimizer.test_convergence():
                converged = self.iteration == 0

                result = optimizer.close()
                factor = Q_(1.0, "a_0").m_as("Å")
                xyzf = np.array(result["final_molecule"]["geometry"]) * factor
                # configuration.atoms.set_coordinates(
                #     xyzf.reshape(configuration.n_atoms, 3)
                # )

                tmp = " ".join([f"{val:9.6f}" for val in xyz.tolist()])
                print(f" xyz: {tmp}")
                tmp = " ".join([f"{val:9.6f}" for val in xyzf.tolist()])
                print(f"xyzf: {tmp}")
                break

        return converged

    def run_old(self, node0, memory, n_threads):
        """Run the accelerated optimization.

        Parameters
        ----------
        node0 : seamm.Node
            The first node in the Psi4 subflowchart

        memory : str
            Amount of memory to allocate Psi4

        n_threads : int
            Number of threads to use in Psi4
        """
        # Create the directory
        directory = Path(self.directory)
        directory.mkdir(parents=True, exist_ok=True)

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

        # Handle the system
        _, configuration = self.get_system_configuration(P)

        # Get the input other than the structure
        self.n_threads = n_threads
        if self._restart_wfn is not None:
            self.prolog, self.epilog = self.input_text(
                node0, memory, restart="old_wfn.npy"
            )
        else:
            self.prolog, self.epilog = self.input_text(node0, memory)

        # Iterate
        max_iterations = P["max geometry steps"]
        if max_iterations == "default":
            max_iterations = 2 * 3 * configuration.n_atoms

        # The initial coordinates as NUMPY array. Psi4 may reorient, so need
        # to run and cache the first results
        self.iteration = 0

        x0 = configuration.atoms.get_coordinates(as_array=True).flatten()
        self._E, self._derivatives = self._step(x0)

        # Set the iteration count back this initial time
        self.iteration = 0

        subdirectory = directory / f"step_{self.iteration:03d}"
        json_file = subdirectory / "final_structure.json"
        if json_file.exists():
            with json_file.open() as fd:
                data = json.load(fd)
        else:
            raise RuntimeError("Structure and gradients from Psi4 not available.")

        self.iteration = 0
        x0 = np.array(data["geom"])
        configuration.atoms.set_coordinates(x0.reshape(configuration.n_atoms, 3))

        # Run MOPAC
        self._triangle = self._run_MOPAC_Hessian(subdirectory / "MOPAC")

        H = self._hessian(None)
        print("MOPAC Hessian")
        for n, row in enumerate(H.tolist()):
            txt = [f"{val:6.3f}" for val in row[: n + 1]]
            print(" ".join(txt))

        print("")

        singular = np.linalg.svd(H, compute_uv=False, hermitian=True).tolist()
        txt = [f"{val:.6f}" for val in singular]
        print("Singular values")
        print(" ".join(txt))
        print(f"{singular[-7]} {singular[-6]}")
        self._rcond = (singular[-7] + singular[-6]) / 2
        print("")

        Hk = np.linalg.pinv(H, rcond=self._rcond, hermitian=True)

        print("Inverse")
        for n, row in enumerate(Hk.tolist()):
            txt = [f"{val:6.3f}" for val in row[: n + 1]]
            print(" ".join(txt))

        print("")

        eigvals = np.linalg.eigvalsh(H).tolist()
        txt = [f"{val:6.3f}" for val in eigvals]
        print("Eigenvalues")
        print(" ".join(txt))
        print("")

        # Set up to restart from the previous wavefunction
        self.prolog, self.epilog = self.input_text(node0, memory, restart="old_wfn.npy")

        result = minimize(self._step, x0, method=self._bfgs, jac=True)
        # result = minimize(self._step, x0, method="BFGS", jac=True)
        # result = minimize(
        #     self._step,
        #     x0,
        #     method="trust-exact",
        #     jac=True,
        #     hess=self._hessian,
        # )

        # for key, value in result.items():
        #     print(key)
        #     print(value)
        print("BFGS inverse Hessian")
        for n, row in enumerate(result["hess_inv"].tolist()):
            txt = [f"{val:6.3f}" for val in row[: n + 1]]
            print(" ".join(txt))

        # Output the summary of the minimization
        if result["success"]:
            text = (
                "\nThe optimization converged in {nit} steps to {fun:.7f} E_h. "
                "This required {nfev} evaluations of the energy and derivatives."
            )
            printer.normal(__(text, **result, indent=4 * " "))
        else:
            text = (
                "\nThe optimization failed: {message} after {nit} steps, which "
                "required {nfev} evaluations of the energy and derivatives. The "
                "final energy was {fun:.7f} E_h."
            )

    def _run_MOPAC_SPE(self, directory):
        """Run the MOPAC flowchart in the given directory."""
        # Get the MOPAC flowchart
        self._counters["surrogate"] += 1
        directory.mkdir(parents=True, exist_ok=True)
        path = Path(pkg_resources.resource_filename(__name__, "data/"))
        flowchart = seamm.Flowchart(
            name="MOPAC",
            parent=self,
            namespace="org.molssi.seamm",
            directory=directory,
            parser_name="MOPAC subflowchart",
        )
        flowchart.read(path / "MOPAC_PM7_SPE.flow")

        # Find the handler for job.out and set the level up
        job_handler = None
        out_handler = None
        for handler in job.handlers:
            if (
                isinstance(handler, logging.FileHandler)
                and "job.out" in handler.baseFilename
            ):
                job_handler = handler
                job_level = job_handler.level
                job_handler.setLevel(printing.JOB)
            elif isinstance(handler, logging.StreamHandler):
                out_handler = handler
                out_level = out_handler.level
                out_handler.setLevel(printing.JOB)

        # Redirect output
        _file_handler = logging.FileHandler(directory / "subflowchart.out")
        _file_handler.setLevel(printing.NORMAL)
        formatter = logging.Formatter(fmt="{message:s}", style="{")
        _file_handler.setFormatter(formatter)
        job.addHandler(_file_handler)

        # Set up the argument parser for this node.
        parser = flowchart.parser

        parser.epilog = textwrap.dedent(
            """
            The plug-ins in this flowchart are listed above.
            Options, if any, for plug-ins are placed after
            the name of the plug-in, e.g.:

               test.flow lammps-step --log-level DEBUG --np 4

            To get help for a plug-in, use --help or -h after the
            plug-in name. E.g.

               test.flow lammps-step --help
            """
        )
        parser.usage = "%(prog)s [options] plug-in [options] plug-in [options] ..."

        # Now traverse the flowchart, setting up the ids and parsers
        flowchart.set_ids()
        flowchart.create_parsers()

        # And handle the command-line arguments and ini file options.
        parser.parse_args([])
        logger.info("Parsed the command-line arguments")

        # Get the command line options
        options = parser.get_options()

        # Set the options in each step
        for node in flowchart:
            step_type = node.step_type
            logger.info(f"    setting options for {step_type}")
            if step_type in options:
                node._options = options[step_type]
            if "SEAMM" in options:
                node._global_options = options["SEAMM"]

        # Write out an initial summary of the flowchart before doing anything
        # Reset the visited flag for traversal
        flowchart.reset_visited()

        # Get the start node
        next_node = flowchart.get_node("1")

        # describe ourselves
        printer.normal("\nDescription of the flowchart\n----------------------------")

        while next_node:
            # and print the description
            try:
                next_node = next_node.describe()
            except Exception:
                message = "Error describing the flowchart\n\n" + traceback.format_exc()
                print(message)
                logger.critical(message)
                raise
            except:  # noqa: E722
                message = (
                    "Unexpected error describing the flowchart\n\n"
                    + traceback.format_exc()
                )
                print(message)
                logger.critical(message)
                raise

        printer.normal("")

        # And actually run it!
        printer.normal(("Running the flowchart\n" "---------------------"))

        try:
            next_node = flowchart.get_node("1")
            while next_node is not None:
                try:
                    next_node = next_node.run()
                except DeprecationWarning as e:
                    print("\nDeprecation warning: " + str(e))
                    traceback.print_exc(file=sys.stderr)
                    traceback.print_exc(file=sys.stdout)
        finally:
            pass

        # Get the derivatives and return them.
        E = self.get_variable("_mopac_energy")
        gradients = self.get_variable("_mopac_gradients")
        dE = np.array(gradients).flatten()

        # Set output back
        _file_handler.close()
        job.removeHandler(_file_handler)
        if job_handler is not None:
            job_handler.setLevel(job_level)
        if out_handler is not None:
            out_handler.setLevel(out_level)

        return E, dE

    def _run_MOPAC_Hessian(self, directory):
        """Run the MOPAC flowchart in the given directory."""
        self._counters["hessian"] += 1
        # Get the MOPAC flowchart
        directory.mkdir(parents=True, exist_ok=True)
        path = Path(pkg_resources.resource_filename(__name__, "data/"))
        flowchart = seamm.Flowchart(
            name="MOPAC",
            parent=self,
            namespace="org.molssi.seamm",
            directory=directory,
            parser_name="MOPAC subflowchart",
        )
        flowchart.read(path / "MOPAC_PM7_Hessian.flow")

        # Find the handler for job.out and set the level up
        job_handler = None
        out_handler = None
        for handler in job.handlers:
            if (
                isinstance(handler, logging.FileHandler)
                and "job.out" in handler.baseFilename
            ):
                job_handler = handler
                job_level = job_handler.level
                job_handler.setLevel(printing.JOB)
            elif isinstance(handler, logging.StreamHandler):
                out_handler = handler
                out_level = out_handler.level
                out_handler.setLevel(printing.JOB)

        # Redirect output
        _file_handler = logging.FileHandler(directory / "subflowchart.out")
        _file_handler.setLevel(printing.NORMAL)
        formatter = logging.Formatter(fmt="{message:s}", style="{")
        _file_handler.setFormatter(formatter)
        job.addHandler(_file_handler)

        # Set up the argument parser for this node.
        parser = flowchart.parser

        parser.epilog = textwrap.dedent(
            """
            The plug-ins in this flowchart are listed above.
            Options, if any, for plug-ins are placed after
            the name of the plug-in, e.g.:

               test.flow lammps-step --log-level DEBUG --np 4

            To get help for a plug-in, use --help or -h after the
            plug-in name. E.g.

               test.flow lammps-step --help
            """
        )
        parser.usage = "%(prog)s [options] plug-in [options] plug-in [options] ..."

        # Now traverse the flowchart, setting up the ids and parsers
        flowchart.set_ids()
        flowchart.create_parsers()

        # And handle the command-line arguments and ini file options.
        parser.parse_args([])
        logger.info("Parsed the command-line arguments")

        # Get the command line options
        options = parser.get_options()

        # Set the options in each step
        for node in flowchart:
            step_type = node.step_type
            logger.info(f"    setting options for {step_type}")
            if step_type in options:
                node._options = options[step_type]
            if "SEAMM" in options:
                node._global_options = options["SEAMM"]

        # Write out an initial summary of the flowchart before doing anything
        # Reset the visited flag for traversal
        flowchart.reset_visited()

        # Get the start node
        next_node = flowchart.get_node("1")

        # describe ourselves
        printer.normal("\nDescription of the flowchart\n----------------------------")

        while next_node:
            # and print the description
            try:
                next_node = next_node.describe()
            except Exception:
                message = "Error describing the flowchart\n\n" + traceback.format_exc()
                print(message)
                logger.critical(message)
                raise
            except:  # noqa: E722
                message = (
                    "Unexpected error describing the flowchart\n\n"
                    + traceback.format_exc()
                )
                print(message)
                logger.critical(message)
                raise

        printer.normal("")

        # And actually run it!
        printer.normal(("Running the flowchart\n" "---------------------"))

        try:
            next_node = flowchart.get_node("1")
            while next_node is not None:
                try:
                    next_node = next_node.run()
                except DeprecationWarning as e:
                    print("\nDeprecation warning: " + str(e))
                    traceback.print_exc(file=sys.stderr)
                    traceback.print_exc(file=sys.stdout)
        finally:
            pass

        # Get the Hessian matrix and return it.
        triangle = []
        with open(directory / "1" / "hessian.dat") as fd:
            lines = fd.read().splitlines()

        if lines[0] != "!molssi hessian 1.0":
            raise IOError("Not a hessian file!")

        for line in lines[1:]:
            if line[0] == "@":
                pass
            else:
                triangle.append([float(val) for val in line.split()])

        # Set output back
        _file_handler.close()
        job.removeHandler(_file_handler)
        if job_handler is not None:
            job_handler.setLevel(job_level)
        if out_handler is not None:
            out_handler.setLevel(out_level)

        return triangle

    def _hessian(self, x, *args):
        """Return the Hessian matrix to SciPY."""
        n = len(self._triangle)
        matrix = np.zeros((n, n))
        i = 0
        for i, row in enumerate(self._triangle):
            for j, val in enumerate(row):
                matrix[i, j] = val
                matrix[j, i] = val

        return matrix

    def _run_psi4(self, directory):
        self._counters["main"] += 1
        directory.mkdir(parents=True, exist_ok=True)
        _, configuration = self.get_system_configuration()

        # Add the structure
        structure = self.parent._convert_structure(name="initial")

        # And run the calculation
        files = {}
        if self._restart_wfn is not None:
            files["old_wfn.npy"] = self._restart_wfn
            self.prolog, self.epilog = self.input_text(
                self._node0, self._memory, restart="old_wfn.npy"
            )
        else:
            self.prolog, self.epilog = self.input_text(self._node0, self._memory)
        files["input.dat"] = self.prolog + structure + self.epilog

        return_files = ["output.dat", "wfn.npy", "final_structure.json"]

        exe_path = Path(self.parent.options["psi4_path"])
        env = {
            "PSIPATH": str(exe_path),
            "PATH": str(exe_path),
        }

        local = seamm.ExecLocal()
        exe = exe_path / "psi4"
        result = local.run(
            cmd=[str(exe), f"-n {self.n_threads}"],
            files=files,
            return_files=return_files,
            env=env,
            in_situ=True,
            directory=directory,
        )

        if result is None:
            self.logger.error("There was an error running Psi4")
            return None

        # Capture the wavefunction for restarts
        self._restart_wfn = result["wfn.npy"]["data"]

        json_file = directory / "final_structure.json"
        if json_file.exists():
            with json_file.open() as fd:
                data = json.load(fd)
        else:
            raise RuntimeError("Structure and gradients from Psi4 not available.")

        self._E = data["energy"]
        self._derivatives = np.array(data["gradient"]).flatten()

        return self._E, self._derivatives

    def _step(self, x):
        _, configuration = self.get_system_configuration()

        x0 = configuration.atoms.get_coordinates(as_array=True).flatten()

        if self.iteration != 0 or self._E is None:
            # Update the coordinates in the configuration
            xyz = x.reshape(configuration.n_atoms, 3)
            configuration.atoms.set_coordinates(xyz)

            # Add the structure
            structure = self.parent._convert_structure(name="initial")

            # Make directory to work in
            directory = Path(self.directory)
            subdirectory = directory / f"step_{self.iteration:03d}"
            subdirectory.mkdir(parents=True, exist_ok=True)

            # And run the calculation
            files = {"input.dat": self.prolog + structure + self.epilog}

            if self._restart_wfn is not None:
                files["old_wfn.npy"] = self._restart_wfn

            return_files = ["output.dat", "wfn.npy", "final_structure.json"]

            exe_path = Path(self.parent.options["psi4_path"])
            env = {
                "PSIPATH": str(exe_path),
                "PATH": str(exe_path),
            }

            local = seamm.ExecLocal()
            exe = exe_path / "psi4"
            result = local.run(
                cmd=[str(exe), f"-n {self.n_threads}"],
                files=files,
                return_files=return_files,
                env=env,
                in_situ=True,
                directory=subdirectory,
            )

            if result is None:
                self.logger.error("There was an error running Psi4")
                return None

            # Capture the wavefunction for restarts
            self._restart_wfn = result["wfn.npy"]["data"]

            json_file = subdirectory / "final_structure.json"
            if json_file.exists():
                with json_file.open() as fd:
                    data = json.load(fd)
            else:
                raise RuntimeError("Structure and gradients from Psi4 not available.")

            self._E = data["energy"]
            self._derivatives = np.array(data["gradient"]).flatten()

            # Run MOPAC
            # self._triangle = self._run_MOPAC_Hessian(subdirectory / "MOPAC")

            if self.iteration == 0:
                text = " It    Energy (Ha)   RMS deriv  Max deriv  RMS xyz    Max xyz"
                printer.normal(__(text, indent=self.indent + 4 * " "))
                text = "--- --------------- ---------- ---------- ---------- ----------"
                printer.normal(__(text, indent=self.indent + 4 * " "))

            rmsd = np.sqrt(np.mean(np.square(self._derivatives)))
            maxd = np.amax(np.abs(self._derivatives))
            if self.iteration > 0:
                rms_step = np.sqrt(np.mean(np.square(x - x0)))
                max_step = np.amax(np.abs(x - x0))
                text = (
                    f"{self.iteration:3d} {self._E:15.7f} {rmsd:10.6f} {maxd:10.6f} "
                    f"{rms_step:10.6f} {max_step:10.6f}"
                )
            else:
                text = f"{self.iteration:3d} {self._E:15.7f} {rmsd:10.6f} {maxd:10.6f}"
            printer.normal(self.indent + 4 * " " + text)

        self.iteration += 1

        return self._E, self._derivatives

    def _bfgs(
        self,
        fun,
        x0,
        args=(),
        jac=None,
        callback=None,
        gtol=1e-5,
        norm=np.Inf,
        eps=_epsilon,
        maxiter=None,
        disp=False,
        return_all=False,
        finite_diff_rel_step=None,
        xrtol=0,
        **unknown_options,
    ):
        """
        Minimization of scalar function of one or more variables using the
        BFGS algorithm.

        Options
        -------
        disp : bool
            Set to True to print convergence messages.
        maxiter : int
            Maximum number of iterations to perform.
        gtol : float
            Terminate successfully if gradient norm is less than `gtol`.
        norm : float
            Order of norm (Inf is max, -Inf is min).
        eps : float or ndarray
            If `jac is None` the absolute step size used for numerical
            approximation of the jacobian via forward differences.
        return_all : bool, optional
            Set to True to return a list of the best solution at each of the
            iterations.
        finite_diff_rel_step : None or array_like, optional
            If `jac in ['2-point', '3-point', 'cs']` the relative step size to
            use for numerical approximation of the jacobian. The absolute step
            size is computed as ``h = rel_step * sign(x) * max(1, abs(x))``,
            possibly adjusted to fit into the bounds. For ``method='3-point'``
            the sign of `h` is ignored. If None (default) then step is selected
            automatically.
        xrtol : float, default: 0
            Relative tolerance for `x`. Terminate successfully if step size is
            less than ``xk * xrtol`` where ``xk`` is the current parameter vector.
        """
        # _check_unknown_options(unknown_options)
        retall = return_all

        x0 = np.asarray(x0).flatten()
        if x0.ndim == 0:
            x0.shape = (1,)
        if maxiter is None:
            maxiter = len(x0) * 200

        sf = _prepare_scalar_function(
            fun,
            x0,
            jac,
            args=args,
            epsilon=eps,
            finite_diff_rel_step=finite_diff_rel_step,
        )

        f = sf.fun
        myfprime = sf.grad

        old_fval = f(x0)
        gfk = myfprime(x0)

        k = 0
        N = len(x0)
        I = np.eye(N, dtype=int)  # noqa: E741

        H = self._hessian(None)
        Hk = np.linalg.pinv(H, rcond=self._rcond, hermitian=True)
        # Hk = I

        # Sets the initial step guess to dx ~ 1
        old_old_fval = old_fval + np.linalg.norm(gfk) / 2
        old_old_fval = None

        xk = x0
        if retall:
            allvecs = [x0]
        warnflag = 0
        gnorm = vecnorm(gfk, ord=norm)
        while (gnorm > gtol) and (k < maxiter):
            pk = -np.dot(Hk, gfk)
            print(f"{old_fval=} {old_old_fval=}")
            try:
                alpha_k, fc, gc, old_fval, old_old_fval, gfkp1 = _line_search_wolfe12(
                    f,
                    myfprime,
                    xk,
                    pk,
                    gfk,
                    old_fval,
                    old_old_fval,
                    amin=1e-100,
                    amax=1e100,
                )
            except _LineSearchError:
                # Line search failed to find a better solution.
                warnflag = 2
                break

            print(f"{old_fval=} {old_old_fval=} {alpha_k=}")
            sk = alpha_k * pk
            xkp1 = xk + sk

            if retall:
                allvecs.append(xkp1)
            xk = xkp1
            if gfkp1 is None:
                gfkp1 = myfprime(xkp1)

            yk = gfkp1 - gfk
            gfk = gfkp1
            if callback is not None:
                callback(xk)
            k += 1
            gnorm = vecnorm(gfk, ord=norm)
            if gnorm <= gtol:
                break

            #  See Chapter 5 in  P.E. Frandsen, K. Jonasson, H.B. Nielsen,
            #  O. Tingleff: "Unconstrained Optimization", IMM, DTU.  1999.
            #  These notes are available here:
            #  http://www2.imm.dtu.dk/documents/ftp/publlec.html
            if alpha_k * vecnorm(pk) <= xrtol * (xrtol + vecnorm(xk)):
                break

            if not np.isfinite(old_fval):
                # We correctly found +-Inf as optimal value, or something went
                # wrong.
                warnflag = 2
                break

            rhok_inv = np.dot(yk, sk)
            # this was handled in numeric, let it remaines for more safety
            # Cryptic comment above is preserved for posterity. Future reader:
            # consider change to condition below proposed in gh-1261/gh-17345.
            if rhok_inv == 0.0:
                rhok = 1000.0
                if disp:
                    print("Divide-by-zero encountered: rhok assumed large")
            else:
                rhok = 1.0 / rhok_inv

            A1 = I - sk[:, np.newaxis] * yk[np.newaxis, :] * rhok
            A2 = I - yk[:, np.newaxis] * sk[np.newaxis, :] * rhok
            Hk = np.dot(A1, np.dot(Hk, A2)) + (
                rhok * sk[:, np.newaxis] * sk[np.newaxis, :]
            )

        fval = old_fval

        if warnflag == 2:
            msg = _status_message["pr_loss"]
        elif k >= maxiter:
            warnflag = 1
            msg = _status_message["maxiter"]
        elif np.isnan(gnorm) or np.isnan(fval) or np.isnan(xk).any():
            warnflag = 3
            msg = _status_message["nan"]
        else:
            msg = _status_message["success"]

        if disp:
            print("%s%s" % ("Warning: " if warnflag != 0 else "", msg))
            print("         Current function value: %f" % fval)
            print("         Iterations: %d" % k)
            print("         Function evaluations: %d" % sf.nfev)
            print("         Gradient evaluations: %d" % sf.ngev)

        result = OptimizeResult(
            fun=fval,
            jac=gfk,
            hess_inv=Hk,
            nfev=sf.nfev,
            njev=sf.ngev,
            status=warnflag,
            success=(warnflag == 0),
            message=msg,
            x=xk,
            nit=k,
        )
        if retall:
            result["allvecs"] = allvecs
        return result

    def runPsi4(self, x, needed=("E", "dE")):
        self._counters["main"] += 1

        directory = Path(self.directory) / f"Psi4_{self._counters['main']}"
        directory.mkdir(parents=True, exist_ok=True)

        _, configuration = self.get_system_configuration()

        x_orig = configuration.atoms.get_coordinates()
        configuration.atoms.set_coordinates(x.reshape(configuration.n_atoms, 3))

        # Add the structure
        structure = self.parent._convert_structure(name="initial")

        configuration.atoms.set_coordinates(x_orig)

        # And run the calculation
        files = {}
        if self._restart_wfn is not None:
            files["old_wfn.npy"] = self._restart_wfn
            self.prolog, self.epilog = self.input_text(
                self._node0, self._memory, restart="old_wfn.npy"
            )
        else:
            self.prolog, self.epilog = self.input_text(self._node0, self._memory)
        files["input.dat"] = self.prolog + structure + self.epilog

        return_files = ["output.dat", "wfn.npy", "final_structure.json"]

        exe_path = Path(self.parent.options["psi4_path"])
        env = {
            "PSIPATH": str(exe_path),
            "PATH": str(exe_path),
        }

        local = seamm.ExecLocal()
        exe = exe_path / "psi4"
        result = local.run(
            cmd=[str(exe), f"-n {self.n_threads}"],
            files=files,
            return_files=return_files,
            env=env,
            in_situ=True,
            directory=directory,
        )

        if result is None:
            self.logger.error("There was an error running Psi4")
            return None

        # Capture the wavefunction for restarts
        self._restart_wfn = result["wfn.npy"]["data"]

        json_file = directory / "final_structure.json"
        if json_file.exists():
            with json_file.open() as fd:
                data = json.load(fd)
        else:
            raise RuntimeError("Structure and gradients from Psi4 not available.")

        E = data["energy"]
        dE = np.array(data["gradient"]).flatten() * Q_("E_h/a_0").m_as("E_h/Å")

        return E, dE

    def runMOPAC(self, x, needed=("E", "dE")):
        """Run the MOPAC flowchart in the given directory."""
        self._counters["surrogate"] += 1

        directory = Path(self.directory) / f"MOPAC_{self._counters['surrogate']}"
        directory.mkdir(parents=True, exist_ok=True)

        # Set the coordinates
        _, configuration = self.get_system_configuration()
        x_orig = configuration.atoms.get_coordinates()
        configuration.atoms.set_coordinates(x.reshape(configuration.n_atoms, 3))

        # Get the MOPAC flowchart
        path = Path(pkg_resources.resource_filename(__name__, "data/"))
        flowchart = seamm.Flowchart(
            name="MOPAC",
            parent=self,
            namespace="org.molssi.seamm",
            directory=directory,
            parser_name="MOPAC subflowchart",
        )
        flowchart.read(path / "MOPAC_PM7_surrogate.flow")

        # Find the handler for job.out and set the level up
        job_handler = None
        out_handler = None
        for handler in job.handlers:
            if (
                isinstance(handler, logging.FileHandler)
                and "job.out" in handler.baseFilename
            ):
                job_handler = handler
                job_level = job_handler.level
                job_handler.setLevel(printing.JOB)
            elif isinstance(handler, logging.StreamHandler):
                out_handler = handler
                out_level = out_handler.level
                out_handler.setLevel(printing.JOB)

        # Redirect output
        _file_handler = logging.FileHandler(directory / "subflowchart.out")
        _file_handler.setLevel(printing.NORMAL)
        formatter = logging.Formatter(fmt="{message:s}", style="{")
        _file_handler.setFormatter(formatter)
        job.addHandler(_file_handler)

        # Set up the argument parser for this node.
        parser = flowchart.parser

        parser.epilog = textwrap.dedent(
            """
            The plug-ins in this flowchart are listed above.
            Options, if any, for plug-ins are placed after
            the name of the plug-in, e.g.:

               test.flow lammps-step --log-level DEBUG --np 4

            To get help for a plug-in, use --help or -h after the
            plug-in name. E.g.

               test.flow lammps-step --help
            """
        )
        parser.usage = "%(prog)s [options] plug-in [options] plug-in [options] ..."

        # Now traverse the flowchart, setting up the ids and parsers
        flowchart.set_ids()
        flowchart.create_parsers()

        # And handle the command-line arguments and ini file options.
        parser.parse_args([])
        logger.info("Parsed the command-line arguments")

        # Get the command line options
        options = parser.get_options()

        # Set the options in each step
        for node in flowchart:
            step_type = node.step_type
            logger.info(f"    setting options for {step_type}")
            if step_type in options:
                node._options = options[step_type]
            if "SEAMM" in options:
                node._global_options = options["SEAMM"]

        # Write out an initial summary of the flowchart before doing anything
        # Reset the visited flag for traversal
        flowchart.reset_visited()

        # Get the start node
        next_node = flowchart.get_node("1")

        # describe ourselves
        printer.normal("\nDescription of the flowchart\n----------------------------")

        while next_node:
            # and print the description
            try:
                next_node = next_node.describe()
            except Exception:
                message = "Error describing the flowchart\n\n" + traceback.format_exc()
                print(message)
                logger.critical(message)
                raise
            except:  # noqa: E722
                message = (
                    "Unexpected error describing the flowchart\n\n"
                    + traceback.format_exc()
                )
                print(message)
                logger.critical(message)
                raise

        printer.normal("")

        # And actually run it!
        printer.normal(("Running the flowchart\n" "---------------------"))

        try:
            next_node = flowchart.get_node("1")
            while next_node is not None:
                try:
                    next_node = next_node.run()
                except DeprecationWarning as e:
                    print("\nDeprecation warning: " + str(e))
                    traceback.print_exc(file=sys.stderr)
                    traceback.print_exc(file=sys.stdout)
        finally:
            pass

        # Set the coordinates back
        configuration.atoms.set_coordinates(x_orig)

        # Get the derivatives and return them.
        E = self.get_variable("_mopac_energy")
        gradients = self.get_variable("_mopac_gradients")
        dE = np.array(gradients).flatten()

        # Get the Hessian matrix and return it.
        triangle = []
        with open(directory / "1" / "hessian.dat") as fd:
            lines = fd.read().splitlines()

        if lines[0] != "!molssi hessian 1.0":
            raise IOError("Not a hessian file!")

        for line in lines[1:]:
            if line[0] == "@":
                pass
            else:
                triangle.append([float(val) for val in line.split()])

        n = len(triangle)
        d2E = np.zeros((n, n))
        i = 0
        for i, row in enumerate(triangle):
            for j, val in enumerate(row):
                d2E[i, j] = val
                d2E[j, i] = val

        # Set output back
        _file_handler.close()
        job.removeHandler(_file_handler)
        if job_handler is not None:
            job_handler.setLevel(job_level)
        if out_handler is not None:
            out_handler.setLevel(out_level)

        if "d2E" in needed:
            return E, dE, d2E
        else:
            return E, dE
