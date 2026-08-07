# -*- coding: utf-8 -*-

"""The Psi4 counterpoise (BSSE) sub-step.

Computes the counterpoise-corrected (Boys--Bernardi) energy and gradient of
an N-fragment complex, with independent per-fragment charge, as a thin
wrapper around Psi4's own native N-fragment counterpoise driver
(``bsse_type=['cp', 'nocp']`` on ``energy()``/``gradient()``). Unlike
``orca_step``'s BSSE sub-step, this does **not** use
``seamm_bsse.generate_job_specs()``/``combine()`` -- Psi4 already does that
2N + 1-equivalent orchestration internally and hands back an
already-full-cluster-indexed corrected energy/gradient in one call
(confirmed empirically; see the ``seamm_bsse`` design doc,
``docs/developer_guide/campaigns/2026-08-03/bsse_architecture.rst``).
``seamm_bsse.Fragment``/``validate_fragments`` are reused here only for the
fragment-definition/charge-validation layer, for parameter/GUI consistency
with the ORCA sub-step -- ``generate_job_specs``/``combine`` go unused.

Unlike the other ``psi4_step`` sub-steps, this one builds its own
fragment-aware, ``--``-separated molecule block: ``Psi4._convert_structure``
reads the single global "current" configuration, with no concept of
fragments/ghosts.
"""

import json
import logging
from pathlib import Path
import textwrap

from tabulate import tabulate

import psi4_step
import seamm
import seamm_bsse
from seamm_util import Q_
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter("psi4")

#: Hartree -> kcal/mol, for reporting the (small) correction in familiar units.
_HARTREE_TO_KCAL = 627.509474

#: Boys & Bernardi counterpoise method -- the correction this step applies.
_BSSE_CITATION = """\
@article{Boys1970,
    author = {Boys, S. F. and Bernardi, F.},
    title = {The calculation of small molecular interactions by the differences
             of separate total energies. Some procedures with reduced errors},
    journal = {Molecular Physics},
    volume = {19},
    number = {4},
    pages = {553--566},
    year = {1970},
    doi = {10.1080/00268977000101561}
}
"""


class BSSE(psi4_step.Energy):
    """A counterpoise (BSSE) correction with Psi4, wrapping its native
    N-fragment ``bsse_type`` driver.

    See Also
    --------
    TkBSSE, BSSEParameters, Energy
    """

    def __init__(self, flowchart=None, title="BSSE", extension=None, logger=logger):
        logger.debug(f"Creating Psi4 BSSE {self}")
        super().__init__(flowchart=flowchart, title=title, extension=extension)

        self._calculation = "bsse"
        self._model = None
        self._metadata = psi4_step.metadata
        self.parameters = psi4_step.BSSEParameters()
        self.description = "A counterpoise (BSSE) correction"

    def description_text(self, P=None, configuration=None):
        if not P:
            P = self.parameters.values_to_dict()
        text = super().description_text(
            P,
            calculation_type="Counterpoise (BSSE) corrected energy",
            configuration=configuration,
        )
        return text

    # ------------------------------------------------------------------
    # Fragments -- the same grammar as orca_step's BSSE ("auto (molecules)"
    # / "specified", "fragment atoms", "fragment charges"), own copy of the
    # parsing since it is small and engine-agnostic (a candidate to move
    # into seamm_bsse itself later, rather than duplicate a third time).
    # ------------------------------------------------------------------
    def _fragments(self, P, configuration):
        """Return the N ``seamm_bsse.Fragment`` for this complex, per the
        'fragments' parameter, with per-fragment charge threaded in from
        'fragment charges' -- or, if that is left empty, from the sum of each
        fragment's atoms' ``formal_charge`` (set on ``configuration.atoms``
        when the structure came from a format that carries per-atom formal
        charges, e.g. an SDF/MOL file's ``M  CHG`` record), falling back to
        all-neutral if neither is available. Validates against the complex's
        own charge/multiplicity before returning."""
        n_atoms = configuration.n_atoms
        mode = P["fragments"]
        if mode == "specified":
            groups = self._parse_fragment_groups(P["fragment atoms"], n_atoms)
        else:
            # "auto (molecules)"
            molecules = configuration.find_molecules(as_indices=True)
            if len(molecules) < 2:
                raise RuntimeError(
                    f"BSSE 'auto' fragments require at least two separate "
                    f"molecules, but the structure has {len(molecules)}. Use "
                    "'specified' to define the fragments by atom, or check "
                    "the bonding."
                )
            groups = [sorted(molecule) for molecule in molecules]

        charges = self._parse_charges(P["fragment charges"])
        if not charges:
            if "formal_charge" in configuration.atoms:
                formal_charges = configuration.atoms["formal_charge"]
                charges = [sum(formal_charges[i] for i in group) for group in groups]
                printer.important(
                    __(
                        "No 'fragment charges' given; using each fragment's "
                        "net formal charge from the input structure.",
                        indent=self.indent + 4 * " ",
                    )
                )
            else:
                charges = [0] * len(groups)
        elif len(charges) != len(groups):
            raise RuntimeError(
                f"BSSE: {len(charges)} 'fragment charges' given but "
                f"{len(groups)} fragments found/specified; give one charge "
                "per fragment (in the same order), or leave 'fragment "
                "charges' empty for all-neutral."
            )

        fragments = [
            seamm_bsse.Fragment(label=str(i + 1), atom_indices=group, charge=charge)
            for i, (group, charge) in enumerate(zip(groups, charges))
        ]
        try:
            seamm_bsse.validate_fragments(
                fragments,
                cluster_charge=configuration.charge,
                cluster_multiplicity=configuration.spin_multiplicity,
                atomic_numbers=configuration.atoms.atomic_numbers,
            )
        except ValueError as e:
            raise RuntimeError(f"BSSE: {e}") from e
        return fragments

    @staticmethod
    def _parse_fragment_groups(text, n_atoms):
        """Parse 'specified'-mode fragment atoms: semicolon-separated
        1-based index/range groups, one per fragment, e.g. ``'1-3; 4-6; 7'``
        for three fragments -- to 0-based index lists."""
        groups = [group.strip() for group in str(text).split(";") if group.strip()]
        if len(groups) < 2:
            raise RuntimeError(
                "BSSE: 'Fragment atoms' (specified mode) needs at least two "
                "semicolon-separated groups of atoms, e.g. '1-3; 4-6' for two "
                "fragments."
            )
        return [BSSE._parse_indices(group, n_atoms) for group in groups]

    @staticmethod
    def _parse_charges(text):
        """Parse 'fragment charges': a comma/space separated list of
        integers, in fragment order. '' (the default, all-neutral) -> []."""
        text = str(text).strip()
        if not text:
            return []
        return [int(token) for token in text.replace(",", " ").split()]

    @staticmethod
    def _parse_indices(text, n_atoms):
        """Parse a 1-based index/range list (``'1-3, 5 7'``) to 0-based
        indices."""
        indices = []
        for token in str(text).replace(",", " ").split():
            if "-" in token[1:]:  # a range like 1-3 (not a leading minus)
                lo, hi = token.split("-", 1)
                indices.extend(range(int(lo), int(hi) + 1))
            elif token:
                indices.append(int(token))
        out = []
        for i in indices:
            j = i - 1
            if j < 0 or j >= n_atoms:
                raise RuntimeError(
                    f"BSSE: atom index {i} is out of range (1..{n_atoms})."
                )
            if j not in out:
                out.append(j)
        return out

    # ------------------------------------------------------------------
    # Input generation
    # ------------------------------------------------------------------
    def _molecule_block(self, fragments, configuration, name="bsse_cluster"):
        """The Psi4 ``--``-separated, per-fragment charge/multiplicity
        molecule block for `fragments`. ``no_com``/``no_reorient`` keep the
        geometry (and so the gradient) in the original frame, matching
        ``Psi4._convert_structure``'s convention -- essential here, since the
        returned gradient must map back onto `configuration`'s own atom
        order."""
        symbols = configuration.atoms.symbols
        xyzs = configuration.atoms.get_coordinates(fractionals=False, in_cell=True)

        lines = [f"molecule {name} {{"]
        for i, fragment in enumerate(fragments):
            if i > 0:
                lines.append("--")
            lines.append(f"{fragment.charge} {fragment.multiplicity}")
            for atom_index in fragment.atom_indices:
                symbol = symbols[atom_index]
                x, y, z = xyzs[atom_index]
                lines.append(f"{symbol:2s} {x:15.8f} {y:15.8f} {z:15.8f}")
        lines.append("")
        lines.append("no_com")
        lines.append("no_reorient")
        lines.append("}")
        lines.append(f"activate({name})")
        return "\n".join(lines)

    def get_input(self):
        """The Psi4 input for the counterpoise correction: a fragment-aware
        molecule block, then a ``bsse_type=['cp', 'nocp']``
        ``energy()``/``gradient()`` call -- Psi4's native N-fragment
        counterpoise driver, not ``seamm_bsse.generate_job_specs()``/
        ``combine()`` (which ``orca_step`` needs, but Psi4 does not: 'nocp'
        alongside 'cp' costs no extra energy evaluations, since both reuse
        the same underlying 1-body/N-body sub-calculations, and gives the
        uncorrected total/interaction energies orca_step also reports)."""
        directory = Path(self.directory)
        directory.mkdir(parents=True, exist_ok=True)

        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )

        _, configuration = self.get_system_configuration(None)
        self.description = []
        self.description.append(
            __(
                self.description_text(P, configuration=configuration),
                indent=self.indent,
            )
        )

        fragments = self._fragments(P, configuration)

        method, functional, extended_functional, method_string = self.get_method()
        psi4_method = extended_functional if method == "dft" else method

        if self.parent.basis is not None:
            basis_set = self.parent.basis
            if method == "dft":
                self.parent.model = f"{functional.upper()}/{basis_set}"
            else:
                self.parent.model = f"{method.upper()}/{basis_set}"

        want_gradient = P["compute gradient"] == "yes"

        multiplicity = configuration.spin_multiplicity
        reference = "rks" if method == "dft" else "rhf"
        if multiplicity != 1:
            reference = "u" + reference[1:]

        lines = []
        lines.append("")
        lines.append("#" * 80)
        lines.append(f"# {self.header}")
        lines.append("#" * 80)
        lines.append("")
        lines.append(self._molecule_block(fragments, configuration))
        lines.append("")
        lines.append(f"set reference {reference}")
        if P.get("freeze-cores", "yes") == "yes":
            lines.append("set freeze_core True")
        else:
            lines.append("set freeze_core False")
        lines.append("set fail_on_maxiter False")  # Psi4 hangs! Cope in the plug-in
        lines.append(f"set maxiter {P['maximum iterations']}")
        if P["density convergence"] != "default":
            lines.append(f"set d_convergence {P['density convergence']}")
        if P["energy convergence"] != "default":
            lines.append(f"set e_convergence {P['energy convergence']}")
        lines.append("")

        call = "gradient" if want_gradient else "energy"
        lines.append(
            f"E, wfn = {call}('{psi4_method}', bsse_type=['cp', 'nocp'], "
            "return_wfn=True)"
        )

        # Absolute path: this node's own directory, not the shared psi4
        # process's cwd (the main Psi4 node's directory), which is where
        # every sub-step's input.dat actually runs.
        result_path = (Path(self.directory) / "bsse.json").as_posix()
        lines.append("""
variables = wfn.scalar_variables()
result = {
    "uncorrected energy": variables["NOCP-CORRECTED TOTAL ENERGY"],
    "energy": variables["CP-CORRECTED TOTAL ENERGY"],
    "bsse correction": (
        variables["CP-CORRECTED TOTAL ENERGY"]
        - variables["NOCP-CORRECTED TOTAL ENERGY"]
    ),
    "interaction energy": variables["CP-CORRECTED INTERACTION ENERGY"],
    "uncorrected interaction energy": variables["NOCP-CORRECTED INTERACTION ENERGY"],
}""")
        if want_gradient:
            lines.append('result["gradients"] = wfn.gradient().to_array().tolist()')
        lines.append(f"""
with open(r"{result_path}", "w") as fd:
    json.dump(result, fd, indent=2)
""")

        return "\n".join(lines)

    def _cite_bsse(self):
        """Cite the counterpoise method (Boys & Bernardi). Best-effort."""
        try:
            self.references.cite(
                raw=_BSSE_CITATION,
                alias="boys-bernardi-1970",
                module="psi4_step",
                level=1,
                note="The counterpoise correction for BSSE.",
            )
        except Exception as e:  # pragma: no cover
            logger.warning(f"Could not cite the counterpoise method: {e}")

    # ------------------------------------------------------------------
    # Analysis
    # ------------------------------------------------------------------
    def analyze(self, indent="", data=None, out=[], table=None):
        """Read the counterpoise-corrected results (written by
        ``get_input()``'s Psi4-side code to ``<self.directory>/bsse.json``,
        since the Psi4 subprocess runs in the shared main node's directory,
        not this sub-step's own) and store/report them."""
        if data is None:
            json_file = Path(self.directory) / "bsse.json"
            if not json_file.exists():
                text = (
                    "\nThere are no results from the Psi4 BSSE correction. "
                    f"Perhaps it failed? Looking for {json_file}."
                )
                printer.normal(__(text, indent=self.indent + 4 * " "))
                raise RuntimeError(text)
            with json_file.open() as fd:
                data = json.load(fd)

        # Interaction energies are more useful in kJ/mol, the conventional
        # unit for a binding/interaction energy (matches orca_step's BSSE).
        props = dict(data)
        for key in ("interaction energy", "uncorrected interaction energy"):
            if key in props:
                props[key] = Q_(props[key], "E_h").m_as("kJ/mol")

        self._cite_bsse()

        try:
            self.store_results(
                data=props,
                create_tables=self.parameters["create tables"].get(),
            )
        except Exception as e:  # pragma: no cover
            logger.warning(f"Could not store results: {e}")

        rows = []
        if "interaction energy" in props:
            rows.append(
                [
                    "Interaction energy (CP-corrected)",
                    f"{props['interaction energy']:.4f}",
                    "kJ/mol",
                ]
            )
        if "uncorrected interaction energy" in props:
            rows.append(
                [
                    "Interaction energy (uncorrected)",
                    f"{props['uncorrected interaction energy']:.4f}",
                    "kJ/mol",
                ]
            )
        if "uncorrected energy" in props:
            rows.append(
                ["Uncorrected energy", f"{props['uncorrected energy']:.8f}", "E_h"]
            )
        if "energy" in props:
            rows.append(["BSSE-corrected energy", f"{props['energy']:.8f}", "E_h"])
        if "bsse correction" in props:
            corr = props["bsse correction"]
            rows.append(["BSSE correction", f"{corr:.8f}", "E_h"])
            rows.append(
                ["BSSE correction", f"{corr * _HARTREE_TO_KCAL:.4f}", "kcal/mol"]
            )
        if rows:
            tmp = tabulate(
                rows,
                headers=["Property", "Value", "Units"],
                tablefmt="rounded_outline",
                colalign=("left", "right", "left"),
                disable_numparse=True,
            )
            printer.normal("")
            printer.normal(textwrap.indent(tmp, self.indent + 7 * " "))

        grad = props.get("gradients")
        if grad is not None:
            printer.normal(
                __(
                    f"Computed the counterpoise-corrected gradient on {len(grad)} "
                    "atoms.",
                    indent=self.indent + 4 * " ",
                )
            )
