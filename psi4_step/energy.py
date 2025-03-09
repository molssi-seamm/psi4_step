# -*- coding: utf-8 -*-

"""Setup and run Psi4"""

from collections import Counter
import json
import logging
from math import isnan
from pathlib import Path
import pkg_resources
import pprint
import textwrap

import cclib
from openbabel import openbabel
import pandas
import numpy as np
from tabulate import tabulate

from molsystem import elements
import psi4_step
import seamm
import seamm.data
from seamm_util import Q_, units_class
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter("psi4")

_subscript = {
    "0": "\N{SUBSCRIPT ZERO}",
    "1": "\N{SUBSCRIPT ONE}",
    "2": "\N{SUBSCRIPT TWO}",
    "3": "\N{SUBSCRIPT THREE}",
    "4": "\N{SUBSCRIPT FOUR}",
    "5": "\N{SUBSCRIPT FIVE}",
    "6": "\N{SUBSCRIPT SIX}",
    "7": "\N{SUBSCRIPT SEVEN}",
    "8": "\N{SUBSCRIPT EIGHT}",
    "9": "\N{SUBSCRIPT NINE}",
}

superscript = {
    "1": "\N{SUPERSCRIPT ONE}",
    "2": "\N{SUPERSCRIPT TWO}",
    "3": "\N{SUPERSCRIPT THREE}",
    "4": "\N{SUPERSCRIPT FOUR}",
    "5": "\N{SUPERSCRIPT FIVE}",
    "6": "\N{SUPERSCRIPT SIX}",
    "7": "\N{SUPERSCRIPT SEVEN}",
    "8": "\N{SUPERSCRIPT EIGHT}",
    "9": "\N{SUPERSCRIPT NINE}",
}


def subscript(n):
    """Return the number using Unicode subscript characters."""
    return "".join([_subscript[c] for c in str(n)])


one_half = "\N{VULGAR FRACTION ONE HALF}"
degree_sign = "\N{DEGREE SIGN}"
standard_state = {
    "H": f"{one_half}H{subscript(2)}(g)",
    "He": "He(g)",
    "Li": "Li(s)",
    "Be": "Be(s)",
    "B": "B(s)",
    "C": "C(s,gr)",
    "N": f"{one_half}N{subscript(2)}(g)",
    "O": f"{one_half}O{subscript(2)}(g)",
    "F": f"{one_half}F{subscript(2)}(g)",
    "Ne": "Ne(g)",
    "Na": "Na(s)",
    "Mg": "Mg(s)",
    "Al": "Al(s)",
    "Si": "Si(s)",
    "P": "P(s)",
    "S": "S(s)",
    "Cl": f"{one_half}Cl{subscript(2)}(g)",
    "Ar": "Ar(g)",
    "K": "K(s)",
    "Ca": "Ca(s)",
    "Sc": "Sc(s)",
    "Ti": "Ti(s)",
    "V": "V(s)",
    "Cr": "Cr(s)",
    "Mn": "Mn(s)",
    "Fe": "Fe(s)",
    "Co": "Co(s)",
    "Ni": "Ni(s)",
    "Cu": "Cu(s)",
    "Zn": "Zn(s)",
    "Ga": "Ga(s)",
    "Ge": "Ge(s)",
    "As": "As(s)",
    "Se": "Se(s)",
    "Br": f"{one_half}Br{subscript(2)}(l)",
    "Kr": "(g)",
}


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

    @property
    def model(self):
        """The model chemistry"""
        return self._model

    @model.setter
    def model(self, value):
        self._model = value

    def calculate_enthalpy_of_formation(self, data):
        """Calculate the enthalpy of formation from the results of a calculation.

        This uses tabulated values of the enthalpy of formation of the atoms for
        the elements and tabulated energies calculated for atoms with the current
        method.

        Parameters
        ----------
        data : dict
            The results of the calculation.
        """

        # Read the tabulated values from either user or data directory
        personal_file = Path("~/.seamm.d/data/atom_energies.csv").expanduser()
        if personal_file.exists():
            personal_table = pandas.read_csv(personal_file, index_col=False)
        else:
            personal_table = None

        path = Path(pkg_resources.resource_filename(__name__, "data/"))
        csv_file = path / "atom_energies.csv"
        table = pandas.read_csv(csv_file, index_col=False)

        self.logger.debug(f"self.parent.model = {self.parent.model}")

        # Check if have the data
        atom_energies = None
        correction_energy = None
        if self.parent.model.startswith("U") or self.parent.model.startswith("R"):
            column = self.parent.model[1:]
        else:
            column = self.parent.model

        self.logger.debug(f"Looking for '{column}'")

        column2 = column + " correction"
        if personal_table is not None and column in personal_table.columns:
            atom_energies = personal_table[column].to_list()
            if column2 in personal_table.columns:
                correction_energy = personal_table[column2].to_list()
        elif column in table.columns:
            atom_energies = table[column].to_list()
            if column2 in table.columns:
                correction_energy = table[column2].to_list()

        if atom_energies is None:
            self.logger.debug("     and didn't find it!")

        DfH0gas = None
        references = None
        term_symbols = None
        if personal_table is not None and "ΔfH°gas" in personal_table.columns:
            DfH0gas = personal_table["ΔfH°gas"].to_list()
            if "Reference" in personal_table.columns:
                references = personal_table["Reference"].to_list()
            if "Term Symbol" in personal_table.columns:
                term_symbols = personal_table["Term Symbols"].to_list()
        elif "ΔfH°gas" in table.columns:
            DfH0gas = table["ΔfH°gas"].to_list()
            if "Reference" in table.columns:
                references = table["Reference"].to_list()
            if "Term Symbol" in table.columns:
                term_symbols = table["Term Symbol"].to_list()
        if references is not None:
            len(references)

        if atom_energies is None:
            return f"There are no tabulated atom energies for {column}"

        # Get the atomic numbers and counts
        _, configuration = self.get_system_configuration(None)
        counts = Counter(configuration.atoms.atomic_numbers)

        # Get the Hill formula as a list
        symbols = sorted(elements.to_symbols(counts.keys()))
        composition = []
        if "C" in symbols:
            composition.append((6, "C", counts[6]))
            symbols.remove("C")
            if "H" in symbols:
                composition.append((1, "H", counts[1]))
                symbols.remove("H")

        for symbol in symbols:
            atno = elements.symbol_to_atno[symbol]
            composition.append((atno, symbol, counts[atno]))

        # And the reactions. First, for atomization energy
        middot = "\N{MIDDLE DOT}"
        lDelta = "\N{GREEK CAPITAL LETTER DELTA}"
        formula = ""
        tmp = []
        for atno, symbol, count in composition:
            if count == 1:
                formula += symbol
                tmp.append(f"{symbol}(g)")
            else:
                formula += f"{symbol}{subscript(count)}"
                tmp.append(f"{count}{middot}{symbol}(g)")
        gas_atoms = " + ".join(tmp)
        tmp = []
        for atno, symbol, count in composition:
            if count == 1:
                tmp.append(standard_state[symbol])
            else:
                tmp.append(f"{count}{middot}{standard_state[symbol]}")
        standard_elements = " + ".join(tmp)

        # The atomization energy is the electronic energy minus the energy of the atoms
        try:
            name = "SMILES: " + configuration.canonical_smiles
            if name is None:
                name = "Formula: " + formula
        except Exception:
            name = "Formula: " + formula
        try:
            name = configuration.PC_iupac_name(fallback=name)
        except Exception:
            pass

        if name is None:
            name = "Formula: " + formula

        text = f"Thermochemistry of {name} with {column}\n\n"
        text += "Atomization Energy\n"
        text += "------------------\n"
        text += textwrap.fill(
            f"The atomization energy,  {lDelta}atE{degree_sign}, is the energy to break"
            " all the bonds in the system, separating the atoms from each other."
        )
        text += f"\n\n    {formula} --> {gas_atoms}\n\n"
        text += textwrap.fill(
            "The following table shows in detail the calculation. The first line is "
            "the system and its calculated energy. The next lines are the energies "
            "of each type of atom in the system. These have been tabulated by running "
            "calculations on each atom, and are included in the SEAMM release. "
            "The last two lines give the formation energy from atoms in atomic units "
            "and as kJ/mol.",
        )
        text += "\n\n"
        table = {
            "System": [],
            "Term": [],
            "Value": [],
            "Units": [],
        }

        E = data["energy"]

        Eatoms = 0.0
        for atno, symbol, count in composition:
            Eatom = atom_energies[atno - 1]
            if isnan(Eatom):
                # Don't have the data for this element
                return f"Do not have tabulated atom energies for {symbol} in {column}"
            Eatoms += count * Eatom
            tmp = Q_(Eatom, "kJ/mol").m_as("E_h")
            table["System"].append(f"{symbol}(g)")
            table["Term"].append(f"{count} * {tmp:.6f}")
            table["Value"].append(f"{count * tmp:.6f}")
            table["Units"].append("")

        table["Units"][0] = "E_h"

        table["System"].append("^")
        table["Term"].append("-")
        table["Value"].append("-")
        table["Units"].append("")

        table["System"].append(formula)
        table["Term"].append(f"{-E:.6f}")
        table["Value"].append(f"{-E:.6f}")
        table["Units"].append("E_h")

        data["E atomization"] = Eatoms - Q_(E, "E_h").m_as("kJ/mol")

        table["System"].append("")
        table["Term"].append("")
        table["Value"].append("=")
        table["Units"].append("")

        result = f'{Q_(data["E atomization"], "kJ/mol").m_as("E_h"):.6f}'
        table["System"].append(f"{lDelta}atE")
        table["Term"].append("")
        table["Value"].append(result)
        table["Units"].append("E_h")

        table["System"].append("")
        table["Term"].append("")
        table["Value"].append(f'{data["E atomization"]:.2f}')
        table["Units"].append("kJ/mol")

        tmp = tabulate(
            table,
            headers="keys",
            tablefmt="rounded_outline",
            colalign=("center", "center", "decimal", "center"),
            disable_numparse=True,
        )
        length = len(tmp.splitlines()[0])
        text_lines = []
        text_lines.append(f"Atomization Energy for {formula}".center(length))
        text_lines.append(tmp)
        text += textwrap.indent("\n".join(text_lines), 4 * " ")

        if "H_tot" not in data:
            text += "\n\n"
            text += "Cannot calculate enthalpy of formation without the enthalpy"
            return text
        if DfH0gas is None:
            text += "\n\n"
            text += "Cannot calculate enthalpy of formation without the tabulated\n"
            text += "atomization enthalpies of the elements."
            return text

        # Atomization enthalpy of the elements, experimental
        table = {
            "System": [],
            "Term": [],
            "Value": [],
            "Units": [],
            "Reference": [],
        }

        E = data["energy"]

        DfH_at = 0.0
        refno = 1
        for atno, symbol, count in composition:
            DfH_atom = DfH0gas[atno - 1]
            DfH_at += count * DfH_atom
            tmp = Q_(DfH_atom, "kJ/mol").m_as("E_h")
            table["System"].append(f"{symbol}(g)")
            if count == 1:
                table["Term"].append(f"{tmp:.6f}")
            else:
                table["Term"].append(f"{count} * {tmp:.6f}")
            table["Value"].append(f"{count * tmp:.6f}")
            table["Units"].append("")
            refno += 1
            table["Reference"].append(refno)

        table["Units"][0] = "E_h"

        table["System"].append("^")
        table["Term"].append("-")
        table["Value"].append("-")
        table["Units"].append("")
        table["Reference"].append("")

        table["System"].append(standard_elements)
        table["Term"].append("")
        table["Value"].append("0.0")
        table["Units"].append("E_h")
        table["Reference"].append("")

        table["System"].append("")
        table["Term"].append("")
        table["Value"].append("=")
        table["Units"].append("")
        table["Reference"].append("")

        result = f'{Q_(DfH_at, "kJ/mol").m_as("E_h"):.6f}'
        table["System"].append(f"{lDelta}atH{degree_sign}")
        table["Term"].append("")
        table["Value"].append(result)
        table["Units"].append("E_h")
        table["Reference"].append("")

        table["System"].append("")
        table["Term"].append("")
        table["Value"].append(f"{DfH_at:.2f}")
        table["Units"].append("kJ/mol")
        table["Reference"].append("")

        tmp = tabulate(
            table,
            headers="keys",
            tablefmt="rounded_outline",
            colalign=("center", "center", "decimal", "center", "center"),
            disable_numparse=True,
        )
        length = len(tmp.splitlines()[0])
        text_lines = []
        text_lines.append(
            "Atomization enthalpy of the elements (experimental)".center(length)
        )
        text_lines.append(tmp)

        text += "\n\n"
        text += "Enthalpy of Formation\n"
        text += "---------------------\n"
        text += textwrap.fill(
            f"The enthalpy of formation, {lDelta}fHº, is the enthalpy of creating the "
            "molecule from the elements in their standard state:"
        )
        text += f"\n\n   {standard_elements} --> {formula} (1)\n\n"
        text += textwrap.fill(
            "The standard state of the element, denoted by the superscript º,"
            " is its form at 298.15 K and 1 atm pressure, e.g. graphite for carbon, "
            "H2 gas for hydrogen, etc."
        )
        text += "\n\n"
        text += textwrap.fill(
            "Since it is not easy to calculate the enthalpy of e.g. graphite we will "
            "use two sequential reactions that are equivalent. First, we will create "
            "gas phase atoms from the elements:"
        )
        text += f"\n\n    {standard_elements} --> {gas_atoms} (2)\n\n"
        text += textwrap.fill(
            "This will use the experimental values of the enthalpy of formation of the "
            "atoms in the gas phase to calculate the enthalpy of this reaction. "
            "Then we react the atoms to get the desired system:"
        )
        text += f"\n\n    {gas_atoms} --> {formula} (3)\n\n"
        text += textwrap.fill(
            "Note that this is reverse of the atomization reaction, so "
            f"{lDelta}H = -{lDelta}atH."
        )
        text += "\n\n"
        text += textwrap.fill(
            "First we calculate the enthalpy of the atomization of the elements in "
            "their standard state, using tabulated experimental values:"
        )
        text += "\n\n"
        text += textwrap.indent("\n".join(text_lines), 4 * " ")

        # And the calculated atomization enthalpy
        table = {
            "System": [],
            "Term": [],
            "Value": [],
            "Units": [],
        }

        Hatoms = 0.0
        dH = Q_(6.197, "kJ/mol").m_as("E_h")
        for atno, symbol, count in composition:
            Eatom = atom_energies[atno - 1]
            # 6.197 is the H298-H0 for an atom
            Hatoms += count * (Eatom + 6.197)
            if correction_energy is not None and not isnan(correction_energy[atno - 1]):
                Hatoms += count * correction_energy[atno - 1]

            tmp = Q_(Eatom, "kJ/mol").m_as("E_h")
            table["System"].append(f"{symbol}(g)")
            if count == 1:
                table["Term"].append(f"{-tmp:.6f} + {dH:.6f}")
            else:
                table["Term"].append(f"{count} * ({-tmp:.6f} + {dH:.6f})")
            table["Value"].append(f"{-count * (tmp + dH):.6f}")
            table["Units"].append("")

        table["System"].append("^")
        table["Term"].append("-")
        table["Value"].append("-")
        table["Units"].append("")

        H = data["H_tot"]

        table["System"].append(formula)
        table["Term"].append(f"{H:.6f}")
        table["Value"].append("")
        table["Units"].append("E_h")

        data["H atomization"] = Hatoms - Q_(H, "E_h").m_as("kJ/mol")
        data["DfH0"] = DfH_at - data["H atomization"]
        table["System"].append("")
        table["Term"].append("")
        table["Value"].append("=")
        table["Units"].append("")

        result = f'{Q_(data["H atomization"], "kJ/mol").m_as("E_h"):.6f}'
        table["System"].append(f"{lDelta}atH{degree_sign}")
        table["Term"].append("")
        table["Value"].append(result)
        table["Units"].append("E_h")

        table["System"].append("")
        table["Term"].append("")
        table["Value"].append(f'{data["H atomization"]:.2f}')
        table["Units"].append("kJ/mol")

        tmp = tabulate(
            table,
            headers="keys",
            tablefmt="rounded_outline",
            colalign=("center", "center", "decimal", "center"),
            disable_numparse=True,
        )
        length = len(tmp.splitlines()[0])
        text_lines = []
        text_lines.append("Atomization Enthalpy (calculated)".center(length))
        text_lines.append(tmp)
        text += "\n\n"

        text += textwrap.fill(
            "Next we calculate the atomization enthalpy of the system. We have the "
            "calculated enthalpy of the system, but need the enthalpy of gas phase "
            f"atoms at the standard state (25{degree_sign}C, 1 atm). The tabulated "
            "energies for the atoms, used above, are identical to H0 for an atom. "
            "We will add H298 - H0 to each atom, which [1] is 5/2RT = 0.002360 E_h"
        )
        text += "\n\n"
        text += textwrap.indent("\n".join(text_lines), 4 * " ")
        text += "\n\n"
        text += textwrap.fill(
            "The enthalpy change for reaction (3) is the negative of this atomization"
            " enthalpy. Putting the two reactions together with the negative for Rxn 3:"
        )
        text += "\n\n"
        text += f"{lDelta}fH{degree_sign} = {lDelta}H(rxn 2) - {lDelta}H(rxn 3)\n"
        text += f"     = {DfH_at:.2f} - {data['H atomization']:.2f}\n"
        text += f"     = {DfH_at - data['H atomization']:.2f} kJ/mol\n"

        text += "\n\n"
        text += "References\n"
        text += "----------\n"
        text += "1. https://en.wikipedia.org/wiki/Monatomic_gas\n"
        refno = 1
        for atno, symbol, count in composition:
            refno += 1
            text += f"{refno}. {lDelta}fH{degree_sign} = {DfH0gas[atno - 1]} kJ/mol"
            if term_symbols is not None:
                text += f" for {term_symbols[atno - 1]} {symbol}"
            else:
                text += f" for {symbol}"
            if references is not None:
                text += f" from {references[atno-1]}\n"

        return text

    def description_text(
        self,
        P=None,
        calculation_type="Single-point energy",
        configuration=None,
    ):
        """Prepare information about what this node will do"""

        if P is None:
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
        method, functional, extended_functional, method_string = self.get_method()

        if self.parent.basis is not None:
            basis_set = self.parent.basis
            if method == "dft":
                self.parent.model = f"{functional.upper()}/{basis_set}"
                self.parent.extended_model = (
                    f"{extended_functional.upper()}/{basis_set}"
                )
            else:
                self.parent.model = f"{method.upper()}/{basis_set}"

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
            if restart is None:
                lines.append(
                    f"Eelec, wfn = {calculation_type}('{extended_functional}', "
                    "return_wfn=True)"
                )
                if calculation_type == "energy":
                    lines.append(f"G = gradient('{extended_functional}', ref_wfn=wfn)")
            else:
                if calculation_type == "gradient":
                    lines.append(
                        f"Eelec, wfn = energy('{extended_functional}', return_wfn=True,"
                        f" restart_file='{restart}')"
                    )
                    lines.append(f"G = gradient('{extended_functional}', ref_wfn=wfn)")
                else:
                    lines.append(
                        f"Eelec, wfn = {calculation_type}('{extended_functional}', "
                        f"return_wfn=True, restart_file='{restart}')"
                    )
        else:
            if restart is None:
                lines.append(
                    f"Eelec, wfn = {calculation_type}('{method}', return_wfn=True)"
                )
                if calculation_type == "energy":
                    lines.append(f"G = gradient('{method}', ref_wfn=wfn)")
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
variables["energy"] = Eelec
try:
    variables["gradients"] = np.array(G).tolist()
except Exception as e:
    print("Problem with gradients!")
    print(e)

variables["_method"] = "{method}"
variables["_method_string"] = "{method_string}"

tmp = fix_multipoles(variables)
with open("{filename}", "w") as fd:
    json.dump(tmp, fd, sort_keys=True, indent=3)
"""
            )

        # Orbital plots
        lines.append(self.plot_input())

        return "\n".join(lines)

    def get_method(self):
        """Get the method and functional to use"""
        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )

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
                extended_functional = functional + "-" + P["dispersion"]
            else:
                extended_functional = functional
        else:
            functional = method
            extended_functional = functional
        return method, functional, extended_functional, method_string

    def analyze(self, indent="", data=None, out=[], table=None):
        """Parse the output and generating the text output and store the
        data in variables for other stages to access
        """
        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )
        system, configuration = self.get_system_configuration()
        method, functional, extended_functional, method_string = self.get_method()

        directory = Path(self.directory)

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

        # Put any requested results into variables or tables
        self.store_results(
            data=data,
            create_tables=self.parameters["create tables"].get(),
        )

        text = ""

        # Calculate the enthalpy of formation, if possible
        tmp_text = self.calculate_enthalpy_of_formation(data)
        if tmp_text != "":
            path = Path(self.directory) / "Thermochemistry.txt"
            path.write_text(tmp_text)

        if table is None:
            table = {
                "Property": [],
                "Value": [],
                "Units": [],
            }

        # Special handling for DfH0 if it exists
        if "DfH0" in data:
            tmp = data["DfH0"]
            table["Property"].append(
                "\N{GREEK CAPITAL LETTER DELTA}fH\N{SUPERSCRIPT ZERO}"
            )
            table["Value"].append(f"{Q_(tmp, 'kJ/mol').m_as('kcal/mol'):.2f}")
            table["Units"].append("kcal/mol")
            table["Property"].append("")
            table["Value"].append(f"{tmp:.2f}")
            table["Units"].append("kJ/mol")

        if "ZPE_corr" in data:
            tmp = data["ZPE_corr"]
            table["Property"].append("ZPE")
            table["Value"].append(f"{Q_(tmp, 'kJ/mol').m_as('kcal/mol'):.2f}")
            table["Units"].append("kcal/mol")
            table["Property"].append("")
            table["Value"].append(f"{tmp:.2f}")
            table["Units"].append("kJ/mol")

        keys = [
            "H atomization",
            "DfE0",
            "E atomization",
            "energy",
        ]
        metadata = psi4_step.metadata["results"]
        for key in keys:
            if key in data:
                tmp = data[key]
                mdata = metadata[key]
                table["Property"].append(key)
                if "format" in mdata:
                    table["Value"].append(f"{tmp:{mdata['format']}}")
                else:
                    table["Value"].append(f"{tmp}")
                if "units" in mdata:
                    table["Units"].append(mdata["units"])
                else:
                    table["Units"].append("")

        keys = (
            ("metadata/symmetry_detected", "Symmetry"),
            ("metadata/symmetry_used", "Symmetry used"),
            ("E(gap)", ""),
            ("E(lumo+1)", "E(LUMO+1)"),
            ("E(lumo)", "E(LUMO)"),
            ("E(homo)", "E(HOMO)"),
            ("E(homo-1)", "E(HOMO-1)"),
            ("dipole_moment_magnitude", "Dipole moment"),
        )
        for key, name in keys:
            if name == "":
                name = key
            if key in data:
                tmp = data[key]
                if key == "state":
                    tmp = superscript[tmp[0]] + tmp[1:]
                mdata = metadata[key]
                table["Property"].append(name)
                table["Value"].append(f"{tmp:{mdata['format']}}")
                if "units" in mdata:
                    table["Units"].append(mdata["units"])
                else:
                    table["Units"].append("")

        tmp = tabulate(
            table,
            headers="keys",
            tablefmt="rounded_outline",
            colalign=("center", "decimal", "left"),
            disable_numparse=True,
        )
        length = len(tmp.splitlines()[0])
        text_lines = []

        multiplicity = configuration.spin_multiplicity
        spin_restricted = P["spin-restricted"]
        spin_text = ""
        if spin_restricted == "default":
            if multiplicity == 1:
                spin_text = "R-"
            else:
                spin_text = "U-"
        elif spin_restricted == "yes":
            if multiplicity == 1:
                spin_text = "R-"
            else:
                spin_text = "RO-"
        else:
            spin_text = "U-"

        text_lines.append(
            f"Results for {spin_text}{self.parent.extended_model}".center(length)
        )
        text_lines.append(method_string.center(length))
        text_lines.append(tmp)

        if text != "":
            text = str(__(text, **data, indent=self.indent + 4 * " "))
            text += "\n\n"
        text += textwrap.indent("\n".join(text_lines), self.indent + 7 * " ")
        printer.normal(text)

        # Write the structure locally for use in density and orbital plots
        obConversion = openbabel.OBConversion()
        obConversion.SetOutFormat("sdf")
        obMol = configuration.to_OBMol(properties="*")
        title = f"SEAMM={system.name}/{configuration.name}"
        obMol.SetTitle(title)
        sdf = obConversion.WriteString(obMol)
        path = directory / "structure.sdf"
        path.write_text(sdf)

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

    def process_data(self, data):
        """Massage the cclib data to a more easily used form."""
        logger.debug(pprint.pformat(data))
        # Convert numpy arrays to Python lists
        new = {}
        for key, value in data.items():
            if isinstance(value, np.ndarray):
                new[key] = value.tolist()
            elif isinstance(value, list):
                if len(value) > 0 and isinstance(value[0], np.ndarray):
                    new[key] = [i.tolist() for i in value]
                else:
                    new[key] = value
            elif isinstance(value, dict):
                for k, v in value.items():
                    newkey = f"{key}/{k}"
                    if isinstance(v, np.ndarray):
                        new[newkey] = v.tolist()
                    else:
                        new[newkey] = v
            else:
                new[key] = value

        for key in ("metadata/cpu_time", "metadata/wall_time"):
            if key in new:
                time = new[key][0]
                for tmp in new[key][1:]:
                    time += tmp
                new[key] = str(time).lstrip("0:")
                if "." in new[key]:
                    new[key] = new[key].rstrip("0")

        # Pull out the HOMO and LUMO energies as scalars
        if "homos" in new and "moenergies" in new:
            homos = new["homos"]
            if len(homos) == 2:
                for i, letter in enumerate(["α", "β"]):
                    Es = new["moenergies"][i]
                    homo = homos[i]
                    new[f"N({letter}-homo)"] = homo + 1
                    new[f"E({letter}-homo)"] = Es[homo]
                    if homo > 0:
                        new[f"E({letter}-homo-1)"] = Es[homo - 1]
                    if homo + 1 < len(Es):
                        new[f"E({letter}-lumo)"] = Es[homo + 1]
                        new[f"E({letter}-gap)"] = Es[homo + 1] - Es[homo]
                    if homo + 2 < len(Es):
                        new[f"E({letter}-lumo+1)"] = Es[homo + 2]
                    if "mosyms" in new:
                        syms = new["mosyms"][i]
                        new[f"Sym({letter}-homo)"] = syms[homo]
                        if homo > 0:
                            new[f"Sym({letter}-homo-1)"] = syms[homo - 1]
                        if homo + 1 < len(syms):
                            new[f"Sym({letter}-lumo)"] = syms[homo + 1]
                        if homo + 2 < len(syms):
                            new[f"Sym({letter}-lumo+1)"] = syms[homo + 2]
            else:
                Es = new["moenergies"][0]
                homo = homos[0]
                new["N(homo)"] = homo + 1
                new["E(homo)"] = Es[homo]
                if homo > 0:
                    new["E(homo-1)"] = Es[homo - 1]
                if homo + 1 < len(Es):
                    new["E(lumo)"] = Es[homo + 1]
                    new["E(gap)"] = Es[homo + 1] - Es[homo]
                if homo + 2 < len(Es):
                    new["E(lumo+1)"] = Es[homo + 2]
                if "mosyms" in new:
                    syms = new["mosyms"][0]
                    new["Sym(homo)"] = syms[homo]
                    if homo > 0:
                        new["Sym(homo-1)"] = syms[homo - 1]
                    if homo + 1 < len(syms):
                        new["Sym(lumo)"] = syms[homo + 1]
                    if homo + 2 < len(syms):
                        new["Sym(lumo+1)"] = syms[homo + 2]

        # moments
        if "moments" in new:
            moments = new["moments"]
            new["multipole_reference"] = moments[0]
            new["dipole_moment"] = moments[1]
            new["dipole_moment_magnitude"] = np.linalg.norm(moments[1])
            if len(moments) > 2:
                new["quadrupole_moment"] = moments[2]
            if len(moments) > 3:
                new["octapole_moment"] = moments[3]
            if len(moments) > 4:
                new["hexadecapole_moment"] = moments[4]
            del new["moments"]

        for key in ("metadata/symmetry_detected", "metadata/symmetry_used"):
            if key in new:
                new[key] = new[key].capitalize()

        return new
