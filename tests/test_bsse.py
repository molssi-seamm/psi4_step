#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the Psi4 BSSE (counterpoise) sub-step.

Pure-Python logic only (fragment parsing, the molecule-block builder) -- no
Psi4 execution. See docs/developer_guide/campaigns/2026-08-04/ for the real
Psi4, end-to-end validation script.
"""

from types import SimpleNamespace

import pytest

import psi4_step


def test_bsse_factory():
    assert psi4_step.BSSEStep().description()["name"] == "BSSE"


def test_bsse_extends_energy():
    assert issubclass(psi4_step.BSSE, psi4_step.Energy)
    node = psi4_step.BSSE()
    assert node._calculation == "bsse"
    P = psi4_step.BSSEParameters()
    assert P["fragments"].value == "auto (molecules)"
    # Inherits the energy parameters.
    assert P["level"].value == "recommended"


def test_bsse_parse_indices():
    assert psi4_step.BSSE._parse_indices("1-3, 5 7", 8) == [0, 1, 2, 4, 6]
    assert psi4_step.BSSE._parse_indices("2 2 3", 5) == [1, 2]  # de-duplicated
    with pytest.raises(RuntimeError):
        psi4_step.BSSE._parse_indices("9", 8)


def test_bsse_parse_charges():
    assert psi4_step.BSSE._parse_charges("") == []
    assert psi4_step.BSSE._parse_charges("  ") == []
    assert psi4_step.BSSE._parse_charges("1, -1") == [1, -1]
    assert psi4_step.BSSE._parse_charges("0 0 0") == [0, 0, 0]


def test_bsse_parse_fragment_groups():
    groups = psi4_step.BSSE._parse_fragment_groups("1-3; 4-6", 6)
    assert groups == [[0, 1, 2], [3, 4, 5]]
    with pytest.raises(RuntimeError, match="at least two"):
        psi4_step.BSSE._parse_fragment_groups("1-4", 4)


def test_bsse_fragments_specified_charged_na_cl():
    """The Na+/Cl- pilot case: 'specified' fragments, per-fragment charge,
    neutral overall complex."""
    node = psi4_step.BSSE()
    configuration = SimpleNamespace(
        n_atoms=2,
        charge=0,
        spin_multiplicity=1,
        atoms=SimpleNamespace(atomic_numbers=[11, 17]),
    )
    P = {
        "fragments": "specified",
        "fragment atoms": "1; 2",
        "fragment charges": "1, -1",
    }
    fragments = node._fragments(P, configuration)
    assert [f.atom_indices for f in fragments] == [(0,), (1,)]
    assert [f.charge for f in fragments] == [1, -1]


def test_bsse_fragments_charge_mismatch_raises():
    node = psi4_step.BSSE()
    configuration = SimpleNamespace(
        n_atoms=2,
        charge=0,
        spin_multiplicity=1,
        atoms=SimpleNamespace(atomic_numbers=[11, 17]),
    )
    P = {
        "fragments": "specified",
        "fragment atoms": "1; 2",
        "fragment charges": "1, 1",  # should be 1, -1 for a neutral complex
    }
    with pytest.raises(RuntimeError, match="sum to"):
        node._fragments(P, configuration)


def test_bsse_fragments_auto_needs_at_least_two_molecules():
    node = psi4_step.BSSE()
    configuration = SimpleNamespace(
        n_atoms=3,
        charge=0,
        spin_multiplicity=1,
        find_molecules=lambda as_indices=True: [[0, 1, 2]],
    )
    P = {"fragments": "auto (molecules)", "fragment atoms": "", "fragment charges": ""}
    with pytest.raises(RuntimeError, match="at least two"):
        node._fragments(P, configuration)


class _FakeAtoms:
    """A minimal stand-in for `configuration.atoms` -- symbols/get_coordinates
    for `_molecule_block`, `atomic_numbers` for the BSSE electron-parity
    check, and `formal_charge` (`in`/`[]`, like the real molsystem atoms
    table) for the 'fragment charges' structure-derived default."""

    def __init__(
        self, symbols=None, coords=None, atomic_numbers=None, formal_charge=None
    ):
        self.symbols = symbols
        self._coords = coords
        self.atomic_numbers = atomic_numbers
        self._formal_charge = formal_charge

    def get_coordinates(self, fractionals=False, in_cell=True):
        return self._coords

    def __contains__(self, key):
        return key == "formal_charge" and self._formal_charge is not None

    def __getitem__(self, key):
        if key == "formal_charge" and self._formal_charge is not None:
            return self._formal_charge
        raise KeyError(key)


def test_bsse_fragments_default_charges_from_formal_charge():
    """No 'fragment charges' given, but the structure (e.g. read from an SDF
    with an 'M  CHG' record) carries per-atom formal charges -- the Na+..H2O
    pilot case: water (O, H, H) neutral, Na+ +1. Each fragment's default
    charge must be the sum of its own atoms' formal_charge, not 0."""
    node = psi4_step.BSSE()
    node._id = ("1",)
    configuration = SimpleNamespace(
        n_atoms=4,
        charge=1,
        spin_multiplicity=1,
        find_molecules=lambda as_indices=True: [[0, 1, 2], [3]],
        atoms=_FakeAtoms(atomic_numbers=[8, 1, 1, 11], formal_charge=[0, 0, 0, 1]),
    )
    P = {"fragments": "auto (molecules)", "fragment atoms": "", "fragment charges": ""}
    fragments = node._fragments(P, configuration)
    assert [f.atom_indices for f in fragments] == [(0, 1, 2), (3,)]
    assert [f.charge for f in fragments] == [0, 1]


def test_bsse_fragments_explicit_charges_override_formal_charge():
    """An explicit 'fragment charges' still wins over the structure's own
    (here, deliberately wrong/absent-looking) formal charges."""
    node = psi4_step.BSSE()
    configuration = SimpleNamespace(
        n_atoms=2,
        charge=0,
        spin_multiplicity=1,
        atoms=_FakeAtoms(atomic_numbers=[11, 17], formal_charge=[0, 0]),
    )
    P = {
        "fragments": "specified",
        "fragment atoms": "1; 2",
        "fragment charges": "1, -1",
    }
    fragments = node._fragments(P, configuration)
    assert [f.charge for f in fragments] == [1, -1]


def test_bsse_molecule_block_two_fragments():
    """The '--'-separated molecule block: per-fragment charge/multiplicity
    header, no leading '--' before the first fragment, no_com/no_reorient,
    and an explicit activate() so this molecule (not whatever the shared
    script's own "initial" molecule is) is the active one."""
    node = psi4_step.BSSE()
    configuration = SimpleNamespace(
        atoms=_FakeAtoms(
            ["Na", "Cl"],
            [(0.0, 0.0, 0.0), (0.0, 0.0, 2.44)],
            atomic_numbers=[11, 17],
        )
    )
    fragments = node._fragments(
        {
            "fragments": "specified",
            "fragment atoms": "1; 2",
            "fragment charges": "1, -1",
        },
        SimpleNamespace(
            n_atoms=2,
            charge=0,
            spin_multiplicity=1,
            atoms=SimpleNamespace(atomic_numbers=[11, 17]),
        ),
    )
    block = node._molecule_block(fragments, configuration, name="bsse_cluster")
    lines = block.splitlines()

    assert lines[0] == "molecule bsse_cluster {"
    assert lines[1] == "1 1"
    assert lines[2].startswith("Na")
    assert lines[3] == "--"
    assert lines[4] == "-1 1"
    assert lines[5].startswith("Cl")
    assert "no_com" in lines
    assert "no_reorient" in lines
    assert lines[-1] == "activate(bsse_cluster)"
    assert lines[-2] == "}"
    # Exactly one '--' -- between the two fragments, not before the first.
    assert block.count("--") == 1
