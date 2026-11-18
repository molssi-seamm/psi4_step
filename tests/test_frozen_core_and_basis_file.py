# -*- coding: utf-8 -*-
"""The frozen-core choices (Psi4's freeze_core policies and an explicit count) and
basis files, e.g. ORCA's cc-pwCVTZ for Ge-Br with explicit fitting sets."""

import re
from types import SimpleNamespace

import pytest

import psi4_step
from psi4_step.derivatives import freeze_core_options, gradient_route
from psi4_step.initialization import (
    BASIS_DIRECTORY,
    basis_file_name,
    basis_file_path,
    shipped_basis_files,
)

CCSD_T = "CCSD with perturbative triples (CCSD(T))"
INNER = "inner core ([Ne] on Ga-Kr, [He] on Al-Ar)"
INNERMOST = "innermost core ([He] on Ga-Kr)"


@pytest.fixture
def variables(monkeypatch):
    import seamm

    monkeypatch.setattr(seamm, "flowchart_variables", seamm.Variables())


# -- frozen core ---------------------------------------------------------------


@pytest.mark.parametrize(
    "choice, lines, frozen",
    [
        ("yes", ["set freeze_core True"], True),
        ("no", ["set freeze_core False"], False),
        (INNER, ["set freeze_core -2"], True),
        (INNERMOST, ["set freeze_core -3"], True),
    ],
)
def test_frozen_core_choices(choice, lines, frozen):
    assert freeze_core_options({"freeze-cores": choice}) == (lines, frozen)


def test_frozen_core_count():
    P = {"freeze-cores": "number of orbitals", "frozen orbitals": 6}
    assert freeze_core_options(P) == (["set num_frozen_docc 6"], True)
    P["frozen orbitals"] = 0
    assert freeze_core_options(P) == (["set num_frozen_docc 0"], False)
    P["frozen orbitals"] = -1
    with pytest.raises(ValueError):
        freeze_core_options(P)


def test_unknown_frozen_core_choice():
    with pytest.raises(ValueError, match="Unknown frozen-core choice"):
        freeze_core_options({"freeze-cores": "-1"})


def test_frozen_orbitals_applies_only_to_a_count():
    P = psi4_step.EnergyParameters()
    values = P.current_values()
    values["method"] = CCSD_T
    assert P.applies("freeze-cores", values)
    assert not P.applies("frozen orbitals", values)
    values["freeze-cores"] = "number of orbitals"
    assert P.applies("frozen orbitals", values)
    values["freeze-cores"] = INNER
    assert not P.applies("frozen orbitals", values)


def test_frozen_core_dfocc_gradient_is_analytic():
    """Level (b): DF-CCSD(T) with [Ne] frozen on Ge has an analytic gradient
    (DFOCC); conventional CCSD(T) with a frozen core does not."""
    _, frozen = freeze_core_options({"freeze-cores": INNER})
    assert gradient_route("ccsd(t)", "rhf", frozen, "DF").analytic
    assert not gradient_route("ccsd(t)", "rhf", frozen, "CONV").analytic


def _energy_input(**values):
    import seamm

    flowchart = seamm.Flowchart()
    flowchart.root_directory = "/tmp"
    node = psi4_step.Energy(flowchart=flowchart)
    node._id = ("2", "2")
    node.parent = SimpleNamespace(
        basis="cc-pwCVTZ-ORCA", model=None, extended_model=None
    )
    node.get_system_configuration = lambda *a, **k: (
        None,
        SimpleNamespace(spin_multiplicity=1, n_atoms=5),
    )
    for key, value in values.items():
        node.parameters[key].value = value
    node.description = []
    return node, node.get_input()


def test_energy_input_inner_core(variables):
    node, text = _energy_input(
        method=CCSD_T,
        **{"freeze-cores": INNER, "correlation integrals": "density fitting"},
    )
    assert "set freeze_core -2" in text.splitlines()
    assert "set qc_module occ" in text.splitlines()
    assert "number of frozen core orbitals" in text  # reported from the wfn
    assert node.parent.model == "DF-CCSD(T)/cc-pwCVTZ-ORCA"
    description = "".join("".join(str(d) for d in node.description).split())
    assert "innercore([Ne]onGa-Kr,[He]onAl-Ar)willbefrozen" in description


def test_energy_input_frozen_count(variables):
    _, text = _energy_input(
        method=CCSD_T,
        **{"freeze-cores": "number of orbitals", "frozen orbitals": 6},
    )
    assert "set num_frozen_docc 6" in text.splitlines()
    assert not any(line.startswith("set freeze_core") for line in text.splitlines())


def test_hf_freezes_nothing(variables):
    _, text = _energy_input(method="Hartree-Fock (HF) self consistent field (SCF)")
    assert "set freeze_core False" in text.splitlines()


# -- basis files ---------------------------------------------------------------


def test_shipped_basis_file():
    assert "cc-pwCVTZ-ORCA" in shipped_basis_files()
    path = basis_file_path("cc-pwcvtz-orca")  # any case
    assert path == BASIS_DIRECTORY / "cc-pwCVTZ-ORCA.basis"
    assert basis_file_name("cc-pwCVTZ-ORCA.basis") == "cc-pwCVTZ-ORCA"
    with pytest.raises(FileNotFoundError, match="cc-pwCVTZ-ORCA"):
        basis_file_path("no-such-basis")


def test_basis_file_by_path(tmp_path):
    path = tmp_path / "my-basis.basis"
    path.write_text("basis {\n   assign cc-pvdz\n}\n")
    assert basis_file_path(str(path)) == path
    assert basis_file_name(str(path)) == "my-basis"


def _sections(text):
    """{(block, section): {element: [angular momenta]}} of a basis file."""
    result = {}
    block = section = None
    shells = None
    for line in text.splitlines():
        line = line.strip()
        m = re.match(r"^(\w+) \{$", line)
        if m:
            block = m.group(1)
            continue
        m = re.match(r"^\[ (\S+) \]$", line)
        if m:
            section = m.group(1)
            result[(block, section)] = {}
            continue
        m = re.match(r"^([A-Z]{1,2})\s+0$", line)
        if m and section is not None:
            shells = result[(block, section)].setdefault(m.group(1), [])
            continue
        m = re.match(r"^([SPDFGHI])\s+(\d+)\s+1\.00$", line)
        if m and shells is not None:
            shells.append("SPDFGHI".index(m.group(1)))
    return result


def _nbf(shells):
    """Spherical basis functions of a list of shells' angular momenta."""
    return sum(2 * am + 1 for am in shells)


def test_cc_pwcvtz_orca_function_counts():
    """The counts ORCA gives for cc-pwCVTZ: GeH4 140, GeO 127 (spherical)."""
    text = basis_file_path("cc-pwCVTZ-ORCA").read_text()
    sections = _sections(text)
    orbital = sections[("basis", "orca-cc-pwcvtz")]
    assert "spherical" in text.split("[ orca-cc-pwcvtz ]")[1].splitlines()[1]
    assert _nbf(orbital["GE"]) == 84  # [8s7p5d3f1g]
    assert _nbf(orbital["GE"]) + 4 * _nbf(orbital["H"]) == 140
    assert _nbf(orbital["GE"]) + _nbf(orbital["O"]) == 127
    for element in ("AS", "SE", "BR"):
        assert _nbf(orbital[element]) == 84
    assert _nbf(orbital["CL"]) == 59  # [7s6p4d2f]
    assert _nbf(orbital["C"]) == 43  # [6s5p3d1f]


def test_cc_pwcvtz_orca_fitting_sets():
    text = basis_file_path("cc-pwCVTZ-ORCA").read_text()
    sections = _sections(text)
    ri = sections[("df_basis_cc", "bse-cc-pwcvtz-rifit")]
    assert {"GE", "AS", "SE", "BR", "CL", "C"} <= set(ri)
    assert set(sections[("df_basis_cc", "bse-cc-pvtz-rifit")]) == {"H", "HE"}
    assert "   assign H bse-cc-pvtz-rifit" in text.splitlines()
    orbital = sections[("basis", "orca-cc-pwcvtz")]
    jk = sections[("df_basis_scf", "bse-def2-universal-jkfit")]
    assert set(jk) == set(orbital)


def test_initialization_input_with_basis_file(variables):
    import seamm

    flowchart = seamm.Flowchart()
    flowchart.root_directory = "/tmp"
    node = psi4_step.Initialization(flowchart=flowchart)
    node._id = ("2", "1")
    node.parent = SimpleNamespace(basis=None)
    node.parameters["basis source"].value = "basis file"
    node.parameters["basis file"].value = "cc-pwCVTZ-ORCA"
    text = node.get_input()
    assert node.parent.basis == "cc-pwCVTZ-ORCA"
    assert "set basis" not in text
    for block in ("basis {", "df_basis_cc {", "df_basis_scf {"):
        assert block in text.splitlines()

    node.parameters["basis source"].value = "library"
    node.parameters["basis"].value = "cc-pVTZ"
    text = node.get_input()
    assert node.parent.basis == "cc-pVTZ"
    assert "set basis cc-pVTZ" in text.splitlines()


def test_basis_applies_by_source():
    P = psi4_step.InitializationParameters()
    values = P.current_values()
    assert P.applies("basis", values) and not P.applies("basis file", values)
    values["basis source"] = "basis file"
    assert P.applies("basis file", values) and not P.applies("basis", values)


# -- with Psi4 itself (skipped where it is not installed) ----------------------

MOLECULES = {
    # name: (geometry, frozen orbitals with [Ne] on Ge/As, [He] on Cl)
    "GeH4": (
        "Ge 0 0 0\nH 0.88 0.88 0.88\nH -0.88 -0.88 0.88\nH -0.88 0.88 -0.88\n"
        "H 0.88 -0.88 -0.88",
        5,
    ),
    "GeH3Cl": (
        "Ge 0 0 0\nCl 0 0 2.15\nH 1.45 0 -0.45\nH -0.72 1.25 -0.45\n"
        "H -0.72 -1.25 -0.45",
        6,
    ),
    "ClAsH2": ("As 0 0 0\nCl 0 0 2.2\nH 1.5 0 -0.3\nH -0.4 1.4 -0.3", 6),
    "CH3Cl": (
        "C 0 0 0\nCl 0 0 1.78\nH 1.03 0 -0.36\nH -0.51 0.89 -0.36\n"
        "H -0.51 -0.89 -0.36",
        1,
    ),
}


@pytest.mark.parametrize("name", MOLECULES)
def test_psi4_inner_core_counts(name):
    """Psi4 freezes [Ne] on Ge/As and [He] on Cl for the 'inner core' choice."""
    psi4 = pytest.importorskip("psi4")
    geometry, expected = MOLECULES[name]
    psi4.core.be_quiet()
    molecule = psi4.geometry(geometry + "\nsymmetry c1\nno_com\nno_reorient")
    basis = psi4.core.BasisSet.build(molecule, "BASIS", "cc-pVTZ")
    (line,), _ = freeze_core_options({"freeze-cores": INNER})
    policy = line.split()[-1]
    assert basis.n_frozen_core(policy, molecule) == expected
    assert basis.n_frozen_core("True", molecule) > expected  # [Ar] on Ge/As
