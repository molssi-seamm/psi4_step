# -*- coding: utf-8 -*-
"""The CC amplitude convergence (r_convergence) and the derivatives choice."""

from types import SimpleNamespace

import pytest

import psi4_step

CCSD_T = "CCSD with perturbative triples (CCSD(T))"
DFT = "Kohn-Sham (KS) density functional theory (DFT)"


@pytest.fixture
def variables(monkeypatch):
    import seamm

    monkeypatch.setattr(seamm, "flowchart_variables", seamm.Variables())


def _input(cls, **values):
    """The Psi4 input lines of interest for a step with these values."""
    import seamm

    flowchart = seamm.Flowchart()
    flowchart.root_directory = "/tmp"
    node = cls(flowchart=flowchart)
    node._id = ("2", "2")
    node.parent = SimpleNamespace(basis="cc-pCVTZ", model=None, extended_model=None)
    node.get_system_configuration = lambda *a, **k: (
        None,
        SimpleNamespace(spin_multiplicity=1, n_atoms=5),
    )
    for key, value in values.items():
        node.parameters[key].value = value
    node.description = []
    text = node.get_input()
    return [
        line
        for line in text.splitlines()
        if not line.startswith("#")
        and (
            line.startswith("set r_convergence")
            or "optimize(" in line
            or "frequency(" in line
            or "_driver.gradient = _clean_gradient" in line
        )
    ]


def test_amplitude_convergence_applies_to_correlated_methods():
    P = psi4_step.EnergyParameters()
    values = P.current_values()
    assert not P.applies("CC amplitude convergence", values)
    assert "amplitude" in P.not_applicable_reason("CC amplitude convergence", values)
    values["method"] = CCSD_T
    assert P.applies("CC amplitude convergence", values)


def test_derivatives_choice_on_optimization_and_thermochemistry():
    for cls in (psi4_step.OptimizationParameters, psi4_step.ThermochemistryParameters):
        P = cls()
        assert P["derivatives"].value == "analytic where available"
        assert "finite differences of energies" in P["derivatives"].enumeration
    assert "derivatives" not in psi4_step.EnergyParameters()


def test_amplitude_convergence_value(variables):
    node = psi4_step.Energy()
    node.parameters["method"].value = CCSD_T
    P = node.parameters.current_values_to_dict()
    assert node.cc_amplitude_convergence(P) is None  # Psi4's own
    node.parameters["CC amplitude convergence"].value = "1e-10"
    P = node.parameters.current_values_to_dict()
    assert float(node.cc_amplitude_convergence(P)) == 1e-10
    node.parameters["method"].value = DFT
    P = node.parameters.current_values_to_dict()
    assert node.cc_amplitude_convergence(P) is None


def test_optimization_input(variables):
    lines = _input(psi4_step.Optimization, method=CCSD_T, **{"freeze-cores": "no"})
    assert lines == [
        "    _driver.gradient = _clean_gradient",
        "Eelec, wfn = optimize('ccsd(t)', return_wfn=True)",
    ]
    lines = _input(
        psi4_step.Optimization,
        method=CCSD_T,
        **{
            "freeze-cores": "no",
            "derivatives": "finite differences of energies",
            "CC amplitude convergence": "1e-10",
        },
    )
    assert lines == [
        "set r_convergence 1e-10",
        "    _driver.gradient = _clean_gradient",
        "Eelec, wfn = optimize('ccsd(t)', return_wfn=True, dertype='energy')",
    ]
    # DFT: neither
    lines = _input(psi4_step.Optimization, method=DFT)
    assert lines == ["Eelec, wfn = optimize('b3lyp-d3bj', return_wfn=True)"]


def test_thermochemistry_input_with_its_own_parameters(variables):
    lines = _input(
        psi4_step.Thermochemistry,
        method=CCSD_T,
        **{
            "use existing parameters": "no",
            "freeze-cores": "no",
            "derivatives": "finite differences of energies",
        },
    )
    assert lines == [
        "Eelec, wfn = frequency('ccsd(t)', return_wfn=True, dertype='energy')",
    ]


def test_thermochemistry_own_parameters_keep_the_correlation_integrals(variables):
    """Bugfix: the rules saw 'use existing parameters' as False, not 'no', so a
    Thermochemistry step with its own parameters ignored the integrals choice."""
    node = psi4_step.Thermochemistry()
    node.parameters["method"].value = CCSD_T
    node.parameters["use existing parameters"].value = "no"
    node.parameters["correlation integrals"].value = "density fitting"
    P = node.parameters.current_values_to_dict()
    assert node.integral_choice(P) == "density fitting"


def test_thermochemistry_describes_the_previous_steps_method(variables):
    """Bugfix: with 'use existing parameters', the description named this step's own
    unused settings (the default B3LYP), not the method that runs (Paul 11-08)."""
    import seamm

    flowchart = seamm.Flowchart()
    optimization = psi4_step.Optimization(flowchart=flowchart)
    optimization.parameters["method"].value = CCSD_T
    optimization.parameters["correlation integrals"].value = "density fitting"
    thermo = psi4_step.Thermochemistry(flowchart=flowchart)
    thermo._id = ("2", "3")
    thermo.previous = lambda: optimization
    text = thermo.description_text()
    assert "previous step" in text
    assert "CCSD(T)" in text and "density fitting integrals, frozen core" in text
    assert "B3LYP" not in text and "Kohn-Sham" not in text

    thermo.parameters["use existing parameters"].value = "no"
    text = thermo.description_text()
    assert "previous step" not in text


def test_clean_scratch_only_for_correlated_optimizations(variables):
    """Psi4's optimize() keeps the previous geometry's scratch files, which made a
    conventional CCSD(T) gradient wrong (5652): correlated optimizations clean them
    before each gradient; SCF/DFT optimizations and Hessians do not need it."""
    from psi4_step.energy import CLEAN_GRADIENTS

    lines = _input(psi4_step.Optimization, method=CCSD_T, **{"freeze-cores": "yes"})
    assert "    _driver.gradient = _clean_gradient" in lines
    lines = _input(psi4_step.Optimization, method=DFT)
    assert "    _driver.gradient = _clean_gradient" not in lines
    lines = _input(
        psi4_step.Thermochemistry,
        method=CCSD_T,
        **{"use existing parameters": "no", "freeze-cores": "no"},
    )
    assert "    _driver.gradient = _clean_gradient" not in lines
    compile(CLEAN_GRADIENTS, "CLEAN_GRADIENTS", "exec")


def test_clean_gradients_fix_in_psi4():
    """With Psi4 installed: in optimize(), every step's AE CCSD(T) gradient equals a
    fresh one (without the fix they differed by 1.7e-5 and 3.2e-5 Eh/bohr)."""
    psi4 = pytest.importorskip("psi4")
    import sys

    import numpy as np

    from psi4_step.energy import CLEAN_GRADIENTS

    psi4.set_memory("1 GB")
    psi4.core.be_quiet()
    exec(CLEAN_GRADIENTS, {"core": psi4.core})
    driver = sys.modules["psi4.driver.driver"]
    patched = driver.gradient
    history = []

    def record(*args, **kwargs):
        G, wfn = patched(*args, **kwargs)
        history.append((wfn.molecule().geometry().np.copy(), np.array(G)))
        return G, wfn

    driver.gradient = record
    try:
        molecule = psi4.geometry(
            "0 1\nO 0 0 0\nH 0 0.80 0.60\nH 0 -0.80 0.60\nno_com\nno_reorient"
        )
        psi4.set_options(
            {
                "basis": "cc-pvdz",
                "freeze_core": False,
                "qc_module": "ccenergy",
                "d_convergence": 1e-10,
                "e_convergence": 1e-10,
                "geom_maxiter": 3,
                "g_convergence": "gau_verytight",
            }
        )
        try:
            psi4.optimize("ccsd(t)", molecule=molecule)
        except psi4.OptimizationConvergenceError:
            pass
    finally:
        driver.gradient = patched
    assert len(history) >= 2
    for geometry, gradient in history:
        psi4.core.clean()
        fresh = psi4.core.Molecule.from_arrays(
            geom=geometry,
            elez=[8, 1, 1],
            units="Bohr",
            fix_com=True,
            fix_orientation=True,
        )
        expected = np.array(psi4.gradient("ccsd(t)", molecule=fresh))
        assert np.abs(gradient - expected).max() < 1e-8
