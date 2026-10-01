#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The rules the dialogs and the flowchart builder share (see seamm.Parameters):
which settings apply, the dispersion corrections a functional offers, and the
values implied by others."""

import pytest

import psi4_step
from seamm.builder import FlowchartBuildError, set_parameters

DFT = "Kohn-Sham (KS) density functional theory (DFT)"
HF = "Hartree-Fock (HF) self consistent field (SCF)"
MP2 = "2nd-order Møller–Plesset perturbation theory (MP2)"
B3LYP = "B3LYP Hyb-GGA Exchange-Correlation Functional"


def functional_with(dispersions):
    return next(
        name
        for name, record in psi4_step.dft_functionals.items()
        if record["dispersion"] == dispersions
    )


def test_method_follows_the_level():
    node = psi4_step.Energy()
    with pytest.raises(FlowchartBuildError, match="it applies when 'level' is"):
        set_parameters(node, advanced_method=MP2)
    set_parameters(node, level="advanced", advanced_method=MP2)
    assert node.parameters["advanced_method"].value == MP2


def test_functional_only_for_dft():
    node = psi4_step.Energy()
    with pytest.raises(FlowchartBuildError) as e:
        set_parameters(node, method=MP2, functional=B3LYP)
    assert "'functional' has no effect" in str(e.value)
    assert f"it applies when 'method' is '{DFT}'" in str(e.value)
    # The step is left as it was
    assert node.parameters["method"].value == DFT


def test_freeze_cores_only_for_correlated_methods():
    node = psi4_step.Energy()
    with pytest.raises(FlowchartBuildError, match="does not freeze core orbitals"):
        set_parameters(node, method=HF, freeze_cores="no")
    set_parameters(node, method=MP2, freeze_cores="no")
    assert node.parameters["freeze-cores"].value == "no"
    # A variable could be any method
    P = node.parameters
    assert P.applies("freeze-cores", {**P.current_values(), "method": "$method"})


def test_dispersion_narrowed_and_implied():
    node = psi4_step.Energy()
    P = node.parameters
    f2 = functional_with(["none", "d3bj"])
    values = {**P.current_values(), "level": "advanced", "advanced_functional": f2}
    assert P.choices("dispersion", values) == ("none", "d3bj")

    # The current correction is kept when the functional has it ...
    set_parameters(node, level="advanced", advanced_functional=f2)
    assert P["dispersion"].value == "d3bj"
    # ... replaced when it does not ...
    P["dispersion"].value = "d3mbj"
    set_parameters(node, advanced_functional=functional_with(["none", "d3bj", "nl"]))
    assert P["dispersion"].value == "d3bj"
    # ... and refused when given with it.
    with pytest.raises(FlowchartBuildError, match="it must be 'd3bj'"):
        set_parameters(node, advanced_functional=f2, dispersion="nl")


def test_dispersion_needs_a_functional_with_corrections():
    node = psi4_step.Energy()
    with pytest.raises(FlowchartBuildError, match="has no dispersion corrections"):
        set_parameters(
            node,
            level="advanced",
            advanced_functional=functional_with(["none"]),
            dispersion="none",
        )
    with pytest.raises(FlowchartBuildError, match="it needs 'functional'"):
        set_parameters(node, method=HF, dispersion="none")
    # A variable for the functional could have any
    set_parameters(node, functional="$functional", dispersion="nl")
    assert node.parameters["dispersion"].value == "nl"


def test_convergence_sub_controls():
    node = psi4_step.Energy()
    with pytest.raises(FlowchartBuildError, match="applies when 'use damping' is"):
        set_parameters(node, damping_percentage=30)
    set_parameters(node, use_damping=True, damping_percentage=30)
    for switch, key in (
        ("use level shift", "level shift"),
        ("use soscf", "soscf convergence"),
        ("orbitals", "selected orbitals"),
    ):
        P = node.parameters
        assert not P.applies(key)
        assert P.applies(key, {**P.current_values(), switch: "yes"})


def test_thermochemistry_uses_the_previous_parameters():
    node = psi4_step.Thermochemistry()
    with pytest.raises(FlowchartBuildError, match="previous step are used"):
        set_parameters(node, method=HF)
    set_parameters(node, use_existing_parameters="no", method=HF)
    assert node.parameters["method"].value == HF
    # Its own rules still hold
    with pytest.raises(FlowchartBuildError, match="'method' is 'Kohn-Sham"):
        set_parameters(node, functional=B3LYP)


def test_optimization_has_no_subsequent_structures():
    node = psi4_step.Optimization()
    with pytest.raises(FlowchartBuildError, match="does not use it"):
        set_parameters(node, subsequent_structure_handling="Discard the structure")


def test_bsse_fragment_atoms():
    node = psi4_step.BSSE()
    with pytest.raises(FlowchartBuildError, match="'fragments' is 'specified'"):
        set_parameters(node, fragment_atoms="1-3; 4-6")
    set_parameters(node, fragments="specified", fragment_atoms="1-3; 4-6")
    assert node.parameters["fragment atoms"].value == "1-3; 4-6"


def test_builder():
    """Build a Psi4 flowchart and write it, with the rules checked."""
    from seamm.builder import FlowchartBuilder

    fb = FlowchartBuilder()
    psi4 = fb.add("Psi4")
    psi4.add("Initialization")
    with pytest.raises(FlowchartBuildError, match="'functional' has no effect"):
        psi4.add("Energy", method=MP2, functional=B3LYP)
    energy = psi4.add(
        "Energy", method=DFT, functional="PBE0 Hyb-GGA Exchange-Correlation Functional"
    )
    assert energy.parameters["dispersion"].value == "d3bj"
    assert "Psi4" in fb.to_text()


@pytest.mark.parametrize(
    "functional",
    ["B1LYP Hyb-GGA Exchange-Correlation Functional", "b1lyp", "B1LYP"],
)
def test_get_method_accepts_the_short_name(functional, monkeypatch):
    """A functional given by its short name (e.g. from a variable) raised KeyError."""
    import seamm
    import psi4_step

    monkeypatch.setattr(seamm, "flowchart_variables", seamm.Variables())
    node = psi4_step.Energy(flowchart=seamm.Flowchart())
    P = node.parameters
    P["method"].value = "Kohn-Sham (KS) density functional theory (DFT)"
    P["functional"].value = functional
    P["dispersion"].value = "d3bj"
    method, name, extended, _ = node.get_method()
    assert (method, name, extended) == ("dft", "b1lyp", "b1lyp-d3bj")


def test_bsse_refuses_the_settings_it_ignores():
    """The counterpoise calculation uses Psi4's own SCF settings and makes no plots,
    so the Energy settings for those are refused, saying why."""
    import seamm

    flowchart = seamm.Flowchart(namespace="org.molssi.seamm.psi4", directory=".")
    node = flowchart.create_node("BSSE")
    flowchart.add_node(node)
    for key, value in (("use damping", "yes"), ("orbitals", "yes")):
        with pytest.raises(FlowchartBuildError, match="does not use it"):
            set_parameters(node, {key: value})
    set_parameters(node, {"maximum iterations": 200, "fragments": "specified"})


def test_thermochemistry_plots_only_with_its_own_settings():
    P = psi4_step.ThermochemistryParameters()
    values = {**P.current_values(), "use existing parameters": "yes"}
    assert not P.applies("orbitals", values)
    assert "previous step" in P.not_applicable_reason("density", values)
    values["use existing parameters"] = "no"
    assert P.applies("orbitals", values)
