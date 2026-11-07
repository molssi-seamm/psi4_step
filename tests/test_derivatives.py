#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""How Psi4 gets gradients: analytic, or by finite differences of energies.

Without an analytic gradient Psi4 quietly differentiates energies, so the step
must decide, set ``qc_module`` where Psi4 needs it, and say so."""

from types import SimpleNamespace

import pytest

import psi4_step
from psi4_step.derivatives import (
    gradient_route,
    integral_options,
    integral_type,
    model_prefix,
)
from seamm.builder import set_parameters

CCSD_T = "CCSD with perturbative triples (CCSD(T))"


@pytest.mark.parametrize(
    "method, reference, frozen, analytic, options",
    [
        ("ccsd(t)", "rhf", False, True, ["set qc_module ccenergy"]),
        ("ccsd(t)", "uhf", False, True, ["set qc_module ccenergy"]),
        ("ccsd(t)", "rhf", True, False, []),
        ("ccsd(t)", "rohf", False, False, []),
        ("ccsd", "rhf", False, True, []),
        ("ccsd", "rohf", False, True, []),
        ("ccsd", "rhf", True, False, []),
        ("lccd", "rhf", True, False, []),
        ("ccd", "rhf", False, False, []),
        ("mp2", "rhf", True, True, []),
        ("mp2", "rohf", False, False, []),
        ("mp3", "uhf", True, True, []),
        ("hf", "rhf", False, True, []),
        ("dft", "rks", False, True, []),
        ("omp2", "rohf", True, True, []),
    ],
)
def test_gradient_route(method, reference, frozen, analytic, options):
    route = gradient_route(method, reference, frozen)
    assert route.analytic is analytic
    assert route.options == options
    assert bool(route.reason) is not analytic


def test_frozen_core_reason_says_what_to_do():
    assert (
        "'Freeze core orbitals' to 'no'"
        in gradient_route("ccsd(t)", "rhf", True).reason
    )


def water():
    return SimpleNamespace(spin_multiplicity=1)


@pytest.mark.parametrize("frozen, analytic", [("no", True), ("yes", False)])
def test_ccsd_t_route_from_the_step(frozen, analytic, monkeypatch):
    import seamm

    monkeypatch.setattr(seamm, "flowchart_variables", seamm.Variables())
    node = psi4_step.Energy(flowchart=seamm.Flowchart())
    set_parameters(node, method=CCSD_T, **{"freeze-cores": frozen})
    route = node.gradient_route(water())
    assert route.analytic is analytic
    assert route.options == (["set qc_module ccenergy"] if analytic else [])


@pytest.mark.parametrize(
    "method, reference, frozen, integrals, analytic",
    [
        # Density fitting: DFOCC (or DFMP2), frozen core allowed, RHF and UHF
        ("ccsd(t)", "rhf", True, "DF", True),
        ("ccsd(t)", "uhf", True, "DF", True),
        ("ccsd(t)", "rohf", False, "DF", False),
        ("ccsd", "rhf", True, "DF", True),
        ("ccd", "rhf", True, "DF", True),
        ("lccd", "uhf", True, "DF", True),
        ("mp2", "rhf", True, "DF", True),
        # Conventional MP2/MP3/LCCD gradients are all-electron (occ)
        ("mp2", "rhf", True, "CONV", False),
        ("mp2", "rhf", False, "CONV", True),
        ("mp3", "uhf", True, "CONV", False),
        ("ccd", "rhf", False, "CONV", False),
        # Cholesky: never analytic for correlated methods
        ("ccsd(t)", "rhf", False, "CD", False),
        ("mp2", "rhf", False, "CD", False),
        ("omp2", "rhf", False, "CD", False),
        ("omp2", "rhf", True, "DF", True),
        # SCF methods ignore the correlation integrals
        ("hf", "rhf", False, "CD", True),
    ],
)
def test_gradient_route_integrals(method, reference, frozen, integrals, analytic):
    route = gradient_route(method, reference, frozen, integrals)
    assert route.analytic is analytic
    assert bool(route.reason) is not analytic
    # The DF route needs no qc_module for the gradient itself
    if integrals != "CONV":
        assert route.options == []


def test_conventional_frozen_core_reason_mentions_density_fitting():
    assert "density-fitted" in gradient_route("ccsd(t)", "rhf", True).reason
    assert "density-fitted" not in gradient_route("eom-ccsd", "rhf", True).reason


def test_cholesky_reason():
    assert "Cholesky" in gradient_route("ccsd(t)", "rhf", False, "CD").reason


@pytest.mark.parametrize(
    "method, choice, expected",
    [
        ("ccsd(t)", "default", "CONV"),
        ("ccsd(t)", "density fitting", "DF"),
        ("mp2", "default", "DF"),
        ("mp3", "default", "DF"),
        ("omp3", "default", "CONV"),
        ("mp2", "conventional", "CONV"),
        ("lccd", "Cholesky decomposition", "CD"),
        ("hf", "density fitting", None),
        ("dft", "default", None),
    ],
)
def test_integral_type(method, choice, expected):
    assert integral_type(method, choice) == expected


def test_model_prefix():
    assert model_prefix("ccsd(t)", "density fitting") == "DF-"
    assert model_prefix("ccsd(t)", "default") == ""
    assert model_prefix("ccsd(t)", "conventional") == ""
    assert model_prefix("mp2", "default") == ""
    assert model_prefix("mp2", "density fitting") == ""
    assert model_prefix("mp2", "conventional") == "CONV-"
    assert model_prefix("mp2", "Cholesky decomposition") == "CD-"
    assert model_prefix("hf", "density fitting") == ""


def test_integral_options():
    assert integral_options("ccsd(t)", "rhf") == []
    assert integral_options("ccsd(t)", "rhf", "density fitting") == [
        "set cc_type df",
        "set qc_module occ",
    ]
    assert integral_options("ccsd(t)", "uhf", "density fitting", "cc-pVTZ-RI") == [
        "set cc_type df",
        "set qc_module occ",
        "set df_basis_cc cc-pVTZ-RI",
    ]
    assert integral_options("mp2", "rhf", "conventional") == ["set mp2_type conv"]
    assert integral_options("lccd", "rhf", "Cholesky decomposition") == [
        "set cc_type cd"
    ]
    assert integral_options("dft", "rks", "density fitting") == []


def _energy_node(monkeypatch):
    import seamm

    monkeypatch.setattr(seamm, "flowchart_variables", seamm.Variables())
    return psi4_step.Energy(flowchart=seamm.Flowchart())


def test_df_ccsd_t_frozen_core_route_from_the_step(monkeypatch):
    node = _energy_node(monkeypatch)
    set_parameters(
        node,
        method=CCSD_T,
        **{"freeze-cores": "yes", "correlation integrals": "density fitting"},
    )
    route = node.gradient_route(water())
    assert route.analytic is True
    assert route.options == []


def test_correlation_integral_rules():
    P = psi4_step.EnergyParameters()
    values = P.current_values()
    # DFT: no correlation integrals
    assert not P.applies("correlation integrals", values)
    values["method"] = CCSD_T
    assert P.applies("correlation integrals", values)
    assert not P.applies("correlation fitting basis", values)
    values["correlation integrals"] = "density fitting"
    assert P.applies("correlation fitting basis", values)
    assert P.problems(values) == []
    # CCD needs DF or CD: the conventional/default choice is replaced
    values = P.current_values()
    values["level"] = "advanced"
    values["advanced_method"] = "coupled cluster doubles (CCD)"
    assert P.applies("correlation integrals", values)
    assert P.choices("correlation integrals", values) == (
        "density fitting",
        "Cholesky decomposition",
    )
    assert P.implied(values)["correlation integrals"] == "density fitting"
