#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""How Psi4 gets gradients: analytic, or by finite differences of energies.

Without an analytic gradient Psi4 quietly differentiates energies, so the step
must decide, set ``qc_module`` where Psi4 needs it, and say so."""

from types import SimpleNamespace

import pytest

import psi4_step
from psi4_step.derivatives import gradient_route
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
