# -*- coding: utf-8 -*-
"""The timing records Psi4 runs write (seamm_exec campaign 2026-10-05)."""

from types import SimpleNamespace

from psi4_step.psi4 import timing_descriptors

OUT = """\
  Number of basis functions: 25
    Nalpha       = 5
    Nbeta        = 5
   @DF-RKS iter   1:   -76.30000000000000   -7.63000e+01   1.00000e-02 DIIS
   @DF-RKS iter   2:   -76.40000000000000   -1.00000e-01   1.00000e-03 DIIS
  Energy and wave function converged.
    Psi4 wall time for execution: 0:00:03.12

*** Psi4 exiting successfully. Buy a developer a beer!
"""


def test_descriptors():
    conf = SimpleNamespace(
        atoms=SimpleNamespace(atomic_numbers=[8, 1, 1]),
        charge=0,
        spin_multiplicity=1,
        periodicity=0,
    )
    control = [
        [["basis", "6-31G**"]],
        [["method", "Kohn-Sham (KS) DFT"], ["functional", "B3LYP Hyb-GGA"]],
    ]
    d = timing_descriptors(control, OUT, conf)
    assert d["n_calculations"] == 2
    assert d["basis"] == "6-31G**" and d["method"] == "Kohn-Sham"
    assert d["nbf"] == 25 and d["n_electrons"] == 10
    assert d["scf_runs"] == 1 and d["scf_iterations"] == 2
    assert abs(d["code_seconds"] - 3.12) < 1e-9
    assert d["terminated_normally"] is True
