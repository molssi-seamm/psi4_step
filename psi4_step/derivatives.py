# -*- coding: utf-8 -*-

"""Which Psi4 calculations have analytic gradients.

Psi4 (1.11) picks the code for a gradient from the method, the reference, the
algorithm (``scf_type``, ``mp2_type``, ``mp_type``, ``cc_type``; this step leaves
them at their defaults: density fitting for SCF, MP2, MP2.5 and MP3, conventional
for the coupled-cluster methods and LCCD) and ``qc_module``. When no code fits, Psi4
does not stop: it quietly differentiates energies instead, which for a correlated
method costs many times more -- every energy of a Hessian, roughly the square of
the number of coordinates. So the step decides here, writes ``qc_module`` where
Psi4 needs it, asks for energies explicitly where there is no analytic gradient,
and says so in the output.

The rules follow ``psi4/driver/procrouting/proc.py`` (the ``select_*_gradient``
functions):

* conventional CCSD(T): ``ccenergy``, RHF or UHF, only with ``qc_module ccenergy``
  set, and all-electron (``ccenergy`` refuses frozen core);
* conventional CCSD and EOM-CCSD: ``ccenergy``, all-electron;
* conventional LCCD: ``occ``, RHF or UHF, all-electron;
* conventional CCD: none;
* DF-MP2, DF-MP2.5, DF-MP3: RHF or UHF, frozen core or not;
* the orbital-optimized methods, HF and DFT: analytic.
"""

from dataclasses import dataclass, field

import psi4_step


@dataclass
class GradientRoute:
    """How Psi4 gets the gradient: analytically or from energies."""

    analytic: bool
    options: list = field(default_factory=list)  # Psi4 'set' lines it needs
    reason: str = ""  # why there is no analytic gradient


ALWAYS = {"hf", "scf", "dft", "efp", "omp2", "omp3", "omp2.5", "olccd", "oremp2"}
DF_RHF_UHF = {"mp2", "mp2.5", "mp3"}  # density fitted by default
CONV_ALL_ELECTRON = {
    # method: references with an analytic all-electron gradient
    "ccsd": {"rhf", "uhf", "rohf"},
    "eom-ccsd": {"rhf", "uhf", "rohf"},
    "lccd": {"rhf", "uhf"},
    "ccsd(t)": {"rhf", "uhf"},
}
NEEDS_CCENERGY = {"ccsd(t)"}  # Psi4 routes it to ccenergy only if asked
NONE = {"ccd": "Psi4 has an analytic CCD gradient only with density fitting"}


def gradient_route(method, reference, frozen_core):
    """The route Psi4 takes to the gradient of ``method``.

    Parameters
    ----------
    method : str
        Psi4's name for the method, e.g. "ccsd(t)" or "dft".
    reference : str
        "rhf", "uhf", "rohf", "rks" or "uks".
    frozen_core : bool
        Whether the core orbitals are frozen.

    Returns
    -------
    GradientRoute
    """
    method = method.lower()
    reference = reference.lower()
    if method in ALWAYS or reference in ("rks", "uks"):
        return GradientRoute(True)
    if method in DF_RHF_UHF:
        if reference in ("rhf", "uhf"):
            return GradientRoute(True)
        return GradientRoute(
            False, reason=f"Psi4 has no analytic {reference.upper()} gradient for it"
        )
    if method in NONE:
        return GradientRoute(False, reason=NONE[method])
    if method in CONV_ALL_ELECTRON:
        if reference not in CONV_ALL_ELECTRON[method]:
            return GradientRoute(
                False,
                reason=f"Psi4 has no analytic {reference.upper()} gradient for it",
            )
        if frozen_core:
            return GradientRoute(
                False,
                reason=(
                    "Psi4's analytic gradient for it is all-electron only; set "
                    "'Freeze core orbitals' to 'no' for analytic gradients"
                ),
            )
        options = ["set qc_module ccenergy"] if method in NEEDS_CCENERGY else []
        return GradientRoute(True, options=options)
    # Anything else: what the metadata says.
    for record in psi4_step.methods.values():
        if record["method"] == method:
            if record.get("gradients") == "analytic":
                return GradientRoute(True)
            break
    return GradientRoute(False, reason="Psi4 has no analytic gradient for it")
