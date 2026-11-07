# -*- coding: utf-8 -*-

"""Which Psi4 calculations have analytic gradients, and how the integrals are done.

Psi4 (1.11) picks the code for a gradient from the method, the reference, the
algorithm of the correlation integrals (``mp2_type``, ``mp_type`` or ``cc_type``:
conventional, density fitted or Cholesky decomposed) and ``qc_module``. When no code
fits, Psi4 does not stop: it quietly differentiates energies instead, which for a
correlated method costs many times more -- every energy of a Hessian, roughly the
square of the number of coordinates. So the step decides here, writes the options
Psi4 needs, asks for energies explicitly where there is no analytic gradient, and
says so in the output.

The rules follow ``psi4/driver/procrouting/proc.py`` (the ``select_*`` functions):

* conventional CCSD(T): ``ccenergy``, RHF or UHF, only with ``qc_module ccenergy``
  set, and all-electron (``ccenergy`` refuses frozen core);
* conventional CCSD and EOM-CCSD: ``ccenergy``, all-electron;
* conventional MP2, MP2.5, MP3 and LCCD: ``occ``, RHF or UHF, all-electron;
* conventional CCD: no gradient (nor energy);
* density-fitted MP2, MP2.5, MP3, LCCD, CCD, CCSD and CCSD(T): RHF or UHF, frozen
  core or not (``dfmp2`` or ``dfocc``). The step sets ``qc_module occ`` for
  DF-CCSD and DF-CCSD(T), so that the energy comes from the same code as the
  gradient (Psi4 would otherwise take an RHF energy from ``fnocc``) and UHF
  DF-CCSD(T) works at all;
* Cholesky-decomposed integrals: no analytic gradient for any correlated method;
* the orbital-optimized methods: analytic, conventional or density fitted;
* HF and DFT: analytic.
"""

from dataclasses import dataclass, field

import psi4_step


@dataclass
class GradientRoute:
    """How Psi4 gets the gradient: analytically or from energies."""

    analytic: bool
    options: list = field(default_factory=list)  # Psi4 'set' lines it needs
    reason: str = ""  # why there is no analytic gradient


ALWAYS = {"hf", "scf", "dft", "efp"}
ORBITAL_OPTIMIZED = {"omp2", "omp3", "omp2.5", "olccd", "oremp2"}
CONV_ALL_ELECTRON = {
    # method: references with an analytic all-electron conventional gradient
    "mp2": {"rhf", "uhf"},
    "mp2.5": {"rhf", "uhf"},
    "mp3": {"rhf", "uhf"},
    "lccd": {"rhf", "uhf"},
    "ccsd": {"rhf", "uhf", "rohf"},
    "eom-ccsd": {"rhf", "uhf", "rohf"},
    "ccsd(t)": {"rhf", "uhf"},
}
DF_RHF_UHF = {"mp2", "mp2.5", "mp3", "lccd", "ccd", "ccsd", "ccsd(t)"}
NEEDS_CCENERGY = {"ccsd(t)"}  # Psi4 routes it to ccenergy only if asked
NEEDS_OCC = {"ccsd", "ccsd(t)"}  # DF/CD: dfocc for the energy too, and for UHF
NONE = {"ccd": "Psi4 has CCD only with density-fitted or Cholesky integrals"}

#: The Psi4 option that sets each method's correlation integrals, and its default.
#: (Psi4 treats MP2.5 and MP3 as density fitted unless mp_type is set.)
TYPE_KEYWORD = {
    "mp2": ("mp2_type", "DF"),
    "mp2.5": ("mp_type", "DF"),
    "mp3": ("mp_type", "DF"),
    "lccd": ("cc_type", "CONV"),
    "ccd": ("cc_type", "CONV"),
    "ccsd": ("cc_type", "CONV"),
    "ccsd(t)": ("cc_type", "CONV"),
    "omp2": ("mp2_type", "DF"),
    "omp2.5": ("mp_type", "CONV"),
    "omp3": ("mp_type", "CONV"),
    "olccd": ("cc_type", "CONV"),
}
#: The step's choices for the correlation integrals, as Psi4's values.
INTEGRAL_CHOICES = {
    "conventional": "CONV",
    "density fitting": "DF",
    "Cholesky decomposition": "CD",
}
#: How the model name shows integrals other than the method's default.
MODEL_PREFIX = {"CONV": "CONV-", "DF": "DF-", "CD": "CD-"}


def integral_type(method, choice="default"):
    """Psi4's correlation-integral type for ``method``: "CONV", "DF" or "CD".

    Parameters
    ----------
    method : str
        Psi4's name for the method.
    choice : str
        The step's setting: "default" (Psi4's default for the method),
        "conventional", "density fitting" or "Cholesky decomposition".

    Returns
    -------
    str or None
        None for methods without correlation integrals (HF, DFT, ...).
    """
    method = method.lower()
    if method not in TYPE_KEYWORD:
        return None
    if choice in INTEGRAL_CHOICES:
        return INTEGRAL_CHOICES[choice]
    return TYPE_KEYWORD[method][1]


def model_prefix(method, choice="default"):
    """The prefix for the model name ("DF-", "CD-", "CONV-") when the integrals are
    not the method's default, so results of different integrals are kept apart."""
    method = method.lower()
    integrals = integral_type(method, choice)
    if integrals is None or integrals == TYPE_KEYWORD[method][1]:
        return ""
    return MODEL_PREFIX[integrals]


def integral_options(method, reference, choice="default", fitting_basis="default"):
    """The Psi4 'set' lines for the correlation integrals of ``method``."""
    method = method.lower()
    integrals = integral_type(method, choice)
    if integrals is None:
        return []
    lines = []
    if choice in INTEGRAL_CHOICES:
        lines.append(f"set {TYPE_KEYWORD[method][0]} {integrals.lower()}")
    if integrals in ("DF", "CD") and method in NEEDS_OCC:
        lines.append("set qc_module occ")
    if integrals == "DF" and fitting_basis not in ("default", "", None):
        lines.append(f"set df_basis_cc {fitting_basis}")
    return lines


def gradient_route(method, reference, frozen_core, integrals=None):
    """The route Psi4 takes to the gradient of ``method``.

    Parameters
    ----------
    method : str
        Psi4's name for the method, e.g. "ccsd(t)" or "dft".
    reference : str
        "rhf", "uhf", "rohf", "rks" or "uks".
    frozen_core : bool
        Whether the core orbitals are frozen.
    integrals : str, optional
        "CONV", "DF" or "CD" (see :func:`integral_type`); by default Psi4's
        default for the method.

    Returns
    -------
    GradientRoute
    """
    method = method.lower()
    reference = reference.lower()
    if method in ALWAYS or reference in ("rks", "uks"):
        return GradientRoute(True)
    if integrals is None:
        integrals = integral_type(method)
    if integrals == "CD":
        return GradientRoute(
            False,
            reason=(
                "Psi4 has no analytic gradients with Cholesky-decomposed "
                "integrals; use density fitting for analytic gradients"
            ),
        )
    if method in ORBITAL_OPTIMIZED:
        return GradientRoute(True)
    if integrals == "DF" and method in DF_RHF_UHF:
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
            extra = (
                ", or density-fitted integrals, which allow a frozen core"
                if method in DF_RHF_UHF
                else ""
            )
            return GradientRoute(
                False,
                reason=(
                    "Psi4's analytic gradient for it with conventional integrals is "
                    "all-electron only; set 'Freeze core orbitals' to 'no' for "
                    f"analytic gradients{extra}"
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
