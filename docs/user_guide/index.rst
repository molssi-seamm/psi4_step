.. _user-guide:

**********
User Guide
**********
The Psi4 plug-in provides access to many different quantum chemical calculations usign
the Psi4 code.

..
   The following sections cover accessing and controlling this functionality.

   .. toctree::
      :maxdepth: 2
      :titlesonly:

Settings that depend on each other
==================================

Which settings apply depends on others: the functional only for DFT, the frozen core
only for correlated methods, a dispersion correction only from the functional's own
list, and the damping, level-shift and second-order SCF controls only when switched on.
The step's dialogs show only the settings that apply with the current choices, and the
same rules are used when a flowchart is built or edited without the editor
(``seamm-flowchart`` or SEAMM's MCP server): a setting that would have no effect is
refused, with the reason, and a value that contradicts another is refused too. See
"Flowcharts without the editor" in SEAMM's user guide.

BSSE uses Psi4's own SCF settings and makes no plots, so it does not offer those Energy
settings. Thermochemistry's plots apply only when it uses its own settings rather than
the previous step's; placed right after Initialization, where there is no previous
calculation, it always uses its own.


Correlation integrals
=====================

For MP2, MP2.5, MP3, LCCD, CCD, CCSD, CCSD(T) and the orbital-optimized methods,
"Correlation integrals" chooses how the integrals of the correlated method are done:

* **default**: Psi4's default for the method -- density fitted for MP2, MP2.5 and
  MP3, conventional for LCCD and the coupled-cluster methods.
* **conventional**: exact four-index integrals.
* **density fitting**: much faster and lighter on memory, and the only way to get
  analytic coupled-cluster gradients with a frozen core (Psi4's DFOCC module). The
  energies differ slightly from conventional ones -- of the order of tens of
  microhartree per atom in total energies with the matching fitting basis, largely
  cancelling in relative energies -- so it is a different model: its name starts with ``DF-`` (e.g. ``DF-CCSD(T)/cc-pVTZ``) and its results are
  stored under that name, apart from conventional ones. "Fitting basis" is the
  auxiliary basis (Psi4's ``df_basis_cc``); ``default`` uses the RI basis matching
  the orbital basis, e.g. cc-pVTZ-RI. The SCF is density fitted too, as it is by
  default.
* **Cholesky decomposition**: like density fitting with an automatically built
  auxiliary basis, but with no analytic gradients.

A model whose integrals are not the method's default gets a prefix: ``DF-``,
``CD-``, or ``CONV-`` (e.g. conventional MP2 is ``CONV-MP2``). CCD exists in Psi4 only
with density fitting or Cholesky decomposition, so the dialog offers only those.

For DF-CCSD and DF-CCSD(T) the step tells Psi4 to use DFOCC for the energy as well as
the gradient, so both come from the same code, and UHF DF-CCSD(T) works.

Gradients and Hessians
======================

Optimizations and Thermochemistry need gradients. Psi4 has analytic gradients for most
methods, but where it has none it differentiates energies instead, which for a
correlated method costs many more calculations: a Hessian needs roughly the square of
the number of coordinates. The step decides which applies, asks Psi4 for it
explicitly, and says so in the output when energies are used:

* **Conventional CCSD(T), CCSD, EOM-CCSD, LCCD, MP2, MP2.5 and MP3** have analytic
  gradients only for **all-electron** calculations, RHF or UHF (CCSD and EOM-CCSD
  also ROHF). Set "Freeze core orbitals" to "no" for analytic gradients, or use
  density fitting; with the core frozen and conventional integrals, gradients and
  Hessians come from energies. For CCSD(T) the step also tells Psi4 to use its
  ``ccenergy`` module, without which Psi4 used energies even for all-electron
  calculations (before version 2026.11.7 of this step).
* **Density-fitted MP2, MP2.5, MP3, LCCD, CCD, CCSD and CCSD(T)** have analytic
  gradients for RHF and UHF, with or without a frozen core. Frozen-core
  DF-CCSD(T) is the affordable route to coupled-cluster Hessians.
* **Cholesky-decomposed** integrals have no analytic gradients for any correlated
  method.
* The orbital-optimized methods (conventional or density fitted), HF and DFT have
  analytic gradients.

Psi4 has analytic Hessians only for HF (and LDA functionals), so Thermochemistry finds
the Hessian by finite differences of gradients, using the molecule's symmetry to reduce
the number of displaced geometries. With "use existing parameters", Thermochemistry
uses the previous step's settings, including its correlation integrals.

**Derivatives** (Optimization and Thermochemistry) chooses between "analytic where
available", the default described above, and "finite differences of energies", which
uses only energies even where Psi4 has an analytic gradient. That costs many more
calculations, but it is a check on an analytic gradient, or a way around one.

**Optimizations with correlated methods** clean Psi4's scratch files before each
gradient. Psi4 1.11's optimizer keeps the previous step's files, and a correlated
gradient computed with them present can be wrong: for the all-electron CCSD(T)
gradient of CH\ :sub:`4` near its minimum, ten times too large, so a tight
optimization never converged. The previous step's orbitals are still used as the
starting guess. Hessians are not affected: Psi4 computes each displaced gradient
afresh.

**CC amplitude convergence** (correlated methods) is how tightly Psi4 converges the
amplitude equations, and for gradients the lambda equations: its ``r_convergence``.
"default" leaves Psi4's own, 10\ :sup:`-7`.

Frozen core
===========

For the correlated methods, "Freeze core orbitals" chooses which core orbitals are left
out of the correlation:

* **yes**: the previous noble-gas core of each atom -- [He] on B-Ne, [Ne] on Al-Ar,
  [Ar] on Ga-Kr (Psi4's ``freeze_core true``).
* **no**: all electrons are correlated.
* **inner core ([Ne] on Ga-Kr, [He] on Al-Ar)**: the noble-gas core before that, so
  the outermost core shell is correlated with the valence -- 3s3p3d on Ga-Kr, 2s2p on
  Al-Ar -- and nothing is frozen on H-Ne (``freeze_core -2``). This goes with
  core-valence basis sets such as cc-pwCVTZ.
* **innermost core ([He] on Ga-Kr)**: one shell further: [He] on Ga-Kr, nothing on
  lighter atoms (``freeze_core -3``).
* **number of orbitals**: freezes the given number of lowest doubly occupied orbitals
  of the whole system (``num_frozen_docc``).

The number of frozen orbitals is printed in the results and stored as "number of
frozen core orbitals", so it can be checked: for GeH\ :sub:`4` with the inner core it
is 5 ([Ne] on Ge), and for GeH\ :sub:`3`\ Cl 6 (5 on Ge, 1 on Cl).

Basis files
===========

The Initialization sub-step takes the basis set either from Psi4's library, by name, or
from a **basis file**: Psi4 input blocks copied into the input -- a ``basis {...}``
block for the orbital basis and, optionally, ``df_basis_cc {...}`` and
``df_basis_scf {...}`` blocks giving the density-fitting sets. A basis file is either
one shipped with the step, by name, or the path of a file on the machine running the
flowchart. The name of the file, without its extension, is the basis set's name in the
model, e.g. ``DF-CCSD(T)/cc-pwCVTZ-ORCA``. A fitting basis set in the Energy step's
"Fitting basis" takes precedence over the file's.

The step ships one basis file:

* **cc-pwCVTZ-ORCA**: cc-pwCVTZ exactly as ORCA 6.1.1 has it (cc-pVTZ on H and He), for
  H, He, B-Ne, Al-Ar and Ge-Br, so that Psi4 calculations use the same basis as ORCA
  ones. Neither Psi4's library nor the Basis Set Exchange has an all-electron cc-pwCVTZ
  for Ge-Br. It gives cc-pwCVTZ-RIFIT (cc-pVTZ-RIFIT on H and He) for the correlation
  and def2-universal-JKFIT for the SCF, so DF-CCSD(T) fits the core-valence
  correlation of the outer core too. Its provenance is in ``data/basis/README.md`` in
  the package.

Index
=====

* :ref:`genindex`
