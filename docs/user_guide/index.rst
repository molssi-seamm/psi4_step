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

Index
=====

* :ref:`genindex`
