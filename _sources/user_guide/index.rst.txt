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


Gradients and Hessians
======================

Optimizations and Thermochemistry need gradients. Psi4 has analytic gradients for most
methods, but where it has none it differentiates energies instead, which for a
correlated method costs many more calculations: a Hessian needs roughly the square of
the number of coordinates. The step decides which applies, asks Psi4 for it
explicitly, and says so in the output when energies are used:

* **CCSD(T), CCSD, EOM-CCSD and LCCD** (conventional, the default algorithm) have
  analytic gradients only for **all-electron** calculations, RHF or UHF (CCSD and
  EOM-CCSD also ROHF). Set "Freeze core orbitals" to "no" for analytic gradients; with
  the core frozen, gradients and Hessians come from energies. For CCSD(T) the step
  also tells Psi4 to use its ``ccenergy`` module, without which Psi4 used energies
  even for all-electron calculations (before version 2026.11.7 of this step).
* **CCD** has no analytic conventional gradient.
* **MP2, MP2.5 and MP3** (density fitted, the default) have analytic gradients for RHF
  and UHF, with or without a frozen core; the orbital-optimized methods, HF and DFT
  have them too.

Psi4 has analytic Hessians only for HF (and LDA functionals), so Thermochemistry finds
the Hessian by finite differences of gradients, using the molecule's symmetry to reduce
the number of displaced geometries.

Index
=====

* :ref:`genindex`
