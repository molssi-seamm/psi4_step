2026-11-07 ALMO-EDA in Psi4 (design)
====================================

**Status:** design, not started. Option C of the ALMO-EDA work (Paul's route,
2026-11-06: skip a PySCF CPU port, build the Interaction EDA plug-in on GPU4PySCF,
then add Psi4 as its CPU backend). Estimate 3--5 weeks for the Python-only scope
below, plus 2--4 weeks if exact ωB97M-V parity is wanted, which needs C++ in Psi4.

Goal
----
Second-generation ALMO energy decomposition analysis (ALMO-EDA2) of intermolecular
interactions on CPUs, with an unmodified Psi4 1.11, giving the same terms as the
GPU4PySCF backend of the Interaction EDA step (``interaction_eda_step``):

* **frozen**, split into **electrostatic** (orthogonal decomposition),
  **classical electrostatic**, **Pauli** and **dispersion**;
* **polarization**, by fragment electric-field response functions (FERF);
* **charge transfer**, which carries the basis-set superposition error;
* the **total**, which must equal the supermolecular DFT interaction energy.

Where it runs: the Interaction EDA step gains a ``Psi4`` backend. Like the GPU4PySCF
one, it writes ``input.json`` and a driver script into the step's directory and runs
the driver with the Python of the Psi4 environment (``psi4.ini``). The driver uses
Psi4 as a library (``import psi4``); it is not a Psi4 input file. This page lives in
the Psi4 step's campaigns because the Psi4-side work is the substance; the step
itself (``psi4_step``) does not change.

The reference is the GPU4PySCF implementation the Interaction EDA step already uses,
validated to 1.4 × 10⁻⁶ kJ/mol against science's standalone runs (see
``interaction_eda_step``'s docs). Its definitions are the specification: where this
page and GPU4PySCF's ``eda.py`` differ, ``eda.py`` wins.

Python only, or C++?
--------------------
Psi4 has no ALMO-EDA, no SCF-MI and no FERF. But the expensive pieces are C++ objects
that Psi4 already exports to Python, and ALMO-EDA2 needs nothing else from C++, with
one exception: the exchange-correlation **response kernel**.

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Needed
     - Psi4 1.11, from Python (checked in the v1.11 source)
   * - Coulomb and exchange, incl. range-separated, for any density
     - ``core.JK.build(basis)``, ``jk.C_left_add/C_right_add``, ``set_do_wK/set_omega``, ``compute``, ``J()/K()/wK()``
   * - XC potential and energy, incl. VV10
     - ``core.VBase.build(basis, functional, "RV")``, ``set_D``, ``compute_V``, ``quadrature_values``
   * - XC response kernel (coupled response)
     - ``VBase.compute_Vx`` -- **LDA and GGA only**: ``RV::compute_Vx_full`` throws for VV10 and for meta-GGAs
   * - Functionals, VV10 on/off, ω
     - ``build_superfunctional``, ``set_do_vv10``, ``x_omega``, ``is_meta``, ``needs_vv10``
   * - Fragment SCF and its CPKS solver
     - ``energy(..., return_wfn=True)``; ``wfn.cphf_solve``, ``cphf_Hx``, ``onel_Hx/twoel_Hx`` (these use Vx, so the same limit applies)
   * - Overlap, dipole and quadrupole integrals
     - ``MintsHelper.ao_overlap``, ``ao_dipole``, ``ao_quadrupole``
   * - D4, if a -D4 functional is used
     - Psi4's empirical dispersion (dftd4)

So the scope splits by functional:

* **ωB97X-V** (a range-separated hybrid GGA plus VV10): VV10 in the SCF and SCF-MI,
  coupled FERF with the GGA kernel and **without the VV10 kernel**. The VV10 response
  term is negligible: dropping it moved energy between polarization and charge
  transfer by at most 0.0014 kJ/mol for a Li⁺ shell, and does not change the total
  (GPU4PySCF, mlff's xc_cmp-15003). **All Python.**
* **ωB97M-V with uncoupled FERF** (orbital-energy denominators, no kernel): **all
  Python.** Uncoupled FERF moves about 1 % of the polarization into charge transfer
  (shell A: +1.49 kJ/mol); frozen and total are unchanged.
* **ωB97M-V with coupled FERF**, exactly as GPU4PySCF: needs a **meta-GGA Vx kernel in
  Psi4's C++** (``RV::compute_Vx_full``, "does not support rotated V builds for
  MGGA's"). 2--4 weeks, best done with the Psi4 developers; who contacts them is
  Paul's decision.

The method
----------
The fragments' orbitals are expanded only in their own atoms' basis functions
(absolutely localized MOs, ALMOs), so the coefficient matrix is block-diagonal and the
MOs of different fragments are non-orthogonal. Every energy below is the full system's
DFT energy functional, evaluated at the density

.. math::

   P = C\,(C^{\mathsf T} S\, C)^{-1} C^{\mathsf T},

built from the block-diagonal :math:`C`; :math:`S` is the full AO overlap. That one
formula, with Psi4's ``JK`` and ``VBase`` for the two-electron and XC parts, is the
whole engine.

1. **Isolated fragments.** An SCF of each fragment in its own basis (monomer basis,
   not the dimer basis: the BSSE goes into charge transfer, as in Q-Chem and
   GPU4PySCF). Energies :math:`E_A`, occupied orbitals :math:`C_A`.
2. **Frozen state.** :math:`C_{\rm frz} = \mathrm{diag}(C_A, C_B, \dots)`; the
   frozen interaction is :math:`E[P_{\rm frz}] - \sum_A E_A`.
3. **Frozen decomposition** (Horn, Mao and Head-Gordon 2016):

   * *classical electrostatic*: the Coulomb interaction of the isolated fragments'
     charge densities (electrons and nuclei), from ``J`` of each fragment density;
   * *electrostatic* (orthogonal decomposition): from the fragment densities
     projected as in the paper;
   * *dispersion*: the frozen energy with the functional minus that with its
     dispersion-free counterpart, exactly as ``eda.py`` defines it for each
     functional (milestone 0 transcribes the definition);
   * *Pauli* = frozen − electrostatic − dispersion; the common (SAPT-comparable)
     Pauli = frozen − classical electrostatic − dispersion.

4. **FERF polarization space.** For each fragment, the response of its occupied
   orbitals to uniform fields (dipole) and field gradients (quadrupole), i.e. 3 + 6
   CPKS solves per fragment with the dipole and quadrupole integrals as
   perturbations:

   * coupled: ``wfn.cphf_solve`` on the fragment's wavefunction, whose Hessian uses
     Vx (GGA; VV10 switched off for the response only);
   * uncoupled: orbital-energy denominators, no kernel.

   The responses are orthogonalized per fragment and truncated by singular value
   (threshold 10⁻⁴, the validated default; the step warns above it for cations).
   Together with the occupied orbitals they span each fragment's polarization
   space.
5. **SCF-MI in the FERF space.** Minimize :math:`E[P(C)]` over block-diagonal
   :math:`C` whose fragment blocks rotate only within their polarization spaces.
   Gradient: the full Fock matrix :math:`F[P]` from ``JK`` + ``VBase``, projected
   onto each fragment's occupied-virtual rotations in the non-orthogonal metric
   (the locally projected equations of Stoll and of Khaliullin et al.). Solver:
   DIIS or a quasi-Newton step in the rotation parameters; the space is small (9 ×
   n\ :sub:`occ` per fragment), so this is cheap besides the Fock builds.
   Polarization = :math:`E_{\rm pol} - E_{\rm frz}`.
6. **Full SCF** of the whole system in the full basis (Psi4's own ``energy()``,
   started from the SCF-MI density). Charge transfer = :math:`E_{\rm full} -
   E_{\rm pol}`. Total = :math:`E_{\rm full} - \sum_A E_A`.

Each Fock build is one ``JK`` + ``VBase`` call on the whole system, so an SCF-MI
iteration costs about one SCF iteration. Expected cost: the fragment SCFs, one full
SCF and an SCF-MI of comparable length, plus 9 CPKS solves per fragment. Psi4 runs
the SCF about 1.7× faster than CPU PySCF (Owl, 32 cores: shell A 489 s vs 826 s).

Checks built into the driver
----------------------------
* frozen + polarization + charge transfer = total, to rounding;
* the total equals Psi4's supermolecular interaction energy (full SCF minus the
  fragment SCFs) -- an independent check of the energy functional at ``P``;
* frozen = electrostatic + Pauli + dispersion = classical electrostatic + common Pauli
  + dispersion;
* SCF-MI converged in its space, and its energy between frozen and full.

Validation against GPU4PySCF
----------------------------
The same structures and settings as the GPU4PySCF validation (``interaction_eda_step``
tests and science's references):

* S22 water dimer; Li⁺–EC at R = 1.20, 1.50, 1.70 Å; shell A (Li⁺, 2 EC, 2 DMC) pairs
  and cluster; Li⁺–PF₆⁻ and PF₆⁻–EC (anions);
* ωB97X-V/def2-TZVPD, coupled FERF without the VV10 kernel, against GPU4PySCF with
  ``nlc`` in the response off;
* ωB97M-V/def2-TZVPD with uncoupled FERF, against GPU4PySCF uncoupled;
* exact integrals in both (no density fitting), the same grids where the codes allow.

Gates: every term within 0.1 kJ/mol, aiming for 0.01; the total within 0.01 kJ/mol of
each code's own supermolecular interaction energy. Differences beyond that are
expected only from the DFT grids, and are measured first on the water dimer by
comparing plain DFT interaction energies.

Milestones
----------
.. list-table::
   :header-rows: 1
   :widths: 6 74 20

   * - M
     - Work
     - Estimate
   * - 0
     - Transcribe ``eda.py``'s definitions (dispersion-free functional, electrostatic projection, FERF space, SVD); reference numbers; a Psi4 environment for development
     - 2 days
   * - 1
     - Fragment SCFs, the :math:`P(C)` energy engine, the frozen energy; checked against Psi4's own SCF energy at a converged density
     - 1 week
   * - 2
     - Frozen decomposition: classical and orthogonal electrostatics, dispersion, Pauli
     - 3--4 days
   * - 3
     - FERF space, uncoupled
     - 3 days
   * - 4
     - SCF-MI in the FERF space; polarization and charge transfer; ωB97M-V uncoupled parity
     - 1 week
   * - 5
     - Coupled FERF (GGA kernel, no VV10 kernel); ωB97X-V parity
     - 4--5 days
   * - 6
     - The ``Psi4`` backend of the Interaction EDA step: driver, input, results, tests with a stand-in like the GPU4PySCF one, docs
     - 3--4 days
   * - 7
     - Validation runs on ARC CPUs, timings, release
     - 3 days
   * - (8)
     - Optional: meta-GGA Vx in Psi4's C++, for coupled ωB97M-V
     - 2--4 weeks

Risks
-----
* **The orthogonal electrostatic decomposition** has details (which projector, which
  densities) that must match ``eda.py`` exactly; milestone 0 exists for this.
* **SCF-MI convergence** for strongly bound ions (Li⁺) in the non-orthogonal metric;
  fall back to a trust-region step if DIIS struggles.
* **Grids:** Psi4 and GPU4PySCF build different DFT grids, so their terms can differ
  at the 0.01--0.1 kJ/mol level. If that exceeds the gate, tighten both grids rather
  than chase it in the EDA.
* **Python-side cost** for big clusters: the dense algebra in the metric is
  :math:`O(N^3)` with small prefactors; the Fock builds dominate. Measured in
  milestone 7.
* **Psi4 internals:** the driver uses exported but lightly documented objects
  (``VBase``, ``cphf_solve``); pin the Psi4 version the backend supports.

Results (2026-11-08/09)
-----------------------
The driver (``interaction_eda_step/data/psi4_eda.py``, released in
interaction-eda-step 2026.11.8.1) was run against the GPU4PySCF references of
TinkerCliffs job 793976 (ωB97X-V/def2-TZVPD, coupled FERF without VV10 in the response,
exact integrals) on TinkerCliffs AMD nodes, 16 threads unless noted. ωB97M-V is not
possible: Psi4 has no coupled meta-GGA kernel, and the backend refuses it.

**Parity with exact integrals** (job 8546834): HF on the water dimer and Li⁺–EC
agree to better than 1e-4 kJ/mol; ωB97X-V on the water dimer, Li⁺–EC and PF₆⁻–EC to
0.0036 kJ/mol at worst. The method is right.

**Two driver fixes** (jobs 8548610, 8548842):

* Keeping the integrals in core (PK) was 4.5 times *slower* than computing them
  directly for the EDA's many Fock builds (HF Li⁺–EC: 1179 s against 260 s, the same
  terms): the driver now uses DIRECT, or MEM_DF with an auxiliary basis.
* The density-fitted path stopped at once ("point group has not been set"):
  ``JK.build`` makes the fitting basis itself only for SCF_TYPE DF, so the driver now
  builds it explicitly.

**Density fitting** (def2-universal-JKFIT) on the whole of shell A (job 8552392; the
EC–EC pair 1+2 alone in job 8551357, with 8 threads in 8551537):

.. list-table:: Time and memory, 16 threads (GPU4PySCF on a V100 for comparison)
   :header-rows: 1
   :widths: 30 10 30 30

   * - Calculation
     - nao
     - Psi4 + DF
     - GPU4PySCF
   * - Li⁺–EC pairs 0+2, 0+3
     - 284
     - 424, 428 s
     - 87, 85 s
   * - Li⁺–EC pairs 0+1, 0+4
     - 302
     - 456, 478 s
     - 93, 94 s
   * - EC–EC pair 2+3
     - 534
     - 903 s
     - 166 s
   * - EC–EC pairs 1+2, 1+3, 2+4, 3+4
     - 552
     - 1120, 1048, 1134, 1011 s
     - 200, 191, 203, 180 s
   * - EC–EC pair 1+4
     - 570
     - 1191 s
     - 211 s
   * - The cluster
     - 1121
     - 7712 s, peak 141 GB
     - 1538 s
   * - All of shell A
     -
     - 4 h 26 min, 70.8 core-h
     - 3046 s

On the EC–EC pair, 8 threads took 2235 s (4.97 core-h) against 956 s (4.27 core-h) for
16, so 16 threads is the more efficient. The time grows as nao\ :sup:`2.8` from the
EC–EC pairs to the cluster.

.. list-table:: Density-fitting error against exact integrals (kJ/mol)
   :header-rows: 1
   :widths: 22 10 10 10 10 10 10 10 10

   * - Set
     - total
     - frozen
     - cls. elec.
     - orth. elec.
     - disp.
     - orth. Pauli
     - pol.
     - CT
   * - Pairs, worst
     - 0.011
     - 0.012
     - 0.009
     - 0.007
     - 0.044
     - 0.049
     - 0.007
     - 0.007
   * - Cluster
     - +0.004
     - −0.000
     - +0.008
     - +0.001
     - −0.081
     - +0.081
     - −0.000
     - +0.005
   * - ≥3-body
     - −0.042
     - −0.043
     - −0.004
     - −0.008
     - +0.086
     - −0.121
     - +0.009
     - −0.009

The error is in the split between dispersion and Pauli, and grows with size: the
modified (SAPT-like) Pauli, frozen − classical electrostatics − dispersion, is off by
+0.073 kJ/mol for the cluster and −0.124 in its ≥3-body part. The total, the classical
electrostatics, polarization and charge transfer stay within 0.01 kJ/mol.

**Cost** of shell A: 70.8 core-h, billed 203.6 SU on a node with 240 GB (46 SU/h, two
thirds of it for the memory), against 52.5 SU for GPU4PySCF on a V100: 3.9 times as
billed, 1.35 times in CPU time alone. For a typical exchange cluster (about 1500 nao,
21 pairs) the same scaling gives about 177 core-h, roughly 4.5 times the GPU as billed
(1.6 times in CPU time), and 250--340 GB of memory, at or beyond a 256 GB node.

**Conclusion** (science, 2026-11-09). The terms the fits read -- polarization, charge
transfer and its ≥3-body part, and the pairs' classical electrostatics and modified
Pauli -- are within the 0.05 kJ/mol gate up to the size of shell A (polarization and
charge transfer within 0.01, pair Pauli within 0.049). The cluster-level
dispersion/Pauli error (0.08 for the cluster, 0.12 in the ≥3-body modified Pauli) is
far below the fits' tolerances but beyond the gate. So Psi4 with density fitting is a
valid fallback for the fitted terms up to shell-A size; its real limits are the
memory (about 1500 nao needs 250--340 GB) and the cost, 3.9--4.5 times the GPU as
billed. Production labels stay ωB97M-V on GPUs, and Psi4 is the cross-check and
fallback.

Decisions for Paul
------------------
Decided (2026-11-07 to 11-09): the CPU backend targets ωB97X-V (coupled) and HF, and
refuses meta-GGAs rather than run them uncoupled; production labels stay ωB97M-V on
GPUs, with Psi4 as the cross-check and fallback (see the results above); milestone 8 is
not started. The questions as first put:

* Which functional the CPU backend targets first: ωB97X-V (coupled, Python only) or
  ωB97M-V (uncoupled, Python only, about 1 % of polarization moved to charge
  transfer).
* Whether to do milestone 8 (C++ meta-GGA kernel), and who contacts the Psi4
  developers about it -- and about contributing the ALMO-EDA driver to Psi4 itself.

References
----------
* R. Z. Khaliullin, E. A. Cobar, R. C. Lochan, A. T. Bell, M. Head-Gordon, J. Phys.
  Chem. A 111, 8753 (2007) -- ALMO-EDA, SCF-MI.
* P. R. Horn, M. Head-Gordon, J. Chem. Phys. 143, 114111 (2015) -- FERF polarization.
* P. R. Horn, Y. Mao, M. Head-Gordon, J. Chem. Phys. 144, 114107 (2016) -- frozen
  decomposition in DFT: electrostatics, Pauli, dispersion.
* Y. Mao, M. Loipersberger, P. R. Horn, et al., Annu. Rev. Phys. Chem. 72, 641 (2021)
  -- ALMO-EDA review.
* GPU4PySCF ``gpu4pyscf/properties/eda.py`` (1.9.1) -- the reference implementation.
