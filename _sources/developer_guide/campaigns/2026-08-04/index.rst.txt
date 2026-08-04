2026-08-04 BSSE (counterpoise) sub-step
========================================

Adds a **BSSE** sub-step, computing the counterpoise-corrected
(Boys--Bernardi) energy and gradient of an N-fragment complex with
independent per-fragment charge -- the Psi4 counterpart of ``orca_step``'s
BSSE sub-step.

The canonical design document -- the physics, the architecture decision, and
the ``orca_step`` implementation -- lives in the sibling ``seamm_bsse``
repository:
``seamm_bsse/docs/developer_guide/campaigns/2026-08-03/bsse_architecture.rst``.
This page is the Psi4-side implementation record.

Architecture: a thin wrapper, not a second ``seamm_bsse.combine()`` consumer
--------------------------------------------------------------------------------

Unlike ``orca_step``'s BSSE sub-step, this one does **not** use
``seamm_bsse.generate_job_specs()``/``combine()``. Psi4 has a native
N-fragment counterpoise driver -- ``bsse_type=['cp', 'nocp']`` on
``energy()``/``gradient()``, given a ``--``-separated, per-fragment-charge
molecule block -- that does the whole 2N + 1-equivalent orchestration
internally and returns an already-full-cluster-indexed corrected
energy/gradient in one call. ORCA has no such driver, which is why
``seamm_bsse`` exists for it; Psi4 does not have that gap. ``seamm_bsse`` is
still used here for ``Fragment``/``validate_fragments`` -- the
fragment-definition/charge-validation layer, shared with ``orca_step`` for
parameter/GUI consistency -- but ``generate_job_specs``/``combine`` go
unused.

Requesting both ``'cp'`` and ``'nocp'`` costs no extra energy evaluations
(both reuse the same 1-body/N-body sub-calculations) and gives the
uncorrected total/interaction energies ``orca_step``'s BSSE also reports, so
the two sub-steps' result shape matches: ``energy``, ``uncorrected energy``,
``bsse correction``, ``interaction energy``, ``uncorrected interaction
energy``, ``gradients``.

A structural difference from ``orca_step``: Psi4 sub-steps do not drive their
own execution (no ``run()``/executor call). ``psi4_step``'s main ``Psi4``
node builds *one* shared, multi-section Psi4 input file by concatenating
every sub-step's ``get_input()`` text and runs it as a single Psi4 process;
each sub-step's ``analyze()`` then parses the shared output afterward. Since
every node's Python-side result-writing code executes with the *main* node's
directory as its current working directory, ``BSSE.get_input()`` embeds an
**absolute path** to its own sub-step directory (``<self.directory>/
bsse.json``) rather than a bare filename, sidestepping the multi-node
shared-directory collision the other sub-steps handle differently (an
``@<node id>+`` filename prefix convention). ``BSSE.get_input()`` also
builds its own fragment-aware, ``no_com``/``no_reorient`` molecule block --
``Psi4._convert_structure`` reads the single global "current" configuration,
with no concept of fragments or ghosts.

Validation
----------

:download:`validate_psi4_bsse.py <validate_psi4_bsse.py>`
   Real Psi4 (1.10, local ``seamm-psi4`` conda env), water dimer
   (``H2O-H2O.sdf``, the same geometry/system as ``orca_step``'s M2
   regression), HF/def2-SVP. Confirms the whole pipeline works end to end:
   fragment auto-detection, the molecule-block builder, the
   ``bsse_type=['cp', 'nocp']`` call, the JSON round trip through the shared
   process, and ``analyze()``'s parsing. Cross-checked against the
   already-validated ORCA HF/def2-SVP result for the identical geometry --
   uncorrected/corrected energies agree to ~2e-4 E\ :sub:`h`, gradient
   components to ~5e-5 E\ :sub:`h`/bohr. This is *not* the SCF-noise-level
   agreement ORCA's own Compound-script regression showed (~1e-9) -- that
   was the same code compared with itself; ~2e-4 E\ :sub:`h` between two
   independently-implemented quantum chemistry codes at the same nominal
   level of theory is the expected, healthy level of cross-code agreement,
   and a genuine independent validation of the counterpoise machinery (unlike
   ``orca_step``'s Compound-script comparison, which validates the new
   *plumbing* against old plumbing of the *same* code).

:download:`validate_psi4_bsse_m3.py <validate_psi4_bsse_m3.py>`
   Charged fragments and N = 3, **through this sub-step's own
   ``get_input()``/``analyze()``** (not hand-written ``psi4`` scripts) --
   the gap the first validation script left open. Same geometries and level
   of theory (B3LYP-D3BJ/def2-TZVP) as ``orca_step``'s M3/N3 legs, so the
   two codes' numbers are directly comparable:

   .. list-table::
      :header-rows: 1

      * - System
        - Psi4
        - ORCA (M3/N3)
      * - Na\ :sup:`+`\ ···Cl\ :sup:`-`
        - -136.352
        - -136.289
      * - Na\ :sup:`+`\ ···H\ :sub:`2`\ O
        - -26.176
        - -26.110
      * - Cl\ :sup:`-`\ ···H\ :sub:`2`\ O
        - -15.699
        - -15.638
      * - Na\ :sup:`+`\ ···Cl\ :sup:`-`\ ···H\ :sub:`2`\ O (N=3)
        - -150.332
        - -150.204

   (kcal/mol, CP-corrected interaction energy.) Agreement is ~0.05-0.13
   kcal/mol throughout -- the same healthy, code-vs-code level the neutral
   water-dimer check showed, now confirmed for per-fragment charge (``"1;
   2"``/``"1, -1"`` in ``specified`` mode -- these hand-built configurations
   have no bond table, so ``"auto (molecules)"`` would see every atom as its
   own fragment, the same gotcha ``orca_step``'s M3/N3 scripts hit) and for
   N = 3 (2N + 1 = 7 equivalent internal sub-calculations, all through
   Psi4's native driver in one call). One bug caught and fixed *in this
   validation script* (not the sub-step): an early run divided the raw
   Hartree value in ``bsse.json`` by 4.184 (the kJ->kcal factor) instead of
   converting Hartree->kcal directly, producing near-zero "interaction
   energies" even though ``bsse.json`` and Psi4's own printed ManyBody table
   already agreed -- a reminder to sanity-check a validation script's own
   arithmetic before suspecting the code under test.

Not done
--------

* Unit tests for the fragment/molecule-block-building logic in isolation
  (no Psi4) -- done (``tests/test_bsse.py``), except real-Psi4 execution
  itself is (deliberately) not unit-tested, only exercised by the
  validation scripts above.
* Energy-of-formation (``DfE0``) support -- ``orca_step``'s BSSE computes
  this via ``seamm_thermochemistry``; ``psi4_step`` uses a different,
  older, self-contained mechanism (``Energy.calculate_enthalpy_of_formation``,
  tabulated atom energies) not yet wired up for BSSE.
* Advanced SCF convergence controls (damping, level shift, SOSCF, stability
  analysis) that ``Energy.get_input()`` supports -- omitted here to keep the
  wrapper thin; can be added if a production case needs them.
