# Basis files shipped with the Psi4 step

A basis file holds Psi4 input blocks that the Initialization step copies into the input
when "Basis set from" is "basis file": a `basis {...}` block for the orbital basis and,
optionally, `df_basis_cc {...}` and `df_basis_scf {...}` blocks for the density-fitting
sets. The file's name without `.basis` is the basis set's name in the model, e.g.
`DF-CCSD(T)/cc-pwCVTZ-ORCA`.

## cc-pwCVTZ-ORCA.basis

cc-pwCVTZ exactly as ORCA 6.1.1 has it internally, so that Psi4 calculations use the same
basis as ORCA ones (e.g. CCSD(T)/cc-pwCVTZ geometries computed with ORCA, re-optimized in
Psi4 with analytic DF-CCSD(T) gradients). It is spherical (puream).

- **Orbital basis:** cc-pwCVTZ for B-Ne, Al-Ar and Ge-Br, cc-pVTZ for H and He.
  - **Ge, As, Se, Br:** from ORCA 6.1.1, which has an all-electron cc-pwCVTZ for them.
    Neither the Basis Set Exchange (0.12) nor Psi4's own library has one.
  - **H, B, C, O, Si, S, Cl:** from ORCA. For each of them, the script checks that ORCA's
    functions span the same space as the Basis Set Exchange's cc-pwCVTZ (cc-pVTZ for H).
    ORCA removes from the contractions the primitives that are also separate shells,
    which changes the coefficients but not the functions.
  - **He, N, F, Ne, Al, P, Ar:** from the Basis Set Exchange 0.12, since they were not
    exported from ORCA.
  - Ga and Kr are not included.
- **Correlation fitting** (`df_basis_cc`): cc-pwCVTZ-RIFIT from the Basis Set Exchange
  0.12, which fits the core-valence correlation (Ge-Br included); cc-pVTZ-RIFIT for H and
  He, which cc-pwCVTZ-RIFIT does not cover.
- **SCF fitting** (`df_basis_scf`): def2-universal-JKFIT from the Basis Set Exchange 0.12.
- **Checks:** GeH4 has 140 basis functions and GeO 127, as in ORCA. The tests check these.

### How it was made

The ORCA exports are one output per element, from the input
`! UHF cc-pwCVTZ PrintBasis NoIter` with ORCA 6.1.1 (2026-11-18). Their sections
"BASIS SET IN INPUT FORMAT" are kept in `devtools/basis/orca_6.1.1_cc-pwCVTZ/<El>.out`
in the repository. To rebuild the file from them:

```
python devtools/basis/make_cc_pwcvtz_orca.py devtools/basis/orca_6.1.1_cc-pwCVTZ \
    psi4_step/data/basis/cc-pwCVTZ-ORCA.basis
```

To add an element, export it from ORCA the same way and rerun the script.
