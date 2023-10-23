=======
History
=======

2023.8.22 -- Enhancement of orbital plots
   * Added structure to the orbital plots
   * Fixed a bug if the default nouber of cores was not 'available'

2023.6.4 -- Enhancements
   * Added thermochemistry substep to compute the vibrational and other corrections for
     thermochemistry.
   * Added ability to create .cube files for plotting the density and orbitals.
     
2023.2.21 -- Added control over SCF convergence
   * Able to set convergence criteria.
   * Options for damping. level shifting and second-order SCF
   * Set default for number of steps for optimization to 6*nAtoms to
     make it more explicit.
     
2023.2.17 -- Checking more thoroughly for errors in Psi4
   * Check and throw errors for various issues when running Psi4, including whenever it
     does not complete successfully according to the output.
     
2023.2.16.1 -- OptKing will not die!!!
   * Yet another part needed to be removed.
     
2023.2.16 -- Thoroughly removing calls to OptKing
   * OptKing is causing compatibility issues, and is not currently being used in the
     Python, so am temporarily removing it.

2023.1.23 (23 January 2023)
---------------------------

* Corrected bugs, including if no dispersion term requested
* Added preliminary version of accelerated optimization using
  MOPAC for Hessians, etc. May not fully work at the moment.
* Updated documentation theme and structure to match new style.

2022.11.18 -- Added total charge and multiplicity to Psi4 input

2021.2.11 (11 February 2021)
----------------------------

* Updated the README file to give a better description.
* Updated the short description in setup.py to work with the new installer.
* Added keywords for better searchability.

2021.2.4 (4 February 2021)
--------------------------

* Updated for compatibility with the new system classes in MolSystem
  2021.2.2 release.

2020.12.5 (5 December 2020)
---------------------------

* Internal: switching CI from TravisCI to GitHub Actions, and in the
  process moving documentation from ReadTheDocs to GitHub Pages where
  it is consolidated with the main SEAMM documentation.

2020.11.2 (2 November 2020)
---------------------------

* Updated to be compatible with the new command-line argument
  handling.

2020.9.28 (28 September 2020)
-----------------------------

* Improved handling of results

  - Selection depends on both calculation and method.
  - Arrays can be stored in variables (but not tables!).
  - Multipoles gathered as arrays rather than a large number of scalars.

2020.8.1 (1 August 2020)
------------------------

* Fixed minor bugs.

2020.7.3 (24 July 2020)
------------------------

* Added control options for number of threads and amount of memory to
  use.

2020.7.2 (23 July 2020)
------------------------

* First pass at specifying number of threads and amount of memory.

2020.7.1 (23 July 2020)
------------------------

* Added substep for optimization.

2020.7.0 (21 July 2020)
------------------------

* First release on PyPI of initial working version.
