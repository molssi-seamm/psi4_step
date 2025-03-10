=======
History
=======
2025.3.9 -- Optimization and thermochemistry enhancements
   * Expanded the output to include gap, HOMO, LUMO, RMSD between initial and final
     structure for optimizations, dipole moment, etc.
   * Added control for updating the structure after optimization.
   * Added enthalpy of formation to the thermochemistry output for known basis sets.
   * Updated the installation to work with recent changes.
     
2024.12.7 -- Minor update due to changes in molsystem
   * The molsystem was updated to improve the handling of properties. This required
     changes in this module to match the new molsystem.
    
2024.10.15 -- Bugfix: error if used in a loop and previous directories deleted.
   * The code crashed if called with a loop in the flowchart, and the last directory of
     a previous loop iteration was deleted before running the next iteration.
     
2024.10.5 -- Enhancements and bug fixes for thermochemistry
   * Improved GUI for thermochemistry so that it automatically recognizes whether it is
     after e.g. an optimization and configures appropriately.
   * Fixed and issue with transferring the multipole moments to the JSON file
     
2024.7.30 -- Enhanced results and fixed problem with psi4.ini
   * Added to the results by using cclib to parse the output.
   * If ~/SEAMM/psi4.ini did not exist the version that was automatically created was
     not complete, causing calculations to fail.
     
2024.5.23.3 -- Added auxiliary codes that help Psi4
   * Psi4 uses codes like DFTD4 and DFTD4 for parts of e.g. dispersion
     calculations. This release adds them to the Conda install SEAMM uses for Psi4.
     
2024.5.23.2 -- Bugfix: incorrect name for gradients
   * There was a typo in the name for the gradients, such that they could not be output
     to Results.json.
   * The units for the energy and gradients in Results.json were incorrect.
     
2024.5.23.1 -- Internal fix for creating Docker image.

2024.5.23 -- Added standard energy and gradients to results
   * Added 'energy' and 'gradients' to optional results to support e.g. Energy Scan
   * Fixed crashing bug in description of the Energy substep.
     
2024.3.17 -- Updated the installer
   * Updated the installer to use the new version of the SEAMM installer.
   * Finalizes installing either with Conda or Docker
     
2024.3.4 -- Allowing short names for method and DFT functionals
   * Added short names for the methods (Hamiltonians)  and DFT functionals.
   * Catch errors in Psi4 calculating properties for e.g. CISD(T) method

2024.2.29 -- Completed support for containers
   * Fixed issues with running amd64 containers on arm64 systems.
     
2024.1.11 -- Changes to allow running in containers.
   * Moved to the new executor and ensured it still runs directly.
   * Fixed bugs in printing the summary output.

2023.8.23 -- Fix for installation of Psi4
   * Psi4 is now available on CondaForge, so install from there if requested.
   * Psi4 crashes if asked to optimize and atom, so change to doing a single point
     energy is optimization is requested for and atom.

2023.8.22 -- Enhancement of orbital plots
   * Added structure to the orbital plots
   * Fixed a bug if the default number of cores was not 'available'

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
