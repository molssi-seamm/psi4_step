# -*- coding: utf-8 -*-

"""
psi4_step
A step for Psi4 in a SEAMM flowchart
"""

# Bring up the classes so that they appear to be directly in
# the psi4_step package.

from psi4_step.psi4_metadata import methods  # noqa: F401
from psi4_step.psi4_metadata import dft_functionals  # noqa: F401
from psi4_step.psi4_metadata import optimization_methods  # noqa: F401
from psi4_step.psi4_metadata import optimization_convergence  # noqa: F401
from psi4_step.psi4_metadata import metadata  # noqa: F401

from psi4_step.psi4 import Psi4  # noqa: F401
from psi4_step.psi4_parameters import Psi4Parameters  # noqa: F401
from psi4_step.psi4_step import Psi4Step  # noqa: F401
from psi4_step.tk_psi4 import TkPsi4  # noqa: F401

from psi4_step.initialization import Initialization  # noqa: F401
from psi4_step.initialization_parameters import InitializationParameters  # noqa: F401
from psi4_step.initialization_step import InitializationStep  # noqa: F401
from psi4_step.tk_initialization import TkInitialization  # noqa: F401

from psi4_step.energy import Energy  # noqa: F401
from psi4_step.energy_parameters import EnergyParameters  # noqa: F401
from psi4_step.energy_step import EnergyStep  # noqa: F401
from psi4_step.tk_energy import TkEnergy  # noqa: F401

from psi4_step.optimization import Optimization  # noqa: F401
from psi4_step.optimization_parameters import OptimizationParameters  # noqa: F401
from psi4_step.optimization_step import OptimizationStep  # noqa: F401
from psi4_step.tk_optimization import TkOptimization  # noqa: F401

# Handle versioneer
from ._version import get_versions

__author__ = """Paul Saxe"""
__email__ = "psaxe@molssi.org"
versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions
