# -*- coding: utf-8 -*-

"""
psi4_step
A step for Psi4 in a SEAMM flowchart
"""

# Bring up the classes so that they appear to be directly in
# the psi4_step package.

from .psi4_metadata import methods  # noqa: F401
from .psi4_metadata import dft_functionals  # noqa: F401
from .psi4_metadata import optimization_methods  # noqa: F401
from .psi4_metadata import optimization_convergence  # noqa: F401
from .psi4_metadata import metadata  # noqa: F401

from .psi4 import Psi4  # noqa: F401
from .psi4_parameters import Psi4Parameters  # noqa: F401
from .psi4_step import Psi4Step  # noqa: F401
from .tk_psi4 import TkPsi4  # noqa: F401

from .initialization import Initialization  # noqa: F401
from .initialization_parameters import InitializationParameters  # noqa: F401
from .initialization_step import InitializationStep  # noqa: F401
from .tk_initialization import TkInitialization  # noqa: F401

from .energy import Energy  # noqa: F401
from .energy_parameters import EnergyParameters  # noqa: F401
from .energy_step import EnergyStep  # noqa: F401
from .tk_energy import TkEnergy  # noqa: F401

from .optimization import Optimization  # noqa: F401
from .optimization_parameters import OptimizationParameters  # noqa: F401
from .optimization_step import OptimizationStep  # noqa: F401
from .tk_optimization import TkOptimization  # noqa: F401

# from .accelerated_optimization import AcceleratedOptimization  # noqa: F401
# from .accelerated_optimization_parameters import (  # noqa: F401
#     AcceleratedOptimizationParameters,
# )
# from .accelerated_optimization_step import (  # noqa: F401
#     AcceleratedOptimizationStep,
# )
# from .tk_accelerated_optimization import (  # noqa: F401
#     TkAcceleratedOptimization,
# )

from .thermochemistry_step import ThermochemistryStep  # noqa: F401
from .thermochemistry import Thermochemistry  # noqa: F401
from .thermochemistry_parameters import ThermochemistryParameters  # noqa: F401
from .tk_thermochemistry import TkThermochemistry  # noqa: F401

# Handle versioneer
from ._version import get_versions

__author__ = """Paul Saxe"""
__email__ = "psaxe@molssi.org"
versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions
