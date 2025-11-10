"""
GCTB: Genome-wide Complex Trait Bayesian Analysis
Python interface with C++ computational core

Authors: Jian Zeng, Luke Lloyd-Jones, Zhili Zheng, Shouye Liu
License: MIT
"""

# Import C++ core
try:
    from ._core import (
        # Data structures
        SnpInfo,
        IndInfo,
        Data,
        # Model and MCMC
        Model,
        MCMC,
        McmcSamples,
        # Factory functions
        build_model,
        build_model_summary,
        run_mcmc,
        # Utilities
        Timer,
        # Metadata
        __version__,
        __author__,
    )
except ImportError as e:
    raise ImportError(
        "Failed to import GCTB C++ core. "
        "Please ensure the package is properly installed with: pip install -e ."
    ) from e

# Optional diagnostics module (requires scipy)
try:
    from . import diagnostics
    _has_diagnostics = True
except ImportError:
    _has_diagnostics = False

# High-level workflows (pure Python)
from . import workflows
from .workflows import compute_credible_sets

__all__ = [
    # Data structures
    'SnpInfo',
    'IndInfo',
    'Data',
    # Model and MCMC
    'Model',
    'MCMC',
    'McmcSamples',
    # Factory functions  
    'build_model',
    'build_model_summary',
    'run_mcmc',
    # Utilities
    'Timer',
    '__version__',
    'workflows',
    'compute_credible_sets',
]

if _has_diagnostics:
    __all__.append('diagnostics')

# Version info
__version__ = "1.2.0"
__author__ = "Jian Zeng, Luke Lloyd-Jones, Zhili Zheng, Shouye Liu"

