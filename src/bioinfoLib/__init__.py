"""
bioinfoLib - A comprehensive Python library for bioinformatics analysis.
"""

__version__ = "0.1.0"

# Import main modules
from . import GP, linalg
from . import imageAnalysis as img
from . import julia as jl
from . import scRNAseq as scr

# Export commonly used components
__all__ = ["scr", "GP", "img", "jl"]

# Add convenience imports for commonly used functions
try:
    from .GP import *
    from .imageAnalysis import *
    from .julia import *
    from .linalg import *
    from .scRNAseq import *
except ImportError as e:
    print(f"Warning: Some modules could not be imported: {e}")
