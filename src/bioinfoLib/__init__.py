"""
bioinfoLib - A comprehensive Python library for bioinformatics analysis.
"""

__version__ = "0.1.0"

# Import main modules
from . import GP
from . import imageAnalysis as img
from . import scRNAseq as scr

# Export commonly used components
__all__ = [
    "scr",
    "GP",
    "img",
]

# Add convenience imports for commonly used functions
try:
    from .GP import *
    from .imageAnalysis import *
    from .scRNAseq import *
except ImportError as e:
    print(f"Warning: Some modules could not be imported: {e}")
