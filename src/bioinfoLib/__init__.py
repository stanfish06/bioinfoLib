"""
bioinfoLib - A comprehensive Python library for bioinformatics analysis.
"""

__version__ = "0.1.0"

# Import main modules
from . import GP
from . import julia as jl
from . import linearAlgebra as la
from . import scRNAseq as scr

# Export commonly used components
__all__ = ["scr", "GP", "jl", "la"]
