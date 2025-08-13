"""
bioinfoLib - A comprehensive Python library for bioinformatics analysis.
"""

__version__ = "0.1.1"

# Import main modules
try:
    from . import GP
except Exception as e:
    print(f"Could not import GP modules {e}")
try:
    from . import julia as jl
except Exception as e:
    print(f"Could not import julia modules {e}")
try:
    from . import linearAlgebra as la
except Exception as e:
    print(f"Could not import linearAlgebra modules {e}")
try:
    from . import RNAseq as rn
except Exception as e:
    print(f"Could not import RNAseq modules {e}")
try:
    from . import scRNAseq as scr
except Exception as e:
    print(f"Could not import scRNAseq modules {e}")
try:
    from . import topology as tp
except Exception as e:
    print(f"Could not import topology modules {e}")

# Export commonly used components
__all__ = ["rn", "scr", "GP", "jl", "la", "tp"]
