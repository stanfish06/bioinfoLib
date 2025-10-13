"""
bioinfoLib - A comprehensive Python library for bioinformatics analysis.
"""

import logging

# Set up environment variables
from dotenv import load_dotenv

from . import utils

__version__ = "0.1.3"
load_dotenv()

logger = logging.getLogger(__name__)

# Import main modules
try:
    from . import GP
except ImportError as e:
    logger.warning(f"Could not import GP modules: {e}")
    GP = None
try:
    from . import julia as jl
except ImportError as e:
    logger.warning(f"Could not import julia modules: {e}")
    jl = None
try:
    from . import linearAlgebra as la
except ImportError as e:
    logger.warning(f"Could not import linearAlgebra modules: {e}")
    la = None
try:
    from . import RNAseq as rn
except ImportError as e:
    logger.warning(f"Could not import RNAseq modules: {e}")
    rn = None
try:
    from . import scRNAseq as scr
except ImportError as e:
    logger.warning(f"Could not import scRNAseq modules: {e}")
    scr = None
try:
    from . import topology as tp
except ImportError as e:
    logger.warning(f"Could not import topology modules: {e}")
    tp = None
try:
    from . import image as im
except ImportError as e:
    logger.warning(f"Could not import image modules: {e}")
    im = None
try:
    from . import R as r
except ImportError as e:
    logger.warning(f"Could not import R modules: {e}")
    r = None
try:
    from . import neuralNet as nn
except ImportError as e:
    logger.warning(f"Could not import neuralNet modules: {e}")
    nn = None

# Export commonly used components (only those that were successfully imported)
__all__ = [
    name
    for name in ["rn", "scr", "GP", "jl", "la", "tp", "im", "r", "nn"]
    if globals().get(name) is not None
]
