import sys

from . import gauss_mod2 as gmod2

if sys.platform == "linux":
    from . import gauss_mod2_linbox as gmod2_linbox
    from . import gauss_mod2_m4ri as gmod2_m4ri
