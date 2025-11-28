# Copyright 2025 Zhiyuan Yu (Heemskerk's lab, University of Michigan)

import sys

if sys.platform == "linux":
    try:
        from . import gauss_mod2_linbox as gmod2_linbox
        from . import gauss_mod2_m4ri as gmod2_m4ri
    except Exception as e:
        print(f"Could not import gauss_mod2_linbox/m4ri {e}")
