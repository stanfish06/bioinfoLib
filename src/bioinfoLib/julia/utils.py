# Copyright 2025 Zhiyuan Yu (Heemskerk's lab, University of Michigan)

from juliacall import Main as julia


# not sure why julia automatically load my script but this seems not required
def ripserer_helper():
    # Initialize required Julia packages
    julia.seval("using Ripserer")
    julia.seval("using Base.Threads")

    # Get the path to the Julia helper file
    import os

    current_dir = os.path.dirname(os.path.abspath(__file__))
    julia_file = os.path.join(current_dir, "ripserer_helper.jl")

    # Load the Julia functions from file
    julia.seval(f'include("{julia_file}")')
