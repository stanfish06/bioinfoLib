from juliacall import Main as julia

from bioinfoLib.julia.utils import ripserer_helper

from . import containers, utils

ripserer_helper()
julia.seval("using Ripserer")
