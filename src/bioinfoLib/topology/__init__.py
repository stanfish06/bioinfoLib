from juliacall import Main as julia

from bioinfoLib.julia.setup import install_packages
from bioinfoLib.julia.utils import ripserer_helper

from . import containers, utils

install_packages()
ripserer_helper()
julia.seval("using Ripserer")
