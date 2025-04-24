from Cython.Build import cythonize
from setuptools import setup

setup(
    name="gauss_mod2",
    ext_modules=cythonize("gauss_mod2.pyx"),  # gdb_debug=True
)
