from setuptools import setup
from Cython.Build import cythonize

setup(
    name = 'test2',
    ext_modules = cythonize("gauss_mod2.pyx", gdb_debug=True),
)
