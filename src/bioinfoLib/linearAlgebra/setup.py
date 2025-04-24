from Cython.Build import cythonize
from setuptools import setup, Extension

modules = [
    Extension(
        "gauss_mod2_m4ri",
        sources = ["gauss_mod2_m4ri.pyx"],
        libraries = ["m4ri"],
        library_dirs = ["/usr/local/lib"],
        include_dirs = ["/usr/local/include"]
    ), 
    Extension(
        "gauss_mod2",
        sources = ["gauss_mod2.pyx"]
    )
]

setup(
    name="gaussian_elimination",
    ext_modules=cythonize(modules),  # gdb_debug=True
)
