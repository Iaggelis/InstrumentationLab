# python setup.py build_ext --inplace
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy


ext_modules = [
    Extension(
        "kalman",
        ["kalman_filter.pyx"],
        extra_compile_args=["-fopenmp","-march=native"],
        extra_link_args=["-fopenmp"],
        include_dirs=[numpy.get_include()],
    )
]

setup(name="kalman-parallel", ext_modules=cythonize(ext_modules))

