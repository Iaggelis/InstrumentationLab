# python setup.py build_ext --inplace
from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(
    name="smoothing",
    ext_modules=cythonize("smooth.pyx"),
    include_dirs=[numpy.get_include()],
)
# setup(name="event display", ext_modules=cythonize("ev_dis.pyx"))

