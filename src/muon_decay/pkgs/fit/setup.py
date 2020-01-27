# python setup.py build_ext --inplace
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize


ext_modules = [
    Extension(
        "fit",
        ["fit.pyx"],
    )
]

setup(name="kalman-parallel", ext_modules=cythonize(ext_modules))
