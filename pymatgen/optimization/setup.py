from distutils.core import setup
from Cython.Build import cythonize
import numpy as np
from distutils.extension import Extension


setup(
            ext_modules = [Extension("linear_assignment_cython", ["linear_assignment_cython.c"], include_dirs=[np.get_include()])]
            )
