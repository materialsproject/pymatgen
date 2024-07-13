"""Pymatgen package configuration."""

from __future__ import annotations

import platform
import sys

import numpy as np
from setuptools import Extension, setup

is_win_64 = sys.platform.startswith("win") and platform.machine().endswith("64")
extra_link_args = ["-Wl,--allow-multiple-definition"] if is_win_64 else []

setup(
    ext_modules=[
        Extension(
            "pymatgen.optimization.linear_assignment",
            ["src/pymatgen/optimization/linear_assignment.pyx"],
            extra_link_args=extra_link_args,
        ),
        Extension(
            "pymatgen.util.coord_cython",
            ["src/pymatgen/util/coord_cython.pyx"],
            extra_link_args=extra_link_args,
        ),
        Extension(
            "pymatgen.optimization.neighbors",
            ["src/pymatgen/optimization/neighbors.pyx"],
            extra_link_args=extra_link_args,
        ),
    ],
    include_dirs=[np.get_include()],
)
