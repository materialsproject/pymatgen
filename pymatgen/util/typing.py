# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module defines convenience types for type hinting purposes.
Type hinting is new to pymatgen, so this module is subject to
change until best practices are established.
"""

from pathlib import Path
from typing import Union, Sequence

import numpy as np

try:
    from numpy.typing import ArrayLike
except ImportError:
    ArrayLike = Union[Sequence[float], Sequence[Sequence[float]], Sequence[np.ndarray], np.ndarray]  # type: ignore

VectorLike = Union[Sequence[float], np.ndarray]
MatrixLike = Union[Sequence[Sequence[float]], Sequence[np.ndarray], np.ndarray]

PathLike = Union[str, Path]
