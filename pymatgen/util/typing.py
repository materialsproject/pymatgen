# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module defines convenience types for type hinting purposes.
Type hinting is new to pymatgen, so this module is subject to
change until best practices are established.
"""

from typing import Union, List
from pathlib import Path
import numpy as np

Vector3Like = Union[List[float], np.ndarray]
PathLike = Union[str, Path]
