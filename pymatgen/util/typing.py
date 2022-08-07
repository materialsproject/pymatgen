# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module defines convenience types for type hinting purposes.
Type hinting is new to pymatgen, so this module is subject to
change until best practices are established.
"""

from __future__ import annotations

from pathlib import Path
from typing import Sequence, Union

import numpy as np

try:
    from numpy.typing import ArrayLike
except ImportError:
    ArrayLike = Union[Sequence[float], Sequence[Sequence[float]], Sequence[np.ndarray], np.ndarray]  # type: ignore
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import DummySpecies, Element, Species

VectorLike = Union[Sequence[float], np.ndarray]
MatrixLike = Union[Sequence[Sequence[float]], Sequence[np.ndarray], np.ndarray]

PathLike = Union[str, Path]

# Things that can be cast into a Species-like object using get_el_sp
SpeciesLike = Union[str, Element, Species, DummySpecies]

# Things that can be cast into a Composition
CompositionLike = Union[str, Element, Species, DummySpecies, dict, Composition]
