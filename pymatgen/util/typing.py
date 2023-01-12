# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module defines convenience types for type hinting purposes.
Type hinting is new to pymatgen, so this module is subject to
change until best practices are established.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, Sequence, Union

import numpy as np

from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import DummySpecies, Element, Species

try:
    from numpy.typing import ArrayLike
except ImportError:
    ArrayLike = Union[Sequence[float], Sequence[Sequence[float]], Sequence[np.ndarray], np.ndarray]  # type: ignore

if TYPE_CHECKING:  # needed to avoid circular imports
    from pymatgen.analysis.cost import CostEntry  # type: ignore[attr-defined]
    from pymatgen.analysis.phase_diagram import (
        GrandPotPDEntry,
        PDEntry,
        TransformedPDEntry,
    )
    from pymatgen.entries import Entry
    from pymatgen.entries.computed_entries import (
        ComputedEntry,
        ComputedStructureEntry,
        GibbsComputedStructureEntry,
    )
    from pymatgen.entries.exp_entries import ExpEntry

VectorLike = Union[Sequence[float], np.ndarray]
MatrixLike = Union[Sequence[Sequence[float]], Sequence[np.ndarray], np.ndarray]

PathLike = Union[str, Path]

# Things that can be cast to a Species-like object using get_el_sp
SpeciesLike = Union[str, Element, Species, DummySpecies]

# Things that can be cast to a Composition
CompositionLike = Union[str, Element, Species, DummySpecies, dict, Composition]

# Entry or any of its subclasses or dicts that can be unpacked into any of them
EntryLike = Union[
    Dict[str, Any],
    "Entry",
    "PDEntry",
    "ComputedEntry",
    "ComputedStructureEntry",
    "ExpEntry",
    "TransformedPDEntry",
    "GrandPotPDEntry",
    "CostEntry",
    "GibbsComputedStructureEntry",
]
