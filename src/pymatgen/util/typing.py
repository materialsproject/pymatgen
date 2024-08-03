"""This module defines convenience types for type hinting purposes.
Type hinting is new to pymatgen, so this module is subject to
change until best practices are established.
"""

from __future__ import annotations

from collections.abc import Sequence
from os import PathLike as OsPathLike
from typing import TYPE_CHECKING, Any, Literal, Union

from numpy.typing import NDArray

from pymatgen.core import Composition, DummySpecies, Element, Species
from pymatgen.electronic_structure.core import Magmom, Spin

if TYPE_CHECKING:  # needed to avoid circular imports
    from pymatgen.analysis.cost import CostEntry  # type: ignore[attr-defined]
    from pymatgen.analysis.phase_diagram import GrandPotPDEntry, PDEntry, TransformedPDEntry
    from pymatgen.entries import Entry
    from pymatgen.entries.computed_entries import ComputedEntry, ComputedStructureEntry, GibbsComputedStructureEntry
    from pymatgen.entries.exp_entries import ExpEntry

# Commonly used composite types
Tuple3Ints = tuple[int, int, int]
Tuple3Floats = tuple[float, float, float]

PathLike = Union[str, OsPathLike]
PbcLike = tuple[bool, bool, bool]

# Things that can be cast to a Spin
SpinLike = Union[Spin, Literal[-1, 1, "up", "down"]]

# Things that can be cast to a magnetic moment
MagMomentLike = Union[float, Sequence[float], NDArray, Magmom]

# Things that can be cast to a Species-like object using get_el_sp
SpeciesLike = Union[str, Element, Species, DummySpecies]

# Things that can be cast to a Composition
CompositionLike = Union[str, Element, Species, DummySpecies, dict, Composition]

# Entry or any of its subclasses or dicts that can be unpacked into any of them
EntryLike = Union[
    dict[str, Any],
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

Vector3D = Tuple3Floats
Matrix3D = tuple[Vector3D, Vector3D, Vector3D]

SitePropsType = Union[list[dict[Any, Sequence[Any]]], dict[Any, Sequence[Any]]]

# Types specific to io.vasp
Kpoint = Union[Tuple3Floats, tuple[int,]]

# Miller index
MillerIndex = Tuple3Ints
