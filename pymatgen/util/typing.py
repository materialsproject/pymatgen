"""This module defines convenience types for type hinting purposes.
Type hinting is new to pymatgen, so this module is subject to
change until best practices are established.
"""

from __future__ import annotations

from collections.abc import Sequence
from os import PathLike as OsPathLike
from typing import TYPE_CHECKING, Any, Union

from pymatgen.core import Composition, DummySpecies, Element, Species

if TYPE_CHECKING:  # needed to avoid circular imports
    from pymatgen.analysis.cost import CostEntry  # type: ignore[attr-defined]
    from pymatgen.analysis.phase_diagram import GrandPotPDEntry, PDEntry, TransformedPDEntry
    from pymatgen.entries import Entry
    from pymatgen.entries.computed_entries import ComputedEntry, ComputedStructureEntry, GibbsComputedStructureEntry
    from pymatgen.entries.exp_entries import ExpEntry


PathLike = Union[str, OsPathLike]
PbcLike = tuple[bool, bool, bool]

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

Vector3D = tuple[float, float, float]
Matrix3D = tuple[Vector3D, Vector3D, Vector3D]

SitePropsType = Union[list[dict[Any, Sequence[Any]]], dict[Any, Sequence[Any]]]

# Types specific to io.vasp
Kpoint = Union[tuple[float, float, float], tuple[int,]]
