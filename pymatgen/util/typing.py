"""This module defines convenience types for type hinting purposes.
Type hinting is new to pymatgen, so this module is subject to
change until best practices are established.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any, Union

from pymatgen.core import Composition, DummySpecies, Element, Species

if TYPE_CHECKING:  # needed to avoid circular imports
    from pymatgen.analysis.cost import CostEntry  # type: ignore[attr-defined]
    from pymatgen.analysis.phase_diagram import GrandPotPDEntry, PDEntry, TransformedPDEntry
    from pymatgen.entries import Entry
    from pymatgen.entries.computed_entries import ComputedEntry, ComputedStructureEntry, GibbsComputedStructureEntry
    from pymatgen.entries.exp_entries import ExpEntry


PathLike = Union[str, Path]

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
