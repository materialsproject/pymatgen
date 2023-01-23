from __future__ import annotations

from typing import Any

# pymatgen.entries needs to be imported before pymatgen.util.typing
# to avoid circular import.
from pymatgen.entries import Entry
from pymatgen.util.typing import (
    CompositionLike,
    EntryLike,
    MatrixLike,
    PathLike,
    SpeciesLike,
    VectorLike,
)

# This module tests types are as expected and can be imported without circular ImportError.

__author__ = "Janosh Riebesell"
__date__ = "2022-10-20"
__email__ = "janosh@lbl.gov"


def _type_str(some_type: Any) -> str:
    return str(some_type).replace("typing.", "").replace("pymatgen.core.periodic_table.", "")


def test_entry_like():
    entries = (
        "Entry",
        "ComputedEntry",
        "ComputedStructureEntry",
        "PDEntry",
        "ExpEntry",
        "TransformedPDEntry",
        "GrandPotPDEntry",
        "CostEntry",
        "GibbsComputedStructureEntry",
    )
    type_str = _type_str(EntryLike)
    for entry in entries:
        assert entry in type_str
    assert Entry.__name__ in str(EntryLike)


def test_species_like():
    assert _type_str(SpeciesLike) == "Union[str, Element, Species, DummySpecies]"


def test_composition_like():
    assert (
        _type_str(CompositionLike)
        == "Union[str, Element, Species, DummySpecies, dict, pymatgen.core.composition.Composition]"
    )


def test_vector_like():
    assert _type_str(VectorLike) == "Union[Sequence[float], numpy.ndarray]"


def test_matrix_like():
    assert _type_str(MatrixLike) == "Union[Sequence[Sequence[float]], Sequence[numpy.ndarray], numpy.ndarray]"


def test_path_like():
    assert _type_str(PathLike) == "Union[str, pathlib.Path]"
