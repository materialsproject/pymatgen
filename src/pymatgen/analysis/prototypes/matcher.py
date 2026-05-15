"""In this module, the PrototypeDatabaseMatcher defaults to the AFLOW LIBRARY OF
CRYSTALLOGRAPHIC PROTOTYPES. If using the default library, please cite their
publication appropriately:

Mehl, M. J., Hicks, D., Toher, C., Levy, O., Hanson, R. M., Hart, G., & Curtarolo, S. (2017).
The AFLOW library of crystallographic prototypes: part 1.
Computational Materials Science, 136, S1-S828.
https://doi.org/10.1016/j.commatsci.2017.01.017
"""

from __future__ import annotations

import os
from typing import Any

from monty.dev import deprecated
from monty.serialization import loadfn

from pymatgen.core import Structure
from pymatgen.core.structure_matcher import StructureMatcher
from pymatgen.util.due import Doi, due

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
AFLOW_PROTOTYPE_LIBRARY = f"{MODULE_DIR}/aflow_prototypes.json.gz"


class PrototypeDatabaseMatcher:
    """
    This class will match structures to crystal prototypes from a database, and will
    attempt to group species together to match structures derived from
    prototypes (e.g. an A_xB_1-x_C from a binary prototype), and will
    give these the names the "-like" suffix.

    By default, this class uses data from the AFLOW LIBRARY OF CRYSTALLOGRAPHIC
    PROTOTYPES. If using the default library, please cite their publication
    appropriately:

    Mehl, M. J., Hicks, D., Toher, C., Levy, O., Hanson, R. M., Hart, G., & Curtarolo, S. (2017).
    The AFLOW library of crystallographic prototypes: part 1.
    Computational Materials Science, 136, S1-S828.
    https://doi.org/10.1016/j.commatsci.2017.01.017
    """

    def __init__(
        self,
        initial_ltol: float = 0.2,
        initial_stol: float = 0.3,
        initial_angle_tol: float = 5,
        prototype_db: Any = None,
    ) -> None:
        """
        Tolerances as defined in StructureMatcher. Tolerances will be
        gradually decreased until only a single match is found (if possible).

        Args:
            initial_ltol (float): fractional length tolerance.
            initial_stol (float): site tolerance.
            initial_angle_tol (float): angle tolerance.
            prototype_db: Path to a monty-loadable prototype database or an
                iterable of prototype entries. Entries must contain either "snl"
                or "structure".
        """
        self.initial_ltol = initial_ltol
        self.initial_stol = initial_stol
        self.initial_angle_tol = initial_angle_tol
        if prototype_db is None:
            prototype_db = AFLOW_PROTOTYPE_LIBRARY

        self._prototype_db: list[tuple[Structure, dict]] = []
        for dct in self._iter_db_entries(prototype_db):
            structure = self._get_entry_structure(dct)
            reduced_structure = self._preprocess_structure(structure)
            self._prototype_db.append((reduced_structure, self._get_entry_data(dct)))

    @staticmethod
    def _iter_db_entries(prototype_db: Any) -> Any:
        if isinstance(prototype_db, (str, os.PathLike)):
            prototype_db = loadfn(prototype_db)

        if isinstance(prototype_db, dict) and {"columns", "data"} <= set(prototype_db):
            columns = prototype_db["columns"]
            yield from (dict(zip(columns, row, strict=True)) for row in prototype_db["data"])
            return

        if hasattr(prototype_db, "iterrows"):
            yield from (row for _, row in prototype_db.iterrows())
        else:
            yield from prototype_db

    @staticmethod
    def _get_entry_structure(dct: dict) -> Structure:
        if "snl" in dct:
            return dct["snl"].structure
        if "structure" in dct:
            structure = dct["structure"]
            return Structure.from_dict(structure) if isinstance(structure, dict) else structure
        raise KeyError("Prototype database entries must contain either 'snl' or 'structure'.")

    @staticmethod
    def _get_entry_data(dct: dict) -> dict:
        if "mineral" in dct and "structure" in dct:
            structure = dct["structure"]
            return {
                "type": dct["mineral"],
                "distance": dct.get("distance", -1),
                "structure": Structure.from_dict(structure) if isinstance(structure, dict) else structure,
            }
        return dict(dct)

    @staticmethod
    def _preprocess_structure(structure: Structure) -> Structure:
        return structure.get_reduced_structure(reduction_algo="niggli").get_primitive_structure()

    def _match_prototype(
        self,
        structure_matcher: StructureMatcher,
        reduced_structure: Structure,
    ) -> list[dict]:
        tags = []
        for aflow_reduced_structure, dct in self._prototype_db:
            match = structure_matcher.fit_anonymous(
                aflow_reduced_structure, reduced_structure, skip_structure_reduction=True
            )
            if match:
                tags.append(dct)
        return tags

    def _match_single_prototype(self, structure: Structure) -> list[dict]:
        sm = StructureMatcher(
            ltol=self.initial_ltol,
            stol=self.initial_stol,
            angle_tol=self.initial_angle_tol,
            primitive_cell=True,
        )
        reduced_structure = self._preprocess_structure(structure)
        tags = self._match_prototype(sm, reduced_structure)
        while len(tags) > 1:
            sm.ltol *= 0.8
            sm.stol *= 0.8
            sm.angle_tol *= 0.8
            tags = self._match_prototype(sm, reduced_structure)
            if sm.ltol < 0.01:
                break
        return tags

    def get_prototypes(self, structure: Structure) -> list[dict] | None:
        """Get prototype(s) structures for a given input structure.

        Args:
            structure (Structure): structure to match

        Returns:
            list[dict] | None: A list of dicts containing matched prototype data.
                This should be a list containing just a single entry, but it is
                possible a material can match multiple prototypes.
        """
        tags: list[dict] = self._match_single_prototype(structure)

        return tags or None


@deprecated(
    PrototypeDatabaseMatcher,
    "AflowPrototypeMatcher is deprecated. Use PrototypeDatabaseMatcher instead.",
    category=DeprecationWarning,
    deadline=(2026, 11, 15),
)
@due.dcite(
    Doi("10.1016/j.commatsci.2017.01.017"),
    description="The AFLOW library of crystallographic prototypes: part 1.",
)
class AflowPrototypeMatcher(PrototypeDatabaseMatcher):
    """Deprecated alias for :class:`PrototypeDatabaseMatcher`."""

    def __init__(self, *args, **kwargs) -> None:
        if "prototype_db" in kwargs:
            raise TypeError("AflowPrototypeMatcher uses the built-in AFLOW prototype database.")
        super().__init__(*args, prototype_db=AFLOW_PROTOTYPE_LIBRARY, **kwargs)
