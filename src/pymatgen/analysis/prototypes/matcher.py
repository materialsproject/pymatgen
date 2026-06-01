"""In this module, the PrototypeDatabaseMatcher defaults to the AFLOW LIBRARY OF
CRYSTALLOGRAPHIC PROTOTYPES. If using the default library, please cite their
publication appropriately:

Mehl, M. J., Hicks, D., Toher, C., Levy, O., Hanson, R. M., Hart, G., & Curtarolo, S. (2017).
The AFLOW library of crystallographic prototypes: part 1.
Computational Materials Science, 136, S1-S828.
https://doi.org/10.1016/j.commatsci.2017.01.017
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import pandas as pd
from monty.dev import deprecated

from pymatgen.analysis.prototypes._data import AFLOW_PROTOTYPE_LIBRARY
from pymatgen.core.structure_matcher import StructureMatcher
from pymatgen.util.due import Doi, due

if TYPE_CHECKING:
    from typing import Any

    from pymatgen.core import Structure


class PrototypeDatabaseMatcher:
    """Match structures against a database of crystal prototypes using anonymous structure matching.

    Each prototype structure is preprocessed into a Niggli-reduced primitive cell at construction
    time. To match a query structure, it is reduced in the same way and then compared against all
    prototypes via :meth:`StructureMatcher.fit_anonymous`, which matches without regard to species
    identity. If multiple prototypes match, the tolerances are tightened iteratively (by 0.8x per
    step) until a single match is found or the length tolerance falls below 0.01.
    """

    def __init__(
        self,
        prototype_db: pd.DataFrame,
        initial_ltol: float = 0.2,
        initial_stol: float = 0.3,
        initial_angle_tol: float = 5,
    ) -> None:
        """
        Tolerances as defined in StructureMatcher. Tolerances will be
        gradually decreased until only a single match is found (if possible).

        Args:
            prototype_db: DataFrame of prototype entries. Rows must contain "snl".
            initial_ltol (float): fractional length tolerance.
            initial_stol (float): site tolerance.
            initial_angle_tol (float): angle tolerance.
        """
        self.initial_ltol = initial_ltol
        self.initial_stol = initial_stol
        self.initial_angle_tol = initial_angle_tol

        self._prototype_db: list[tuple[Structure, dict[str, Any]]] = []
        for _, row in prototype_db.iterrows():
            structure = self._get_structure(row)
            reduced_structure = self._preprocess_structure(structure)
            self._prototype_db.append((reduced_structure, self._get_entry_data(row)))

    @staticmethod
    def _get_structure(row: pd.Series) -> Structure:
        """Override this method to parse a pmg structure for a given row."""
        return row["snl"].structure

    @staticmethod
    def _get_entry_data(row: pd.Series) -> dict[str, Any]:
        """Override this method to parse tags/metadata to return for a given row."""
        return dict(row)

    @staticmethod
    def _preprocess_structure(structure: Structure) -> Structure:
        return structure.get_reduced_structure(reduction_algo="niggli").get_primitive_structure()

    def _match_prototype(
        self,
        structure_matcher: StructureMatcher,
        reduced_structure: Structure,
    ) -> list[dict[str, Any]]:
        tags = []
        for aflow_reduced_structure, dct in self._prototype_db:
            match = structure_matcher.fit_anonymous(
                aflow_reduced_structure, reduced_structure, skip_structure_reduction=True
            )
            if match:
                tags.append(dct)
        return tags

    def _match_single_prototype(self, structure: Structure) -> list[dict[str, Any]]:
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

    def get_prototypes(self, structure: Structure) -> list[dict[str, Any]] | None:
        """Get prototype(s) structures for a given input structure.

        Args:
            structure (Structure): structure to match

        Returns:
            A list of dicts containing matched prototype data. This should be a
            list containing just a single entry, but it is possible a material
            can match multiple prototypes. The schema of each dict follows
            :meth:`_get_entry_data`.
        """
        tags = self._match_single_prototype(structure)

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
    """Deprecated alias for :class:`PrototypeDatabaseMatcher`.

    This class uses data from the AFLOW LIBRARY OF CRYSTALLOGRAPHIC
    PROTOTYPES. If using this class please cite their publication
    appropriately:

    Mehl, M. J., Hicks, D., Toher, C., Levy, O., Hanson, R. M., Hart, G., & Curtarolo, S. (2017).
    The AFLOW library of crystallographic prototypes: part 1.
    Computational Materials Science, 136, S1-S828.
    https://doi.org/10.1016/j.commatsci.2017.01.017
    """

    def __init__(self, initial_ltol: float = 0.2, initial_stol: float = 0.3, initial_angle_tol: float = 5) -> None:
        super().__init__(
            pd.DataFrame(AFLOW_PROTOTYPE_LIBRARY),
            initial_ltol=initial_ltol,
            initial_stol=initial_stol,
            initial_angle_tol=initial_angle_tol,
        )
