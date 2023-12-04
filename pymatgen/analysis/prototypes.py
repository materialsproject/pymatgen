"""
This module is intended to match crystal structures against known crystallographic "prototype"
structures.

In this module, the AflowPrototypeMatcher uses the AFLOW LIBRARY OF CRYSTALLOGRAPHIC PROTOTYPES.
If using this particular class, please cite their publication appropriately:

Mehl, M. J., Hicks, D., Toher, C., Levy, O., Hanson, R. M., Hart, G., & Curtarolo, S. (2017).
The AFLOW library of crystallographic prototypes: part 1.
Computational Materials Science, 136, S1-S828.
https://doi.org/10.1016/j.commatsci.2017.01.017
"""

from __future__ import annotations

import os
from typing import TYPE_CHECKING

from monty.serialization import loadfn

from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.util.due import Doi, due

if TYPE_CHECKING:
    from pymatgen.core import Structure

module_dir = os.path.dirname(os.path.abspath(__file__))
AFLOW_PROTOTYPE_LIBRARY = loadfn(f"{module_dir}/aflow_prototypes.json")


@due.dcite(
    Doi("10.1016/j.commatsci.2017.01.017"),
    description="The AFLOW library of crystallographic prototypes: part 1.",
)
class AflowPrototypeMatcher:
    """
    This class will match structures to their crystal prototypes, and will
    attempt to group species together to match structures derived from
    prototypes (e.g. an A_xB_1-x_C from a binary prototype), and will
    give these the names the "-like" suffix.

    This class uses data from the AFLOW LIBRARY OF CRYSTALLOGRAPHIC PROTOTYPES.
    If using this class, please cite their publication appropriately:

    Mehl, M. J., Hicks, D., Toher, C., Levy, O., Hanson, R. M., Hart, G., & Curtarolo, S. (2017).
    The AFLOW library of crystallographic prototypes: part 1.
    Computational Materials Science, 136, S1-S828.
    https://doi.org/10.1016/j.commatsci.2017.01.017
    """

    def __init__(self, initial_ltol=0.2, initial_stol=0.3, initial_angle_tol=5):
        """
        Tolerances as defined in StructureMatcher. Tolerances will be
        gradually decreased until only a single match is found (if possible).

        Args:
            initial_ltol: fractional length tolerance
            initial_stol: site tolerance
            initial_angle_tol: angle tolerance
        """
        self.initial_ltol = initial_ltol
        self.initial_stol = initial_stol
        self.initial_angle_tol = initial_angle_tol

    @staticmethod
    def _match_prototype(structure_matcher, structure):
        tags = []
        for d in AFLOW_PROTOTYPE_LIBRARY:
            struct = d["snl"].structure
            match = structure_matcher.fit_anonymous(struct, structure)
            if match:
                tags.append(d)
        return tags

    def _match_single_prototype(self, structure: Structure):
        sm = StructureMatcher(
            ltol=self.initial_ltol,
            stol=self.initial_stol,
            angle_tol=self.initial_angle_tol,
        )
        tags = self._match_prototype(sm, structure)
        while len(tags) > 1:
            sm.ltol *= 0.8
            sm.stol *= 0.8
            sm.angle_tol *= 0.8
            tags = self._match_prototype(sm, structure)
            if sm.ltol < 0.01:
                break
        return tags

    def get_prototypes(self, structure: Structure) -> list | None:
        """
        Get prototype(s) structures for a given input structure. If you use this method in
        your work, please cite the appropriate AFLOW publication:

        Mehl, M. J., Hicks, D., Toher, C., Levy, O., Hanson, R. M., Hart, G., & Curtarolo,
        S. (2017). The AFLOW library of crystallographic prototypes: part 1. Computational
        Materials Science, 136, S1-S828. https://doi.org/10.1016/j.commatsci.2017.01.017

        Args:
            structure: structure to match

        Returns (list): A list of dicts with keys 'snl' for the matched prototype and
            'tags', a dict of tags ('mineral', 'strukturbericht' and 'aflow') of that
            prototype. This should be a list containing just a single entry, but it is
            possible a material can match multiple prototypes.
        """
        tags = self._match_single_prototype(structure)

        return tags or None
