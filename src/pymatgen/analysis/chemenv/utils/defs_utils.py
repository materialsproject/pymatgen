"""This module contains the definition of some objects used in the chemenv package."""

from __future__ import annotations

from typing import TYPE_CHECKING

from pymatgen.analysis.chemenv.utils.coordination_geometry_utils import is_anion_cation_bond
from pymatgen.util.due import Doi, due

if TYPE_CHECKING:
    from typing import ClassVar

    from pymatgen.core import Structure


__author__ = "David Waroquiers"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Geoffroy Hautier"
__version__ = "2.0"
__maintainer__ = "David Waroquiers"
__email__ = "david.waroquiers@gmail.com"
__date__ = "Feb 20, 2016"

"""D. Waroquiers, X. Gonze, G.-M. Rignanese, C. Welker-Nieuwoudt, F. Rosowski,
M. Goebel, S. Schenk, P. Degelmann, R. Andre, R. Glaum, and G. Hautier,
"Statistical analysis of coordination environments in oxides",
Chem. Mater., 2017, 29 (19), pp 8346-8360,
DOI: 10.1021/acs.chemmater.7b02766

D. Waroquiers, J. George, M. Horton, S. Schenk, K. A. Persson, G.-M. Rignanese, X. Gonze, G. Hautier,
"ChemEnv: a fast and robust coordination environment identification tool",
Acta Cryst. B 2020, 76, pp 683-695.
DOI: 10.1107/S2052520620007994
"""

due.cite(
    Doi("10.1021/acs.chemmater.7b02766"),
    description="Statistical analysis of coordination environments in oxides",
)
due.cite(
    Doi("10.1107/S2052520620007994"),
    description="ChemEnv: a fast and robust coordination environment identification tool",
)


class AdditionalConditions:
    """Additional conditions that can be used to filter coordination environments."""

    NO_ADDITIONAL_CONDITION = 0
    ONLY_ANION_CATION_BONDS = 1
    NO_ELEMENT_TO_SAME_ELEMENT_BONDS = 2
    ONLY_ANION_CATION_BONDS_AND_NO_ELEMENT_TO_SAME_ELEMENT_BONDS = 3
    ONLY_ELEMENT_TO_OXYGEN_BONDS = 4
    # Short versions
    NONE = NO_AC = NO_ADDITIONAL_CONDITION
    ONLY_ACB = ONLY_ANION_CATION_BONDS
    NO_E2SEB = NO_ELEMENT_TO_SAME_ELEMENT_BONDS
    ONLY_ACB_AND_NO_E2SEB = ONLY_ANION_CATION_BONDS_AND_NO_ELEMENT_TO_SAME_ELEMENT_BONDS
    ONLY_E2OB = ONLY_ELEMENT_TO_OXYGEN_BONDS
    # Dictionary mapping of integer for the condition and its "description"
    CONDITION_DESCRIPTION: ClassVar = {
        NO_ADDITIONAL_CONDITION: "No additional condition",
        ONLY_ANION_CATION_BONDS: "Only anion-cation bonds",
        NO_ELEMENT_TO_SAME_ELEMENT_BONDS: "No element-element bonds (same elements)",
        ONLY_ANION_CATION_BONDS_AND_NO_ELEMENT_TO_SAME_ELEMENT_BONDS: "Only anion-cation bonds and"
        " no element-element bonds (same elements)",
        ONLY_ELEMENT_TO_OXYGEN_BONDS: "Only element-oxygen bonds",
    }

    ALL = (NONE, ONLY_ACB, NO_E2SEB, ONLY_ACB_AND_NO_E2SEB, ONLY_E2OB)

    def check_condition(self, condition, structure: Structure, parameters):
        """
        Args:
            condition:
            structure:
            parameters:
        """
        if condition == self.NONE:
            return True
        if condition == self.ONLY_ACB:
            valences = parameters["valences"]
            ii = parameters["site_index"]
            jj = parameters["neighbor_index"]
            return is_anion_cation_bond(valences, ii, jj)
        if condition == self.NO_E2SEB:
            ii = parameters["site_index"]
            jj = parameters["neighbor_index"]
            elems_ii = [sp.symbol for sp in structure[ii].species]
            elems_jj = [sp.symbol for sp in structure[jj].species]
            return len(set(elems_ii) & set(elems_jj)) == 0
        if condition == self.ONLY_ACB_AND_NO_E2SEB:
            valences = parameters["valences"]
            ii = parameters["site_index"]
            jj = parameters["neighbor_index"]
            elems_ii = [sp.symbol for sp in structure[ii].species]
            elems_jj = [sp.symbol for sp in structure[jj].species]
            return len(set(elems_ii) & set(elems_jj)) == 0 and is_anion_cation_bond(valences, ii, jj)
        if condition == self.ONLY_E2OB:
            ii = parameters["site_index"]
            jj = parameters["neighbor_index"]
            elems_ii = [sp.symbol for sp in structure[ii].species]
            elems_jj = [sp.symbol for sp in structure[jj].species]
            return ("O" in elems_jj and "O" not in elems_ii) or ("O" in elems_ii and "O" not in elems_jj)
        return None
