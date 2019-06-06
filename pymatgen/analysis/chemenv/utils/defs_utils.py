# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


"""
This module contains the definition of some objects used in the chemenv package.
"""

__author__ = "David Waroquiers"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Geoffroy Hautier"
__version__ = "2.0"
__maintainer__ = "David Waroquiers"
__email__ = "david.waroquiers@gmail.com"
__date__ = "Feb 20, 2016"


from pymatgen.analysis.chemenv.utils.coordination_geometry_utils import is_anion_cation_bond


STATS_ENV_PAPER = 'David Waroquiers, Xavier Gonze, Gian-Marco Rignanese, Cathrin Welker-Nieuwoudt, Frank Rosowski,\n' \
                  'Michael Goebel, Stephan Schenk, Peter Degelmann, Rute Andre, Robert Glaum, and Geoffroy Hautier,\n' \
                  '"Statistical analysis of coordination environments in oxides",\n' \
                  'Chem. Mater., 2017, 29 (19), pp 8346-8360,\n' \
                  'DOI: 10.1021/acs.chemmater.7b02766\n'



def chemenv_citations():
    out = ''
    out += '\nIf you use the ChemEnv tool for your research, please consider citing the following reference(s) :\n'
    out += '==================================================================================================\n'
    out += STATS_ENV_PAPER
    return out


class AdditionalConditions():
    NO_ADDITIONAL_CONDITION = 0
    ONLY_ANION_CATION_BONDS = 1
    NO_ELEMENT_TO_SAME_ELEMENT_BONDS = 2
    ONLY_ANION_CATION_BONDS_AND_NO_ELEMENT_TO_SAME_ELEMENT_BONDS = 3
    ONLY_ELEMENT_TO_OXYGEN_BONDS = 4
    #Short versions
    NONE = NO_ADDITIONAL_CONDITION
    NO_AC = NO_ADDITIONAL_CONDITION
    ONLY_ACB = ONLY_ANION_CATION_BONDS
    NO_E2SEB = NO_ELEMENT_TO_SAME_ELEMENT_BONDS
    ONLY_ACB_AND_NO_E2SEB = ONLY_ANION_CATION_BONDS_AND_NO_ELEMENT_TO_SAME_ELEMENT_BONDS
    ONLY_E2OB = ONLY_ELEMENT_TO_OXYGEN_BONDS
    #Dictionary mapping of integer for the condition and its "description"
    CONDITION_DESCRIPTION = {NO_ADDITIONAL_CONDITION: 'No additional condition',
                             ONLY_ANION_CATION_BONDS: 'Only anion-cation bonds',
                             NO_ELEMENT_TO_SAME_ELEMENT_BONDS: 'No element-element bonds (same elements)',
                             ONLY_ANION_CATION_BONDS_AND_NO_ELEMENT_TO_SAME_ELEMENT_BONDS: 'Only anion-cation bonds and'
                                                                                           ' no element-element bonds'
                                                                                           ' (same elements)',
                             ONLY_ELEMENT_TO_OXYGEN_BONDS: 'Only element-oxygen bonds'}

    ALL = [NONE, ONLY_ACB, NO_E2SEB, ONLY_ACB_AND_NO_E2SEB, ONLY_E2OB]

    def check_condition(self, condition, structure, parameters):
        if condition == self.NONE:
            return True
        elif condition == self.ONLY_ACB:
            valences = parameters['valences']
            ii = parameters['site_index']
            jj = parameters['neighbor_index']
            return is_anion_cation_bond(valences, ii, jj)
        elif condition == self.NO_E2SEB:
            ii = parameters['site_index']
            jj = parameters['neighbor_index']
            elmts_ii = [sp.symbol for sp in structure[ii].species]
            elmts_jj = [sp.symbol for sp in structure[jj].species]
            return len(set(elmts_ii) & set(elmts_jj)) == 0
        elif condition == self.ONLY_ACB_AND_NO_E2SEB:
            valences = parameters['valences']
            ii = parameters['site_index']
            jj = parameters['neighbor_index']
            elmts_ii = [sp.symbol for sp in structure[ii].species]
            elmts_jj = [sp.symbol for sp in structure[jj].species]
            return len(set(elmts_ii) & set(elmts_jj)) == 0 and is_anion_cation_bond(valences, ii, jj)
        elif condition == self.ONLY_E2OB:
            ii = parameters['site_index']
            jj = parameters['neighbor_index']
            elmts_ii = [sp.symbol for sp in structure[ii].species]
            elmts_jj = [sp.symbol for sp in structure[jj].species]
            return ('O' in elmts_jj and 'O' not in elmts_ii) or ('O' in elmts_ii and 'O' not in elmts_jj)
