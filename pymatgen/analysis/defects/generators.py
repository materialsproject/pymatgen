# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals
"""
This module defines classes to generate point defect structures
"""

__author__ = "Danny Broberg, Shyam Dwaraknath, Bharat Medasani, Nils E. R. Zimmermann, Geoffroy Hautier"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyam Dwaraknath"
__email__ = "shyamd@lbl.gov"
__status__ = "Development"
__date__ = "Mar 15, 2018"

import six
import logging
from abc import ABCMeta, abstractmethod

from monty.json import MSONable

from pymatgen.analysis.bond_valence import BVAnalyzer
from pymatgen.analysis.defects.core import Vacancy, Interstitial, Substitution
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

logger = logging.getLogger(__name__)


class DefectGenerator(six.with_metaclass(ABCMeta, MSONable)):
    """
    Abstract class for point defects
    Implements generator pattern
    """

    def __iter__(self):
        """
        Return self as this should be an iterator
        """
        return self

    @abstractmethod
    def __next__(self):
        """
        Abstract method to return defects
        """
        return


class VacancyGenerator(DefectGenerator):
    """
    Simple generator for vacancies based on periodically
    equivalent sites
    """

    def __init__(self, structure, include_bv_charge=False):
        """
        Initializes a Vacancy Generator
        Args:
            structure: pymatgen.core.structure.Structure
        """
        self.structure = structure
        self.include_bv_charge = include_bv_charge

        # Find equivalent site list
        sga = SpacegroupAnalyzer(self.structure)
        self.symm_structure = sga.get_symmetrized_structure()
        self.equiv_site_seq = list(self.symm_structure.equivalent_sites)

        self.struct_valences = None
        if self.include_bv_charge:
            bv = BVAnalyzer()
            self.struct_valences = bv.get_valences(self.structure)

    def __next__(self):
        """
        Returns the next letter in the sequence or
        raises StopIteration
        """
        if len(self.equiv_site_seq) > 0:
            vac_site = self.equiv_site_seq.pop(0)
            charge = 0.0
            if self.struct_valences:
                site_index = self.structure.get_sites_in_sphere(vac_site[0].coords, 0.1, include_index=True)[0][2]
                charge = -1 * self.struct_valences[site_index]

            return Vacancy(self.structure, vac_site[0], multiplicity=len(vac_site), charge=charge)
        else:
            raise StopIteration
