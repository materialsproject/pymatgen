# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals
"""
This module defines classes to define point defect objects
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
from monty.functools import lru_cache

from pymatgen.core.composition import Composition

logger = logging.getLogger(__name__)


class Defect(six.with_metaclass(ABCMeta, MSONable)):
    """
    Abstract class for a single point defect
    """

    def __init__(self, structure, defect_site, charge=0., multiplicity=1):
        """
        Initializes an abstract defect

        Args:
            structure: Pymatgen Structure without any defects
            charge: (int or float) defect charge
                default is zero, meaning no change to NELECT after defect is created in the structure
        """
        self._structure = structure
        self._charge = charge
        self._defect_site = defect_site
        self._multiplicity = multiplicity

    @property
    def structure(self):
        """
        Returns the structure without any defects.
        """
        return self._structure

    @property
    def charge(self):
        """
        Returns the charge of a defect
        """
        return self._charge

    @property
    def defect_site(self):
        """
        Returns the defect position as a site object
        """
        return self._defect_site

    @property
    def multiplicity(self):
        """
        Returns the multiplicity of a defect site within the structure (needed for concentration analysis)
        """
        return self._multiplicity

    @property
    @abstractmethod
    def defect_composition(self):
        """
        Returns the defect composition as a Composition object
        """
        return

    @abstractmethod
    def generate_defect_structure(self, supercell=(1, 1, 1)):
        """
        Given structure and defect_site (and type of defect) should return a defect_structure that is charged
        Args:
            supercell (int, [3x1], or [[]] (3x3)): supercell integer, vector, or scaling matrix
        """
        return


class Vacancy(Defect):
    """
    Subclass of Defect to capture essential information for a single Vacancy defect structure.
    """

    @property
    def defect_composition(self):
        temp_comp = self.structure.composition.as_dict()
        temp_comp[str(self.defect_site.specie)] -= 1
        return Composition(temp_comp)

    def generate_defect_structure(self, supercell=(1, 1, 1)):
        """
        Returns Defective Vacancy structure, decorated with charge
        Args:
            supercell (int, [3x1], or [[]] (3x3)): supercell integer, vector, or scaling matrix
        """
        defect_structure = self.structure.copy()
        defect_structure.make_supercell(supercell)
        poss_deflist = sorted(
            defect_structure.get_sites_in_sphere(self.defect_site.coords, 2, include_index=True), key=lambda x: x[1])
        defindex = poss_deflist[0][2]
        defect_structure.remove_sites([defindex])
        defect_structure.set_charge(self.charge)
        return defect_structure


class Substitution(Defect):
    """
    Subclass of Defect to capture essential information for a single Substitution defect structure.
    """

    @property
    @lru_cache(1)
    def defect_composition(self):
        poss_deflist = sorted(
            self.structure.get_sites_in_sphere(self.defect_site.coords, 2, include_index=True), key=lambda x: x[1])
        defindex = poss_deflist[0][2]

        temp_comp = self.structure.composition.as_dict()
        temp_comp[str(self.defect_site.specie)] += 1
        temp_comp[str(self.structure[defindex].specie)] -= 1
        return Composition(temp_comp)

    def generate_defect_structure(self, supercell=(1, 1, 1)):
        """
        Returns Defective Substitution structure, decorated with charge
        Args:
            supercell (int, [3x1], or [[]] (3x3)): supercell integer, vector, or scaling matrix
        """
        defect_structure = self.structure.copy()
        defect_structure.make_supercell(supercell)
        poss_deflist = sorted(
            defect_structure.get_sites_in_sphere(self.defect_site.coords, 2, include_index=True), key=lambda x: x[1])
        defindex = poss_deflist[0][2]

        subsite = defect_structure.pop(defindex)
        defect_structure.append(self.defect_site.specie.symbol, subsite.coords, coords_are_cartesian=True)
        defect_structure.set_charge(self.charge)
        return defect_structure


class Interstitial(Defect):
    """
    Subclass of Defect to capture essential information for a single Interstitial defect structure.
    """

    @property
    def defect_composition(self):
        temp_comp = self.structure.composition.as_dict()
        temp_comp[str(self.defect_site.specie)] += 1
        return Composition(temp_comp)

    def generate_defect_structure(self, supercell=(1, 1, 1)):
        """
        Returns Defective Interstitial structure, decorated with charge
        Args:
            supercell (int, [3x1], or [[]] (3x3)): supercell integer, vector, or scaling matrix
        """
        defect_structure = self.structure.copy()
        defect_structure.make_supercell(supercell)
        defect_structure.append(self.defect_site.specie.symbol, self.defect_site.coords, coords_are_cartesian=True)
        defect_structure.set_charge(self.charge)
        return defect_structure