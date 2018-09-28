# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import six
import logging
from abc import ABCMeta, abstractmethod

from monty.json import MSONable

from pymatgen.core import PeriodicSite
from pymatgen.analysis.bond_valence import BVAnalyzer
from pymatgen.analysis.defects.core import Vacancy, Interstitial, Substitution
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.defects.utils import StructureMotifInterstitial, TopographyAnalyzer

__author__ = "Danny Broberg, Shyam Dwaraknath"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyam Dwaraknath"
__email__ = "shyamd@lbl.gov"
__status__ = "Development"
__date__ = "Mar 15, 2018"
"""
This module defines classes to generate point defect structures
"""

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
            structure(Structure): pymatgen structure object
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
        Returns the next vacancy in the sequence or
        raises StopIteration
        """
        if len(self.equiv_site_seq) > 0:
            vac_site = self.equiv_site_seq.pop(0)
            charge = 0.0
            if self.struct_valences:
                site_index = self.structure.get_sites_in_sphere(vac_site[0].coords, 0.1, include_index=True)[0][2]
                charge = -1 * self.struct_valences[site_index]

            return Vacancy(self.structure, vac_site[0], charge=charge)
        else:
            raise StopIteration


class SubstitutionGenerator(DefectGenerator):
    """
    Simple generator for substitution based on periodically
    equivalent sites
    """

    def __init__(self, structure, element):
        """
        Initializes a Substitution Generator
        note: an Antisite is considered a type of substitution
        Args:
            structure(Structure): pymatgen structure object
            element (str or Element or Specie): element for the substitution
        """
        self.structure = structure
        self.element = element

        # Find equivalent site list
        sga = SpacegroupAnalyzer(self.structure)
        self.symm_structure = sga.get_symmetrized_structure()

        self.equiv_sub = []
        for equiv_site_set in list(self.symm_structure.equivalent_sites):
            vac_site = equiv_site_set[0]
            if isinstance(element, str):  # make sure you compare with specie symbol or Element type
                vac_specie = vac_site.specie.symbol
            else:
                vac_specie = vac_site.specie
            if element != vac_specie:
                defect_site = PeriodicSite(element, vac_site.coords, structure.lattice, coords_are_cartesian=True)
                sub = Substitution(structure, defect_site)
                self.equiv_sub.append(sub)

    def __next__(self):
        """
        Returns the next Substitution in the sequence or
        raises StopIteration
        """
        if len(self.equiv_sub) > 0:
            return self.equiv_sub.pop(0)
        else:
            raise StopIteration


class InterstitialGenerator(DefectGenerator):
    """
    Generator for interstitials at positions
    where the interstitialcy is coordinated by nearest neighbors
    in a way that resembles basic structure motifs
    (e.g., tetrahedra, octahedra).  The algorithm is called InFiT
    (Interstitialcy Finding Tool), it was introducted by
    Nils E. R. Zimmermann, Matthew K. Horton, Anubhav Jain,
    and Maciej Haranczyk (Front. Mater., 4, 34, 2017),
    and it is used by the Python Charged Defect Toolkit
    (PyCDT: D. Broberg et al., Comput. Phys. Commun., in press, 2018).
    """

    def __init__(self, structure, element):
        """
        Initializes an Interstitial generator using structure motifs
        Args:
            structure (Structure): pymatgen structure object
            element (str or Element or Specie): element for the interstitial
        """
        self.structure = structure
        self.element = element
        interstitial_finder = StructureMotifInterstitial(self.structure, self.element)
        self.defect_sites = list(interstitial_finder.enumerate_defectsites())

        # for multiplicity, neccessary to get prim_structure
        spa = SpacegroupAnalyzer(self.structure, symprec=1e-2)
        prim_struct = spa.get_primitive_standard_structure()
        conv_prim_rat = int(self.structure.num_sites / prim_struct.num_sites)
        self.multiplicities = [
            int(interstitial_finder.get_defectsite_multiplicity(def_ind) / conv_prim_rat)
            for def_ind in range(len(self.defect_sites))
        ]
        self.count_def = 0  # for counting the index of the generated defect

    def __next__(self):
        """
        Returns the next interstitial or
        raises StopIteration
        """
        if len(self.defect_sites) > 0:
            int_site = self.defect_sites.pop(0)
            mult = self.multiplicities.pop(0)
            self.count_def += 1
            site_name = 'InFiT' + str(self.count_def)
            return Interstitial(self.structure, int_site, site_name=site_name, multiplicity=mult)
        else:
            raise StopIteration


class VoronoiInterstitialGenerator(DefectGenerator):
    """
    Generator for interstitials based on a simple Voronoi analysis
    """

    def __init__(self, structure, element):
        """
        Initializes an Interstitial generator using Voronoi sites
        Args:
            structure (Structure): pymatgen structure object
            element (str or Element or Specie): element for the interstitial
        """
        self.structure = structure
        self.element = element

        framework = list(self.structure.symbol_set)
        get_voronoi = TopographyAnalyzer(self.structure, framework, [], check_volume=False)
        get_voronoi.cluster_nodes()
        get_voronoi.remove_collisions()

        # trim equivalent nodes with symmetry analysis
        struct_to_trim = self.structure.copy()
        for poss_inter in get_voronoi.vnodes:
            struct_to_trim.append(self.element, poss_inter.frac_coords, coords_are_cartesian=False)

        symmetry_finder = SpacegroupAnalyzer(struct_to_trim, symprec=1e-1)
        equiv_sites_list = symmetry_finder.get_symmetrized_structure().equivalent_sites

        self.equiv_site_seq = []
        for poss_site_list in equiv_sites_list:
            if poss_site_list[0] not in self.structure:
                self.equiv_site_seq.append(poss_site_list)

        self.count_def = 0  # for counting the index of the generated defect

    def __next__(self):
        """
        Returns the next interstitial or
        raises StopIteration
        """
        if len(self.equiv_site_seq) > 0:
            inter_site_list = self.equiv_site_seq.pop(0)
            self.count_def += 1
            site_name = 'Voronoi' + str(self.count_def)
            return Interstitial(
                self.structure, inter_site_list[0], site_name=site_name, multiplicity=len(inter_site_list))
        else:
            raise StopIteration


class SimpleChargeGenerator(DefectGenerator):
    """
    Does an extremely simple/limited charge generation scheme (only one charge generated)

    for vacancies: use bond valence method to assign oxidation states and consider
                    negative of the vacant site's oxidation state as single charge to try
    for antisites and subs: use bond valence method to assign oxidation states and consider
                    negative of the vacant site's oxidation state as single charge to try +
                    added to likely charge of substitutional site (closest to zero)
    for interstitial: charge zero
    """

    def __init__(self, defect):
        """
        Args:
            defect(Defect): pymatgen Defect object
        """
        self.defect = defect

        try:
            bv = BVAnalyzer()
            struct_valences = bv.get_valences(self.defect.bulk_structure)
            site_index = self.defect.bulk_structure.get_sites_in_sphere(
                self.defect.site.coords, 0.1, include_index=True)[0][2]
            def_site_valence = struct_valences[site_index]
        except Exception:  # sometimes valences cant be assigned
            def_site_valence = 0

        if isinstance(defect, Vacancy):
            self.charges = [-1 * def_site_valence]
        elif isinstance(defect, Substitution):
            #(minimize difference with host site specie)
            probable_chgs = [ox - def_site_valence for ox in self.defect.site.specie.oxidation_states]
            self.charges = [min(probable_chgs, key=abs)]
        elif isinstance(defect, Interstitial):
            self.charges = [0]
        else:
            raise ValueError("Defect Type not recognized.")

    def __next__(self):
        """
        Returns the next defect type with the correct charge appended
        raises StopIteration
        """
        if len(self.charges) > 0:
            charge = self.charges.pop(0)
            defect = self.defect.copy()
            defect.set_charge(charge)
            return defect
        else:
            raise StopIteration
