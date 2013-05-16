#!/usr/bin/env python

"""
Abstract class for defects
"""
from __future__ import division

__author__ = "Bharat K. Medasani"

import abc 

from pymatgen.core.sites import PeriodicSite
from pymatgen.symmetry.finder import SymmetryFinder

#from pymatgen.core.structure import PeriodicSize
#from pymatgen.core.structure_modifier import SuperCellMaker

class Defect:
    """
    Abstract class for point defects
    """
    __metaclass__ = abc.ABCMeta
    
    @abc.abstractmethod
    def enumerate_uniq_defectsites(self):
        """
        Enumerates all the unique defects
        """
        print 'Not implemented'
    
    @abc.abstractmethod
    def make_supercells_with_defects(self, scaling_matrix):
        """
        Generate the supercell with input multipliers and create the defect
        """ 
        print 'Not implemented'
    
    
class Vacancy(Defect):
    """
    Subclass of Defect to generate vacancies
    """
    def __init__(self, inp_structure):
        """
        Given a structure, generate the unique vacancy sites
        
        Args:
            inp_structure:
                pymatgen.core.Structure
        """
        
        self._structure = inp_structure
        symm_finder = SymmetryFinder(self._structure)
        symm_structure = symm_finder.get_symmetrized_structure()
        equiv_site_seq = symm_structure.equivalent_sites
        self._uniq_vac_sites = []
        for equiv_sites in equiv_site_seq:
            self._uniq_vac_sites.append(equiv_sites[0])
        
        
    def enumerate_uniq_defectsites(self):
        """
        Enumerate the unique defect sites
        
        Returns:
            List of unique defect sites

        """
        return self._uniq_vac_sites
    
    def _supercell_with_defect(self, scaling_matrix, defect_site):
        sc = self._structure.copy()
        sc.make_supercell(scaling_matrix)
        oldf_coords = defect_site.frac_coords
        coords = defect_site.lattice.get_cartesian_coords(oldf_coords)
        newf_coords = sc.lattice.get_fractional_coords(coords)
        sc_defect_site = PeriodicSite(defect_site.species_and_occu, newf_coords,
                                      sc.lattice,
                                      properties=defect_site.properties)
        for i in range(len(sc.sites)):
            if sc_defect_site == sc.sites[i]:
                sc.remove(i)
                return sc
        
    def make_supercells_with_defects(self, scaling_matrix):
        """
        Returns sequence of supercells containing unique defects
        """
        sc_with_vac = []
        for uniq_defect_site in self.enumerate():
            sc_with_vac.append(self._supercell_with_defect(scaling_matrix, 
                                                           uniq_defect_site))
        return sc_with_vac
             
        
class Interstitial(Defect):
    """
    Subclass of Defect to generate interstitials
    """
    def __init__(self, inp_structure ):
        """
        Given a structure, generate unique interstitial sites
        
        Args:
            inp_structure:
                pymatgen.core.Structure
        """
        
        self._structure = inp_structure
        self._uniq_interstitial_sites = []
        #symm_finder = SymmetryFinder(self._structure)
        #Use Zeo++ to analyze the void space and identify the void space
    
    def enumerate_uniq_defectsites(self):
        return self._uniq_interstitial_sites

    def add_interstitial(self, site, rad

    def make_supercells_with_defects(self, scaling_matrix, element):
        pass
