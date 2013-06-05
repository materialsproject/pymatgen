#!/usr/bin/env python

"""
This module defines classes for point defects
"""
from __future__ import division
import abc 

from pymatgen.core.sites import PeriodicSite
from pymatgen.symmetry.finder import SymmetryFinder
from pymatgen.io.zeoio import get_voronoi_nodes

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
        To create unit cell with defect pass unit matrix.
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
        Returns sequence of supercells in pymatgen.core.structure.Structure
        format, with each supercell containing unique vacancy
        """
        sc_with_vac = []
        for uniq_defect_site in self.enumerate_uniq_defectsites():
            sc_with_vac.append(self._supercell_with_defect(scaling_matrix, 
                                                           uniq_defect_site))
        return sc_with_vac
             
        
class Interstitial(Defect):
    """
    Subclass of Defect to generate interstitials
    """
    def __init__(self, structure):
        """
        Given a structure, generate symmetrically unique interstitial sites.
        
        Args:
            structure:
                pymatgen.core.Structure
        """
        
        self._structure = structure
        self._uniq_interstitial_sites = []

        #Use Zeo++ to obtain the voronoi nodes. The unique voronoi nodes
        #are possible candidates for interstitial sites
        possible_interstitial_sites = _get_uniq_voronoi_nodes(structure)
        
        #Do futher processing on possibleInterstitialSites to obtain 
        #interstitial sites
        self._uniq_interstitial_sites = possible_interstitial_sites
        
        
    def enumerate_uniq_defectsites(self):
        """
        Enumerate all the defect sites unique w.r.t. symmetry.
        The defect site has "Al" as occupied specie.
        """
        return self._uniq_interstitial_sites

            
    def _supercell_with_defect(self, scaling_matrix, defect_site, element):
        sc = self._structure.copy()
        sc.make_supercell(scaling_matrix)
        oldf_coords = defect_site.frac_coords
        coords = defect_site.lattice.get_cartesian_coords(oldf_coords)
        newf_coords = sc.lattice.get_fractional_coords(coords)
        #sc_defect_site = PeriodicSite(element, newf_coords,
        #                              sc.lattice)
        try:
            sc.append(element, newf_coords, coords_are_cartesian=False, 
                      validate_proximity=True)
        except ValueError:      
            sc = None
        finally:
            return sc
           
    def make_supercells_with_defects(self, scaling_matrix, element):
        """
        Returns sequence of supercells in pymatgen.core.structure.Structure
        format, with each supercell containing unique interstitial
        """
        sc_list_with_interstitial = []
        for defect_site in self.enumerate_uniq_defectsites():
            sc_with_inter = self._supercell_with_defect(
                    scaling_matrix, defect_site, element
                    )
            if sc_with_inter:
                sc_list_with_interstitial.append(sc_with_inter)
        return sc_list_with_interstitial

def _get_uniq_voronoi_nodes(structure):
    """
    Obtain symmetrically unique voronoi nodes using Zeo++ and 
    pymatgen.symmetry.finder.SymmetryFinder
    Args:
        strucutre:
            pymatgen Structure object

    Returns:
        Symmetrically unique voronoi nodes as pymatgen Strucutre
    """
    vor_node_struct = get_voronoi_nodes(structure)
    vor_symmetry_finder = SymmetryFinder(vor_node_struct)
    vor_symm_struct = vor_symmetry_finder.get_symmetrized_structure()
    equiv_sites = vor_symm_struct.equivalent_sites
    uniq_sites = []
    for equiv_site in vor_symm_struct.equivalent_sites:
        uniq_sites.append(equiv_site[0])

    #lat = structure.lattice
    #sp = [site.specie for site in uniq_sites]   # "Al" because to Zeo++
    #coords = [site.coords for site in uniq_sites]
    #vor_node_radii = [site.properties['voronoi_radius'] for site in uniq_sites]
    #uniq_vor_node_struct = Structure(lat, sp, coords, 
    #        coords_are_cartesian=True, 
    #        site_properties={'voronoi_radius':vor_node_radii}
    #        )
    return uniq_sites
    
    
