#!/usr/bin/env python

"""
This module defines classes for point defects
"""
from __future__ import division
import abc 
from sets import Set

from pymatgen.core.sites import PeriodicSite
from pymatgen.symmetry.finder import SymmetryFinder
from pymatgen.io.zeoio import get_voronoi_nodes
from pymatgen.analysis.structure_analyzer import VoronoiCoordFinder
from pymatgen.analysis.bond_valence import BVAnalyzer

#from pymatgen.core.structure import PeriodicSize
#from pymatgen.core.structure_modifier import SuperCellMaker

class Defect:
    """
    Abstract class for point defects
    """
    __metaclass__ = abc.ABCMeta
    
    @abc.abstractmethod
    def enumerate_defectsites(self):
        """
        Enumerates all the symmetrically distinct defects.
        """
        print 'Not implemented'

    @abc.abstractproperty
    def structure(self):
        """
        Returns the structure without any defects
        Useful for Mott-Littleton calculations.
        """
        print 'Not implemented'
    
    @abc.abstractmethod
    def make_supercells_with_defects(self, scaling_matrix):
        """
        Generate the supercell with input multipliers and create the defect.
        First supercell has no defects.
        To create unit cell with defect pass unit matrix.
        """ 
        print 'Not implemented'
    
    
class Vacancy(Defect):
    """
    Subclass of Defect to generate vacancies and their analysis
    """
    def __init__(self, structure):
        """

        Args:
            inp_structure:
                pymatgen.core.Structure
        """
        
        self._structure = structure
        # Store symmetrically distinct sites, their coordination numbers
        # coordinated_sites, effective charge
        symm_finder = SymmetryFinder(self._structure)
        symm_structure = symm_finder.get_symmetrized_structure()
        equiv_site_seq = symm_structure.equivalent_sites
        self._vac_sites = []
        for equiv_sites in equiv_site_seq:
            self._vac_sites.append(equiv_sites[0])
        self._vac_site_indices = []
        for site in self._vac_sites:
            for i in range(len(self._structure.sites)):
                if site == self._structure[i]:
                    self._vac_site_indices.append(i)
        coord_finder = VoronoiCoordFinder(self._structure)
        self._vac_coordination_numbers = []
        self._vac_coordinated_sites = []
        for i in self._vac_site_indices:
            self._vac_coordination_numbers.append(
                    coord_finder.get_coordination_number(i))
            self._vac_coordinated_sites.append(
                    coord_finder.get_coordinated_sites(i))
        # Effective charge (In Kroger-Vink notation, cation vacancy has 
        # effectively -ve charge and anion vacancy has +ve charge.) Inverse
        # the BVAnalyzer.get_valences result.
        bv = BVAnalyzer()
        valences = bv.get_valences(self._structure)
        self._vac_eff_charges = []
        for i in self._vac_site_indices:
            self._vac_eff_charges.append(-valences[i])


    def defectsite_count(self):
        """
        Returns the number of symmetrically distinct vacancy sites
        """
        return len(self._vac_sites)
        
    def enumerate_defectsites(self):
        """
        Returns symmetrically distinct vacancy sites
        
        """
        return self._vac_sites
    
    def get_defectsite_indices(self):
        """
        Returns indices of symmetrically distinct vacancy sites
        
        """
        return self._vac_site_indices

    def get_defectsite(self,n):
        return self._vac_sites[i]

    def get_defectsite_index(self, n):
        """
        index of the vacacy site in the structure.sites list
        Args:
            n
                Index of vacancy list
        """
        return self._vac_site_indices[n]

    def get_defectsite_coordination_number(self, n):
        """
        Returns the coordination number of vacancy site.
        Args:
            n
                Index of vacancy list
        """
        return self._vac_coordination_numbers[n]

    def get_defectsite_coordinated_sites(self, n):
        """
        Returns the sites surrounding the vacancy site.
        Args:
            n
                Index of vacancy list
        """
        return self._vac_coordinated_sites[n]

    def get_defectsite_coordinated_elements(self, n):
        """
        Returns the elements of sites surrounding the vacancy site.
        Args:
            n
                Index of vacancy list
        """
        coordinated_species = []
        for site in self._vac_coordinated_sites[n]:
            coordinated_species.append(site.species_and_occu)
        return set(coordinated_species)

    def get_defectsite_effective_charge(self, n):
        """
        Effective charge (In Kroger-Vink notation, cation vacancy has 
        effectively -ve charge and anion vacancy has +ve charge.) 
        Args:
            n
                Index of vacancy list

        Returns:
            Effective charnge of defect site
        """
        return self._vac_eff_charges[n]

    @property
    def structure(self):
        """
        Returns the structure without any defects
        """
        return self._structure

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
        Generate sequence of supercells in pymatgen.core.structure.Structure
        format, with each supercell containing one vacancy.

        Args:
            scaling_matrix:
                super cell scale parameters in matrix forms

        Returns:
            Supercells with vacancies. First supercell has no defects.
        """
        sc_with_vac = []
        sc = self._structure.copy()
        sc.make_supercell(scaling_matrix)
        sc_with_vac.append(sc)
        for defect_site in self.enumerate_defectsites():
            sc_with_vac.append(self._supercell_with_defect(scaling_matrix, 
                                                           defect_site))
        return sc_with_vac
             
        
class Interstitial(Defect):
    """
    Subclass of Defect to generate interstitials
    """
    def __init__(self, structure):
        """
        Given a structure, generate symmetrically distinct interstitial sites.
        
        Args:
            structure:
                pymatgen.core.Structure
        """
        
        self._structure = structure
        self._interstitial_sites = []

        #Use Zeo++ to obtain the voronoi nodes. Apply symmetry reduction and 
        #the symmetry reduced voronoi nodes
        #are possible candidates for interstitial sites
        possible_interstitial_sites = _symmetry_reduced_voronoi_nodes(structure)
        
        #Do futher processing on possibleInterstitialSites to obtain 
        #interstitial sites
        self._interstitial_sites = possible_interstitial_sites
        
        
    def enumerate_defectsites(self):
        """
        Enumerate all the symmetrically distinct interstitial sites.
        The defect site has "Al" as occupied specie.
        """
        return self._interstitial_sites
    
    @property
    def structure(self):
        """
        Returns the structure without any defects
        Useful for Mott-Littleton calculations.
        """
        return self._structure
    
            
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
        format, with each supercell containing an interstitial.
        First supercell has no defects.
        """
        sc_list_with_interstitial = []
        sc = self._structure.copy()
        sc.make_supercell(scaling_matrix)
        sc_list_with_interstitial.append(sc)
        for defect_site in self.enumerate_defectsites():
            sc_with_inter = self._supercell_with_defect(
                    scaling_matrix, defect_site, element
                    )
            if sc_with_inter:
                sc_list_with_interstitial.append(sc_with_inter)
        return sc_list_with_interstitial

def _symmetry_reduced_voronoi_nodes(structure):
    """
    Obtain symmetry reduced voronoi nodes using Zeo++ and 
    pymatgen.symmetry.finder.SymmetryFinder
    Args:
        strucutre:
            pymatgen Structure object

    Returns:
        Symmetrically distinct voronoi nodes as pymatgen Strucutre
    """
    vor_node_struct = get_voronoi_nodes(structure)
    vor_symmetry_finder = SymmetryFinder(vor_node_struct)
    vor_symm_struct = vor_symmetry_finder.get_symmetrized_structure()
    equiv_sites = vor_symm_struct.equivalent_sites
    sites = []
    for equiv_site in vor_symm_struct.equivalent_sites:
        sites.append(equiv_site[0])

    #lat = structure.lattice
    #sp = [site.specie for site in sites]   # "Al" because to Zeo++
    #coords = [site.coords for site in sites]
    #vor_node_radii = [site.properties['voronoi_radius'] for site in sites]
    #vor_node_struct = Structure(lat, sp, coords, 
    #        coords_are_cartesian=True, 
    #        site_properties={'voronoi_radius':vor_node_radii}
    #        )
    return sites
    
    
