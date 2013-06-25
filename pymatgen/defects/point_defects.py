#!/usr/bin/env python

"""
This module defines classes for point defects
"""
from __future__ import division
import abc 
from sets import Set
import sys

from pymatgen.core.sites import PeriodicSite
from pymatgen.symmetry.finder import SymmetryFinder
from pymatgen.io.zeoio import get_voronoi_nodes, get_void_volume_surfarea
from pymatgen.analysis.structure_analyzer import VoronoiCoordFinder
from pymatgen.analysis.bond_valence import BVAnalyzer
from pymatgen.core import Specie

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

        bv = BVAnalyzer()
        try:
            valences = bv.get_valences(self._structure)
        except:
            valences = [0]*len(self._structure.sites)
            err_str = "BVAnalyzer failed. The defect effective charge, and"
            err_str += " volume and surface area may not work"
            print err_str
            raise ValueError()
        # Store the ionic radii for the elements in the structure
        # (Used to  computing the surface are and volume)
        # Computed based on valence of each element
        self._valence_dict = {}
        sites = self._structure.sites
        for i in range(len(sites)):
            if sites[i].species_string in self._valence_dict.keys():
                continue
            else:
                self._valence_dict[sites[i].species_string] = valences[i]
                if len(self._valence_dict) == len(self._structure.composition):
                    break
        #print self._valence_dict

        self._rad_dict = {}
        for el in self._valence_dict.keys():
            val = self._valence_dict[el]
            self._rad_dict[el] = Specie(el, val).ionic_radius
            if not self._rad_dict[el]: #get covalent radii
                pass
        #print self._rad_dict
        assert len(self._rad_dict) == len(self._structure.composition)

        self._vac_eff_charges = None
        self._vol = None
        self._sa = None

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
        return self._vac_sites[n]

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
        Coordination number of vacancy site.
        Args:
            n
                Index of vacancy list
        """
        return self._vac_coordination_numbers[n]

    def get_coordinated_sites(self, n):
        """
        The sites surrounding the vacancy site.
        Args:
            n
                Index of vacancy list
        """
        return self._vac_coordinated_sites[n]

    def get_coordinated_elements(self, n):
        """
        Elements of sites surrounding the vacancy site.
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
        # Effective charge (In Kroger-Vink notation, cation vacancy has 
        # effectively -ve charge and anion vacancy has +ve charge.) Inverse
        # the BVAnalyzer.get_valences result.
        if not self._vac_eff_charges:
            self._vac_eff_charges = []
            for site in self.enumerate_defectsites():
                specie = site.species_string
                self._vac_eff_charges.append(-self._valence_dict[specie])
        return self._vac_eff_charges[n]

    def get_coordsites_min_max_charge(self,n):
        """
        Minimum and maximum charge of sites surrounding the vacancy site.
        Args:
            n
                Index of vacancy list
        """
        bv = BVAnalyzer()
        struct_valences = bv.get_valences(self._structure)
        coordinated_site_valences = []

        def _get_index(site):
            for i in range(len(self._structure.sites)):
                if site.is_periodic_image(self._structure.sites[i]):
                    return i
            raise ValueError("Site not found")

        for site in self._vac_coordinated_sites[n]:
            ind = _get_index(site)
            coordinated_site_valences.append(struct_valences[ind])
        coordinated_site_valences.sort()
        return coordinated_site_valences[0], coordinated_site_valences[-1]

    def get_volume(self, n):
        """
        Volume of the nth vacancy

        Args:
            n
                Index of symmetrically distinct vacancies list
        Returns:
            floating number representing volume of vacancy
        """
        if not self._vol:
            self._vol = []
            self._sa = []
            um = [[1,0,0],[0,1,0],[0,0,1]]
            sc = self.make_supercells_with_defects(um)[1:]
            rad_dict = self.ionic_radii
            for i in range(len(sc)):
                site_radi = rad_dict[self._vac_sites[i].species_string]
                vol,sa = get_void_volume_surfarea(sc[i], rad_dict)
                self._vol.append(vol)
                self._sa.append(sa)

        return self._vol[n]


    def get_surface_area(self, n):
        """
        Surface area of the nth vacancy

        Args:
            n
                Index of symmetrically distinct vacancies list
        Returns:
            floating number representing volume of vacancy
        """
        if not self._sa:
            self._vol = []
            self._sa = []
            um = [[1,0,0],[0,1,0],[0,0,1]]
            supercells = self.make_supercells_with_defects(um)[1:]
            rad_dict = self.ionic_radii
            for sc in supercells:
                vol, sa = get_void_volume_surfarea(sc, rad_dict)
                self._vol.append(vol)
                self._sa.append(sa)

        return self._sa[n]

    @property
    def structure(self):
        """
        Structure without any defects
        """
        return self._structure

    @property
    def ionic_radii(self):
        """
        Ionic radii of elements in the structure
        """
        return self._rad_dict

    @property
    def valences(self):
        """
        Valence of elements in the structure
        """
        return self._valence_dict

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

        self._interstitial_coord_no = []
        self._interstitial_coord_sites = []
        self._valence_dict = {}
        self._rad_dict = {}
        self._radius = []

        for site in self._interstitial_sites:
            struct = self._structure.copy()
            struct.append(site.species_string, site.frac_coords)
            coord_finder = VoronoiCoordFinder(struct)
            self._interstitial_coord_no.append(
                    coord_finder.get_coordination_number(-1))
            self._interstitial_coord_sites.append(
                    coord_finder.get_coordinated_sites(-1))

        bv = BVAnalyzer()
        try:
            valences = bv.get_valences(self._structure)
        except:
            valences = [0]*len(self._structure.sites)
            err_str = "BVAnalyzer failed. The defect effective charge, and"
            err_str += " volume and surface area may not work"
            print err_str
            raise ValueError()

        sites = self._structure.sites
        for i in range(len(sites)):
            if sites[i].species_string in self._valence_dict.keys():
                continue
            else:
                self._valence_dict[sites[i].species_string] = valences[i]
                if len(self._valence_dict) == len(self._structure.composition):
                    break

        self._rad_dict = {}
        for el in self._valence_dict.keys():
            val = self._valence_dict[el]
            self._rad_dict[el] = Specie(el, val).ionic_radius
            if not self._rad_dict[el]: #get covalent radii
                pass
        assert len(self._rad_dict) == len(self._structure.composition)

        for site in self._interstitial_sites:
            self._radius.append(float(site.properties['voronoi_radius']))
    
    def defectsite_count(self):
        """
        Returns the number of symmetrically distinct interstitial sites
        """
        return len(self._interstitial_sites)
    
    def enumerate_defectsites(self):
        """
        Enumerate all the symmetrically distinct interstitial sites.
        The defect site has "Al" as occupied specie.
        """
        return self._interstitial_sites

    def get_defectsite(self,n):
        """
        Returns the Site position of the interstitial.
        """
        return self._interstitial_sites[n]

    def get_defectsite_coordination_number(self, n):
        """
        Coordination number of interstitial site.
        Args:
            n
                Index of interstitial list
        """
        return self._interstitial_coord_no[n]

    def get_coordinated_sites(self, n):
        """
        The sites surrounding the interstitial site.
        Args:
            n
                Index of interstitial list
        """
        return self._interstitial_coord_sites[n]

    def get_coordinated_elements(self, n=None):
        """
        Elements of sites surrounding the interstitial site.
        Args:
            n
                Index of interstitial list
        """
        coordinated_species = []
        for site in self._interstitial_coord_sites[n]:
            coordinated_species.append(site.species_string)
        return set(coordinated_species)

    def get_coordsites_charge_sum(self,n):
        """
        Minimum and maximum charge of sites surrounding the interstitial site.
        Args:
            n
                Index of interstitial list
        """
        coordsite_valence_sum = 0.0

        for site in self._interstitial_coord_sites[n]:
            coordsite_valence_sum += self._valence_dict[site.species_string]
        return coordsite_valence_sum

    def get_coordsites_min_max_charge(self,n):
        """
        Minimum and maximum charge of sites surrounding the interstitial site.
        Args:
            n
                Index of interstitial list
        """
        coord_site_valences = []

        for site in self._vac_coordinated_sites[n]:
            coord_site_valences.append(self._valence_dict[site.species_string])
        coord_site_valences.sort()
        return coord_site_valences[0], coord_site_valences[-1]

    def get_radius(self, n):
        """
        Volume of the nth interstitial

        Args:
            n
                Index of symmetrically distinct vacancies list
        Returns:
            floating number representing radius of interstitial sphere
        """
        return self._radius[n]

    def get_surface_area(self, n):
        """
        Surface area of the nth interstitial

        Args:
            n
                Index of symmetrically distinct vacancies list
        Returns:
            floating number representing volume of interstitial
        """
        pass

    @property
    def structure(self):
        """
        Returns the structure without any defects
        Useful for Mott-Littleton calculations.
        """
        return self._structure

    @property
    def ionic_radii(self):
        """
        Ionic radii of elements in the structure
        """
        return self._rad_dict
    
    @property
    def valences(self):
        """
        Valence of elements in the structure
        """
        return self._valence_dict
            
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
    
    
