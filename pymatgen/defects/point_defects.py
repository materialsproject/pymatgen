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
from pymatgen.command_line.gulp_caller import GulpCaller

#from pymatgen.core.structure import PeriodicSize
#from pymatgen.core.structure_modifier import SuperCellMaker

class StructWithValenceIonicRadius:
    """
    Structure with site valences and ionic radii computed
    """
    def __init__(self, structure):
        self._structure = structure
        self._valences = self._get_valences()
        self._ionic_radii = self._get_ionic_radii()

    @property
    def radii(self):
        return self._ionic_radii

    @property
    def valences(self):
        return self._valences

    @property
    def structure(self):
        return self._structure

    def _get_ionic_radii(self):
        rad_dict = {}
        for el in self._valences.keys():
            val = self._valences[el]
            rad_dict[el] = Specie(el, val).ionic_radius
            if not rad_dict[el]: #get covalent radii later
                raise LookupError()
        return rad_dict

    def _get_valences(self):
        bv = BVAnalyzer()
        try:
            valences = bv.get_valences(self._structure)
        except:
            err_str = "BVAnalyzer failed. The defect effective charge, and"
            err_str += " volume and surface area may not work"
            print err_str
            raise LookupError()

        valence_dict = {}
        sites = self._structure.sites
        for i in range(len(sites)):
            if sites[i].species_string in valence_dict.keys():
                continue
            else:
                valence_dict[sites[i].species_string] = valences[i]
                if len(valence_dict) == len(self._structure.composition):
                    break

        return valence_dict

    
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
        raise NotImpementedError()

    @property
    def structure(self):
        """
        Returns the structure without any defects
        Useful for Mott-Littleton calculations.
        """
        return self._structure

    @property
    def struct_radii(self):
        """
        Radii of elements in the structure
        """
        return self._rad_dict
    
    @property
    def struct_valences(self):
        """
        Valence of elements in the structure
        """
        return self._valence_dict

    def defectsite_count(self):
        """
        Returns the number of symmetrically distinct defect sites
        """
        return len(self._defect_sites)

    def get_defectsite(self, n):
        """
        Returns the defect site at the index.
        """
        return self._defect_sites[n]

    def get_defectsite_coordination_number(self, n):
        """
        Coordination number of interstitial site.
        Args:
            n
                Index of interstitial list
        """
        return self._defectsite_coord_no[n]

    def get_coordinated_sites(self, n):
        """
        The sites in structure surrounding the defect site.
        Args:
            n
                Index of defects list
        """
        return self._defect_coord_sites[n]

    def get_coordinated_elements(self, n):
        """
        Elements of sites in structure surrounding the defect site.
        Args:
            n
                Index of defect list
        """
        coordinated_species = []
        for site in self._defect_coord_sites[n]:
            coordinated_species.append(site.species_string)
        return list(set(coordinated_species))

    @abc.abstractmethod
    def make_supercells_with_defects(self, scaling_matrix):
        """
        Generate the supercell with input multipliers and create the defect.
        First supercell has no defects.
        To create unit cell with defect pass unit matrix.
        """ 
        print 'Not implemented'
        raise NotImpementedError()
    
class Vacancy(Defect):
    """
    Subclass of Defect to generate vacancies and their analysis
    """
    def __init__(self, structure_with_val_radii):
        """

        Args:
            structure_with_val_radii:
                pymatgen.defects.point_defects.StructureWithValenceIonicRadius
        """
        
        self._structure = structure_with_val_radii.structure
        self._valence_dict = structure_with_val_radii.valences
        self._rad_dict = structure_with_val_radii.radii
        # Store symmetrically distinct sites, their coordination numbers
        # coordinated_sites, effective charge
        symm_finder = SymmetryFinder(self._structure)
        symm_structure = symm_finder.get_symmetrized_structure()
        equiv_site_seq = symm_structure.equivalent_sites

        self._defect_sites = []
        for equiv_sites in equiv_site_seq:
            self._defect_sites.append(equiv_sites[0])

        self._vac_site_indices = []
        for site in self._defect_sites:
            for i in range(len(self._structure.sites)):
                if site == self._structure[i]:
                    self._vac_site_indices.append(i)

        coord_finder = VoronoiCoordFinder(self._structure)
        self._defectsite_coord_no = []
        self._defect_coord_sites = []
        for i in self._vac_site_indices:
            self._defectsite_coord_no.append(
                    coord_finder.get_coordination_number(i)
                    )
            self._defect_coord_sites.append(
                    coord_finder.get_coordinated_sites(i)
                    )

        # Store the ionic radii for the elements in the structure
        # (Used to  computing the surface are and volume)
        # Computed based on valence of each element

        self._vac_eff_charges = None
        self._vol = None
        self._sa = None

    def enumerate_defectsites(self):
        """
        Returns symmetrically distinct vacancy sites
        
        """
        return self._defect_sites
    
    def get_defectsite_structure_indices(self):
        """
        Returns indices of symmetrically distinct vacancy sites
        
        """
        return self._vac_site_indices

    def get_defectsite_structure_index(self, n):
        """
        index of the vacacy site in the structure.sites list
        Args:
            n
                Index of vacancy list
        """
        return self._vac_site_indices[n]

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

        el = self.get_defectsite(n).species_string
        return -self._valence_dict[el]
        #if not self._vac_eff_charges:
        #    self._vac_eff_charges = []
        #    for site in self.enumerate_defectsites():
        #        specie = site.species_string
        #        self._vac_eff_charges.append(-self._valence_dict[specie])
        #return self._vac_eff_charges[n]

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

        for site in self._defect_coord_sites[n]:
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
            rad_dict = self.struct_radii
            for i in range(len(sc)):
                site_radi = rad_dict[self._defect_sites[i].species_string]
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
            rad_dict = self.struct_radii
            for sc in supercells:
                vol, sa = get_void_volume_surfarea(sc, rad_dict)
                self._vol.append(vol)
                self._sa.append(sa)

        return self._sa[n]

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

class VacancyFormationEnergy:
    """
    
    """
    def __init__(structure):
        struct_with_val_rad = StructWithValenceIonicRadius(structure)
        self._vacancy = Vacancy(struct_with_val_rad)
    def get_energy(self, n):
        """
        Formation Energy for nth symmetrically distinct vacancy.
        GULP is used to obtain the energy.
        """
        #generate defect free structure energy


        #generate defect structure energy




             
class Interstitial(Defect):
    """
    Subclass of Defect to generate interstitials
    """
    def __init__(self, structure_with_val_radii):
        """
        Given a structure, generate symmetrically distinct interstitial sites.
        
        Args:
            structure_with_val_radii:
                pymatgen.defects.point_defects.StructureWithValenceIonicRadius
        """
        
        self._structure = structure_with_val_radii.structure
        self._valence_dict = structure_with_val_radii.valences
        self._rad_dict = structure_with_val_radii.radii

        #Use Zeo++ to obtain the voronoi nodes. Apply symmetry reduction and 
        #the symmetry reduced voronoi nodes
        #are possible candidates for interstitial sites
        try:
            possible_interstitial_sites = _symmetry_reduced_voronoi_nodes(
                self._structure, self._rad_dict
                )
        except:
            raise ValueError("Symmetry_reduced_voronoi_nodes failed")
        
        #Do futher processing on possibleInterstitialSites to obtain 
        #interstitial sites
        self._defect_sites = possible_interstitial_sites 
        self._defectsite_coord_no = []
        self._defect_coord_sites = []
        self._radii = []

        for site in self._defect_sites:
            coord_no, coord_sites = self._get_coord_no_sites(site)
            self._defectsite_coord_no.append(coord_no)
            self._defect_coord_sites.append(coord_sites)

        for site in self._defect_sites:
            self._radii.append(float(site.properties['voronoi_radius']))
    
    def _get_coord_no_sites(self, site):
        struct = self._structure.copy()
        struct.append(site.species_string, site.frac_coords)
        coord_finder = VoronoiCoordFinder(struct)
        coord_no = coord_finder.get_coordination_number(-1)
        coord_sites = coord_finder.get_coordinated_sites(-1)

        # In some cases coordination sites to interstitials include 
        # interstitials also. 
        sites_to_be_deleted = []
        for i in range(len(coord_sites)):
            # In the future replace voronoi node place holder "H" with 
            # some thing else
            if coord_sites[i].species_string == 'H':
                sites_to_be_deleted.append(i)
        sites_to_be_deleted.reverse() # So index won't change in between
        for ind in sites_to_be_deleted:
            del coord_sites[ind]

        return coord_no, coord_sites

    def enumerate_defectsites(self):
        """
        Enumerate all the symmetrically distinct interstitial sites.
        The defect site has "Al" as occupied specie.
        """
        return self._defect_sites

    def append_defectsite(self, site):
        """
        Append a site to list of possible interstitials
        """
        raise NotImplementedError()

    def delete_defectsite(self, n):
        del self._defect_sites[n]

    def get_coordsites_charge_sum(self,n):
        """
        Minimum and maximum charge of sites surrounding the interstitial site.
        Args:
            n
                Index of interstitial list
        """
        coordsite_valence_sum = 0.0

        for site in self._defect_coord_sites[n]:
            try:
                coordsite_valence_sum += self._valence_dict[site.species_string]
            except:
                print len(self._structure.sites)
                print '--------'
                print len(self._defect_coord_sites[n])
                print '--------'
                print "Structure "
                print self._structure
                print '--------'
                print "Interstitial coord sites"
                for site1 in self._defect_coord_sites[n]:
                    print site1.species_string, site1.frac_coords
                print '--------'
                print "Interstitial site"
                print self._defect_sites[n].frac_coords
                print '--------'
                print "Interstitial sites"
                for site1 in self._defect_sites:
                    print site1
                raise ValueError("Interstitial site as coordination site")
        return coordsite_valence_sum

    def get_coordsites_min_max_charge(self,n):
        """
        Minimum and maximum charge of sites surrounding the interstitial site.
        Args:
            n
                Index of interstitial list
        """
        coord_site_valences = []

        for site in self._defect_coord_sites[n]:
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
        return self._radii[n]

    def get_radii(self):
        return self._radii

    def reduce_defectsites(self):
        """
        If multiple defect sites have same voronoi radius, only one is kept.
        Useful if the symmetry based reduction of initial sites returned
        from Zeo++ is not working properly due to deviation in ideal lattice
        coordinates.
        """
        distinct_radii = list(set(self._radii))
        no_dstnt_radii = len(distinct_radii)
        flag = [False]*no_dstnt_radii
        for rad in distinct_radii:
            ind = self._radii.index(rad) # Index of first site with 'rad' 
            for i in range(len(self._radii)-1,ind,-1): 
                #Backward search for remaining sites so index is not changed
                if self._radii[i] == rad:
                    self._defect_sites.pop(i)
                    self._defectsite_coord_no.pop(i) 
                    self._defect_coord_sites.pop(i)
                    self._radii.pop(i)

    def radius_prune_defectsites(self, radius):
        """
        Remove all the defect sites with voronoi radius less than input radius
        """
        for i in range(len(self._radii)-1,0,-1):
            if self._radii[i] < radius:
                self._defect_sites.pop(i)
                self._defectsite_coord_no.pop(i) 
                self._defect_coord_sites.pop(i)
                self._radii.pop(i)

    def prune_defectsites(self, el="Li", oxi_state=1):
        """
        Prune all the defect sites which can't acoomodate the input elment 
        with the input oxidation state.
        """
        rad = Specie(el, oxi_state).ionic_radius
        self.radius_prune_defectsites(rad)

    def prune_close_defectsites(self, dist=None):
        """
        Prune the sites that are very close.
        """
        if not dist:
            dist = Specie("H",1).ionic_radius# * 2
        ind = 0
        while(ind != len(self._defect_sites)):
            for i in range(ind+1,len(self._defect_sites)):
                d = self._defect_sites[ind].distance(self._defect_sites[i])
                if d < dist:
                    self._defect_sites.pop(i)
                    self._defectsite_coord_no.pop(i) 
                    self._defect_coord_sites.pop(i)
                    self._radii.pop(i)
            ind += 1




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

def _symmetry_reduced_voronoi_nodes(structure, rad_dict):
    """
    Obtain symmetry reduced voronoi nodes using Zeo++ and 
    pymatgen.symmetry.finder.SymmetryFinder
    Args:
        strucutre:
            pymatgen Structure object

    Returns:
        Symmetrically distinct voronoi nodes as pymatgen Strucutre
    """
    vor_node_struct = get_voronoi_nodes(structure, rad_dict)
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
    
    
