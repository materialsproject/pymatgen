#!/usr/bin/env python

"""
This module defines classes for point defects
"""
from __future__ import division
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
        Given a structure, generate unique interstitial sites
        
        Args:
            structure:
                pymatgen.core.Structure
        """
        
        self._structure = structure
        self._uniq_interstitial_sites = []

        #Use Zeo++ to obtain the voronoi nodes. The unique voronoi nodes
        #are possible candidates for interstitial sites
        possibleInterstitialSites = _get_zeo_uniq_voronoi_nodes(self._structure)
        
        #Do futher processing on possibleInterstitialSites to obtain 
        #interstitial sites
        self._uniq_interstitial_sites = possibleInterstitialSites
        
        
    def enumerate_uniq_defectsites(self):
        return self._uniq_interstitial_sites

            
    def _supercell_with_defect(self, scaling_matrix, defect_site, element):
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
                sc.replace(i, element)
                return sc
           

    def make_supercells_with_defects(self, scaling_matrix, element):
        """
        Returns sequence of supercells in pymatgen.core.structure.Structure
        format, with each supercell containing unique interstitial
        """
        sc_with_interstitial = []
        for uniq_defect_site in self.enumerate_uniq_defectsites():
            sc_with_interstitial.append(self._supercell_with_defect(
                                    scaling_matrix, uniq_defect_site, element))
        return sc_with_interstitial

def _get_zeo_uniq_voronoi_nodes(structure):
    """
    Obtain uniq voronoi nodes using Zeo++
    """
    from pymatgen.io.cifio import CifWriter
    from zeo.netstorage import AtomNetwork, VoronoiNetwork
    from voronoi_node_symmetry import uniq_voronoi_nodes
    #cssrStruct = Cssr(structure)
    cifwriteobj = CifWriter(structure)
    name = "tmp"
    #tmpCssrFile = name+".cssr"
    tmpCifFile = name+".cif"
    tmpVoroXyzFile = name+"_voro.xyz"
    #cssrStruct.write_file(tmpCssrFile)
    cifwriteobj.write_file(tmpCifFile)
    # Use pymatgen radii for elements in future
    atmnet = AtomNetwork.readfromCif(tmpCifFile)
    vornet = atmnet.perform_voronoi_decomposition()
    tmpProbeRad = 0.1
    vornet.analyze_writeto_xyz(name, tmpProbeRad, atmnet)
    return uniq_voronoi_nodes(tmpCifFile, tmpVoroXyzFile)
    
    
