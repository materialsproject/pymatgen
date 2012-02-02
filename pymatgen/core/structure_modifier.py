#!/usr/bin/env python

"""
This module provides classes used to modify structures.
"""

from __future__ import division

__author__="Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__status__ = "Production"
__date__ ="Sep 23, 2011"

import abc
import itertools

import numpy as np
from pymatgen.core.periodic_table import Specie, Element
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure, PeriodicSite


class StructureModifier(object):
    """
    Abstract class definition for all classes that modify structures.
    """
    __metaclass__ = abc.ABCMeta
    
    @abc.abstractproperty
    def modified_structure(self):
        """
        Returns the modified structure.
        """
        return 
    
    @abc.abstractproperty
    def original_structure(self):
        '''
        Returns the original structure.
        '''
        return 
    
class StructureEditor(StructureModifier):   
    """
    Editor for adding, removing and changing sites from a structure
    """
    DISTANCE_TOLERANCE = 0.01
    
    def __init__(self, structure):
        """
        Arguments:
            structure:
                pymatgen.core.structure Structure object.
        """
        self._original_structure = structure
        self._lattice = structure.lattice
        self._sites = list(structure.sites)
    
    def replace_species(self, species_mapping):
        """
        Swap species in a structure.
        
        Arguments:
            species_mapping:
                dict of species to swap. Species can be elements too.
                e.g., {Element("Li"): Element("Na")} performs a Li for Na substitution.
        """
        
        def mod_site(site):
            new_atom_occu = dict()
            for sp, amt in site.species_and_occu.items():
                if sp in species_mapping:
                    new_atom_occu[species_mapping[sp]] = amt
                else:
                    new_atom_occu[sp] = amt 
            return PeriodicSite(new_atom_occu, self._lattice.get_fractional_coords(site.coords), self._lattice)    
        self._sites = map(mod_site, self._sites)
    
    def replace_single_site(self, index, species = None, atoms_n_occu = None):
        """
        Replace a single site. Takes either a species or a dict of occus
        
        Arguments:
            species: a species object
            index: the index of the site in the _sites list
                
        """
        if atoms_n_occu == None:
            atoms_n_occu = dict()
            atoms_n_occu[species] = 1
        
        self._sites[index] = PeriodicSite(atoms_n_occu, self._lattice.get_fractional_coords(self._sites[index].coords), self._lattice)
    
    def remove_species(self, species):
        """
        Remove all occurrences of a species from a structure.
        
        Args:
            species:
                species to remove
        """    
        new_sites = []
        for site in self._sites:
            new_sp_occu = {sp:amt for sp, amt in site.species_and_occu.items() if sp not in species}
            if len(new_sp_occu) > 0:
                new_sites.append(PeriodicSite(new_sp_occu, site.frac_coords, self._lattice))
        self._sites = new_sites
    
    def append_site(self, species, coords, coords_are_cartesian = False, validate_proximity = True): 
        """
        Append a site to the structure at the end.
        
        Args:
            species:
                species of inserted site
            coords:
                coordinates of inserted site
            fractional_coord:
                Whether coordinates are cartesian. Defaults to False.
            validate_proximity:
                Whether to check if inserted site is too close to an existing site. Defaults to True.
        
        """       
        self.insert_site(len(self._sites), species, coords, coords_are_cartesian, validate_proximity)
        
    def insert_site(self, i, species, coords, coords_are_cartesian = False, validate_proximity = True):
        """
        Insert a site to the structure.
        
        Args:
            i:
                index to insert site
            species:
                species of inserted site
            coords:
                coordinates of inserted site
            coords_are_cartesian:
                Whether coordinates are cartesian. Defaults to False.
            validate_proximity:
                Whether to check if inserted site is too close to an existing site. Defaults to True.
        """
        if not coords_are_cartesian:
            new_site = PeriodicSite(species,coords,self._lattice)
        else:
            new_site = PeriodicSite(species, self._lattice.get_fractional_coords(coords), self._lattice)
            
        if validate_proximity:
            for site in self._sites:
                if site.distance(new_site) < self.DISTANCE_TOLERANCE:
                    raise ValueError("New site is too close to an existing site!")
        
        self._sites.insert(i, new_site)
        
    def delete_site(self, i):
        """
        Delete site at index i.
        
        Args:
            i:
                index of site to delete.
        """
        del(self._sites[i])
        
    def delete_sites(self, indices):
        """
        Delete sites with at indices.
        
        Args:
            indices:
                sequence of indices of sites to delete.
        """
        self._sites = [self._sites[i] for i in xrange(len(self._sites)) if i not in indices]
    
    def apply_operation(self, symmop):
        """
        Apply a symmetry operation to the structure and return the new structure.
        The lattice is operated by the rotation matrix only.
        Coords are operated in full and then transformed to the new lattice.
        
        Args:
            symmop:
                Symmetry operation to apply.
        """        
        self._lattice = Lattice([symmop.apply_rotation_only(row) for row in self._lattice.matrix])
        def operate_site(site):
            new_cart = symmop.operate(site.coords)
            return PeriodicSite(site.species_and_occu, self._lattice.get_fractional_coords(new_cart), self._lattice) 
        self._sites = map(operate_site, self._sites)
    
    def modify_lattice(self, new_lattice):
        """
        Modify the lattice of the structure.  Mainly used for changing the basis.
        
        Args:
            new_lattice:
                New lattice
        """
        self._lattice = new_lattice
        new_sites = []
        for site in self._sites:
            new_sites.append(PeriodicSite(site.species_and_occu, self._lattice.get_fractional_coords(site.coords), self._lattice))
        self._sites = new_sites


    def apply_strain_transformation(self, F=np.identity(3)):
		"""
		Apply a  deformation gradient tensor F to a lattice. Defaults to identity tensor. 
		Note: fractional ionic coordinates should not change upon applying this transformation, but
		Cartesian coordinates do!

		Args:
			F:
				deformation gradient tensor (3x3 numpy matrix)
		"""
		self._lattice = np.matrix(self._lattice._matrix)*F
	

    @property
    def original_structure(self):
        return self._original_structure
    
    @property
    def modified_structure(self):
        return Structure(self._lattice, [site.species_and_occu for site in self._sites], [site.frac_coords for site in self._sites], False)
    
class SupercellMaker(StructureModifier):
    """
    Makes a supercell
    """
    
    def __init__(self, structure, scaling_matrix = [[1,0,0],[0,1,0],[0,0,1]]):
        """
        Create a supercell.
        
        Arguments:
            structure:
                pymatgen.core.structure Structure object.
            scaling_matrix:
                a matrix of transforming the lattice vectors. Defaults to the identity matrix.
                Has to be all integers. e.g., [[2,1,0],[0,3,0],[0,0,1]] generates a new structure
                with lattice vectors a' = 2a + b, b' = 3b, c' = c where a, b, and c are the lattice 
                vectors of the original structure. 
        """
        self._original_structure = structure
        old_lattice = structure.lattice
        scale_matrix = np.array(scaling_matrix)
        new_lattice = Lattice(np.dot(scale_matrix, old_lattice.matrix))
        new_species = []
        new_fcoords = []
        def range_vec(i):
            return xrange(max(scale_matrix[:][:,i]) - min(scale_matrix[:][:,i]))
        for site in structure.sites:
            for (i, j, k) in itertools.product(range_vec(0), range_vec(1), range_vec(2)):
                new_species.append(site.species_and_occu)
                fcoords = site.frac_coords
                coords = old_lattice.get_cartesian_coords(fcoords + np.array([i, j, k]))
                new_fcoords.append(new_lattice.get_fractional_coords(coords))
        self._modified_structure = Structure(new_lattice, new_species, new_fcoords, False)
                        
    @property
    def original_structure(self):
        return self._original_structure
    
    @property
    def modified_structure(self):
        return self._modified_structure

class OxidationStateDecorator(StructureModifier):
    """
    Given a dictionary of oxidation states, decorate a structure 
    by replacing each Element at a site with a Specie with an oxidation state.
    Useful for higher level functions. 
    """
    
    def __init__(self, structure, oxidation_states):
        """
        Decorates a structure with oxidation states.
        
        Args:
            structure:
                pymatgen.core.structure Structure object.
            oxidation_states:
                dict of oxidation states. e.g., {"Li":1, "Fe":2, "P":5, "O": -2} 
        """
        self._original_structure = structure
        try:
            new_species = [{Specie(el.symbol, oxidation_states[el.symbol]) : occu for el, occu in site.species_and_occu.items()} for site in structure]
        except KeyError as ex:
            raise ValueError("Oxidation state of all elements must be specified in the dictionary.")
        self._modified_structure = Structure(structure.lattice, new_species, structure.frac_coords, False)
                        
    @property
    def original_structure(self):
        return self._original_structure
    
    @property
    def modified_structure(self):
        return self._modified_structure
    
class OxidationStateRemover(StructureModifier):
    """
    Replace each Specie at a site with an element.
    Useful for doing structure comparisons after applying higher level functions
    """
    
    def __init__(self, structure):
        """
        Removes oxidation states from a structure
        
        Args:
            structure:
                pymatgen.core.structure Structure object.
        """
        self._original_structure = structure
        new_species = [{Element(el.symbol) : occu for el, occu in site.species_and_occu.items()} for site in structure]
        self._modified_structure = Structure(structure.lattice, new_species, structure.frac_coords, False)
                        
    @property
    def original_structure(self):
        return self._original_structure
    
    @property
    def modified_structure(self):
        return self._modified_structure

class BasisChange(StructureModifier):
    """
    Given a new basis, we express the structure in this new basis
    """
    
    def __init__(self, structure, new_lattice):
        """
        Express a given structure in a new basis.
        
        Arguments:
            structure:
                pymatgen.core.structure Structure object.
            new_lattice:
                a pymatgen.core.Lattice object
        """
        self._original_structure = structure
        sp=[site.species_and_occu for site in structure._sites]
        coords=[site.coords for site in structure._sites]
        self._modified_structure = Structure(new_lattice, sp, coords, validate_proximity = False, to_unit_cell = True, coords_are_cartesian = True)
                        
    @property
    def original_structure(self):
        return self._original_structure
    
    @property
    def modified_structure(self):
        return self._modified_structure    

    
if __name__ == "__main__":
    import doctest
    doctest.testmod()
