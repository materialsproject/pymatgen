#!/usr/bin/env python

'''
This module defines standard transformations which transforms a structure into another structure.
All transformations should inherit the AbstractTransformation ABC.
'''

from __future__ import division

__author__="Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Sep 23, 2011"

import json
import itertools
import warnings

from pymatgen.transformations.transformation_abc import AbstractTransformation
from pymatgen.core.structure import Structure
from pymatgen.core.operations import SymmOp
from pymatgen.core.structure_modifier import StructureEditor
from pymatgen.core.periodic_table import smart_element_or_specie
from pymatgen.analysis.ewald import EwaldSummation

class IdentityTransformation(AbstractTransformation):
    """
    This is a demo transformation which does nothing, i.e. just return the same structure.
    """
    
    def __init__(self):
        pass
    
    def apply_transformation(self, structure):
        return Structure(structure.lattice, structure.species_and_occu, structure.frac_coords)
    
    def __str__(self):
        return "Identity Transformation"
    
    def __repr__(self):
        return self.__str__()
    
    @property
    def inverse(self):
        return self
    
    @property
    def to_dict(self):
        output = {'name' : self.__class__.__name__, 'init_args': {} }
        return output
      
class RotationTransformation(AbstractTransformation):
    """
    The RotationTransformation applies a rotation to a structure.
    """
    
    def __init__(self, axis, angle, angle_in_radians = False):
        """
        Arguments:
            axis - Axis of rotation, e.g., [1, 0, 0]
            angle - angle to rotate
            angle_in_radians - Set to True if angle is supplied in radians. Else degrees are assumed.
        """
        self._axis = axis
        self._angle = angle
        self._angle_in_radians = angle_in_radians
        self._symmop = SymmOp.from_axis_angle_and_translation(self._axis, self._angle, self._angle_in_radians)
    
    def apply_transformation(self, structure):
        editor = StructureEditor(structure)
        editor.apply_operation(self._symmop)
        return editor.modified_structure
    
    def __str__(self):
        return "Rotation Transformation about axis %s with angle = %.4f %s" % (str(self._axis), self._angle, "radians" if self._angle_in_radians else "degrees")
    
    def __repr__(self):
        return self.__str__()
    
    @property
    def inverse(self):
        return RotationTransformation(self._axis, - self._angle, self._angle_in_radians)
    
    @property
    def to_dict(self):
        output = {'name' : self.__class__.__name__}
        output['init_args'] = {'axis': self._axis, 'angle':self._angle, 'angle_in_radians':self._angle_in_radians}
        return output
        
    
    
class SubstitutionTransformation(AbstractTransformation):
    """
    This transformation substitutes species for one another.
    """
    def __init__(self, species_map):
        """
        Arguments:
            species_map - a dict containing the species mapping in string-string pairs. E.g., { "Li":"Na"} or {"Fe2+","Mn2+"}. Multiple substitutions can be done.
        """
        self._species_map = species_map
    
    def apply_transformation(self, structure):
        species_map = {smart_element_or_specie(k): smart_element_or_specie(v) for k, v in self._species_map.items()}
        editor = StructureEditor(structure)
        editor.replace_species(species_map)
        return editor.modified_structure
    
    def __str__(self):
        return "Substitution Transformation :" + ", ".join([k+"->"+v for k, v in self._species_map.items()])
    
    def __repr__(self):
        return self.__str__()
    
    @property
    def inverse(self):
        return SubstitutionTransformation({v:k for k, v in self._species_map.items()})
    
    @property
    def to_dict(self):
        output = {'name' : self.__class__.__name__}
        output['init_args'] = {'species_map': self._species_map}
        return output
        
     
class RemoveSpeciesTransformation(AbstractTransformation):
    """
    Remove all occurrences of some species from a structure.
    """
    def __init__(self, species_to_remove):
        """
        Arguments:
            species_to_remove - List of species to remove. E.g., ["Li", "Mn"] 
        """
        self._species = species_to_remove
    
    def apply_transformation(self, structure):
        editor = StructureEditor(structure)
        map(editor.remove_species, [[smart_element_or_specie(sp)] for sp in self._species])
        return editor.modified_structure
    
    def __str__(self):
        return "Remove Species Transformation :" + ", ".join(self._species)
    
    def __repr__(self):
        return self.__str__()
    
    @property
    def inverse(self):
        return None
    
    @property
    def to_dict(self):
        output = {'name' : self.__class__.__name__}
        output['init_args'] = {'species_to_remove': self._species}
        return output
        
class PartialRemoveSpecieTransformation(AbstractTransformation):
    """
    Remove fraction of specie from a structure. Requires an oxidation state decorated structure for ewald sum to be computed.
    """
    def __init__(self, specie_to_remove, fraction_to_remove):
        """
        Arguments:
            specie_to_remove - Specie to remove. Must have oxidation state E.g., "Li1+"
            fraction_to_remove - Fraction of specie to remove. E.g., 0.5
        """
        self._specie = specie_to_remove
        self._frac = fraction_to_remove
    
    def apply_transformation(self, structure):
        sp = smart_element_or_specie(self._specie)
        num_to_remove = structure.composition[sp] * self._frac
        if abs(num_to_remove - int(num_to_remove)) > 1e-8:
            raise ValueError("Fraction to remove must be consistent with integer amounts in structure.")
        else:
            num_to_remove = int(round(num_to_remove))
        
        specie_indices = [i for i in xrange(len(structure)) if structure[i].specie == sp]
                
        lowestewald = 1e100
        lowestenergy_s = None
        all_structures = list()
        for indices in itertools.combinations(specie_indices, num_to_remove):
            mod = StructureEditor(structure)
            mod.delete_sites(indices)
            s_new = mod.modified_structure.get_sorted_structure()
            all_structures.append(s_new)
            ewaldsum = EwaldSummation(s_new)
            if ewaldsum.total_energy < lowestewald:
                lowestewald = ewaldsum.total_energy
                lowestenergy_s = s_new
        
        self.all_structures = all_structures        
        
        return lowestenergy_s
    
    def __str__(self):
        return "Remove Species Transformation :" + ", ".join(self._specie)
    
    def __repr__(self):
        return self.__str__()
    
    @property
    def inverse(self):
        return None
    
    @property
    def to_dict(self):
        output = {'name' : self.__class__.__name__}
        output['init_args'] = {'specie_to_remove': self._specie, 'fraction_to_remove': self._frac}
        return output

class OrderDisorderedStructureTransformation(AbstractTransformation):
    """
    Order a disordered structure. The disordered structure must be oxidation state decorated for ewald sum to be computed.
    Please note that the current form uses a "dumb" algorithm of completely enumerating all possible combinations of
    partially occupied sites.  No attempt is made to perform symmetry determination to reduce the number of combinations.
    Hence, attempting to performing ordering on a large number of disordered sites may be extremely expensive.  Also, simple
    rounding of the occupancies are performed, with no attempt made to achieve a target composition.  This is usually not a problem
    for most ordering problems, but there can be times where rounding errors may result in structures that do not have the desired composition.
    This second step will be implemented in the next iteration of the code.
    
    USE WITH CARE.
    """
    def __init__(self):
        pass
    
    def apply_transformation(self, structure, max_iterations = 100):
        """
        For this transformation, the apply_transformation method will return only the ordered
        structure with the lowest Ewald energy, to be consistent with the method signature of the other transformations.  
        However, all structures are stored in the  all_structures attribute in the transformation object for easy access.
        Arguments:
            structure - Oxidation state decorated disordered structure to order
            max_iterations - Maximum number of structures to consider.  Defaults to 100. This is useful if there are a large number of sites 
                             and there are too many orderings to enumerate.
        """
        ordered_sites = []
        
        sites_to_order = {}
        for site in structure:
            species_and_occu = site.species_and_occu
            if sum(species_and_occu.values()) == 1 and len(species_and_occu) == 1:
                ordered_sites.append(site)
            else:
                spec = tuple([(sp, occu) for sp, occu in species_and_occu.items()])
                if spec not in sites_to_order:
                    sites_to_order[spec] = [site]
                else:
                    sites_to_order[spec].append(site)

        allselections = []
        species = []
        for spec, sites in sites_to_order.items():
            total_sites = len(sites)
            for (sp, fraction) in spec:
                num_to_select = int(round(fraction * total_sites))
                if num_to_select == 0:
                    raise ValueError("Fraction not consistent with selection of at least a single site.  Make a supercell before proceeding further.")
                allselections.append(itertools.combinations(sites, num_to_select))
                species.append(sp)
        
        all_ordered_s = {}
        count = 0
        
        def in_coords(allcoords, coord):
            for test_coord in allcoords:
                if all(coord == test_coord):
                    return True
            return False
        
        for selection in itertools.product(*allselections):
            all_species = [site.species_and_occu for site in ordered_sites]
            all_coords = [site.frac_coords for site in ordered_sites]
            
            contains_dupes = False
            for i in xrange(len(selection)):
                subsel = selection[i]
                sp = species[i]
                for site in subsel:
                    if not in_coords(all_coords, site.frac_coords):
                        all_species.append(sp)
                        all_coords.append(site.frac_coords)
                    else:
                        contains_dupes = True
                        break
                if contains_dupes:
                    break
            
            if not contains_dupes:
                s = Structure(structure.lattice, all_species, all_coords, False).get_sorted_structure()
                ewaldsum = EwaldSummation(s)
                ewald_energy = ewaldsum.total_energy
                all_ordered_s[s] = ewald_energy
                count += 1
                if count == max_iterations:
                    warnings.warn("Maximum number of iterations reached.  Structures will be ordered based on " + str(max_iterations) + " structures.")
                    break

        self.all_structures = all_ordered_s
        sorted_structures = sorted(all_ordered_s.keys(), key=lambda a: all_ordered_s[a])
                
        return sorted_structures[0]
    
    def __str__(self):
        return "Order disordered structure transformation"
     
    def __repr__(self):
        return self.__str__()
    
    @property
    def inverse(self):
        return None
    
    @property
    def to_dict(self):
        output = {'name' : self.__class__.__name__}
        output['init_args'] = {}
        return output


def transformation_from_json(json_string):
    """
    A helper function that can simply get a transformation from a json representation.
    """
    jsonobj = json.loads(json_string)
    trans = globals()[jsonobj['name']]
    return trans(**jsonobj['init_args'])


    