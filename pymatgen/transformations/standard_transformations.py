#!/usr/bin/env python

'''
Created on Sep 23, 2011
'''

__author__="Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Sep 23, 2011"

import json
import itertools

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
        
        for indices in itertools.combinations(specie_indices, num_to_remove):
            mod = StructureEditor(structure)
            mod.delete_sites(indices)
            s_new = mod.modified_structure
            ewaldsum = EwaldSummation(s_new)
            if ewaldsum.total_energy < lowestewald:
                lowestewald = ewaldsum.total_energy
                lowestenergy_s = s_new
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


def transformation_from_json(json_string):
    """
    A helper function that can simply get a transformation from a json representation.
    """
    jsonobj = json.loads(json_string)
    trans = globals()[jsonobj['name']]
    return trans(**jsonobj['init_args'])

