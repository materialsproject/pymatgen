#!/usr/bin/env python

'''
This module defines site transformations which transforms a structure into 
another structure. Site transformations differ from standard transformations 
in that they operate in a site-specific manner.
All transformations should inherit the AbstractTransformation ABC.
'''


from __future__ import division

__author__ = "Shyue Ping Ong, Will Richards"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Sep 23, 2011"


from pymatgen.transformations.transformation_abc import AbstractTransformation
from pymatgen.core.structure_modifier import StructureEditor


class ReplaceSiteSpeciesTransformation(AbstractTransformation):
    """
    This transformation substitutes certain sites with certain species.
    """
    def __init__(self, indices_species_map):
        """
        Args:
            indices_species_map:
                A dict containing the species mapping in int-string pairs. 
                E.g., { 1:"Na"} or {2,"Mn2+"}. Multiple substitutions can 
                be done. Overloaded to accept sp_and_occu dictionary
                E.g. {'Si: {'Ge':0.75, 'C':0.25} }, which substitutes a single
                species with multiple species to generate a disordered structure.
        """
        self._indices_species_map = indices_species_map

    def apply_transformation(self, structure):
        editor = StructureEditor(structure)
        for i, sp in self._indices_species_map.items():
            editor.replace_site(int(i), sp)
        return editor.modified_structure

    def __str__(self):
        return "ReplaceSiteSpeciesTransformationTransformation :" + ", ".join([k + "->" + v for k, v in self._species_map.items()])

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        return None

    @property
    def to_dict(self):
        output = {'name' : self.__class__.__name__, 'version': __version__}
        output['init_args'] = {'indices_species_map': self._indices_species_map}
        return output


class RemoveSitesTransformation(AbstractTransformation):
    """
    Remove certain sites in a structure.
    """
    def __init__(self, indices_to_remove):
        """
        Args:
            indices_to_remove:
                List of indices to remove. E.g., [0, 1, 2] 
        """
        self._indices = indices_to_remove

    def apply_transformation(self, structure):
        editor = StructureEditor(structure)
        editor.delete_sites(self._indices)
        return editor.modified_structure

    def __str__(self):
        return "RemoveSitesTransformation :" + ", ".join(self._indices)

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        return None

    @property
    def to_dict(self):
        output = {'name' : self.__class__.__name__, 'version': __version__}
        output['init_args'] = {'indices_to_remove': self._indices}
        return output


class TranslateSitesTransformation(AbstractTransformation):
    """
    This class translates a set of sites by a certain vector.
    """
    def __init__(self, indices_to_move, translation_vector, vector_in_frac_coords = True):
        self._indices = indices_to_move
        self._vector = translation_vector
        self._frac = vector_in_frac_coords

    def apply_transformation(self, structure):
        editor = StructureEditor(structure)
        editor.translate_sites(self._indices, self._vector, self._frac)
        return editor.modified_structure

    def __str__(self):
        return "TranslateSitesTransformation for indices {}, vector {} and vector_in_frac_coords = {}".format(self._indices, self._translation_vector, self._frac)

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        return TranslateSitesTransformation(self._indices, [-c for c in self._vector], self._frac)

    @property
    def to_dict(self):
        output = {'name' : self.__class__.__name__, 'version': __version__}
        output['init_args'] = {'indices_to_move': self._indices,
                               'translation_vector': self._vector,
                               'vector_in_frac_coords': self._frac}
        return output
