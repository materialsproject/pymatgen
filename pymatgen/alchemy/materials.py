#!/usr/bin/env python

"""
This module provides various representations of transformed structures. A
TransformedStructure is a structure that has been modified by undergoing a
series of transformations.
"""

from __future__ import division

__author__ = "Shyue Ping Ong, Will Richards"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Mar 2, 2012"

import os
import re
import json
import datetime
from copy import deepcopy

from pymatgen.core.structure import Structure
from pymatgen.io.cifio import CifParser
from pymatgen.io.vaspio.vasp_input import Poscar
from pymatgen.serializers.json_coders import MSONable, PMGJSONDecoder
from pymatgen.matproj.snl import StructureNL

from warnings import warn

dec = PMGJSONDecoder()

class TransformedStructure(MSONable):
    """
    Container object for new structures that include history of
    transformations.

    Each transformed structure is made up of a sequence of structures with
    associated transformation history.
    """
    
    def __init__(self, structure, transformations=None, history=None,
                 other_parameters=None):
        self.final_structure = structure
        self.history = history or []
        self.other_parameters = other_parameters or {}
        self._undone = []
        
        transformations = transformations or []
        for t in transformations:
            self.append_transformation(t)
            
    def undo_last_change(self):
        """
        Undo the last change in the TransformedStructure.

        Raises:
            IndexError if already at the oldest change.
        """
        if len(self.history) == 0:
            raise IndexError("Can't undo. Already at oldest change.")
        if 'input_structure' not in self.history[-1]:
            raise IndexError("Can't undo. Latest history has no "
                             "input_structure")
        h = self.history.pop()
        self._undone.append((h, self.final_structure))
        s = h["input_structure"]
        if isinstance(s, dict):
            s = Structure.from_dict(s)
        self.final_structure = s

    def redo_next_change(self):
        """
        Redo the last undone change in the TransformedStructure.

        Raises:
            IndexError if already at the latest change.
        """
        if len(self._undone) == 0:
            raise IndexError("Can't redo. Already at latest change.")
        h, s = self._undone.pop()
        self.history.append(h)
        self.final_structure = s
        
    def __getattr__(self, name):
        s = object.__getattribute__(self, 'final_structure')
        return getattr(s, name)
    
    def __len__(self):
        return len(self.history)
    
    def append_transformation(self, transformation, return_alternatives=False,
                              clear_redo=True):
        """
        Appends a transformation to the TransformedStructure.

        Args:
            transformation:
                Transformation to append
            return_alternatives:
                Whether to return alternative TransformedStructures for
                one-to-many transformations. return_alternatives can be a
                number, which stipulates the total number of structures to
                return.
            clear_redo:
                Boolean indicating whether to clear the redo list. By default,
                this is True, meaning any appends clears the history of
                undoing. However, when using append_transformation to do a
                redo, the redo list should not be cleared to allow multiple
                redos.
        """
        if clear_redo:
            self._undone = []
            
        if return_alternatives and transformation.is_one_to_many:
            ranked_list = transformation.apply_transformation(
                self.final_structure, return_ranked_list=return_alternatives)
            
            input_structure = self.final_structure.to_dict
            alts = []
            for x in ranked_list[1:]:
                s = x.pop("structure")
                actual_transformation = x.pop("transformation", transformation)
                hdict = actual_transformation.to_dict
                hdict["input_structure"] = input_structure
                hdict["output_parameters"] = x
                self.final_structure = s
                d = self.to_dict
                d['history'].append(hdict)
                d['final_structure'] = s.to_dict
                alts.append(TransformedStructure.from_dict(d))
                
            x = ranked_list[0]
            s = x.pop("structure")
            actual_transformation = x.pop("transformation", transformation)
            hdict = actual_transformation.to_dict
            hdict["input_structure"] = self.final_structure.to_dict
            hdict["output_parameters"] = x
            self.history.append(hdict)
            self.final_structure = s
            return alts
        else:
            s = transformation.apply_transformation(self.final_structure)
            hdict = transformation.to_dict
            hdict["input_structure"] = self.final_structure.to_dict
            hdict["output_parameters"] = {}
            self.history.append(hdict)
            self.final_structure = s
            
    def append_filter(self, structure_filter):
        """
        Adds a transformation parameter to the last transformation.
        """
        hdict = structure_filter.to_dict
        hdict["input_structure"] = self.final_structure.to_dict
        self.history.append(hdict)
        
    def extend_transformations(self, transformations):
        """
        Extends a sequence of transformations to the TransformedStructure.

        Args:
            transformations:
                Sequence of Transformations
        """
        for t in transformations:
            self.append_transformation(t)
            
    def get_vasp_input(self, vasp_input_set, generate_potcar=True):
        """
        Returns VASP input as a dict of vaspio objects.

        Args:
            vasp_input_set:
                pymatgen.io.vaspio_set.VaspInputSet like object that creates
                vasp input files from structures
            generate_potcar:
                Set to False to generate a POTCAR.spec file instead of a
                POTCAR, which contains the POTCAR labels but not the actual
                POTCAR. Defaults to True.
        """
        d = vasp_input_set.get_all_vasp_input(self.final_structure,
                                              generate_potcar)
        d["transformations.json"] = json.dumps(self.to_dict)
        return d

    def write_vasp_input(self, vasp_input_set, output_dir,
                         create_directory=True):
        """
        Writes VASP input to an output_dir.

        Args:
            vasp_input_set:
                pymatgen.io.vaspio_set.VaspInputSet like object that creates
                vasp input files from structures
            output_dir:
                Directory to output files
            create_directory:
                Create the directory if not present. Defaults to True.
        """
        vasp_input_set.write_input(self.final_structure, output_dir,
                                   make_dir_if_not_present=create_directory)
        with open(os.path.join(output_dir, "transformations.json"), "w") as fp:
            json.dump(self.to_dict, fp)

    def __str__(self):
        output = ["Current structure", "------------",
                  str(self.final_structure),
                  "\nHistory",
                  "------------"]
        for h in self.history:
            h.pop('input_structure', None)
            output.append(str(h))
        output.append("\nOther parameters")
        output.append("------------")
        output.append(str(self.other_parameters))
        return "\n".join(output)
    
    def set_parameter(self, key, value):
        self.other_parameters[key] = value
    
    @property
    def was_modified(self):
        """
        Boolean describing whether the last transformation on the structure
        made any alterations to it one example of when this would return false
        is in the case of performing a substitution transformation on the
        structure when the specie to replace isn't in the structure.
        """
        return not self.final_structure == self.structures[-2]
    
    @property
    def structures(self):
        """
        Returns a copy of all structures in the TransformedStructure. A
        structure is stored after every single transformation.
        """
        hstructs = [Structure.from_dict(s['input_structure'])
                    for s in self.history if 'input_structure' in s]
        return hstructs + [self.final_structure]
            
    @staticmethod
    def from_cif_string(cif_string, transformations=None, primitive=True,
                        occupancy_tolerance=1.):
        """
        Generates TransformedStructure from a cif string.

        Args:
            cif_string:
                Input cif string. Should contain only one structure. For cifs
                containing multiple structures, please use CifTransmuter.
            transformations:
                Sequence of transformations to be applied to the input
                structure.
            primitive:
                Option to set if the primitive cell should be extracted.
                Defaults to True. However, there are certain instances where
                you might want to use a non-primitive cell, e.g., if you are
                trying to generate all possible orderings of partial removals
                or order a disordered structure.
            occupancy_tolerance:
                If total occupancy of a site is between 1 and
                occupancy_tolerance, the occupancies will be scaled down to 1.
        """
        parser = CifParser.from_string(cif_string, occupancy_tolerance)
        raw_string = re.sub("'", "\"", cif_string)
        cif_dict = parser.to_dict
        cif_keys = cif_dict.keys()
        s = parser.get_structures(primitive)[0]
        partial_cif = cif_dict[cif_keys[0]]
        if "_database_code_ICSD" in partial_cif:
            source = partial_cif["_database_code_ICSD"] + "-ICSD"
        else:
            source = "uploaded cif"
        source_info = {"source": source,
                       "datetime": str(datetime.datetime.now()),
                       "original_file": raw_string,
                       "cif_data": cif_dict[cif_keys[0]]}
        return TransformedStructure(s, transformations, history=[source_info])

    @staticmethod
    def from_poscar_string(poscar_string, transformations=None):
        """
        Generates TransformedStructure from a poscar string.

        Args:
            poscar_string:
                Input POSCAR string.
            transformations:
                Sequence of transformations to be applied to the input
                structure.
        """
        p = Poscar.from_string(poscar_string)
        if not p.true_names:
            raise ValueError("Transformation can be craeted only from POSCAR "
                             "strings with proper VASP5 element symbols.")
        raw_string = re.sub("'", "\"", poscar_string)
        s = p.structure
        source_info = {"source": "POSCAR",
                       "datetime": str(datetime.datetime.now()),
                       "original_file": raw_string}
        return TransformedStructure(s, transformations, history=[source_info])
    
    @property
    def to_dict(self):
        """
        Returns a dict representation of the TransformedStructure.
        """
        d = self.final_structure.to_dict
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        d["history"] = deepcopy(self.history)
        d["version"] = __version__
        d["last_modified"] = str(datetime.datetime.utcnow())
        d["other_parameters"] = deepcopy(self.other_parameters)
        return d
    
    @classmethod
    def from_dict(cls, d):
        """
        Creates a TransformedStructure from a dict.
        """
        s = Structure.from_dict(d)
        return cls(s, history=d["history"], 
                   other_parameters=d.get("other_parameters", None))
        
    def to_snl(self, authors, projects=None, references='', remarks=None,
               data=None, created_at=None):
        if self.other_parameters:
            warn('Data in TransformedStructure.other_parameters discarded '
                 'during type conversion to SNL')
        hist = []
        for h in self.history:
            snl_metadata = h.pop('_snl', {})
            hist.append({'name' : snl_metadata.pop('name', 'pymatgen'),
                         'url' : snl_metadata.pop('url', 
                                    'http://pypi.python.org/pypi/pymatgen'),
                         'description' : h})
        return StructureNL(self.final_structure, authors, projects, references,
                           remarks, data, hist, created_at)
        
    @classmethod
    def from_snl(cls, snl):
        """
        Create TransformedStructure from SNL.

        Args:
            snl:
                Starting snl
            copy_metadata:
                update the authors, projects, references, and remarks
                in the last history node with the metadata at the
                root of the SNL object
            copy_data:
                copy the contents of snl.data into the last history
                node

        Returns:
            TransformedStructure.
        """
        hist = []
        for h in snl.history:
            d = h.description
            d['_snl'] = {'url' : h.url, 'name' : h.name}
            hist.append(d)
        return cls(snl.structure, history=hist)

