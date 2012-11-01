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
__email__ = "shyue@mit.edu"
__date__ = "Mar 2, 2012"

import os
import re
import json
import datetime
import warnings
from copy import deepcopy

from pymatgen.core.structure import Structure
from pymatgen.transformations.transformation_abc import AbstractTransformation

from pymatgen.io.cifio import CifParser
from pymatgen.io.vaspio.vasp_input import Poscar
from pymatgen.serializers.json_coders import MSONable


class TransformedStructure(MSONable):
    """
    Container object for new structures that include history of
    transformations.

    Each transformed structure is made up of a sequence of structures with
    associated transformation history.
    """

    def __init__(self, structure, transformations, history=None,
                 other_parameters=None):
        """
        Standard constructor for a TransformedStructure.

        Args:
            structure:
                input structure
            transformations:
                Sequence of transformations to be applied to the input
                structure.
            history:
                optional history for the input structure, which provides a way
                to track structures having undergone multiple series of
                transformations.
            other_parameters:
                optional parameters to store along with the
                TransformedStructure. This can include tags (a list) or author
                which will be parsed.
        """
        history = [] if history == None else history
        self._source = {}
        self._structures = []
        self._changes = []
        self._change_parameters = []
        self._redo_trans = []
        self._other_parameters = {} if other_parameters is None \
            else deepcopy(other_parameters)
        if len(history) > 0:
            self._source = history[0]
            for i in xrange(1, len(history)):
                struct = Structure.from_dict(history[i]["input_structure"])
                trans = AbstractTransformation.from_dict(history[i])
                param = history[i].get("output_parameters", {})
                self._structures.append(struct)
                self._changes.append(trans)
                self._change_parameters.append(param)

        self._structures.append(structure)
        for t in transformations:
            self.append_transformation(t)

    def undo_last_transformation(self):
        """
        .. deprecated:: v2.2.2
        """
        warnings.warn("Deprecated. Use undo_last_change.", DeprecationWarning)
        self.undo_last_change()

    def redo_next_transformation(self):
        """
        .. deprecated:: v2.2.2
        """
        warnings.warn("Deprecated. Use redo_last_change.", DeprecationWarning)
        self.redo_next_change()

    def undo_last_change(self):
        """
        Undo the last change in the TransformedStructure.

        Raises:
            IndexError if already at the oldest change.
        """
        if len(self._changes) == 0:
            raise IndexError("Can't undo. Already at oldest change.")
        self._structures.pop()
        self._change_parameters.pop()
        self._redo_trans.append(self._changes.pop())

    def redo_next_change(self):
        """
        Redo the last undone change in the TransformedStructure.

        Raises:
            IndexError if already at the latest change.
        """
        if len(self._redo_trans) == 0:
            raise IndexError("Can't undo. Already at latest change.")
        t = self._redo_trans.pop()
        if hasattr(t, 'apply_transformation'):
            self.append_transformation(t, clear_redo=False)
        else:
            self.append_filter(t)

    def __getitem__(self, index):
        return (self._structures[index], self._changes[0:index])

    def __getattr__(self, name):
        return getattr(self._structures[-1], name)

    def __len__(self):
        return len(self._structures)

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
            self._redo_trans = []

        if return_alternatives and transformation.is_one_to_many:
            starting_struct = self._structures[-1]
            ranked_list = transformation.apply_transformation(starting_struct,
                                        return_ranked_list=return_alternatives)
            #generate the alternative structures
            alts = []
            for x in ranked_list[1:]:
                struct = x.pop("structure")
                other_paras = self._other_parameters.copy()
                hist = self.history
                actual_transformation = x.pop("transformation", transformation)
                tdict = actual_transformation.to_dict
                tdict["input_structure"] = starting_struct.to_dict
                tdict["output_parameters"] = x
                hist.append(tdict)
                alts.append(TransformedStructure(struct, [], history=hist,
                                                 other_parameters=other_paras))
            #use the first item in the ranked_list and apply it to this
            #transformed_structure
            x = ranked_list[0]
            struct = x.pop("structure")
            actual_transformation = x.pop("transformation", transformation)
            self._structures.append(struct)
            self._changes.append(actual_transformation)
            self._change_parameters.append(x)
            return alts
        else:
            new_s = transformation.apply_transformation(self._structures[-1])
            self._structures.append(new_s)
            self._change_parameters.append({})
            self._changes.append(transformation)

    def append_filter(self, structure_filter):
        """
        Adds a transformation parameter to the last transformation.
        """
        self._structures.append(self._structures[-1])
        self._change_parameters.append({})
        self._changes.append(structure_filter)

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
        d = vasp_input_set.get_all_vasp_input(self._structures[-1],
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
        vasp_input_set.write_input(self._structures[-1], output_dir,
                                   make_dir_if_not_present=create_directory)
        with open(os.path.join(output_dir, "transformations.json"), "w") as fp:
            json.dump(self.to_dict, fp)

    def __str__(self):
        output = ["Current structure"]
        output.append("------------")
        output.append(str(self._structures[-1]))
        output.append("\nSource")
        output.append("------------")
        output.append(str(self._source))
        output.append("\nTransformation history")
        output.append("------------")
        for i, t in enumerate(self._changes):
            output.append("{} {}".format(t.to_dict,
                                         self._change_parameters[i]))
        output.append("\nOther parameters")
        output.append("------------")
        output.append(str(self._other_parameters))
        return "\n".join(output)

    def set_parameter(self, key, value):
        self._other_parameters[key] = value

    @property
    def other_parameters(self):
        return self._other_parameters

    @property
    def was_modified(self):
        """
        Boolean describing whether the last transformation on the structure
        made any alterations to it one example of when this would return false
        is in the case of performing a substitution transformation on the
        structure when the specie to replace isn't in the structure.
        """
        return not self._structures[-1] == self._structures[-2]

    @property
    def structures(self):
        """
        Returns a copy of all structures in the TransformedStructure. A
        structure is stored after every single transformation.
        """
        return [s for s in self._structures]

    @property
    def transformations(self):
        """
        Returns a copy of all transformations in the TransformedStructure.
        """
        return [t for t in self._changes]

    @property
    def final_structure(self):
        """
        Returns the final structure in the TransformedStructure.
        """
        return self._structures[-1]

    @staticmethod
    def from_dict(d):
        """
        Creates a TransformedStructure from a dict.
        """
        s = Structure.from_dict(d)
        return TransformedStructure(s, [], d["history"],
                                    d.get("other_parameters", None))

    @property
    def history(self):
        history = [self._source]
        for i, t in enumerate(self._changes):
            tdict = t.to_dict
            tdict["input_structure"] = self._structures[i].to_dict
            tdict["output_parameters"] = self._change_parameters[i]
            history.append(tdict)
        return history

    @property
    def to_dict(self):
        """
        Returns a dict representation of the TransformedStructure.
        """
        d = self._structures[-1].to_dict
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        d["history"] = self.history
        d["version"] = __version__
        d["last_modified"] = str(datetime.datetime.utcnow())
        d["other_parameters"] = self._other_parameters
        return d

    @staticmethod
    def from_cif_string(cif_string, transformations=[], primitive=True,
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
        return TransformedStructure(s, transformations, [source_info])

    @staticmethod
    def from_poscar_string(poscar_string, transformations=[]):
        """
        Generates TransformedStructure from a poscar string.

        Args:
            poscar_string:
                Input POSCAR string.
        """
        p = Poscar.from_string(poscar_string)
        if not p.true_names:
            raise ValueError("Transformation can be craeted only from POSCAR "
                             "strings with proper VASP5 element symbols.")
        raw_string = re.sub("'", "\"", poscar_string)
        s = p.structure
        source_info = {"source": "uploaded POSCAR",
                       "datetime": str(datetime.datetime.now()),
                       "original_file": raw_string}
        return TransformedStructure(s, transformations, [source_info])
