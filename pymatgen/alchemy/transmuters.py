#!/usr/bin/env python

'''
This module implements various transmuter classes.
Transmuters are essentially classes that generate TransformedStructures from
various data sources. They enable the high-throughput generation of new
structures and input files.

It also includes the helper function, batch_write_vasp_input to generate an entire
directory of vasp input files for running.
'''

from __future__ import division

__author__ = "Shyue Ping Ong, Will Richards"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Mar 4, 2012"

import os
import re
import collections

from pymatgen.alchemy.materials import TransformedStructure

from copy import deepcopy

class TransformedStructureTransmuter(object):

    def __init__(self, transformed_structures, transformations = [], extend_collection = False):
        self._transformed_structures = transformed_structures
        for trans in transformations:
            self.append_transformation(trans, extend_collection = extend_collection)

    def __getitem__(self, index):
        return self._transformed_structures[index]

    def __getattr__(self, name):
        return [getattr(x, name) for x in self._transformed_structures]

    def undo_last_transformation(self):
        """
        Undo the last transformation in the TransformedStructure.
        
        Raises:
            IndexError if already at the oldest change.
        """
        for x in self._transformed_structures:
            x.undo_last_transformation()

    def redo_next_transformation(self):
        """
        Redo the last undone transformation in the TransformedStructure.
        
        Raises:
            IndexError if already at the latest change.
        """
        for x in self._transformed_structures:
            x.redo_next_transformation()

    def __len__(self):
        return len(self._transformed_structures)

    def append_transformation(self, transformation, extend_collection = False, clear_redo = True):
        """
        TODO: clean this up a lot
        
        Appends a transformation to the TransformedStructure.
        
        Arguments:
            transformation:
                Transformation to append
            clear_redo:
                Boolean indicating whether to clear the redo list. By default,
                this is True, meaning any appends clears the history of undoing.
                However, when using append_transformation to do a redo, the redo
                list should not be cleared to allow multiple redos.
                
        Returns:
            list of booleans corresponding to initial transformed structures
            each boolean describes whether the transformation altered the structure
        """
        new_structures = []

        for x in self._transformed_structures:
            new = x.append_transformation(transformation, return_alternatives = extend_collection, clear_redo = clear_redo)
            if new is not None:
                new_structures.extend(new)
        output = [x.was_modified for x in self._transformed_structures]
        self._transformed_structures.extend(new_structures)
        return output

    def branch_collection(self, transformations, retention_level = 1, extend_collection = False, clear_redo = True):
        '''
        Copies the structures collection, applying one transformation to each copy.
        
        Args:
            transformations:
                List of transformations to apply (each structure gets one of 
                these transformations. To append multiple transformations to 
                each structure use extend_transformations)
            retention_level:
                Specifies which structures will be kept and which will be thrown out
                0 - throws out all structures that weren't modified by any of 
                the transformations
                1 - keeps structures that weren't modified by anything.
                2 - keeps all structures, including the untransformed ones. Note 
                that this may cause issues with undoing transformations
                since they will have different transformation histories
                
                e.g if you start with 2 structures and apply 2 transformations, 
                and one structure isn't modified by either of them but the other 
                structure is modified by both, for retention_level = 0, you will have 
                2 structures left, for retention_level = 1 you will have 3 structures, 
                and for retention_level 2, you will have 4. In most cases,
                retention_level = 1 will provide the desired functionality.
        '''
        any_modification = [False] * len(self._transformed_structures)
        old_transformed_structures = self._transformed_structures

        new_trans_structures = []
        for transformation in transformations:
            self._transformed_structures = deepcopy(old_transformed_structures)
            modified = self.append_transformation(transformation, extend_collection, clear_redo)
            for structure in self._transformed_structures:
                if structure.was_modified:
                    new_trans_structures.append(structure)
            any_modification = map(lambda x:x[0] or x[1], zip(modified, any_modification))

        self._transformed_structures = new_trans_structures

        if retention_level == 1:
            for i, structure in enumerate(old_transformed_structures):
                if not any_modification[i]:
                    self._transformed_structures.append(structure)
        elif retention_level == 2:
            self._transformed_structures.extend(old_transformed_structures)


    def extend_transformations(self, transformations):
        """
        Extends a sequence of transformations to the TransformedStructure.
        
        Args:
            transformations:
                Sequence of Transformations
        """
        for t in transformations:
            self.append_transformation(t)

    def write_vasp_input(self, vasp_input_set, output_dir, create_directory = True, subfolder = None):
        """
        Batch write vasp input for a sequence of transformed structures to output_dir,
        following the format output_dir/{formula}_{number}.
        
        Args:
            vasp_input_set:
                pymatgen.io.vaspio_set.VaspInputSet like object that creates
                vasp input files from structures
            output_dir:
                Directory to output files
            create_directory:
                Create the directory if not present. Defaults to True.
            subfolder:
                function to create subdirectory name from transformed_structure.
                eg. lambda x: x.other_parameters['tags'][0] to use the first tag
        """
        batch_write_vasp_input(self._transformed_structures, vasp_input_set, output_dir, create_directory, subfolder)

    def set_parameter(self, key, value):
        for x in self._transformed_structures:
            x.set_parameter(key, value)

    def __str__(self):
        output = ["Current structures"]
        output.append("------------")
        for x in self._transformed_structures:
            output.append(str(x._structures[-1]))
        return "\n".join(output)

    def remove_duplicates(self):
        '''
        TODO: write this method
        '''
        pass

    @staticmethod
    def from_cif_string(cif_string, transformations = [], primitive = True, extend_collection = False):
        '''
        Generates a TransformedStructureCollection from a cif string, possibly
        containing multiple structures.
        
        Args:
            cif_filenames:
                List of strings of the cif files
        '''
        transformed_structures = []
        lines = cif_string.split("\n")
        structure_data = []
        read_data = False
        for line in lines:
            if re.match("^\s*data", line):
                structure_data.append([])
                read_data = True
            if read_data:
                structure_data[-1].append(line)
        transformed_structures.extend([TransformedStructure.from_cif_string("".join(data), [], primitive) for data in structure_data])
        return TransformedStructureTransmuter(transformed_structures, transformations, extend_collection)


    @staticmethod
    def from_cifs(cif_filenames, transformations = [], primitive = True):
        '''
        Generates a TransformedStructureCollection from a cif, possibly
        containing multiple structures.
        
        Args:
            cif_filenames:
                List of strings of the cif files
        '''
        transformed_structures = []
        for filename in cif_filenames:
            with open(filename, "r") as f:
                structure_data = []
                read_data = False
                for line in f:
                    if re.match("^\s*data", line):
                        structure_data.append([])
                        read_data = True
                    if read_data:
                        structure_data[-1].append(line)
                transformed_structures.extend([TransformedStructure.from_cif_string("".join(data), transformations, primitive) for data in structure_data])
        return TransformedStructureTransmuter(transformed_structures, [])

    @staticmethod
    def from_poscars(poscar_filenames, transformations = []):
        transformed_structures = []
        for filename in poscar_filenames:
            with open(filename, "r") as f:
                transformed_structures.append(TransformedStructure.from_poscar_string(f.read(), transformations))
        return TransformedStructureTransmuter(transformed_structures, [])


def batch_write_vasp_input(transformed_structures, vasp_input_set, output_dir, create_directory = True, subfolder = None):
    """
    Batch write vasp input for a sequence of transformed structures to output_dir,
    following the format output_dir/{group}/{formula}_{number}.
    
    Args:
        transformed_structures:
            Sequence of TransformedStructures.
        vasp_input_set:
            pymatgen.io.vaspio_set.VaspInputSet like object that creates
            vasp input files from structures
        output_dir:
            Directory to output files
        create_directory:
            Create the directory if not present. Defaults to True.
        subfolder:
            function to create subdirectory name from transformed_structure.
            eg. lambda x: x.other_parameters['tags'][0] to use the first tag
    """
    dnames_count = collections.defaultdict(int)
    for s in transformed_structures:
        formula = re.sub("\s+", "", s.final_structure.formula)
        if subfolder is not None:
            subdir = subfolder(s)
            dirname = os.path.join(output_dir, subdir, '{}_{}'.format(formula, dnames_count[subdir+formula] + 1))
        else:
            dirname = os.path.join(output_dir, '{}_{}'.format(formula, dnames_count[formula] + 1))
        s.write_vasp_input(vasp_input_set, dirname, create_directory = True)
        dnames_count[formula] += 1
