#!/usr/bin/env python

'''
This module implements various transmuter classes.
Transmuters are essentially classes that generate TransformedStructures from
various data sources. They enable the high-throughput generation of new
structures and input files.

It also includes the helper function, batch_write_vasp_input to generate an 
entire directory of vasp input files for running.
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


class StandardTransmuter(object):
    """
    An example of a Transmuter object, which performs a sequence of
    transformations on many structures to generate TransformedStructures.
    """

    def __init__(self, transformed_structures, transformations=[],
                 extend_collection=0):
        """
        Args:
            transformed_structures:
                Input transformed structures
            transformations:
                New transformations to be applied to all structures
            extend_collection:
                Whether to use more than one output structure from one-to-many
                transformations. extend_collection can be a number, which
                determines the maximum branching for each transformation.
        """
        self._transformed_structures = transformed_structures
        for trans in transformations:
            self.append_transformation(trans, extend_collection=extend_collection)

    def get_transformed_structures(self):
        """
        Returns all TransformedStructures.
        """
        return self._transformed_structures

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

    def append_transformation(self, transformation, extend_collection=False,
                              clear_redo=True):
        """
        Appends a transformation to all TransformedStructures.
        
        Args:
            transformation:
                Transformation to append
            extend_collection:
                Whether to use more than one output structure from one-to-many
                transformations. extend_collection can be a number, which
                determines the maximum branching for each transformation.
            clear_redo:
                Boolean indicating whether to clear the redo list. By default,
                this is True, meaning any appends clears the history of undoing.
                However, when using append_transformation to do a redo, the 
                redo list should not be cleared to allow multiple redos.
                
        Returns:
            List of booleans corresponding to initial transformed structures
            each boolean describes whether the transformation altered the 
            structure
        """
        new_structures = []

        for x in self._transformed_structures:
            new = x.append_transformation(transformation, return_alternatives=extend_collection)
            if new is not None:
                new_structures.extend(new)

        self._transformed_structures.extend(new_structures)

    def extend_transformations(self, transformations):
        """
        Extends a sequence of transformations to the TransformedStructure.
        
        Args:
            transformations:
                Sequence of Transformations
        """
        for t in transformations:
            self.append_transformation(t)

    def write_vasp_input(self, vasp_input_set, output_dir,
                         create_directory=True, subfolder=None):
        """
        Batch write vasp input for a sequence of transformed structures to 
        output_dir, following the format output_dir/{formula}_{number}.
        
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

    def append_transformed_structures(self, tstructs_or_transmuter):
        '''
        Method is overloaded to accept either a list of transformed structures
        or transmuter, it which case it appends the 2nd transmuter's structures.
        
        Args:
            tstructs_or_transmuter:
                A list of transformed structures or a transmuter.
        '''
        if isinstance(tstructs_or_transmuter, self.__class__):
            self._transformed_structures.extend(tstructs_or_transmuter._transformed_structures)
        else:
            for ts in tstructs_or_transmuter:
                assert isinstance(ts, TransformedStructure)
            self._transformed_structures.extend(tstructs_or_transmuter)


class CifTransmuter(StandardTransmuter):

    def __init__(self, cif_string, transformations=[], primitive=True,
                  extend_collection=False):
        '''
        Generates a Transmuter from a cif string, possibly
        containing multiple structures.
        
        Args:
            cif_string:
                A string containing a cif or a series of cifs
            transformations:
                New transformations to be applied to all structures
            primitive:
                Whether to generate the primitive cell from the cif.
            extend_collection:
                Whether to use more than one output structure from one-to-many
                transformations.
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
        transformed_structures.extend([TransformedStructure.from_cif_string("\n".join(data), [], primitive) for data in structure_data])
        StandardTransmuter.__init__(self, transformed_structures, transformations, extend_collection)

    @staticmethod
    def from_filenames(filenames, transformations=[], primitive=True,
                  extend_collection=False):
        '''
        Generates a TransformedStructureCollection from a cif, possibly
        containing multiple structures.
        
        Args:
            filenames:
                List of strings of the cif files
            transformations:
                New transformations to be applied to all structures
            primitive:
                Whether to generate the primitive cell from the cif.
            extend_collection:
                Whether to use more than one output structure from one-to-many
                transformations.
        '''

        allcifs = []
        for fname in filenames:
            with open(fname, "r") as f:
                allcifs.append(f.read())
        return CifTransmuter("\n".join(allcifs), transformations, primitive=primitive, extend_collection=extend_collection)


class PoscarTransmuter(StandardTransmuter):

    def __init__(self, poscar_string, transformations=[],
                 extend_collection=False):
        """
        Generates a transmuter from a sequence of POSCARs.
        
        Args:
            poscar_string:
                List of POSCAR strings
            transformations:
                New transformations to be applied to all structures.
            primitive:
                Whether to generate the primitive cell from the cif.
            extend_collection:
                Whether to use more than one output structure from one-to-many
                transformations.
        """
        transformed_structures = []
        transformed_structures.append(TransformedStructure.from_poscar_string(poscar_string, []))
        StandardTransmuter.__init__(self, transformed_structures, transformations, extend_collection=extend_collection)

    @staticmethod
    def from_filenames(poscar_filenames, transformations=[],
                       extend_collection=False):
        transformed_structures = []
        for filename in poscar_filenames:
            with open(filename, "r") as f:
                transformed_structures.append(TransformedStructure.from_poscar_string(f.read(), []))
        return StandardTransmuter(transformed_structures, transformations, extend_collection=extend_collection)


def batch_write_vasp_input(transformed_structures, vasp_input_set, output_dir, create_directory=True, subfolder=None):
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
            dirname = os.path.join(output_dir, subdir, '{}_{}'.format(formula, dnames_count[subdir + formula] + 1))
        else:
            dirname = os.path.join(output_dir, '{}_{}'.format(formula, dnames_count[formula] + 1))
        s.write_vasp_input(vasp_input_set, dirname, create_directory=True)
        dnames_count[formula] += 1
