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

__author__="Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Mar 4, 2012"

import os
import re
import collections

from pymatgen.alchemy.materials import CifTransformedStructure, PoscarTransformedStructure


class CifTransmuter(object):
    """
    The CifTransmuter generates new structures from cifs that may contain
    multiple structures.
    """
    
    def __init__(self, transformations):
        """
        Args:
            transformations:
                List of standard_transformations objects.
        """
        self._transformations = transformations
        
    def transmute(self, cif_filenames):
        """
        Args:
            cif_filenames:
                Sequence of cif filenames. Each cif may or may not contain multiple 
                structures.
        """
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
                transformed_structures.extend([CifTransformedStructure("".join(data), self._transformations) for data in structure_data])
        return transformed_structures

class PoscarTransmuter(object):
    """
    The PoscarTransmuter generates new structures from a sequence of Poscar file names.
    """
    
    def __init__(self, transformations):
        """
        Args:
            transformations:
                List of standard_transformations objects.
        """
        self._transformations = transformations
        
    def transmute(self, poscar_filenames):
        """
        Args:
            poscar_filenames:
                Sequence of poscar filenames.
        """
        transformed_structures = []
        for filename in poscar_filenames:
            with open(filename, "r") as f:
                transformed_structures.append(PoscarTransformedStructure(f.read(), self._transformations))
        return transformed_structures
        

def batch_write_vasp_input(transformed_structures, vasp_input_set, output_dir, create_directory = True):
    """
    Batch write vasp input for a sequence of transformed structures to output_dir,
    following the format output_dir/{formula}_{number}.
    
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
    """
    dnames_count = collections.defaultdict(int)
    for s in transformed_structures:
        formula = re.sub("\s+", "", s.final_structure.formula)
        dirname = os.path.join(output_dir, '{}_{}'.format(formula, dnames_count[formula] + 1))
        s.write_vasp_input(vasp_input_set, dirname, create_directory = True)
        dnames_count[formula] += 1
        