#!/usr/bin/env python

"""
This module sets up VASP calculations for a set of strained structures.
"""

from __future__ import division
import warnings
import sys
import unittest
import pymatgen
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.cif import CifWriter
from pymatgen.io.cif import CifParser
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.structure import IStructure
from pymatgen.transformations.standard_transformations import StrainStructureTransformation
from pymatgen.io.vasp.sets import MPStaticVaspInputSet
from pymatgen.io.vasp.sets import AbstractVaspInputSet
import numpy as np
import os
import subprocess

__author__ = "Cormac Toher"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "1.1"
__maintainer__ = "Cormac Toher"
__email__ = "cormac.toher@duke.edu"
__date__ = "June 7, 2013"
    


class ModifiedVolumeStructureSet(object):
    """
    class that generates a set of modified volume structures that
    can be used to calculate the energy-volume curve
    """

    def __init__(self, rlxd_str, nstructs=28, stepsize=0.01):
        """
        constructs the modified volumes of a structure.  Generates
        nstructs different volume structures, isotropically compressing
        and expanding lattice vectors in steps of stepsize.

        Args:
            rlxd_str (structure): structure to modify volume,
                should be a geometry optimized structure
            nstructs (int): number of different volume structures 
            stepsize: size of step by which lattice vectors should be extended/compressed
        """
        
        if nstructs % 2 != 0:
            mid = (nstructs + 1) / 2
        else:
            mid = nstructs / 2
  
        self.undeformed_structure = rlxd_str
        self.strainfactors = []
        self.modvol_structs = []

        for i in range(nstructs):
            strain = (i - mid) * stepsize
            strainfactor = strain + 1.0
            t = StrainStructureTransformation(strain)
            modified_structure = t.apply_transformation(rlxd_str)
            self.strainfactors.append(strainfactor)
            self.modvol_structs.append(modified_structure)


    def __iter__(self):
        return iter(self.modvol_structs)


