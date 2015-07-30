#!/usr/bin/env python

"""
This module sets up VASP calculations for a set of strained structures.
"""

from __future__ import division
import warnings
import sys
import unittest
import pymatgen
from pymatgen.io.vaspio import Poscar
from pymatgen.io.vaspio import Vasprun
from pymatgen.io.cifio import CifWriter
from pymatgen.io.cifio import CifParser
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.structure import IStructure
from pymatgen.transformations.standard_transformations import StrainStructureTransformation
from pymatgen.io.vaspio_set import MPStaticVaspInputSet
from pymatgen.io.vaspio_set import AbstractVaspInputSet
import numpy as np
import os
import subprocess

__author__ = "Cormac Toher"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "1.1"
__maintainer__ = "Cormac Toher"
__email__ = "cormac.toher@duke.edu"
__date__ = "June 7, 2013"
    
def Vasp_setup(inpfilename, nstructs=28, kfactor=1.0, stepsize=0.01):
    strfilename = str(inpfilename)
    initstruct.from_file(inpfilename)
    natoms = len(initstruct.sites)
    i = 0
    sysname = ""
    for i in range(natoms):
        sysname = sysname + initstruct.species[i].symbol
    print(sysname)
    modified_structs = ModifyVolume(initstruct, nstructs, stepsize)
    strdata = MPStaticVaspInputSet()
    for i, strainfactor in enumerate(modified_structs.keys()):
        vaspdir = './Vasp_run/AGL_Structure_' + str(strainfactor) 
        mod_struct = modified_structs[strainfactor]
        AbstractVaspInputSet.write_input(strdata, mod_struct, vaspdir, True)        
    for i, strainfactor in enumerate(modified_structs.keys()):
        vaspdir = './Vasp_run/AGL_Structure_' + str(strainfactor) 
        subprocess.call("sed '6 d' < " + vaspdir + '/POSCAR > ' + vaspdir + '/POSCAR.temp', shell=True)
        subprocess.call("rm -rf " + vaspdir + "/POSCAR", shell=True)
        subprocess.call("mv " + vaspdir + '/POSCAR.temp ' + vaspdir + '/POSCAR', shell=True)
        i = i + 1



def ModifyVolume(initstruct, nstructs=28, stepsize=0.01):
    i = 0
    mid = nstructs / 2
    mod_struct = {}
    while i < nstructs:
        strain = (i - mid) * stepsize
        strainfactor = strain + 1.0
        t = StrainStructureTransformation(strain)
        modified_structure = t.apply_transformation(initstruct)
        mod_struct[strainfactor] = modified_structure
        i = i + 1
    return mod_struct
        

if __name__ == "__main__":
    Vasp_setup(str(sys.argv[1]))
