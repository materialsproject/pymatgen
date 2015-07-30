#!/usr/bin/env python

"""
This module reads the (E, V) data from the VASP output files.
"""

from __future__ import division
import warnings
import sys
import subprocess
import unittest
import pymatgen
from pymatgen.io.vaspio import Poscar
from pymatgen.io.vaspio import Vasprun
from pymatgen.io.cifio import CifWriter
from pymatgen.io.cifio import CifParser
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.structure import IStructure
from pymatgen.io.vaspio_set import MITVaspInputSet
from pymatgen.io.vaspio_set import AbstractVaspInputSet
from pymatgen.agl_thermal.agl_vasp_setup import ModifyVolume
import numpy as np
import os


__author__ = "Cormac Toher"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Cormac Toher"
__email__ = "cormac.toher@duke.edu"
__date__ = "June 7, 2013"
    
def Vasp_output_read(agl_data, initstruct, nstructs, stepsize):
    ang32bohr3 = 6.74833303710415
    modified_structs = ModifyVolume(initstruct, nstructs, stepsize)
    for i, strainfactor in enumerate(modified_structs.keys()):
        vaspdir = './Vasp_run/AGL_Structure_' + str(strainfactor) 
        agl_data.vol_inp.append(ang32bohr3 * modified_structs[strainfactor].volume)
        vaspfile = vaspdir + '/vasprun.xml'
        vaspres = Vasprun(vaspfile, None, False, False, False)
        agl_data.energ_inp.append(agl_data.ev2hartree * vaspres.final_energy)
    return


if __name__ == "__main__":
    Vasp_output_read()
