from __future__ import division
import warnings
import sys
sys.path.append('/home/MDEJONG1/pythonplayground/pymatgen/pymatgen_repo/pymatgen_repo') # (If one does not want to change $PYTHONPATH)
import unittest
import pymatgen
from pymatgen.io.vaspio import Poscar
from pymatgen.io.vaspio import Vasprun
from pymatgen.io.cifio import CifWriter
from pymatgen.io.cifio import CifParser
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.transformations.standard_transformations import *
from pymatgen.core.structure_modifier import StructureEditor
import numpy as np
import os

struct = CifParser('/home/MDEJONG1/pythonplayground/pymatgen/pymatgen_repo/pymatgen_repo/pymatgen/phonons/aluminum.cif').get_structures()[0]

#print struct.__dict__.keys()

#print len(struct._sites)
#print struct._sites[0]
#print struct._sites[1]

#print struct._lattice.__dict__.keys()

#print struct._lattice._mc2d
#print struct._lattice._matrix
#print struct._lattice._md2c

####print type(struct)
#print struct._lattice._matrix
#print struct._lattice._mc2d
####

s = StructureEditor(struct)
F = np.identity(3)
F[0,0] = 1.05
s.apply_strain_transformation(F) 

####print type(s.modified_structure)
#print s.modified_structure._lattice._matrix
#print s.modified_structure._lattice._mc2d
####



