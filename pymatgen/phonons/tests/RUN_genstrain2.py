from __future__ import division
import warnings
import sys
#sys.path.append('/home/MDEJONG1/pythonplayground/pymatgen/pymatgen_repo/pymatgen_repo') # (If one does not want to change $PYTHONPATH)
#sys.path.append('/home/MDEJONG1/pythonplayground/pymatgen/pymatgen_repo/pymatgen/')
sys.path.append('/home/MDEJONG1/pythonplayground/pymatgen/pymatgen_repo')
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

myCIF = CifParser('aluminum.cif').get_structures()[0]
		
#print type(myCIF)

#print myCIF


print myCIF.__dict__.keys()

print myCIF._sites
print myCIF._lattice

print myCIF._lattice.__dict__.keys()

print myCIF._lattice._matrix

print  myCIF._lattice._mc2d



#apply_strain_transformation()

#print myCIF._sites.__dict__.keys()

#s = StructureEditor(myCIF)
F = np.identity(3)
F[0,0] = 5
F = np.array(F)
#F = np.matrix(F)
#s.apply_strain_transformation(F)
#s.SupercellMaker(F)
#s =  SupercellMaker([[2,0,0],[0,3,0],[0,0,1]])


#print type(s.modified_structure)
#print s.__dict__.keys()
#print s._lattice
#print s._sites




