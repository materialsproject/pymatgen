import warnings
import sys
sys.path.append('/home/MDEJONG1/pythonplayground/pymatgen/pymatgen_repo/') # (If one does not want to change $PYTHONPATH)
import unittest
import pymatgen
from pymatgen.io.vaspio import Poscar
from pymatgen.io.vaspio import Poscar
from pymatgen.io.vaspio import Vasprun
from pymatgen.io.cifio import CifWriter
from pymatgen.io.cifio import CifParser
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.transformations.standard_transformations import *
from pymatgen.core.structure_modifier import StructureEditor
import numpy as np
import random
from pymatgen.phonons.strain import Strain
from pymatgen.phonons.strain import IndependentStrain


for i in range(0, 24):

    path = '/home/MDEJONG1/Re/elasticity/Cij_redo/Re_Ta/M' + str(i) + '/vasprun.xml'
    A = Vasprun(path)
    sigma =  A.ionic_steps[-1]['stress']
    print sigma[0,0], sigma[0,1], sigma[0,2]
    print sigma[1,0], sigma[1,1], sigma[1,2]
    print sigma[2,0], sigma[2,1], sigma[2,2]


