import warnings
import sys
sys.path.append('/home/MDEJONG1/pythonplayground/pymatgen/pymatgen_repo/pymatgen')
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

class TestStrain(unittest.TestCase):

    def setUp(self):
        self.F1 = np.matrix([[ 1.,   0.,   0.      ], [ 0.,   1.,   0.      ], [ 0.,   0.01,  1.      ]])
        self.E1 = 0.5*(np.transpose(self.F1)*self.F1 - np.identity(3))
        self.s1 = Strain(self.F1)
        
        self.F2 = np.matrix([[ 1.06,   0.,   0.      ], [ 0.,   1.,   0.      ], [ 0.,   0.,  1.      ]])
        self.E2 = 0.5*(np.transpose(self.F2)*self.F2 - np.identity(3))
        self.s2 = Strain(self.F2)

    def test_return_F1(self):
        self.assertEqual(np.matrix([[1,2],[3,4]]).tolist(), np.matrix([[1,2],[3,4]]).tolist())
        # print type(self.F1)
        #print type(self.s1.deformation_matrix)
        #self.assertEqual(self.F1,  np.matrix(self.s1.deformation_matrix))

    def test_return_E1(self):
        self.assertEqual(hash(tuple(self.E1)),  hash(tuple(self.s1.strain)))
     
    def test_return_F2(self):
        self.assertEqual(hash(tuple(self.F2)),  hash(tuple(self.s2.deformation_matrix)))

    def test_return_E2(self):
        self.assertEqual(hash(tuple(self.E2)),  hash(tuple(self.s2.strain)))
        
    def test_independent_strain(self):
        self.assertEqual(IndependentStrain(self.F1)._j, 1)
        self.assertEqual(IndependentStrain(self.F2)._i, 0)
        
        







        
if __name__ == '__main__':
    unittest.main()



