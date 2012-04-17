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

mat = np.eye(3)
mat[2,1] = 1.000001
s = Strain(mat)
#print s.deformation_matrix


class TestStrain(unittest.TestCase):

    def setUp(self):
        self.matrix1 = np.matrix([[ 1.,   0.,   0.      ], [ 0.,   1.,   0.      ], [ 0.,   1.000001,  1.      ]])
        self.s = Strain(self.matrix1)

    def testreturn(self):
        self.assertEqual(self.matrix1.all(), self.s.deformation_matrix.all())

        
if __name__ == '__main__':
    unittest.main()



"""
class TestSequenceFunctions(unittest.TestCase):

    def setUp(self):
        self.seq = range(10)

    def testshuffle(self):
    # make sure the shuffled sequence does not lose any elements
        random.shuffle(self.seq)
        self.seq.sort()
        self.assertEqual(self.seq, range(10))

    def testchoice(self):
        element = random.choice(self.seq)
        self.assert_(element in self.seq)

    def testsample(self):
        self.assertRaises(ValueError, random.sample, self.seq, 20)
        for element in random.sample(self.seq, 5):
            self.assert_(element in self.seq)

"""


