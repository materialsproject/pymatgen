import warnings
import sys
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

class StrainTest(unittest.TestCase):

    def setUp(self):
        self.strain = Strain([[ 1., 0., 0.], 
                              [ 0., 1., 0.], 
                              [ 0., 1.000001, 1.]])

    def testreturn(self):
        self.assertEqual(self.matrix1.all(), self.s.deformation_matrix.all())
        
if __name__ == '__main__':
    unittest.main()
