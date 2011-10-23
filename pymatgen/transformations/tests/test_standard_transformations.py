#!/usr/bin/env python

'''
Created on Sep 23, 2011
'''

__author__="Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Sep 23, 2011"

import unittest

from pymatgen.transformations.standard_transformations import transformation_from_json, IdentityTransformation, RotationTransformation
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure

class TransformationsTest(unittest.TestCase):
    
    def setUp(self):
        coords = list()
        coords.append([0,0,0])
        coords.append([0.75,0.5,0.75])
        lattice = Lattice([[ 3.8401979337, 0.00, 0.00],[1.9200989668, 3.3257101909, 0.00],[0.00,-2.2171384943,3.1355090603]])
        self.struct = Structure(lattice,["Si"] *2 ,coords)

    def test_identity_transformation(self):
        t = IdentityTransformation()
        self.assertEqual(self.struct, t.apply_transformation(self.struct))
    
    def test_rotation_transformation(self):
        t = RotationTransformation([0,1,0], 30, False)
        s2 = t.apply_transformation(self.struct)
        s1 = t.inverse.apply_transformation(s2)
        self.assertTrue((abs(s1.lattice.matrix - self.struct.lattice.matrix) < 1e-8).all())

class TransformationJsonTest(unittest.TestCase):
    
    def test_from_json(self):
        self.assertIsInstance(transformation_from_json('{"name": "IdentityTransformation", "init_args": {}}'), IdentityTransformation)
        self.assertIsInstance(transformation_from_json('{"name": "RotationTransformation", "init_args": {"angle": 30, "angle_in_radians": false, "axis": [0, 1, 0]}}'), RotationTransformation)
        
if __name__ == "__main__":
    unittest.main()