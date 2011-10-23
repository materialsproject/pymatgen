#!/usr/bin/python

import unittest
from pymatgen.core.operations import SymmOp
import numpy as np
     
class  SymmOpTestCase(unittest.TestCase):
    
    def setUp(self):
        self.op = SymmOp.from_axis_angle_and_translation([0, 0, 1], 30, False, np.array([0,0,1]))
    
    def test_properties(self):
        rot = self.op.rotation_matrix
        vec = self.op.translation_vector
        self.assertTrue((abs(rot - np.array([[ 0.8660254, -0.5, 0.], [ 0.5, 0.8660254,  0.], [ 0., 0., 1.]])) < 0.01).all())
        self.assertTrue((abs(vec - np.array([0 ,0 ,1])) < 0.01).all())
    
    def test_operate(self):
        point = np.array([1,2,3])
        newcoord = self.op.operate(point)
        self.assertTrue((abs(newcoord - np.array([-0.1339746, 2.23205081, 4.])) < 0.01).all())
    
    def test_inverse(self): 
        point = np.random.rand(3)
        newcoord = self.op.operate(point)
        self.assertTrue((abs(self.op.inverse.operate(newcoord)- point) < 0.01).all())

    def test_apply_rotation_only(self):
        point = np.random.rand(3)
        newcoord = self.op.operate(point)
        rotate_only = self.op.apply_rotation_only(point)
        self.assertTrue((abs(rotate_only - newcoord + self.op.translation_vector) < 0.01).all())

    def test_are_symmetrically_related(self):
        point = np.random.rand(3)
        newcoord = self.op.operate(point)
        self.assertTrue(self.op.are_symmetrically_related(point, newcoord))
        self.assertTrue(self.op.are_symmetrically_related(newcoord,point))

    
if __name__ == '__main__':
    unittest.main()
