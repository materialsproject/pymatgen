#!/usr/bin/python

from pymatgen.util.testing import PymatgenTest
from pymatgen.core.operations import SymmOp
import numpy as np


class SymmOpTestCase(PymatgenTest):

    def setUp(self):
        self.op = SymmOp.from_axis_angle_and_translation([0, 0, 1], 30, False,
                                                         [0, 0, 1])

    def test_properties(self):
        rot = self.op.rotation_matrix
        vec = self.op.translation_vector
        self.assertArrayAlmostEqual(rot, [[0.8660254, -0.5, 0.],
                                          [0.5, 0.8660254, 0.],
                                          [0., 0., 1.]], 2)
        self.assertArrayAlmostEqual(vec, [0, 0, 1], 2)

    def test_operate(self):
        point = np.array([1, 2, 3])
        newcoord = self.op.operate(point)
        self.assertArrayAlmostEqual(newcoord, [-0.1339746, 2.23205081, 4.], 2)

    def test_inverse(self):
        point = np.random.rand(3)
        newcoord = self.op.operate(point)
        self.assertArrayAlmostEqual(self.op.inverse.operate(newcoord),
                                    point, 2)

    def test_reflection(self):
        normal = np.random.rand(3)
        origin = np.random.rand(3)
        refl = SymmOp.reflection(normal, origin)
        point = np.random.rand(3)
        newcoord = refl.operate(point)
        #Distance to the plane should be negatives of each other.
        self.assertAlmostEqual(np.dot(newcoord - origin, normal),
                               -np.dot(point - origin, normal))


    def test_apply_rotation_only(self):
        point = np.random.rand(3)
        newcoord = self.op.operate(point)
        rotate_only = self.op.apply_rotation_only(point)
        self.assertArrayAlmostEqual(
            rotate_only + self.op.translation_vector, newcoord, 2)

    def test_are_symmetrically_related(self):
        point = np.random.rand(3)
        newcoord = self.op.operate(point)
        self.assertTrue(self.op.are_symmetrically_related(point, newcoord))
        self.assertTrue(self.op.are_symmetrically_related(newcoord, point))

    def test_to_from_dict(self):
        d = self.op.to_dict
        op = SymmOp.from_dict(d)
        point = np.random.rand(3)
        newcoord = self.op.operate(point)
        self.assertTrue(op.are_symmetrically_related(point, newcoord))

    def test_inversion(self):
        origin = np.random.rand(3)
        op = SymmOp.inversion(origin)
        pt = np.random.rand(3)
        inv_pt = op.operate(pt)
        self.assertArrayAlmostEqual(pt - origin, origin - inv_pt)

if __name__ == '__main__':
    import unittest
    unittest.main()
