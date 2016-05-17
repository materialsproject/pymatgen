# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

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

    def test_operate_multi(self):
        point = np.array([1, 2, 3])
        newcoords = self.op.operate_multi([point, point])
        self.assertArrayAlmostEqual(newcoords, [[-0.1339746, 2.23205081, 4.]]*2, 2)
        newcoords = self.op.operate_multi([[point, point]]*2)
        self.assertArrayAlmostEqual(newcoords, [[[-0.1339746, 2.23205081, 4.]]*2]*2, 2)

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
        d = self.op.as_dict()
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

    def test_xyz(self):
        op = SymmOp([[1, -1, 0, 0], [0, -1, 0, 0],
                     [0, 0, -1, 0], [0, 0, 0, 1]])
        s = op.as_xyz_string()
        self.assertEqual(s, 'x-y, -y, -z')
        self.assertEqual(op, SymmOp.from_xyz_string(s))

        op2 = SymmOp([[0, -1, 0, 0.5], [1, 0, 0, 0.5],
                      [0, 0, 1, 0.5+1e-7], [0, 0, 0, 1]])
        s2 = op2.as_xyz_string()
        self.assertEqual(s2, '-y+1/2, x+1/2, z+1/2')
        self.assertEqual(op2, SymmOp.from_xyz_string(s2))

        op2 = SymmOp([[3, -2, -1, 0.5], [-1, 0, 0, 12./13],
                      [0, 0, 1, 0.5+1e-7], [0, 0, 0, 1]])
        s2 = op2.as_xyz_string()
        self.assertEqual(s2, '3x-2y-z+1/2, -x+12/13, z+1/2')
        self.assertEqual(op2, SymmOp.from_xyz_string(s2))

        op3 = SymmOp.from_xyz_string('3x - 2y - z+1 /2 , -x+12/ 13, z+1/2')
        self.assertEqual(op2, op3)

        # Ensure strings can be read in any order
        op4 = SymmOp.from_xyz_string('1 /2 + 3X - 2y - z , 12/ 13-x, z+1/2')
        op5 = SymmOp.from_xyz_string('+1 /2 + 3x - 2y - z , 12/ 13-x, +1/2+z')
        self.assertEqual(op4, op3)
        self.assertEqual(op4, op5)
        self.assertEqual(op3, op5)

        self.assertRaises(ValueError, self.op.as_xyz_string)

        o = SymmOp.from_xyz_string('0.5+x, 0.25+y, 0.75+z')
        self.assertArrayAlmostEqual(o.translation_vector, [0.5, 0.25, 0.75])
        o = SymmOp.from_xyz_string('x + 0.5, y + 0.25, z + 0.75')
        self.assertArrayAlmostEqual(o.translation_vector, [0.5, 0.25, 0.75])

if __name__ == '__main__':
    import unittest
    unittest.main()
