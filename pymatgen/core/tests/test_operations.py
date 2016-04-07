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
        self.assertArrayAlmostEqual(
            newcoords, [[-0.1339746, 2.23205081, 4.]] * 2, 2)
        newcoords = self.op.operate_multi([[point, point]] * 2)
        self.assertArrayAlmostEqual(
            newcoords, [[[-0.1339746, 2.23205081, 4.]] * 2] * 2, 2)

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
        # Distance to the plane should be negatives of each other.
        self.assertAlmostEqual(np.dot(newcoord - origin, normal),
                               -np.dot(point - origin, normal))

    def test_apply_rotation_only(self):
        point = np.random.rand(3)
        newcoord = self.op.operate(point)
        rotate_only = self.op.apply_rotation_only(point)
        self.assertArrayAlmostEqual(
            rotate_only + self.op.translation_vector, newcoord, 2)

    def test_transform_r2_tensor(self):
        tensor = np.arange(0, 9.).reshape(3, 3)
        new_tensor = self.op.transform_r2_tensor(tensor)
        self.assertArrayAlmostEqual(new_tensor,
                                    [[2.7320508,  1.73205079,  4.2320508],
                                     [3.73205078,  1.26794917,  3.330127],
                                     [8.6961524,  3.0621778,  8.]], 5)

    def test_transform_r3_tensor(self):
        tensor=np.arange(0, 27.).reshape(3, 3, 3)
        new_tensor=self.op.transform_r3_tensor(tensor)
        print(tensor)
        print(new_tensor)
        self.assertArrayAlmostEqual(new_tensor,
                                    [[[12.12916506,   4.61602535,  11.92820319],
                                     [7.34807613,   2.33493645,   6.19615237],
                                     [18.02627936,   5.83012695,  15.4282032]],

                                    [[15.54422848,   4.53108883,  12.19615233],
                                     [5.26313963,   1.50833021,   4.07179671],
                                     [13.8301269,   3.9737205,  10.7224318]],

                                    [[36.32050788,  10.73205068,  28.820508],
                                     [12.73205066,   3.67949186,   9.9185842],
                                     [33.2846096,   9.650635,  26.]]], 3)

    def test_transform_r4_tensor(self):
        tensor=np.arange(0, 81).reshape(3, 3, 3, 3)
        new_tensor=self.op.transform_r4_tensor(tensor)
        self.assertArrayAlmostEqual(new_tensor,
                                    [[[[50.98076169,   15.5262792,   41.48557134],
                                    [19.25832996,    5.66025391,   15.21410143],
                                    [49.81569829,   14.71410142,   39.51666035]],

                                    [[30.45448225,    8.66025385,   23.41025377],
                                    [9.66025383,    2.72243179,    7.37083474],
                                    [25.64230454,    7.23686014,   19.58845709]],

                                    [[74.80607912,   21.41025373,   57.81088887],
                                    [24.14230451,    6.83493635,    18.49038085],
                                    [63.90896504,   18.12435544,   49.0166604]]],


                                    [[[64.04293911,   17.66025368,   47.99871081],
                                    [18.66025366,    5.13397445,   13.9592919],
                                    [50.23076158,   13.8253173,   37.58845697]],

                                    [[21.6602536,    5.93782201,   16.15544428],
                                    [6.20577119,    1.69872975,    4.62306684],
                                    [16.75352048,    4.58716846,   12.48333931]],

                                    [[56.9269139,   15.61954589,   42.49038069],
                                    [16.35159668,    4.4794733,   12.18911067],
                                    [44.12435527,   12.09103446,   32.8993462]]],


                                    [[[149.77722162,   41.49871066,  112.69357443],
                                    [44.23076145,   12.21762212,   33.19615201],
                                    [118.79165061,   32.8301266,   89.1935748]],

                                    [[52.42691379,   14.4137745,   39.19615197],
                                    [15.1458253,    4.15638784,   11.30642475],
                                    [40.83012655,   11.20834854,   30.4878034]],

                                    [[137.08587913,   37.73205032,  102.5858796],
                                    [39.73205031,   10.9141199,   29.6839558],
                                    [107.0499812,   29.4160066,   80.]]]], 5)

    def test_are_symmetrically_related(self):
        point=np.random.rand(3)
        newcoord=self.op.operate(point)
        self.assertTrue(self.op.are_symmetrically_related(point, newcoord))
        self.assertTrue(self.op.are_symmetrically_related(newcoord, point))

    def test_to_from_dict(self):
        d=self.op.as_dict()
        op=SymmOp.from_dict(d)
        point=np.random.rand(3)
        newcoord=self.op.operate(point)
        self.assertTrue(op.are_symmetrically_related(point, newcoord))

    def test_inversion(self):
        origin=np.random.rand(3)
        op=SymmOp.inversion(origin)
        pt=np.random.rand(3)
        inv_pt=op.operate(pt)
        self.assertArrayAlmostEqual(pt - origin, origin - inv_pt)

    def test_xyz(self):
        op=SymmOp([[1, -1, 0, 0], [0, -1, 0, 0],
                     [0, 0, -1, 0], [0, 0, 0, 1]])
        s=op.as_xyz_string()
        self.assertEqual(s, 'x-y, -y, -z')
        self.assertEqual(op, SymmOp.from_xyz_string(s))

        op2=SymmOp([[0, -1, 0, 0.5], [1, 0, 0, 0.5],
                      [0, 0, 1, 0.5 + 1e-7], [0, 0, 0, 1]])
        s2=op2.as_xyz_string()
        self.assertEqual(s2, '-y+1/2, x+1/2, z+1/2')
        self.assertEqual(op2, SymmOp.from_xyz_string(s2))

        op2=SymmOp([[3, -2, -1, 0.5], [-1, 0, 0, 12. / 13],
                      [0, 0, 1, 0.5 + 1e-7], [0, 0, 0, 1]])
        s2=op2.as_xyz_string()
        self.assertEqual(s2, '3x-2y-z+1/2, -x+12/13, z+1/2')
        self.assertEqual(op2, SymmOp.from_xyz_string(s2))

        op3=SymmOp.from_xyz_string('3x - 2y - z+1 /2 , -x+12/ 13, z+1/2')
        self.assertEqual(op2, op3)

        # Ensure strings can be read in any order
        op4=SymmOp.from_xyz_string('1 /2 + 3X - 2y - z , 12/ 13-x, z+1/2')
        op5=SymmOp.from_xyz_string('+1 /2 + 3x - 2y - z , 12/ 13-x, +1/2+z')
        self.assertEqual(op4, op3)
        self.assertEqual(op4, op5)
        self.assertEqual(op3, op5)

        self.assertRaises(ValueError, self.op.as_xyz_string)

        o=SymmOp.from_xyz_string('0.5+x, 0.25+y, 0.75+z')
        self.assertArrayAlmostEqual(o.translation_vector, [0.5, 0.25, 0.75])
        o=SymmOp.from_xyz_string('x + 0.5, y + 0.25, z + 0.75')
        self.assertArrayAlmostEqual(o.translation_vector, [0.5, 0.25, 0.75])

if __name__ == '__main__':
    import unittest
    unittest.main()
