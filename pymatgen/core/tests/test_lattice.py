#!/usr/bin/python

from __future__ import division

import unittest
from pymatgen.core.lattice import Lattice
import numpy as np


class  LatticeTestCase(unittest.TestCase):

    def setUp(self):
        self.lattice = Lattice.cubic(10.0)
        self.tetragonal = Lattice.tetragonal(10, 20)

    def test_init(self):
        a = 9.026
        lattice = Lattice.cubic(a)
        self.assertIsNotNone(lattice, "Initialization from new_cubic failed")
        lattice2 = Lattice([[a, 0, 0], [0, a, 0], [0, 0, a]])
        for i in range(0, 3):
            for j in range(0, 3):
                self.assertAlmostEqual(lattice.matrix[i][j],
                                       lattice2.matrix[i][j], 5,
                                       "Inconsistent matrix from two inits!")
        primlatt = lattice.get_primitive_lattice('I')
        (lengths, angles) = primlatt.lengths_and_angles
        for i in range(0, 3):
            self.assertAlmostEqual(lengths[i], 7.81674529, 5,
                                   "Wrong primitive lattice obtained!")
            self.assertAlmostEqual(angles[i], 109.47122063, 5,
                                   "Wrong primitive lattice obtained!")
        coord = lattice.get_cartesian_coords(np.array([0.5, 0.5, 0.5]))
        prim_frac = primlatt.get_fractional_coords(coord)
        for i in range(0, 3):
            self.assertAlmostEqual(coord[i], 4.513, 5, "Wrong coord!")
            self.assertAlmostEqual(prim_frac[i], 1, 5,
                                   "Wrong primitive frac coord!")

    def test_get_cartesian_or_frac_coord(self):
        coord = self.lattice.get_cartesian_coords(np.array([0.15, 0.3, 0.4]))
        self.assertTrue(np.allclose(coord, [1.5, 3., 4.]))
        self.assertTrue(np.allclose(self.tetragonal
                                    .get_fractional_coords([12.12312, 45.2134,
                                                            1.3434]),
                                    [1.212312, 4.52134, 0.06717]))

        #Random testing that get_cart and get_frac coords reverses each other.
        rand_coord = np.random.random_sample(3)
        coord = self.tetragonal.get_cartesian_coords(rand_coord)
        fcoord = self.tetragonal.get_fractional_coords(coord)
        self.assertTrue(np.allclose(fcoord, rand_coord))

    def test_reciprocal_lattice(self):
        self.assertTrue(np.allclose(self.lattice.reciprocal_lattice.matrix,
                                    0.628319 * np.eye(3)))
        self.assertTrue(np.allclose(self.tetragonal.reciprocal_lattice.matrix,
                                    [[0.628319, 0., 0.], [0., 0.628319, 0],
                                     [0., 0., 0.3141590]]))

    def test_static_methods(self):
        lengths_c = [3.840198, 3.84019885, 3.8401976]
        angles_c = [119.99998575, 90, 60.00000728]
        mat_c = [[3.840198, 0.000000, 0.0000], [1.920099, 3.325710, 0.000000],
                 [0.000000, -2.217138, 3.135509]]
        #should give the lengths and angles above
        newlatt = Lattice(mat_c)
        (lengths, angles) = newlatt.lengths_and_angles
        for i in range(0, 3):
            self.assertAlmostEqual(lengths[i], lengths_c[i], 5,
                                   "Lengths incorrect!")
            self.assertAlmostEqual(angles[i], angles_c[i], 5,
                                   "Angles incorrect!")
        (lengths, angles) = \
            Lattice.from_lengths_and_angles(lengths, angles).lengths_and_angles
        for i in range(0, 3):
            self.assertAlmostEqual(lengths[i], lengths_c[i], 5,
                                   "Lengths incorrect!")
            self.assertAlmostEqual(angles[i], angles_c[i], 5,
                                   "Angles incorrect!")

    def test_attributes(self):
        """docstring for test_attributes"""
        lattice = Lattice.cubic(10.0)
        self.assertEqual(lattice.a, 10.0)
        self.assertEqual(lattice.b, 10.0)
        self.assertEqual(lattice.c, 10.0)
        self.assertAlmostEqual(lattice.volume, 1000.0)
        xyz = lattice.get_cartesian_coords([0.25, 0.35, 0.45])
        self.assertEqual(xyz[0], 2.5)
        self.assertEqual(xyz[1], 3.5)
        self.assertEqual(xyz[2], 4.5)

    def test_consistency(self):
        '''
        when only lengths and angles are given for constructors, the
        internal matrix representation is ambiguous since the lattice rotation
        is not specified.
        This test makes sure that a consistent definition is specified for the
        lattice rotation when using different constructors from lengths angles
        '''
        l = [3.840198, 3.84019885, 3.8401976]
        a = [119.99998575, 90, 60.00000728]
        mat1 = Lattice.from_lengths_and_angles(l, a).matrix
        mat2 = Lattice.from_parameters(l[0], l[1], l[2],
                                       a[0], a[1], a[2]).matrix
        for i in range(0, 3):
            for j in range(0, 3):
                self.assertAlmostEqual(mat1[i][j], mat2[i][j], 5)

    def test_get_lll_reduced_lattice(self):
        lattice = Lattice([1.0, 1, 1, -1.0, 0, 2, 3.0, 5, 6])
        reduced_latt = lattice.get_lll_reduced_lattice()

        expected_ans = np.array([0.0, 1.0, 0.0,
                                 1.0, 0.0, 1.0,
                                 - 2.0, 0.0, 1.0]).reshape((3, 3))
        self.assertTrue(np.allclose(reduced_latt.matrix, expected_ans))
        self.assertAlmostEqual(reduced_latt.volume, lattice.volume)
        latt = [7.164750, 2.481942, 0.000000,
                - 4.298850, 2.481942, 0.000000,
                0.000000, 0.000000, 14.253000]
        expected_ans = np.array([-4.298850, 2.481942, 0.000000,
                                 2.865900, 4.963884, 0.000000,
                                 0.000000, 0.000000, 14.253000])
        expected_ans = expected_ans.reshape((3, 3))
        reduced_latt = Lattice(latt).get_lll_reduced_lattice()
        self.assertTrue(np.allclose(reduced_latt.matrix,
                                    expected_ans))
        self.assertAlmostEqual(reduced_latt.volume, Lattice(latt).volume)

        expected_ans = np.array([0.0, 10.0, 10.0,
                                 10.0, 10.0, 0.0,
                                 30.0, -30.0, 40.0]).reshape((3, 3))

        lattice = np.array([100., 0., 10., 10., 10., 20., 10., 10., 10.])
        lattice = lattice.reshape(3, 3)
        lattice = Lattice(lattice.T)
        reduced_latt = lattice.get_lll_reduced_lattice()
        self.assertTrue(np.allclose(reduced_latt.matrix, expected_ans))
        self.assertAlmostEqual(reduced_latt.volume, lattice.volume)

        random_latt = Lattice(np.random.random((3, 3)))
        if np.linalg.det(random_latt.matrix) > 1e-8:
            reduced_random_latt = random_latt.get_lll_reduced_lattice()
            self.assertAlmostEqual(reduced_random_latt.volume,
                                   random_latt.volume)

    def test_get_niggli_reduced_lattice(self):
        latt = Lattice.from_parameters(3, 5.196, 2, 103 + 55 / 60,
                                       109 + 28 / 60,
                                       134 + 53 / 60)
        reduced_cell = latt.get_niggli_reduced_lattice()
        abc, angles = reduced_cell.lengths_and_angles
        self.assertAlmostEqual(abc[0], 2, 3)
        self.assertAlmostEqual(abc[1], 3, 3)
        self.assertAlmostEqual(abc[2], 3, 3)
        self.assertAlmostEqual(angles[0], 116.382855225, 3)
        self.assertAlmostEqual(angles[1], 94.769790287999996, 3)
        self.assertAlmostEqual(angles[2], 109.466666667, 3)

        mat = [[5.0, 0, 0], [0, 5.0, 0], [5.0, 0, 5.0]]
        latt = Lattice(np.dot([[1, 1, 1], [1, 1, 0], [0, 1, 1]], mat))
        reduced_cell = latt.get_niggli_reduced_lattice()
        abc, angles = reduced_cell.lengths_and_angles
        for l in abc:
            self.assertAlmostEqual(l, 5, 3)
        for a in angles:
            self.assertAlmostEqual(a, 90, 3)

        latt = Lattice([1.432950, 0.827314, 4.751000, -1.432950, 0.827314,
                        4.751000, 0.0, -1.654628, 4.751000])
        ans = [[-1.432950, -2.481942, 0.0],
               [-2.8659, 0.0, 0.0],
               [-1.432950, -0.827314, -4.751000]]
        self.assertTrue(np.allclose(latt.get_niggli_reduced_lattice().matrix,
                                    ans))

    def test_to_from_dict(self):
        d = self.tetragonal.to_dict
        t = Lattice.from_dict(d)
        for i in range(3):
            self.assertEqual(t.abc[i], self.tetragonal.abc[i])
            self.assertEqual(t.angles[i], self.tetragonal.angles[i])
        #Make sure old style dicts work.
        del d["matrix"]
        t = Lattice.from_dict(d)
        for i in range(3):
            self.assertEqual(t.abc[i], self.tetragonal.abc[i])
            self.assertEqual(t.angles[i], self.tetragonal.angles[i])

if __name__ == '__main__':
    unittest.main()
