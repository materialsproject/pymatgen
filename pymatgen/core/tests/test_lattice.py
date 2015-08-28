# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import itertools
from pymatgen.core.lattice import Lattice
import numpy as np
from pymatgen.util.testing import PymatgenTest
from pymatgen.core.operations import SymmOp


class LatticeTestCase(PymatgenTest):

    def setUp(self):
        self.lattice = Lattice.cubic(10.0)
        self.cubic = self.lattice
        self.tetragonal = Lattice.tetragonal(10, 20)
        self.orthorhombic = Lattice.orthorhombic(10, 20, 30)
        self.monoclinic = Lattice.monoclinic(10, 20, 30, 66)
        self.hexagonal = Lattice.hexagonal(10, 20)
        self.rhombohedral = Lattice.rhombohedral(10, 77)

        family_names = ["cubic", "tetragonal", "orthorhombic", "monoclinic",
                        "hexagonal", "rhombohedral"]

        self.families = {}
        for name in family_names:
            self.families[name] = getattr(self, name)

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

    def test_copy(self):
        cubic_copy = self.cubic.copy()
        self.assertTrue(cubic_copy == self.cubic)
        self.assertFalse(cubic_copy._matrix is self.cubic._matrix)

    def test_get_cartesian_or_frac_coord(self):
        coord = self.lattice.get_cartesian_coords([0.15, 0.3, 0.4])
        self.assertArrayAlmostEqual(coord, [1.5, 3., 4.])
        self.assertArrayAlmostEqual(
            self.tetragonal.get_fractional_coords([12.12312, 45.2134,
                                                   1.3434]),
            [1.212312, 4.52134, 0.06717])

        #Random testing that get_cart and get_frac coords reverses each other.
        rand_coord = np.random.random_sample(3)
        coord = self.tetragonal.get_cartesian_coords(rand_coord)
        fcoord = self.tetragonal.get_fractional_coords(coord)
        self.assertArrayAlmostEqual(fcoord, rand_coord)

    def test_reciprocal_lattice(self):
        recip_latt = self.lattice.reciprocal_lattice
        self.assertArrayAlmostEqual(recip_latt.matrix,
                                    0.628319 * np.eye(3), 5)
        self.assertArrayAlmostEqual(self.tetragonal.reciprocal_lattice.matrix,
                                    [[0.628319, 0., 0.], [0., 0.628319, 0],
                                     [0., 0., 0.3141590]], 5)

        #Test the crystallographic version.
        recip_latt_xtal = self.lattice.reciprocal_lattice_crystallographic
        self.assertArrayAlmostEqual(recip_latt.matrix,
                                    recip_latt_xtal.matrix * 2 * np.pi,
                                    5)


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
        """
        when only lengths and angles are given for constructors, the
        internal matrix representation is ambiguous since the lattice rotation
        is not specified.
        This test makes sure that a consistent definition is specified for the
        lattice rotation when using different constructors from lengths angles
        """
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

        expected_ans = Lattice(np.array(
            [0.0, 1.0, 0.0, 1.0, 0.0, 1.0, -2.0, 0.0, 1.0]).reshape((3, 3)))
        self.assertAlmostEqual(
            np.linalg.det(np.linalg.solve(expected_ans.matrix,
                                          reduced_latt.matrix)),
            1)
        self.assertArrayAlmostEqual(
            sorted(reduced_latt.abc), sorted(expected_ans.abc))
        self.assertAlmostEqual(reduced_latt.volume, lattice.volume)
        latt = [7.164750, 2.481942, 0.000000,
                - 4.298850, 2.481942, 0.000000,
                0.000000, 0.000000, 14.253000]
        expected_ans = Lattice(np.array(
            [-4.298850, 2.481942, 0.000000, 2.865900, 4.963884, 0.000000,
             0.000000, 0.000000, 14.253000]))
        reduced_latt = Lattice(latt).get_lll_reduced_lattice()
        self.assertAlmostEqual(
            np.linalg.det(np.linalg.solve(expected_ans.matrix,
                                          reduced_latt.matrix)),
            1)
        self.assertArrayAlmostEqual(
            sorted(reduced_latt.abc), sorted(expected_ans.abc))

        expected_ans = Lattice([0.0, 10.0, 10.0,
                                10.0, 10.0, 0.0,
                                30.0, -30.0, 40.0])

        lattice = np.array([100., 0., 10., 10., 10., 20., 10., 10., 10.])
        lattice = lattice.reshape(3, 3)
        lattice = Lattice(lattice.T)
        reduced_latt = lattice.get_lll_reduced_lattice()
        self.assertAlmostEqual(
            np.linalg.det(np.linalg.solve(expected_ans.matrix,
                                          reduced_latt.matrix)),
            1)
        self.assertArrayAlmostEqual(
            sorted(reduced_latt.abc), sorted(expected_ans.abc))

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
        self.assertArrayAlmostEqual(latt.get_niggli_reduced_lattice().matrix,
                                    ans)

        latt = Lattice.from_parameters(7.365450, 6.199506, 5.353878,
                                       75.542191, 81.181757, 156.396627)
        ans = [[2.578932, 0.826965, 0.000000],
               [-0.831059, 2.067413, 1.547813],
               [-0.458407, -2.480895, 1.129126]]
        self.assertArrayAlmostEqual(latt.get_niggli_reduced_lattice().matrix,
                                    np.array(ans), 5)

    def test_find_mapping(self):
        m = np.array([[0.1, 0.2, 0.3], [-0.1, 0.2, 0.7], [0.6, 0.9, 0.2]])
        latt = Lattice(m)

        op = SymmOp.from_origin_axis_angle([0, 0, 0], [2, 3, 3], 35)
        rot = op.rotation_matrix
        scale = np.array([[1, 1, 0], [0, 1, 0], [0, 0, 1]])

        latt2 = Lattice(np.dot(rot, np.dot(scale, m).T).T)
        (aligned_out, rot_out, scale_out) = latt2.find_mapping(latt)
        self.assertAlmostEqual(abs(np.linalg.det(rot)), 1)

        rotated = SymmOp.from_rotation_and_translation(rot_out).operate_multi(latt.matrix)

        self.assertArrayAlmostEqual(rotated, aligned_out.matrix)
        self.assertArrayAlmostEqual(np.dot(scale_out, latt2.matrix), aligned_out.matrix)
        self.assertArrayAlmostEqual(aligned_out.lengths_and_angles, latt.lengths_and_angles)
        self.assertFalse(np.allclose(aligned_out.lengths_and_angles,
                                     latt2.lengths_and_angles))

    def test_find_all_mappings(self):
        m = np.array([[0.1, 0.2, 0.3], [-0.1, 0.2, 0.7], [0.6, 0.9, 0.2]])
        latt = Lattice(m)

        op = SymmOp.from_origin_axis_angle([0, 0, 0], [2, -1, 3], 40)
        rot = op.rotation_matrix
        scale = np.array([[0, 2, 0], [1, 1, 0], [0,0,1]])

        latt2 = Lattice(np.dot(rot, np.dot(scale, m).T).T)

        for (aligned_out, rot_out, scale_out) in latt.find_all_mappings(latt2):
            self.assertArrayAlmostEqual(np.inner(latt2.matrix, rot_out),
                                        aligned_out.matrix, 5)
            self.assertArrayAlmostEqual(np.dot(scale_out, latt.matrix),
                                        aligned_out.matrix)
            self.assertArrayAlmostEqual(aligned_out.lengths_and_angles, latt2.lengths_and_angles)
            self.assertFalse(np.allclose(aligned_out.lengths_and_angles,
                                         latt.lengths_and_angles))

        latt = Lattice.orthorhombic(9, 9, 5)
        self.assertEqual(len(list(latt.find_all_mappings(latt))), 16)

        #catch the singular matrix error
        latt = Lattice.from_lengths_and_angles([1,1,1], [10,10,10])
        for l, _, _ in latt.find_all_mappings(latt, ltol=0.05, atol=11):
            self.assertTrue(isinstance(l, Lattice))

    def test_to_from_dict(self):
        d = self.tetragonal.as_dict()
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

    def test_scale(self):
        new_volume = 10
        for (family_name, lattice) in self.families.items():
            new_lattice = lattice.scale(new_volume)
            self.assertAlmostEqual(new_lattice.volume, new_volume)
            self.assertEqual(new_lattice.angles, lattice.angles)

    def test_get_wigner_seitz_cell(self):
        ws_cell = Lattice([[10, 0, 0], [0, 5, 0], [0, 0, 1]])\
            .get_wigner_seitz_cell()
        self.assertEqual(6, len(ws_cell))
        for l in ws_cell[3]:
            self.assertEqual([abs(i) for i in l], [5.0, 2.5, 0.5])

    def test_dot_and_norm(self):
        frac_basis = [[1,0,0], [0,1,0], [0,0,1]]

        for family_name, lattice in self.families.items():
            #print(family_name)
            self.assert_equal(lattice.norm(lattice.matrix, frac_coords=False), lattice.abc)
            self.assert_equal(lattice.norm(frac_basis), lattice.abc)
            for (i, vec) in enumerate(frac_basis):
                length = lattice.norm(vec)
                self.assert_equal(length[0], lattice.abc[i])
                # We always get a ndarray.
                self.assertTrue(hasattr(length, "shape"))

        # Passing complex arrays should raise TypeError
        with self.assertRaises(TypeError):
            lattice.norm(np.zeros(3, dtype=np.complex))

        # Cannot reshape the second argument.
        with self.assertRaises(ValueError):
            lattice.dot(np.zeros(6), np.zeros(8))

        # Passing vectors of different length is invalid.
        with self.assertRaises(ValueError):
            lattice.dot(np.zeros(3), np.zeros(6))

    def test_get_points_in_sphere(self):
        latt = Lattice.cubic(1)
        pts = []
        for a, b, c in itertools.product(range(10), range(10), range(10)):
            pts.append([a / 10, b / 10, c / 10])

        self.assertEqual(len(latt.get_points_in_sphere(
            pts, [0, 0, 0], 0.1)), 7)
        self.assertEqual(len(latt.get_points_in_sphere(
            pts, [0.5, 0.5, 0.5], 0.5)), 515)

    def test_get_all_distances(self):
        fcoords = np.array([[0.3, 0.3, 0.5],
                            [0.1, 0.1, 0.3],
                            [0.9, 0.9, 0.8],
                            [0.1, 0.0, 0.5],
                            [0.9, 0.7, 0.0]])
        lattice = Lattice.from_lengths_and_angles([8, 8, 4],
                                                  [90, 76, 58])
        expected = np.array([[0.000, 3.015, 4.072, 3.519, 3.245],
                             [3.015, 0.000, 3.207, 1.131, 4.453],
                             [4.072, 3.207, 0.000, 2.251, 1.788],
                             [3.519, 1.131, 2.251, 0.000, 3.852],
                             [3.245, 4.453, 1.788, 3.852, 0.000]])
        output = lattice.get_all_distances(fcoords, fcoords)
        self.assertArrayAlmostEqual(output, expected, 3)
        #test just one input point
        output2 = lattice.get_all_distances(fcoords[0], fcoords)
        self.assertArrayAlmostEqual(output2, [expected[0]], 2)
        #test distance when initial points are not in unit cell
        f1 = [0, 0, 17]
        f2 = [0, 0, 10]
        self.assertEqual(lattice.get_all_distances(f1, f2)[0, 0], 0)

    def test_monoclinic(self):
        lengths, angles = self.monoclinic.lengths_and_angles
        self.assertNotAlmostEqual(angles[1], 90)
        self.assertAlmostEqual(angles[0], 90)
        self.assertAlmostEqual(angles[2], 90)

    def test_is_hexagonal(self):
        self.assertFalse(self.cubic.is_hexagonal())
        self.assertFalse(self.tetragonal.is_hexagonal())
        self.assertFalse(self.orthorhombic.is_hexagonal())
        self.assertFalse(self.monoclinic.is_hexagonal())
        self.assertFalse(self.rhombohedral.is_hexagonal())
        self.assertTrue(self.hexagonal.is_hexagonal())

    def test_get_distance_and_image(self):
        dist, image = self.cubic.get_distance_and_image([0, 0, 0.1], [0, 0.,
                                                                     0.9])
        self.assertAlmostEqual(dist, 2)
        self.assertArrayAlmostEqual(image, [0, 0, -1])

    def test_get_all_distance_and_image(self):
        r = self.cubic.get_all_distance_and_image([0, 0, 0.1],
                                                  [0, 0., 0.9])
        self.assertEqual(len(r), 8)
        dist, image = min(r, key=lambda x: x[0])
        self.assertAlmostEqual(dist, 2)
        self.assertArrayAlmostEqual(image, [0, 0, -1])
        dist, image = max(r, key=lambda x: x[0])
        self.assertAlmostEqual(dist, 16.24807680927192)
        self.assertArrayAlmostEqual(image, [1, 1, 0])


if __name__ == '__main__':
    import unittest
    unittest.main()
