#!/usr/bin/env python

"""
Created on Apr 25, 2012
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Apr 25, 2012"

import itertools

import numpy as np
from pymatgen.core.lattice import Lattice
from pymatgen.util.coord_utils import get_linear_interpolated_value,\
    in_coord_list, pbc_diff, in_coord_list_pbc, get_points_in_sphere_pbc,\
    find_in_coord_list, find_in_coord_list_pbc, pbc_all_distances,\
    barycentric_coords, pbc_shortest_vectors
from pymatgen.util.testing import PymatgenTest


class CoordUtilsTest(PymatgenTest):

    def test_get_linear_interpolated_value(self):
        xvals = [0, 1, 2, 3, 4, 5]
        yvals = [3, 6, 7, 8, 10, 12]
        self.assertEqual(get_linear_interpolated_value(xvals, yvals, 3.6), 9.2)
        self.assertRaises(ValueError, get_linear_interpolated_value, xvals,
                          yvals, 6)

    def test_in_coord_list(self):
        coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
        test_coord = [0.1, 0.1, 0.1]
        self.assertFalse(in_coord_list(coords, test_coord))
        self.assertTrue(in_coord_list(coords, test_coord, atol=0.15))
        self.assertFalse(in_coord_list([0.99, 0.99, 0.99], test_coord,
                                       atol=0.15))

    def test_find_in_coord_list(self):
        coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
        test_coord = [0.1, 0.1, 0.1]
        self.assertFalse(find_in_coord_list(coords, test_coord))
        self.assertEqual(find_in_coord_list(coords, test_coord, atol=0.15)[0],
                         0)
        self.assertFalse(find_in_coord_list([0.99, 0.99, 0.99], test_coord,
                                            atol=0.15))
        coords = [[0, 0, 0], [0.5, 0.5, 0.5], [0.1, 0.1, 0.1]]
        self.assertArrayEqual(find_in_coord_list(coords, test_coord,
                                                 atol=0.15), [0, 2])

    def test_pbc_diff(self):
        self.assertArrayAlmostEqual(pbc_diff([0.1, 0.1, 0.1], [0.3, 0.5, 0.9]),
                                    [-0.2, -0.4, 0.2])
        self.assertArrayAlmostEqual(pbc_diff([0.9, 0.1, 1.01],
                                             [0.3, 0.5, 0.9]),
                                    [-0.4, -0.4, 0.11])
        self.assertArrayAlmostEqual(pbc_diff([0.1, 0.6, 1.01],
                                             [0.6, 0.1, 0.9]),
                                    [-0.5, 0.5, 0.11])
        self.assertArrayAlmostEqual(pbc_diff([100.1, 0.2, 0.3],
                                             [0123123.4, 0.5, 502312.6]),
                                    [-0.3, -0.3, -0.3])

    def test_pbc_all_distances(self):
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
        output = pbc_all_distances(lattice, fcoords, fcoords)
        self.assertArrayAlmostEqual(output, expected, 3)
        #test just one input point
        output2 = pbc_all_distances(lattice, fcoords[0], fcoords)
        self.assertArrayAlmostEqual(output2, [expected[0]], 2)
        #test distance when initial points are not in unit cell
        f1 = [0, 0, 17]
        f2 = [0, 0, 10]
        self.assertEqual(pbc_all_distances(lattice, f1, f2)[0, 0], 0)

    def test_in_coord_list_pbc(self):
        coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
        test_coord = [0.1, 0.1, 0.1]
        self.assertFalse(in_coord_list_pbc(coords, test_coord))
        self.assertTrue(in_coord_list_pbc(coords, test_coord, atol=0.15))
        test_coord = [0.99, 0.99, 0.99]
        self.assertFalse(in_coord_list_pbc(coords, test_coord, atol=0.01))

    def test_find_in_coord_list_pbc(self):
        coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
        test_coord = [0.1, 0.1, 0.1]
        self.assertFalse(find_in_coord_list_pbc(coords, test_coord))
        self.assertEqual(find_in_coord_list_pbc(coords, test_coord,
                                                atol=0.15)[0], 0)
        test_coord = [0.99, 0.99, 0.99]
        self.assertEqual(
            find_in_coord_list_pbc(coords, test_coord, atol=0.02)[0], 0)
        test_coord = [-0.499, -0.499, -0.499]
        self.assertEqual(
            find_in_coord_list_pbc(coords, test_coord, atol=0.01)[0], 1)

    def test_get_points_in_sphere_pbc(self):
        latt = Lattice.cubic(1)
        pts = []
        for a, b, c in itertools.product(xrange(10), xrange(10), xrange(10)):
            pts.append([a / 10, b / 10, c / 10])

        self.assertEqual(len(get_points_in_sphere_pbc(latt, pts, [0, 0, 0],
                                                      0.1)), 7)
        self.assertEqual(len(get_points_in_sphere_pbc(latt, pts,
                                                      [0.5, 0.5, 0.5],
                                                      0.5)), 515)

    def test_barycentric(self):
        #2d test
        simplex1 = np.array([[0.3, 0.1], [0.2, -1.2], [1.3, 2.3]])
        pts1 = np.array([[0.6, 0.1], [1.3, 2.3], [0.5, 0.5], [.7, 1]])
        output1 = barycentric_coords(pts1, simplex1)
        #do back conversion to cartesian
        o_dot_s = np.sum(output1[:, :, None] * simplex1[None, :, :], axis=1)
        self.assertTrue(np.allclose(pts1, o_dot_s))

        #do 3d tests
        simplex2 = np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0], [0, 0, 0]])
        pts2 = np.array([[0, 0, 1], [0, 0.5, 0.5], [1./3, 1./3, 1./3]])
        output2 = barycentric_coords(pts2, simplex2)
        self.assertTrue(np.allclose(output2[1], [0.5, 0.5, 0, 0]))
        #do back conversion to cartesian
        o_dot_s = np.sum(output2[:, :, None] * simplex2[None, :, :], axis=1)
        self.assertTrue(np.allclose(pts2, o_dot_s))
        #test single point
        self.assertTrue(np.allclose(output2[2],
                                    barycentric_coords(pts2[2], simplex2)))
        
    def test_pbc_shortest_vectors(self):
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
                             [3.519, 1.131, 2.251, 0.000, 3.852]])
        vectors = pbc_shortest_vectors(lattice, fcoords[:-1], fcoords)
        dists = np.sum(vectors**2, axis = -1)**0.5
        self.assertArrayAlmostEqual(dists, expected, 3)
        

if __name__ == "__main__":
    import unittest
    unittest.main()
