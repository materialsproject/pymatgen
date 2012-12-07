#!/usr/bin/env python

"""
Created on Apr 25, 2012
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Apr 25, 2012"

import unittest
import itertools

import numpy as np
from pymatgen.core.lattice import Lattice
from pymatgen.util.coord_utils import get_linear_interpolated_value,\
    in_coord_list, pbc_diff, in_coord_list_pbc, get_points_in_sphere_pbc,\
    find_in_coord_list, find_in_coord_list_pbc


class CoordUtilsTest(unittest.TestCase):

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
        self.assertTrue(all(np.equal(find_in_coord_list(coords, test_coord,
                                                        atol=0.15), [0, 2])))

    def test_pbc_diff(self):
        self.assertTrue(np.allclose(pbc_diff([0.1, 0.1, 0.1], [0.3, 0.5, 0.9]),
                                    [-0.2, -0.4, 0.2]))
        self.assertTrue(np.allclose(pbc_diff([0.9, 0.1, 1.01],
                                             [0.3, 0.5, 0.9]),
                                    [-0.4, -0.4, 0.11]))
        self.assertTrue(np.allclose(pbc_diff([0.1, 0.6, 1.01],
                                             [0.6, 0.1, 0.9]),
                                    [-0.5, 0.5, 0.11]))
        self.assertTrue(np.allclose(pbc_diff([100.1, 0.2, 0.3],
                                             [0123123.4, 0.5, 502312.6]),
                                    [-0.3, -0.3, -0.3]))

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

if __name__ == "__main__":
    unittest.main()
