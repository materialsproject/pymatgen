#!/usr/bin/env python

'''
Created on Apr 25, 2012
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Apr 25, 2012"

import unittest

from pymatgen.util.coord_utils import get_linear_interpolated_value, in_coord_list, get_convex_hull

class CoordUtilsTest(unittest.TestCase):

    def test_get_linear_interpolated_value(self):
        xvals = [0, 1, 2, 3, 4, 5]
        yvals = [3, 6, 7, 8, 10, 12]
        self.assertEqual(get_linear_interpolated_value(xvals, yvals, 3.6), 9.2)
        self.assertRaises(ValueError, get_linear_interpolated_value, xvals, yvals, 6)

    def test_in_coord_list(self):
        coords = [[0, 0, 0], [1, 1, 0.3]]
        test_coord = [0.1, 0.1, 0.1]
        self.assertFalse(in_coord_list(coords, test_coord))
        self.assertTrue(in_coord_list(coords, test_coord, atol=0.15))

    def test_get_convex_hull(self):
        coords = [[0, 0, 0], [0, 0, 1], [0, 1, 0], [1, 0, 0], [0.2, 0.1, 0.2]]
        hull1 = get_convex_hull(coords)
        hull2 = get_convex_hull(coords, True)
        hull1 = [set(l) for l in hull1]
        hull2 = [set(l) for l in hull2]
        self.assertEqual(hull1, hull2)
        for l in hull1:
            self.assertNotIn(4, l)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
