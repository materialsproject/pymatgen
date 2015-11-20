# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

__author__ = 'Shyue Ping Ong'
__copyright__ = 'Copyright 2013, The Materials Project'
__version__ = '0.1'
__maintainer__ = 'Shyue Ping Ong'
__email__ = 'ongsp@ucsd.edu'
__date__ = '9/25/14'

import unittest

import random

from pymatgen.util.num_utils import abs_cap, min_max_indexes


class FuncTestCase(unittest.TestCase):

    def test_abs_cap(self):
        self.assertEqual(abs_cap(1.000000001), 1.0)
        self.assertEqual(abs_cap(-1.000000001), -1.0)

        v = random.uniform(-1, 1)
        self.assertEqual(abs_cap(v), v)

        self.assertEqual(abs_cap(1.000000001, 2), 1.000000001)
        self.assertEqual(abs_cap(-2.000000001, 2), -2.0)

    def test_min_max_indexes(self):
        val = ['b', 'a', 'm', 'z', 'y']
        min_ind, max_ind = min_max_indexes(val)
        self.assertEqual(min_ind, 1)
        self.assertEqual(max_ind, 3)

if __name__ == '__main__':
    unittest.main()
