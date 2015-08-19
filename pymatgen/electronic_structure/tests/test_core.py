# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest

from pymatgen.electronic_structure.core import Orbital, Spin


class SpinTest(unittest.TestCase):

    def test_init(self):
        self.assertEqual(int(Spin.up), 1)
        self.assertEqual(int(Spin.down), -1)

    def test_from_int(self):
        self.assertEqual(Spin.from_int(1), Spin.up)
        self.assertEqual(Spin.from_int(-1), Spin.down)
        self.assertRaises(ValueError, Spin.from_int, 0)

    def test_cached(self):
        self.assertEqual(id(Spin.from_int(1)), id(Spin.up))


class OrbitalTest(unittest.TestCase):

    def test_init(self):
        for i, orb in enumerate(Orbital.all_orbitals):
            self.assertEqual(Orbital.from_vasp_index(i), orb)
        self.assertRaises(IndexError, Orbital.from_vasp_index, 100)

    def test_cached(self):
        self.assertEqual(id(Orbital.from_vasp_index(0)), id(Orbital.s))

if __name__ == '__main__':
    unittest.main()
