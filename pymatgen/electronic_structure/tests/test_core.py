#!/usr/bin/python

import unittest

from pymatgen.electronic_structure.core import Orbital, Spin


class SpinTest(unittest.TestCase):

    def test_init(self):
        self.assertEquals(int(Spin.up), 1)
        self.assertEquals(int(Spin.down), -1)

    def test_from_int(self):
        self.assertEquals(Spin.from_int(1), Spin.up)
        self.assertEquals(Spin.from_int(-1), Spin.down)
        self.assertRaises(ValueError, Spin.from_int, 0)

    def test_cached(self):
        self.assertEquals(id(Spin.from_int(1)), id(Spin.up))


class OrbitalTest(unittest.TestCase):

    def test_init(self):
        for i, orb in enumerate(Orbital.all_orbitals):
            self.assertEqual(Orbital.from_vasp_index(i), orb)
        self.assertRaises(IndexError, Orbital.from_vasp_index, 100)

    def test_cached(self):
        self.assertEquals(id(Orbital.from_vasp_index(0)), id(Orbital.s))

if __name__ == '__main__':
    unittest.main()
