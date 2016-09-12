# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest2 as unittest

from pymatgen.electronic_structure.core import Orbital, Spin


class SpinTest(unittest.TestCase):

    def test_init(self):
        self.assertEqual(int(Spin.up), 1)
        self.assertEqual(int(Spin.down), -1)

    def test_from_int(self):
        self.assertEqual(Spin(1), Spin.up)
        self.assertEqual(Spin(-1), Spin.down)
        self.assertRaises(ValueError, Spin, 0)

    def test_cached(self):
        self.assertEqual(id(Spin(1)), id(Spin.up))


class OrbitalTest(unittest.TestCase):

    def test_init(self):
        for orb in Orbital:
            self.assertEqual(Orbital(orb.value), orb)
        self.assertRaises(ValueError, Orbital, 100)

    def test_cached(self):
        self.assertEqual(id(Orbital(0)), id(Orbital.s))

if __name__ == '__main__':
    unittest.main()
