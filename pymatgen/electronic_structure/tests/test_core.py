#!/usr/bin/python

import unittest
import os

from pymatgen.electronic_structure.core import Orbital, Spin

import pymatgen

test_dir = os.path.join(os.path.dirname(os.path.abspath(pymatgen.__file__)), '..', 'test_files')

class SpinTest(unittest.TestCase):

    def test_init(self):
        self.assertEquals(int(Spin.up), 1)
        self.assertEquals(int(Spin.down), -1)

class OrbitalTest(unittest.TestCase):

    def test_init(self):
        self.assertEqual(Orbital.from_vasp_index(1), Orbital.py)


if __name__ == '__main__':
    unittest.main()

