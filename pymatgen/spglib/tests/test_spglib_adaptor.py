#!/usr/bin/env python

'''
Created on Mar 9, 2012
'''

from __future__ import division

__author__="Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Mar 9, 2012"

import unittest
import os

import pymatgen.io
from pymatgen.io.vaspio import Poscar
from pymatgen.spglib.spglib_adaptor import SymmetryFinder

module_dir = os.path.dirname(os.path.abspath(pymatgen.io.__file__))

class SymmetryFinderTest(unittest.TestCase):

    def setUp(self):
        p = Poscar.from_file(os.path.join(module_dir, 'tests','vasp_testfiles', 'POSCAR'))
        self.sg = SymmetryFinder(p.struct)
        
    def test_get_space_group(self):
        self.assertEqual(self.sg.get_spacegroup(), "Pnma       (62)")

    def test_get_space_symbol(self):
        self.assertEqual(self.sg.get_spacegroup_symbol(), "Pnma")
    
    def test_get_space_number(self):
        self.assertEqual(self.sg.get_spacegroup_number(), 62)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()