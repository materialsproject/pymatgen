#!/usr/bin/env python

'''
Created on Jul 30, 2012
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Jul 30, 2012"

import unittest
import os

from pymatgen import __file__
from pymatgen.io.smartio import read_structure
from pymatgen.core.structure import Structure

test_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..',
                            'test_files')


class MethodsTest(unittest.TestCase):

    def test_read_structure(self):
        for fname in ("Li2O.cif", "Li2O2.cif", "vasprun.xml",
                      "vasprun_Si_bands.xml", "Si.cssr"):
            filename = os.path.join(test_dir, fname)
            struct = read_structure(filename)
            self.assertIsInstance(struct, Structure)
            print struct


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
