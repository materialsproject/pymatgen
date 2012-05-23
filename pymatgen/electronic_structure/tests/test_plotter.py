#!/usr/bin/env python

'''
Created on May 1, 2012
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "May 1, 2012"

import unittest
import os
import json

from pymatgen.electronic_structure.dos import CompleteDos
from pymatgen.electronic_structure.plotter import DosPlotter

import pymatgen

test_dir = os.path.join(os.path.dirname(os.path.abspath(pymatgen.__file__)), '..', 'test_files')


class DosPlotterTest(unittest.TestCase):

    def setUp(self):
        with open(os.path.join(test_dir, "complete_dos.json"), "r") as f:
            self.dos = CompleteDos.from_dict(json.load(f))
            self.plotter = DosPlotter(sigma=0.2, stack=True)

    def test_add_dos_dict(self):
        d = self.plotter.get_dos_dict()
        self.assertEqual(len(d), 0)
        self.plotter.add_dos_dict(self.dos.get_element_dos(), key_sort_func=lambda x:x.X)
        d = self.plotter.get_dos_dict()
        self.assertEqual(len(d), 4)

    def test_get_dos_dict(self):
        self.plotter.add_dos_dict(self.dos.get_element_dos(), key_sort_func=lambda x:x.X)
        d = self.plotter.get_dos_dict()
        for el in ["Li", "Fe", "P", "O"]:
            self.assertIn(el, d)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
