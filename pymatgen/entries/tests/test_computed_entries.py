#!/usr/bin/env python

'''
Created on Mar 18, 2012
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Mar 18, 2012"

import unittest
import os

from pymatgen.io.vaspio import Vasprun
from pymatgen.entries.computed_entries import ComputedEntry, ComputedStructureEntry
import pymatgen

test_dir = os.path.join(os.path.dirname(os.path.abspath(pymatgen.__file__)), '..', 'test_files')
filepath = os.path.join(test_dir, 'vasprun.xml')
vasprun = Vasprun(filepath)

class ComputedEntryTest(unittest.TestCase):

    def setUp(self):
        self.entry = ComputedEntry(vasprun.final_structure.composition,
                                   vasprun.final_energy, parameters = vasprun.incar)

    def test_energy(self):
        self.assertAlmostEqual(self.entry.energy, -269.38319884)
        self.entry.correction = 1.0
        self.assertAlmostEqual(self.entry.energy, -268.38319884)

    def test_composition(self):
        self.assertEqual(self.entry.composition.reduced_formula, "LiFe4(PO4)4")

    def test_to_from_dict(self):
        d = self.entry.to_dict
        e = ComputedEntry.from_dict(d)
        self.assertAlmostEqual(e.energy, -269.38319884)

class ComputedStructureEntryTest(unittest.TestCase):

    def setUp(self):
        self.entry = ComputedStructureEntry(vasprun.final_structure,
                                   vasprun.final_energy, parameters = vasprun.incar)

    def test_energy(self):
        self.assertAlmostEqual(self.entry.energy, -269.38319884)
        self.entry.correction = 1.0
        self.assertAlmostEqual(self.entry.energy, -268.38319884)

    def test_composition(self):
        self.assertEqual(self.entry.composition.reduced_formula, "LiFe4(PO4)4")

    def test_to_from_dict(self):
        d = self.entry.to_dict
        e = ComputedStructureEntry.from_dict(d)
        self.assertAlmostEqual(e.energy, -269.38319884)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
