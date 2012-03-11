#!/usr/bin/env python

'''
Created on Mar 5, 2012
'''

from __future__ import division

__author__="Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Mar 5, 2012"

import unittest
import os

from pymatgen.transformations.standard_transformations import SubstitutionTransformation
from pymatgen.alchemy.transmuters import CifTransmuter, PoscarTransmuter

import pymatgen

test_dir = os.path.join(os.path.dirname(os.path.abspath(pymatgen.__file__)), '..', 'test_files')


class CifTransmuterTest(unittest.TestCase):

    def setUp(self):
        trans = []
        trans.append(SubstitutionTransformation({"Fe":"Mn", "Fe2+":"Mn2+"}))
        self.qep = CifTransmuter(trans)
        
    def test_transmute(self):
        trans_structures = self.qep.transmute([os.path.join(test_dir, "MultiStructure.cif")])
        self.assertEqual(len(trans_structures), 2)
        expected_ans = set(["Mn", "O", "Li", "P"])
        for s in trans_structures:
            els = set([el.symbol for el in s.final_structure.composition.elements])
            self.assertEqual(expected_ans, els)

class PoscarTransmuterTest(unittest.TestCase):

    def setUp(self):
        trans = []
        trans.append(SubstitutionTransformation({"Fe":"Mn"}))
        self.qep = PoscarTransmuter(trans)
        
    def test_transmute(self):
        trans_structures = self.qep.transmute([os.path.join(test_dir, "POSCAR"), os.path.join(test_dir, "POSCAR")])
        self.assertEqual(len(trans_structures), 2)
        expected_ans = set(["Mn", "O", "P"])
        for s in trans_structures:
            els = set([el.symbol for el in s.final_structure.composition.elements])
            self.assertEqual(expected_ans, els)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()