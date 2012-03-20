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
from pymatgen.alchemy.materials import TransformedStructure
from pymatgen.alchemy.transmuters import TransformedStructureTransmuter
from pymatgen.transformations.standard_transformations import SubstitutionTransformation, RemoveSpeciesTransformation, OrderDisorderedStructureTransformation

import pymatgen
test_dir = os.path.join(os.path.dirname(os.path.abspath(pymatgen.__file__)), '..', 'test_files')


class TransformedStructureTransmuterTest(unittest.TestCase):
        
    def test_cif_transmute(self):
        trans = []
        trans.append(SubstitutionTransformation({"Fe":"Mn", "Fe2+":"Mn2+"}))
        tsc = TransformedStructureTransmuter.from_cifs([os.path.join(test_dir, "MultiStructure.cif")], trans)
        self.assertEqual(len(tsc), 2)
        expected_ans = set(["Mn", "O", "Li", "P"])
        for s in tsc:
            els = set([el.symbol for el in s.final_structure.composition.elements])
            self.assertEqual(expected_ans, els)
          
    def test_poscar_transmute(self):
        trans = []
        trans.append(SubstitutionTransformation({"Fe":"Mn"}))
        tsc = TransformedStructureTransmuter.from_poscars([os.path.join(test_dir, "POSCAR"), os.path.join(test_dir, "POSCAR")], trans)
        self.assertEqual(len(tsc), 2)
        expected_ans = set(["Mn", "O", "P"])
        for s in tsc:
            els = set([el.symbol for el in s.final_structure.composition.elements])
            self.assertEqual(expected_ans, els)
            
    def test_transmuter(self):
        tsc = TransformedStructureTransmuter.from_poscars([os.path.join(test_dir, "POSCAR")])
        tsc.append_transformation(RemoveSpeciesTransformation('O'))
        self.assertEqual(len(tsc[0].final_structure), 8)
        
        tsc.append_transformation(SubstitutionTransformation({"Fe":{"Fe2+":.25, "Mn3+":.75}, "P":"P5+"}))
        tsc.append_transformation(OrderDisorderedStructureTransformation(num_structures = 50), extend_collection = True)
        self.assertEqual(len(tsc), 4)
        
        tsc.branch_collection([SubstitutionTransformation({"Fe2+":"Mg2+"}),SubstitutionTransformation({"Fe2+":"Zn2+"})], retention_level = 1)
        self.assertEqual(len(tsc), 8)
        
        tsc.branch_collection([SubstitutionTransformation({"Zn2+":"Fe2+"})], retention_level = 0)
        self.assertEqual(len(tsc), 4)

        self.assertEqual(len(tsc[0]), 6)
        #5 transformations plus initial structure
        
        tsc = TransformedStructureTransmuter([tsc[0]],[SubstitutionTransformation({"Mn3+":{"Fe2+":0.333}}), OrderDisorderedStructureTransformation(num_structures = 50)], extend_collection = False)
        self.assertEqual(len(tsc), 1)
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()