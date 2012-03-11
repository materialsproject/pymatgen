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

import numpy as np

from pymatgen.core.structure import PeriodicSite
from pymatgen.io.vaspio import Poscar
from pymatgen.symmetry.spglib_adaptor import SymmetryFinder, get_pointgroup
from pymatgen.io.cifio import CifParser

import pymatgen

test_dir = os.path.join(os.path.dirname(os.path.abspath(pymatgen.__file__)), '..', 'test_files')

class SymmetryFinderTest(unittest.TestCase):

    def setUp(self):
        p = Poscar.from_file(os.path.join(test_dir, 'POSCAR'))
        self.structure = p.struct
        self.sg = SymmetryFinder(self.structure, 0.1)
        
    def test_get_space_symbol(self):
        self.assertEqual(self.sg.get_spacegroup_symbol(), "Pnma")
    
    def test_get_space_number(self):
        self.assertEqual(self.sg.get_spacegroup_number(), 62)
        
    def test_get_symmetry_dataset(self):
        ds = self.sg.get_symmetry_dataset()
        self.assertEqual(ds['international'], 'Pnma')
    
    def test_get_symmetry_operations(self):
        fracsymmops = self.sg.get_symmetry_operations()
        symmops = self.sg.get_symmetry_operations(True)
        self.assertEqual(len(symmops), 8)
        latt = self.structure.lattice
        for fop, op in zip(fracsymmops, symmops):
            for site in self.structure:
                newfrac = fop.operate(site.frac_coords)
                newcart = op.operate(site.coords)
                self.assertTrue(np.allclose(latt.get_fractional_coords(newcart), newfrac))
                found = False
                newsite = PeriodicSite(site.species_and_occu, newcart, latt, coords_are_cartesian = True)
                for testsite in self.structure:
                    if newsite.is_periodic_image(testsite, 1e-3):
                        found = True
                        break
                self.assertTrue(found)
                    
    def test_get_refined_structure(self):
        for a in self.sg.get_refined_structure().lattice.angles:
            self.assertEqual(a, 90)
            
    def test_get_symmetrized_structure(self):
        symm_struct = self.sg.get_symmetrized_structure()
        for a in symm_struct.lattice.angles:
            self.assertEqual(a, 90)
        self.assertEqual(len(symm_struct.equivalent_sites), 5)
        
    def test_get_primitive(self):
        #Pnma spacegroup has no primitive cell
        self.assertIsNone(self.sg.find_primitive())
        parser = CifParser(os.path.join(test_dir, 'Li2O.cif'))
        s = SymmetryFinder(parser.get_structures(False)[0])
        self.assertEqual(s.find_primitive().formula, "Li2 O1")

class HelperFunctionsTest(unittest.TestCase):

    def setUp(self):
        p = Poscar.from_file(os.path.join(test_dir, 'POSCAR'))
        self.sg = SymmetryFinder(p.struct, 0.1)
        
    def test_get_pointgroup(self):
        (rots, trans) = self.sg.get_symmetry()
        pg = get_pointgroup(rots)
        self.assertEqual(pg[0].strip(), "mmm")

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()