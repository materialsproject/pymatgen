#!/usr/bin/env python

'''
Created on Mar 8, 2012
'''

from __future__ import division

__author__="Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Mar 8, 2012"

import unittest
import os

from pymatgen.io.vaspio import Poscar
import pymatgen.io.aseio as aio

module_dir = os.path.dirname(os.path.abspath(__file__))

class AseAtomsAdaptorTest(unittest.TestCase):

    def test_get_atoms(self):
        p = Poscar.from_file(os.path.join(module_dir, 'vasp_testfiles', 'POSCAR'))
        atoms = aio.AseAtomsAdaptor.get_atoms(p.struct)
        self.assertEqual(atoms.get_name(), "P4Fe4O16")

    def test_get_structure(self):
        p = Poscar.from_file(os.path.join(module_dir, 'vasp_testfiles', 'POSCAR'))
        atoms = aio.AseAtomsAdaptor.get_atoms(p.struct)
        self.assertEqual(aio.AseAtomsAdaptor.get_structure(atoms).formula, "Fe4 P4 O16")
    
    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    if aio.ase_loaded:
        unittest.main()
    else:
        print "Skipping tests"