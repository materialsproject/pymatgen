#!/usr/bin/env python

'''
Created on Mar 12, 2012
'''

from __future__ import division

__author__="Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Mar 12, 2012"

import unittest
import os

from pymatgen.core.structure import PeriodicSite
from pymatgen.symmetry.spacegroup import Spacegroup
from pymatgen.io.vaspio import Poscar
from pymatgen.symmetry.spglib_adaptor import SymmetryFinder

import pymatgen

test_dir = os.path.join(os.path.dirname(os.path.abspath(pymatgen.__file__)), '..', 'test_files')


class SpacegroupTest(unittest.TestCase):

    def setUp(self):
        p = Poscar.from_file(os.path.join(test_dir, 'POSCAR'))
        self.structure = p.struct
        self.sg1 = SymmetryFinder(self.structure, 0.001).get_spacegroup()
        self.sg2 = Spacegroup.from_spacegroup_number(62)

    def test_are_symmetrically_equivalent(self):
        sites1 = [self.structure[i] for i in [0,1]]
        sites2 = [self.structure[i] for i in [2,3]]
        self.assertTrue(self.sg1.are_symmetrically_equivalent(sites1, sites2, 1e-3))
        self.assertTrue(self.sg2.are_symmetrically_equivalent(sites1, sites2, 1e-3))
        
        sites1 = [self.structure[i] for i in [0,1]]
        sites2 = [self.structure[i] for i in [0,2]]
        self.assertFalse(self.sg1.are_symmetrically_equivalent(sites1, sites2, 1e-3))
        self.assertFalse(self.sg2.are_symmetrically_equivalent(sites1, sites2, 1e-3))
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()