#!/usr/bin/env python

'''
Created on Apr 28, 2012
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Apr 28, 2012"

import unittest
import os

from nose.exc import SkipTest

from pymatgen.core.structure import Molecule
import pymatgen.io.babelio as babelio
from pymatgen.io.babelio import BabelMolAdaptor

import pymatgen

test_dir = os.path.join(os.path.dirname(os.path.abspath(pymatgen.__file__)),
                        '..', 'test_files')


class BabelMolAdaptorTest(unittest.TestCase):

    def test_init(self):
        if not babelio.babel_loaded:
            raise SkipTest("OpenBabel not present. Skipping...")
        coords = [[0.000000, 0.000000, 0.000000],
                  [0.000000, 0.000000, 1.089000],
                  [1.026719, 0.000000, -0.363000],
                  [-0.513360, -0.889165, -0.363000],
                  [-0.513360, 0.889165, -0.363000]]
        mol = Molecule(["C", "H", "H", "H", "H"], coords)
        adaptor = BabelMolAdaptor(mol)
        obmol = adaptor.openbabel_mol
        self.assertEqual(obmol.NumAtoms(), 5)

        adaptor = BabelMolAdaptor(adaptor.openbabel_mol)
        self.assertEqual(adaptor.pymatgen_mol.formula, "H4 C1")

    def test_from_file(self):
        if not babelio.babel_loaded:
            raise SkipTest("OpenBabel not present. Skipping...")
        adaptor = BabelMolAdaptor.from_file(os.path.join(test_dir,
                                                         "Ethane_e.pdb"),
                                            "pdb")
        mol = adaptor.pymatgen_mol
        self.assertEqual(mol.formula, "H6 C2")

if __name__ == "__main__":
    if babelio.babel_loaded:
        unittest.main()
    else:
        print "openbabel not installed. Skipping test."
