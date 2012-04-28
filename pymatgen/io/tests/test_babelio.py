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

from pymatgen.core.structure import Molecule
import pymatgen.io.babelio as babelio
from pymatgen.io.babelio import BabelMolAdaptor


class BabelMolAdaptorTest(unittest.TestCase):

    def test_init(self):
        coords = [[0.000000, 0.000000, 0.000000],
                  [0.000000, 0.000000, 1.089000],
                  [1.026719, 0.000000, -0.363000],
                  [-0.513360, -0.889165, -0.363000],
                  [-0.513360 , 0.889165 , -0.363000]]
        mol = Molecule(["C", "H", "H", "H", "H"], coords)
        adaptor = BabelMolAdaptor(mol)
        obmol = adaptor.openbabel_mol
        self.assertEqual(obmol.NumAtoms(), 5)

        adaptor = BabelMolAdaptor(adaptor.openbabel_mol)
        self.assertEqual(adaptor.pymatgen_mol.formula, "H4 C1")

if __name__ == "__main__":
    if babelio.babel_loaded:
        unittest.main()
    else:
        print "openbabel not installed. Skipping test."
