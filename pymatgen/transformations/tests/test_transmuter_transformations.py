#!/usr/bin/env python
from __future__ import division

'''
Created on Sep 23, 2011
'''

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Sep 23, 2011"

import os
import unittest
import random

from pymatgen.transformations.transmuter_transformations import *
from pymatgen.alchemy.transmuters import TransformedStructureTransmuter
from pymatgen.alchemy.materials import TransformedStructure
from pymatgen.core.structure import Lattice, Structure


class MultipleSubstitutionTransformationTest(unittest.TestCase):

    def test_apply_transformation(self):
        sub_dict = {1: ["Na", "K"]}
        t = MultipleSubstitutionTransformation("Li+", 0.5, sub_dict, None)
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.75, 0.75])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.25, 0.25, 0.25])
        lattice = Lattice([[ 3.8401979337, 0.00, 0.00], [1.9200989668, 3.3257101909, 0.00], [0.00, -2.2171384943, 3.1355090603]])
        struct = Structure(lattice, ["Li+", "Li+", "O2-", "O2-"], coords)
        tst = TransformedStructureTransmuter([TransformedStructure(struct, [])])
        t.apply_transformation(tst)
        self.assertEqual(len(tst), 2)


if __name__ == "__main__":
    unittest.main()
