# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import numpy as np
import unittest
import os

from pymatgen.analysis.find_dimension import find_dimension
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Xdatcar
from pymatgen import Element, Structure, Lattice
from pymatgen.util.testing import PymatgenTest

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')

class FindDimensionTest(PymatgenTest):

    def test_get_dimensionality(self):
        s = self.get_structure('LiFePO4')
        self.assertEqual(find_dimension(s), 'intercalated ion')

        s = self.get_structure('Graphite')
        self.assertEqual(find_dimension(s), '2D')

    def test_get_dimensionality_with_bonds(self):
        s = self.get_structure('CsCl')
        self.assertEqual(find_dimension(s), 'intercalated ion')
        self.assertEqual(find_dimension(s, ldict={"Cs":3.7, "Cl": 3}), '3D')


if __name__ == '__main__':
    unittest.main()
