#!/usr/bin/env python

"""
TODO: Modify unittest doc.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "5/22/14"

import unittest
import os
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.io.smartio import read_structure


test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        'test_files')


class XRDCalculatorTest(unittest.TestCase):

    def test_get_xrd_data(self):
        a = 4.209
        latt = Lattice.cubic(a)
        structure = Structure(latt, ["Cs", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        c = XRDCalculator()
        data = c.get_xrd_data(structure)
        #Check the first two peaks
        self.assertAlmostEqual(data[0][0], 21.107738329639844)
        self.assertAlmostEqual(data[0][1], 36.483184003748946)
        self.assertAlmostEqual(data[1][0], 30.024695921112777)
        self.assertAlmostEqual(data[1][1], 100)
        s = read_structure(os.path.join(test_dir, "LiFePO4.cif"))
        data = c.get_xrd_data(s)
        self.assertAlmostEqual(data[1][0], 17.03504233621785)
        self.assertAlmostEqual(data[1][1], 50.400928948337075)


if __name__ == '__main__':
    unittest.main()
