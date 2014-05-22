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

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.analysis.diffraction.xrd import XRDCalculator


class XRDCalculatorTest(unittest.TestCase):

    def test_get_xrd_data(self):
        a = 4.209 #Angstrom
        latt = Lattice.cubic(a)
        structure = Structure(latt, ["Cs", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        c = XRDCalculator()
        data = c.get_xrd_data(structure)
        #Check the first two peaks
        self.assertAlmostEqual(data[0][0], 21.109953842244817)
        self.assertAlmostEqual(data[0][1][0], 36.483541952310695)
        self.assertAlmostEqual(data[1][0], 30.027884973250128)
        self.assertAlmostEqual(data[1][1][0], 100)


if __name__ == '__main__':
    unittest.main()
