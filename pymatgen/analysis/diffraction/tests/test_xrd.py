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
        data = c.get_xrd_data(structure, two_theta_range=(0, 90))
        #Check the first two peaks
        self.assertAlmostEqual(data[0][0], 21.107738329639844)
        self.assertAlmostEqual(data[0][1], 36.483184003748946)
        self.assertEqual(data[0][2], {(1, 0, 0): 6})
        self.assertAlmostEqual(data[0][3], 4.2089999999999996)
        self.assertAlmostEqual(data[1][0], 30.024695921112777)
        self.assertAlmostEqual(data[1][1], 100)
        self.assertEqual(data[1][2], {(1, 1, 0): 12})
        self.assertAlmostEqual(data[1][3], 2.976212442014178)

        s = read_structure(os.path.join(test_dir, "LiFePO4.cif"))
        data = c.get_xrd_data(s, two_theta_range=(0, 90))
        self.assertAlmostEqual(data[1][0], 17.03504233621785)
        self.assertAlmostEqual(data[1][1], 50.400928948337075)

        s = read_structure(os.path.join(test_dir, "Li10GeP2S12.cif"))
        data = c.get_xrd_data(s, two_theta_range=(0, 90))
        self.assertAlmostEqual(data[1][0], 14.058274883353876)
        self.assertAlmostEqual(data[1][1], 4.4111123641667671)

        # Test a hexagonal structure.
        s = read_structure(os.path.join(test_dir, "Graphite.cif"),
                           primitive=False)
        data = c.get_xrd_data(s, two_theta_range=(0, 90))
        self.assertAlmostEqual(data[0][0], 7.929279053132362)
        self.assertAlmostEqual(data[0][1], 100)
        self.assertAlmostEqual(len(list(data[0][2].keys())[0]), 4)

        #Add test case with different lengths of coefficients.
        #Also test d_hkl.
        coords = [[0.25, 0.25, 0.173], [0.75, 0.75, 0.827], [0.75, 0.25, 0],
                  [0.25, 0.75, 0], [0.25, 0.25, 0.676], [0.75, 0.75, 0.324]]
        sp = ["Si", "Si", "Ru", "Ru", "Pr", "Pr"]
        s = Structure(Lattice.tetragonal(4.192, 6.88), sp, coords)
        data = c.get_xrd_data(s)
        self.assertAlmostEqual(data[0][0], 12.86727341476735)
        self.assertAlmostEqual(data[0][1], 31.448239816769796)
        self.assertAlmostEqual(data[0][3], 6.88)
        self.assertEqual(len(data), 42)
        data = c.get_xrd_data(s, two_theta_range=[0, 60])
        self.assertEqual(len(data), 18)

        #Test with and without Debye-Waller factor
        tungsten = Structure(Lattice.cubic(3.1653), ["W"] * 2,
                             [[0, 0, 0], [0.5, 0.5, 0.5]])
        data = c.get_xrd_data(tungsten, scaled=False)
        self.assertAlmostEqual(data[0][0], 40.294828554672264)
        self.assertAlmostEqual(data[0][1], 2414237.5633093244)
        self.assertAlmostEqual(data[0][3], 2.2382050944897789)
        c = XRDCalculator(debye_waller_factors={"W": 0.1526})
        data = c.get_xrd_data(tungsten, scaled=False)
        self.assertAlmostEqual(data[0][0], 40.294828554672264)
        self.assertAlmostEqual(data[0][1], 2377745.2296686019)
        self.assertAlmostEqual(data[0][3], 2.2382050944897789)


if __name__ == '__main__':
    unittest.main()
