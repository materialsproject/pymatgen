# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals
import unittest
import os

from pymatgen.command_line.bader_caller import *
from monty.os.path import which

"""
TODO: Change the module doc.
"""


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Jul 22, 2012"


@unittest.skipIf(not which('bader'), "bader executable not present.")
class BaderAnalysisTest(unittest.TestCase):

    def test_init(self):
        test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                'test_files')

        # test with reference file
        analysis = BaderAnalysis(os.path.join(test_dir, "CHGCAR.Fe3O4"),
                                 os.path.join(test_dir, "POTCAR.Fe3O4"),
                                 os.path.join(test_dir, "CHGCAR.Fe3O4_ref"))
        self.assertEqual(len(analysis.data), 14)
        self.assertAlmostEqual(analysis.data[0]["charge"], 6.6136782, 3)
        self.assertAlmostEqual(analysis.nelectrons, 96)
        self.assertAlmostEqual(analysis.vacuum_charge, 0)
        ans = [-1.3863218, -1.3812175, -1.3812175, -1.2615902, -1.3812175, -1.3862971, 1.021523, 1.024357, 1.021523, 1.021523, 1.021523, 1.021523, 1.021523, 1.024357]
        for i in range(14):
            self.assertAlmostEqual(ans[i], analysis.get_charge_transfer(i), 3)
        s = analysis.get_oxidation_state_decorated_structure()
        self.assertAlmostEqual(s[0].specie.oxi_state, 1.3863218, 3)

        # make sure bader still runs without reference file
        analysis = BaderAnalysis(os.path.join(test_dir, "CHGCAR.Fe3O4"))
        self.assertEqual(len(analysis.data), 14)

    def test_automatic_runner(self):
        test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                'test_files/bader')

        bader_result = BaderResult.from_path(test_dir)

        summary_ref = {
            'magmom': [4.298761, 4.221997, 4.221997, 3.816685, 4.221997, 4.298763, 0.36292,
                       0.370516, 0.36292, 0.36292, 0.36292, 0.36292, 0.36292, 0.370516],
            'min_dist': [0.835789, 0.92947, 0.92947, 0.973007, 0.92947, 0.835789, 0.94067,
                         0.817381, 0.94067, 0.94067, 0.94067, 0.94067, 0.94067, 0.817381],
            'vacuum_charge': 0.0,
            'vacuum_volume': 0.0,
            'volume': [9.922887, 8.175158, 8.175158, 9.265802, 8.175158, 9.923233, 12.382546,
                       12.566972, 12.382546, 12.382546, 12.382546, 12.382546, 12.382546, 12.566972],
            'charge': [12.248132, 12.26177, 12.26177, 12.600596, 12.26177, 12.248143, 7.267303,
                       7.256998, 7.267303, 7.267303, 7.267303, 7.267303, 7.267303, 7.256998],
            'bader_version': 1.0,
            'reference_used': True
        }

        self.assertEqual(bader_result.__dict__, summary_ref)

        # test as/from dict
        bader_result_dict = bader_result.as_dict()
        bader_result = BaderResult.from_dict(bader_result_dict)

if __name__ == '__main__':
    unittest.main()
