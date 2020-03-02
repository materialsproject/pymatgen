# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import unittest

from pymatgen.command_line.bader_caller import *
from monty.os.path import which
import numpy as np


@unittest.skipIf(not which('bader'), "bader executable not present.")
class BaderAnalysisTest(unittest.TestCase):
    _multiprocess_shared_ = True

    def setUp(self):
        warnings.catch_warnings()
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

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
        ans = [-1.3863218, -1.3812175, -1.3812175, -1.2615902, -1.3812175, -1.3862971,
               1.021523, 1.024357, 1.021523, 1.021523, 1.021523, 1.021523, 1.021523, 1.024357]
        for i in range(14):
            self.assertAlmostEqual(ans[i], analysis.get_charge_transfer(i), 3)
        s = analysis.get_oxidation_state_decorated_structure()
        self.assertAlmostEqual(s[0].specie.oxi_state, 1.3863218, 3)

        # make sure bader still runs without reference file
        analysis = BaderAnalysis(os.path.join(test_dir, "CHGCAR.Fe3O4"))
        self.assertEqual(len(analysis.data), 14)

    def test_from_path(self):
        test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                "test_files", "bader")
        analysis = BaderAnalysis.from_path(test_dir)
        chgcar = os.path.join(test_dir, "CHGCAR.gz")
        chgref = os.path.join(test_dir, "_CHGCAR_sum.gz")
        analysis0 = BaderAnalysis(chgcar_filename=chgcar,
                                  chgref_filename=chgref)
        charge = np.array(analysis.summary["charge"])
        charge0 = np.array(analysis0.summary["charge"])
        self.assertTrue(np.allclose(charge, charge0))
        if os.path.exists("CHGREF"):
            os.remove("CHGREF")

    def test_automatic_runner(self):
        test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                'test_files/bader')

        summary = bader_analysis_from_path(test_dir)

        '''
        Reference summary dict (with bader 1.0)
        summary_ref = {
            'magmom': [4.298761, 4.221997, 4.221997, 3.816685, 4.221997, 4.298763, 0.36292,
                       0.370516, 0.36292, 0.36292, 0.36292, 0.36292, 0.36292, 0.370516],
            'min_dist': [0.835789, 0.92947, 0.92947, 0.973007, 0.92947, 0.835789, 0.94067,
                         0.817381, 0.94067, 0.94067, 0.94067, 0.94067, 0.94067, 0.817381],
            'vacuum_charge': 0.0,
            'vacuum_volume': 0.0,
            'atomic_volume': [9.922887, 8.175158, 8.175158, 9.265802, 8.175158, 9.923233, 12.382546,
                              12.566972, 12.382546, 12.382546, 12.382546, 12.382546, 12.382546, 12.566972],
            'charge': [12.248132, 12.26177, 12.26177, 12.600596, 12.26177, 12.248143, 7.267303,
                       7.256998, 7.267303, 7.267303, 7.267303, 7.267303, 7.267303, 7.256998],
            'bader_version': 1.0,
            'reference_used': True
        }
        '''

        self.assertEqual(set(summary.keys()), {'magmom', 'min_dist', 'vacuum_charge', 'vacuum_volume',
                                               'atomic_volume', 'charge', 'bader_version', 'reference_used'})
        self.assertTrue(summary['reference_used'])
        self.assertAlmostEqual(sum(summary['magmom']), 28, places=1)

    def test_atom_parsing(self):
        test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                'test_files')

        # test with reference file
        analysis = BaderAnalysis(os.path.join(test_dir, "CHGCAR.Fe3O4"),
                                 os.path.join(test_dir, "POTCAR.Fe3O4"),
                                 os.path.join(test_dir, "CHGCAR.Fe3O4_ref"),
                                 parse_atomic_densities=True)

        self.assertEqual(len(analysis.atomic_densities), len(analysis.chgcar.structure))

        self.assertAlmostEqual(np.sum(analysis.chgcar.data['total']),
                               np.sum([np.sum(d['data']) for d in analysis.atomic_densities]))


if __name__ == '__main__':
    unittest.main()
