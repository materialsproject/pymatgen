# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
TODO: Change the module doc.
"""


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Jul 22, 2012"

import unittest2 as unittest
import os

from pymatgen.command_line.bader_caller import BaderAnalysis
from monty.os.path import which


@unittest.skipIf(not which('bader'), "bader executable not present.")
class BaderAnalysisTest(unittest.TestCase):

    def test_init(self):
        test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                'test_files')
        analysis = BaderAnalysis(os.path.join(test_dir, "CHGCAR.Fe3O4"),
                                 os.path.join(test_dir, "POTCAR.Fe3O4"))
        self.assertEqual(len(analysis.data), 14)
        self.assertAlmostEqual(analysis.data[0]["charge"], 6.1485)
        self.assertAlmostEqual(analysis.nelectrons, 96)
        self.assertAlmostEqual(analysis.vacuum_charge, 0)
        ans = [-1.8515, -1.8132, -1.8132, -1.5562, -1.8132, -1.8515, 1.3369,
               1.3378, 1.3373, 1.3374, 1.3374, 1.3369, 1.3373, 1.3378]
        for i in range(14):
            self.assertAlmostEqual(ans[i], analysis.get_charge_transfer(i))
        s = analysis.get_oxidation_state_decorated_structure()
        self.assertAlmostEqual(s[0].specie.oxi_state, -1.8515)


if __name__ == '__main__':
    unittest.main()
