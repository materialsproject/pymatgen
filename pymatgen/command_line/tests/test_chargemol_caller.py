# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import os
import unittest
from monty.os.path import which

from pymatgen.command_line.chargemol_caller import ChargemolAnalysis
from pymatgen.util.testing import PymatgenTest


class ChargemolAnalysisTest(unittest.TestCase):
    def test_parse_chargemol(self):
        ca = ChargemolAnalysis(path=os.path.join(PymatgenTest.TEST_FILES_DIR, "chargemol"), run_chargemol=False)
        self.assertEqual(ca.ddec_charges, [0.8432, -0.8432])
        self.assertEqual(ca.get_partial_charge(0), 0.8432)
        self.assertEqual(ca.get_charge_transfer(0), -0.8432)
        self.assertEqual(ca.dipoles, [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        self.assertEqual(ca.ddec_spin_moments, None)
        self.assertEqual(ca.bond_order_sums, [0.53992, 0.901058])
        self.assertEqual(ca.ddec_rsquared_moments, [8.261378, 34.237274])
        self.assertEqual(ca.ddec_rcubed_moments, [14.496002, 88.169236])
        self.assertEqual(ca.ddec_rfourth_moments, [37.648248, 277.371929])
        self.assertEqual(ca.cm5_charges, [0.420172, -0.420172])
        # needs tests for ca.summary, ca.bond_order_dict

    @unittest.skipIf(
        not (
            which("Chargemol_09_26_2017_linux_parallel")
            or which("Chargemol_09_26_2017_linux_serial")
            or which("chargemol")
        )
    )
    def test_run_chargemol(self):
        ca = ChargemolAnalysis(path=os.path.join(PymatgenTest.TEST_FILES_DIR, "chargemol"))
        self.assertTrue(os.path.exists(os.path.join(PymatgenTest.TEST_FILES_DIR, "chargemol", "job_control.txt")))


if __name__ == "__main__":
    unittest.main()
