# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import annotations

import os
import unittest

from pymatgen.command_line.chargemol_caller import ChargemolAnalysis
from pymatgen.util.testing import PymatgenTest


class ChargemolAnalysisTest(unittest.TestCase):
    def test_parse_chargemol(self):
        ca = ChargemolAnalysis(
            path=os.path.join(PymatgenTest.TEST_FILES_DIR, "chargemol", "spin_unpolarized"), run_chargemol=False
        )
        self.assertEqual(ca.ddec_charges, [0.8432, -0.8432])
        self.assertEqual(ca.get_partial_charge(0), 0.8432)
        self.assertEqual(ca.get_partial_charge(0, charge_type="cm5"), 0.420172)
        self.assertEqual(ca.get_charge_transfer(0), -0.8432)
        self.assertEqual(ca.get_charge_transfer(0, charge_type="cm5"), -0.420172)
        self.assertEqual(ca.get_charge(0, nelect=1), 1 - 0.8432)
        self.assertEqual(ca.get_charge(0, nelect=1, charge_type="cm5"), 1 - 0.420172)
        self.assertEqual(ca.dipoles, [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        self.assertEqual(ca.ddec_spin_moments, None)
        self.assertEqual(ca.bond_order_sums, [0.53992, 0.901058])
        self.assertEqual(ca.ddec_rsquared_moments, [8.261378, 34.237274])
        self.assertEqual(ca.ddec_rcubed_moments, [14.496002, 88.169236])
        self.assertEqual(ca.ddec_rfourth_moments, [37.648248, 277.371929])
        self.assertEqual(ca.cm5_charges, [0.420172, -0.420172])
        self.assertEqual(ca.summary["ddec"]["partial_charges"], ca.ddec_charges)
        self.assertEqual(ca.summary["ddec"]["dipoles"], ca.dipoles)
        self.assertEqual(ca.summary["ddec"]["bond_order_sums"], ca.bond_order_sums)
        self.assertEqual(ca.summary["ddec"]["rsquared_moments"], ca.ddec_rsquared_moments)
        self.assertEqual(ca.summary["ddec"]["rcubed_moments"], ca.ddec_rcubed_moments)
        self.assertEqual(ca.summary["ddec"]["rfourth_moments"], ca.ddec_rfourth_moments)
        self.assertEqual(ca.summary["cm5"]["partial_charges"], ca.cm5_charges)
        self.assertEqual(ca.summary["ddec"]["bond_order_dict"], ca.bond_order_dict)
        self.assertEqual(ca.summary["ddec"].get("spin_moments"), None)
        self.assertEqual(ca.natoms, [1, 1])
        self.assertTrue(ca.structure is not None)
        self.assertEqual(len(ca.bond_order_dict), 2)
        self.assertEqual(ca.bond_order_dict[0]["bonded_to"][0]["spin_polarization"], 0.0)
        self.assertEqual(ca.bond_order_dict[0]["bonded_to"][0]["index"], 1)
        self.assertEqual(ca.bond_order_dict[0]["bonded_to"][0]["direction"], (-1, -1, 0))
        self.assertEqual(ca.bond_order_dict[0]["bonded_to"][0]["bond_order"], 0.0882)
        self.assertEqual(ca.bond_order_dict[1]["bonded_to"][-1]["direction"], (-1, 0, 0))
        self.assertEqual(ca.get_property_decorated_structure().site_properties["partial_charge_ddec6"], ca.ddec_charges)

    def test_parse_chargemol2(self):
        ca = ChargemolAnalysis(
            path=os.path.join(PymatgenTest.TEST_FILES_DIR, "chargemol", "spin_polarized"), run_chargemol=False
        )
        self.assertEqual(ca.ddec_spin_moments, [0.201595, 0.399203, 0.399203])
        self.assertEqual(ca.summary["ddec"]["bond_order_dict"][0]["bonded_to"][0]["spin_polarization"], 0.0490)
        self.assertEqual(ca.summary["ddec"]["spin_moments"], ca.ddec_spin_moments)
        self.assertTrue(ca.natoms is None)
        self.assertTrue(ca.structure is None)


if __name__ == "__main__":
    unittest.main()
