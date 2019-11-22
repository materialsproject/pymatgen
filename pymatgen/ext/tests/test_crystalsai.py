# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import unittest
from pymatgen.ext.crystalsai import CrystalAIRester
from pymatgen.util.testing import PymatgenTest


class CrystalAIResterTest(PymatgenTest):

    def test_query(self):
        m = CrystalAIRester()
        models = m.get_available_models()
        self.assertIn("formation_energy", models)
        self.assertIn("band_gap", models)
        self.assertIn("log10K", models)

        self.assertAlmostEqual(m.predict_mp("formation_energy", "mp-1143"),
                               -3.446291923522949, 4)
        s = PymatgenTest.get_structure("Li2O")
        self.assertAlmostEqual(m.predict_structure("formation_energy", s),
                               -2.015296220779419, 4)


if __name__ == "__main__":
    unittest.main()
