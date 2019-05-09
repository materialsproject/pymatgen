# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import unittest
from pymatgen.ext.crystalsai import CrystalAIRester
from pymatgen.util.testing import PymatgenTest


class CrystalAIResterTest(PymatgenTest):

    def test_query(self):
        m = CrystalAIRester()
        self.assertAlmostEqual(m.predict_mp("formation_energy", "mp-1143"),
                               -3.446291923522949)
        s = PymatgenTest.get_structure("Li2O")
        self.assertAlmostEqual(m.predict_structure("formation_energy", s),
                               -2.015296220779419)


if __name__ == "__main__":
    unittest.main()