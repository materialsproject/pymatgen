# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import unittest

import requests

from pymatgen.ext.crystalsai import CrystalAIRester
from pymatgen.util.testing import PymatgenTest


website_is_up = requests.get("http://megnet.crystals.ai").status_code == 200


@unittest.skipIf(not website_is_up, "megnet.crystals.ai is down.")
class CrystalAIResterTest(PymatgenTest):
    def test_query(self):
        with CrystalAIRester() as m:
            models = m.get_available_models()
            self.assertIn("formation_energy", models)
            self.assertIn("band_gap", models)
            self.assertIn("log10K", models)

            self.assertAlmostEqual(m.predict_mp("formation_energy", "mp-1143"), -3.446291923522949, 4)
            s = PymatgenTest.get_structure("Li2O")
            self.assertAlmostEqual(m.predict_structure("formation_energy", s), -2.015296220779419, 4)


if __name__ == "__main__":
    unittest.main()
