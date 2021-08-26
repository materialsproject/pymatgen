# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import unittest
import warnings

import requests

from monty.os.path import which

from pymatgen.ext.cod import COD

website_is_up = requests.get("https://www.crystallography.net").status_code == 200


@unittest.skipIf(not website_is_up, "www.crystallography.net is down.")
class CODTest(unittest.TestCase):
    _multiprocess_shared_ = True

    def setUp(self):
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    @unittest.skipIf(not which("mysql"), "No mysql.")
    def test_get_cod_ids(self):
        ids = COD().get_cod_ids("Li2O")
        self.assertTrue(len(ids) > 15)

    @unittest.skipIf(not which("mysql"), "No mysql.")
    def test_get_structure_by_formula(self):
        data = COD().get_structure_by_formula("Li2O")
        self.assertTrue(len(data) > 15)
        self.assertEqual(data[0]["structure"].composition.reduced_formula, "Li2O")

    def test_get_structure_by_id(self):
        s = COD().get_structure_by_id(2002926)
        self.assertEqual(s.formula, "Be8 H64 N16 F32")


if __name__ == "__main__":
    unittest.main()
