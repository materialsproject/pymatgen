from __future__ import annotations

import unittest
from shutil import which

import requests

from pymatgen.ext.cod import COD

try:
    website_down = requests.get("https://www.crystallography.net").status_code != 200
except requests.exceptions.ConnectionError:
    website_down = True


@unittest.skipIf(website_down, "www.crystallography.net is down.")
class TestCOD(unittest.TestCase):
    @unittest.skipIf(not which("mysql"), "No mysql.")
    def test_get_cod_ids(self):
        ids = COD().get_cod_ids("Li2O")
        assert len(ids) > 15

    @unittest.skipIf(not which("mysql"), "No mysql.")
    def test_get_structure_by_formula(self):
        data = COD().get_structure_by_formula("Li2O")
        assert len(data) > 15
        assert data[0]["structure"].composition.reduced_formula == "Li2O"

    def test_get_structure_by_id(self):
        struct = COD().get_structure_by_id(2002926)
        assert struct.formula == "Be8 H64 N16 F32"
