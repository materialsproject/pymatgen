from __future__ import annotations

from shutil import which
from unittest import TestCase

import pytest
import requests

from pymatgen.ext.cod import COD

try:
    website_down = requests.get("https://www.crystallography.net").status_code != 200
except requests.exceptions.ConnectionError:
    website_down = True


@pytest.mark.skipif(website_down, reason="www.crystallography.net is down.")
class TestCOD(TestCase):
    @pytest.mark.skipif(not which("mysql"), reason="No mysql.")
    def test_get_cod_ids(self):
        ids = COD().get_cod_ids("Li2O")
        assert len(ids) > 15

    @pytest.mark.skipif(not which("mysql"), reason="No mysql.")
    def test_get_structure_by_formula(self):
        data = COD().get_structure_by_formula("Li2O")
        assert len(data) > 15
        assert data[0]["structure"].reduced_formula == "Li2O"

    def test_get_structure_by_id(self):
        struct = COD().get_structure_by_id(2002926)
        assert struct.formula == "Be8 H64 N16 F32"
