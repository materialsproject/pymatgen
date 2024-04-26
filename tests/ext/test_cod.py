from __future__ import annotations

import os
from shutil import which
from unittest import TestCase

import pytest
import requests

from pymatgen.ext.cod import COD

if "CI" in os.environ:  # test is slow and flaky, skip in CI. see
    # https://github.com/materialsproject/pymatgen/pull/3777#issuecomment-2071217785
    pytest.skip(allow_module_level=True)

try:
    website_down = requests.get("https://www.crystallography.net", timeout=600).status_code != 200
except requests.exceptions.ConnectionError:
    website_down = True


@pytest.mark.skipif(website_down, reason="www.crystallography.net is down")
class TestCOD(TestCase):
    @pytest.mark.skipif(not which("mysql"), reason="No mysql")
    def test_get_cod_ids(self):
        ids = COD().get_cod_ids("Li2O")
        assert len(ids) > 15

    @pytest.mark.skipif(not which("mysql"), reason="No mysql")
    def test_get_structure_by_formula(self):
        data = COD().get_structure_by_formula("Li2O")
        assert len(data) > 15
        assert data[0]["structure"].reduced_formula == "Li2O"

    def test_get_structure_by_id(self):
        struct = COD().get_structure_by_id(2_002_926)
        assert struct.formula == "Be8 H64 N16 F32"
