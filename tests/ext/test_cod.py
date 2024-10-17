from __future__ import annotations

from unittest import TestCase

import pytest
import requests
import urllib3

from pymatgen.ext.cod import COD

TIMEOUT = 10  # use a tighter timeout in CI
# TODO: suppress timeout error in CI

try:
    WEBSITE_DOWN = requests.get("https://www.crystallography.net", timeout=TIMEOUT).status_code != 200
except (requests.exceptions.ConnectionError, urllib3.exceptions.ConnectTimeoutError):
    WEBSITE_DOWN = True

if WEBSITE_DOWN:
    pytest.skip(reason="www.crystallography.net is down", allow_module_level=True)


class TestCOD(TestCase):
    def test_get_cod_ids(self):
        ids = COD(timeout=TIMEOUT).get_cod_ids("Li2O")
        assert len(ids) > 15
        assert set(ids).issuperset({1010064, 1011372})

    def test_get_structure_by_formula(self):
        data = COD(timeout=TIMEOUT).get_structure_by_formula("Li2O")
        assert len(data) > 15
        assert data[0]["structure"].reduced_formula == "Li2O"

    def test_get_structure_by_id(self):
        struct = COD(timeout=TIMEOUT).get_structure_by_id(2_002_926)
        assert struct.formula == "Be8 H64 N16 F32"
