from __future__ import annotations

import os

import pytest
import requests

from pymatgen.ext.cod import COD

TIMEOUT = 10  # use a tighter timeout in CI

try:
    WEBSITE_DOWN = requests.get("https://www.crystallography.net", timeout=TIMEOUT).status_code != 200
except (requests.exceptions.ConnectionError, requests.exceptions.Timeout, requests.exceptions.ReadTimeout):
    WEBSITE_DOWN = True

if WEBSITE_DOWN:
    pytest.skip(reason="www.crystallography.net is down", allow_module_level=True)


class TestCOD:
    def test_get_cod_ids(self):
        try:
            ids = COD(timeout=TIMEOUT).get_cod_ids("Li2O")
        except (requests.exceptions.ConnectionError, requests.exceptions.Timeout, requests.exceptions.ReadTimeout):
            if "CI" in os.environ:
                pytest.skip("request timeout")
            else:
                raise

        assert len(ids) > 15
        assert set(ids).issuperset({1010064, 1011372})

    def test_get_structure_by_formula(self):
        try:
            data = COD(timeout=TIMEOUT).get_structure_by_formula("Li2O")
        except (requests.exceptions.ConnectionError, requests.exceptions.Timeout, requests.exceptions.ReadTimeout):
            if "CI" in os.environ:
                pytest.skip("request timeout")
            else:
                raise

        assert len(data) > 15
        assert data[0]["structure"].reduced_formula == "Li2O"

    def test_get_structure_by_id(self):
        try:
            struct = COD(timeout=TIMEOUT).get_structure_by_id(2_002_926)
        except (requests.exceptions.ConnectionError, requests.exceptions.Timeout, requests.exceptions.ReadTimeout):
            if "CI" in os.environ:
                pytest.skip("request timeout")
            else:
                raise

        assert struct.formula == "Be8 H64 N16 F32"
