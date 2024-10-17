from __future__ import annotations

import os
from functools import wraps

import pytest
import requests

from pymatgen.ext.cod import COD

# Set a tighter timeout in CI
TIMEOUT = 10 if os.getenv("CI") else 60


try:
    WEBSITE_DOWN = requests.get("https://www.crystallography.net", timeout=TIMEOUT).status_code != 200
except (requests.exceptions.ConnectionError, requests.exceptions.Timeout, requests.exceptions.ReadTimeout):
    WEBSITE_DOWN = True

if WEBSITE_DOWN:
    pytest.skip(reason="www.crystallography.net is down", allow_module_level=True)


def skip_on_timeout(func):
    """Skip test in CI when time out."""

    @wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except (requests.exceptions.ConnectionError, requests.exceptions.Timeout, requests.exceptions.ReadTimeout):
            if os.getenv("CI"):
                pytest.skip("Request timeout in CI environment")
            else:
                raise

    return wrapper


class TestCOD:
    @skip_on_timeout
    def test_get_cod_ids(self):
        ids = COD(timeout=TIMEOUT).get_cod_ids("Li2O")
        assert len(ids) > 15
        assert set(ids).issuperset({1010064, 1011372})

    @skip_on_timeout
    def test_get_structure_by_formula(self):
        data = COD(timeout=TIMEOUT).get_structure_by_formula("Li2O")
        assert len(data) > 15
        assert data[0]["structure"].reduced_formula == "Li2O"

    @skip_on_timeout
    def test_get_structure_by_id(self):
        struct = COD(timeout=TIMEOUT).get_structure_by_id(2_002_926)
        assert struct.formula == "Be8 H64 N16 F32"
