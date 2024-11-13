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
except (
    requests.exceptions.ConnectionError,
    requests.exceptions.Timeout,
    requests.exceptions.ReadTimeout,
):
    WEBSITE_DOWN = True

if WEBSITE_DOWN:
    pytest.skip(reason="www.crystallography.net is down", allow_module_level=True)


def skip_on_timeout(func):
    """Skip test in CI when time out."""

    @wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except (
            requests.exceptions.ConnectionError,
            requests.exceptions.Timeout,
            requests.exceptions.ReadTimeout,
        ):
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
        # This formula has only one match (as of 2024-10-17) therefore
        # the runtime is shorter (~ 2s for each match)
        data = COD(timeout=TIMEOUT).get_structure_by_formula("C3 H18 F6 Fe N9")
        assert len(data) >= 1
        assert data[0]["structure"].reduced_formula == "FeH18C3(N3F2)3"

    @skip_on_timeout
    def test_get_structure_by_id(self):
        struct = COD(timeout=TIMEOUT).get_structure_by_id(2_002_926)
        assert struct.formula == "Be8 H64 N16 F32"
