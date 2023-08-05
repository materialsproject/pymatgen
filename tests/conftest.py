"""
Test fixes for pymatgen
"""

from __future__ import annotations

from pathlib import Path

import pytest

from pymatgen.core import SETTINGS

TEST_DIR = Path(__file__).absolute().parent


@pytest.fixture(scope="session", autouse=True)
def TEST_FILES_DIR():
    return Path(SETTINGS.get("PMG_TEST_FILES_DIR", TEST_DIR / "files"))
