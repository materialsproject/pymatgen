from __future__ import annotations

import pytest

from pymatgen.core import SETTINGS
from pymatgen.util.testing import TEST_FILES_DIR


@pytest.fixture(autouse=True)
def _set_aims_species_dir_env_var(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("AIMS_SPECIES_DIR", f"{TEST_FILES_DIR}/io/aims/species_directory")


@pytest.fixture(autouse=True)
def _set_aims_species_dir_settings(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setitem(SETTINGS, "AIMS_SPECIES_DIR", f"{TEST_FILES_DIR}/io/aims/species_directory")
