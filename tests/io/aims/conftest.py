from __future__ import annotations

from pathlib import Path

import pytest

module_dir = Path(__file__).resolve().parent
species_dir = module_dir / "species_directory"


@pytest.fixture(autouse=True)
def _aims_species_dir(monkeypatch):
    monkeypatch.setenv("AIMS_SPECIES_DIR", str(species_dir / "light"))
