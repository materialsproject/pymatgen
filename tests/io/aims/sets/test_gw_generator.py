"""Tests the GW input set generator"""

from __future__ import annotations

from pymatgen.io.aims.sets.bs import GWSetGenerator
from pymatgen.util.testing import TEST_FILES_DIR

from ..conftest import Si, comp_system  # noqa: TID252

SPECIES_DIR = TEST_FILES_DIR / "io/aims/species_directory"
REF_PATH = TEST_FILES_DIR / "io/aims/aims_input_generator_ref"


def test_si_gw(tmp_path):
    parameters = {
        "species_dir": str(SPECIES_DIR / "light"),
        "k_grid": [2, 2, 2],
        "k_point_density": 10,
    }
    comp_system(Si, parameters, "static-si-gw", tmp_path, REF_PATH, GWSetGenerator)
