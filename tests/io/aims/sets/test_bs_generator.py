"""Tests the band structure input set generator"""

from __future__ import annotations

from pymatgen.io.aims.sets.bs import BandStructureSetGenerator
from pymatgen.util.testing import TEST_FILES_DIR

from ..conftest import Si, comp_system  # noqa: TID252

SPECIES_DIR = TEST_FILES_DIR / "io/aims/species_directory"
REF_PATH = TEST_FILES_DIR / "io/aims/aims_input_generator_ref"


def test_si_bs(tmp_path):
    parameters = {
        "species_dir": str(SPECIES_DIR / "light"),
        "k_grid": [8, 8, 8],
    }
    comp_system(Si, parameters, "static-si-bs", tmp_path, REF_PATH, BandStructureSetGenerator)


def test_si_bs_output(tmp_path):
    parameters = {
        "species_dir": str(SPECIES_DIR / "light"),
        "k_grid": [8, 8, 8],
        "output": ["json_log"],
    }
    comp_system(Si, parameters, "static-si-bs-output", tmp_path, REF_PATH, BandStructureSetGenerator)


def test_si_bs_density(tmp_path):
    parameters = {
        "species_dir": str(SPECIES_DIR / "light"),
        "k_grid": [8, 8, 8],
        "k_point_density": 40,
    }
    comp_system(Si, parameters, "static-si-bs-density", tmp_path, REF_PATH, BandStructureSetGenerator)
