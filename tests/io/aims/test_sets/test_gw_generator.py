"""Tests the GW input set generator"""

from __future__ import annotations

from pathlib import Path

from pymatgen.io.aims.sets.bs import GWSetGenerator
from pymatgen.util.testing.aims import Si, comp_system

MODULE_DIR = Path(__file__).resolve().parents[1]
SPECIES_DIR = MODULE_DIR / "species_directory"
REF_PATH = (MODULE_DIR / "aims_input_generator_ref").resolve()


def test_si_gw(tmp_path):
    parameters = {
        "species_dir": str(SPECIES_DIR / "light"),
        "k_grid": [2, 2, 2],
        "k_point_density": 10,
    }
    comp_system(Si, parameters, "static-si-gw", tmp_path, REF_PATH, GWSetGenerator)
