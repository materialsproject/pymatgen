"""Tests the GW input set generator"""

from __future__ import annotations

from pathlib import Path

from pymatgen.io.aims.sets.bs import GWSetGenerator
from pymatgen.util.testing.aims import Si, comp_system

module_dir = Path(__file__).resolve().parents[1]
species_dir = module_dir / "species_directory"
ref_path = (module_dir / "aims_input_generator_ref").resolve()


def test_si_gw(tmp_path):
    parameters = {
        "species_dir": str(species_dir / "light"),
        "k_grid": [2, 2, 2],
        "k_point_density": 10,
    }
    comp_system(Si, parameters, "static-si-gw", tmp_path, ref_path, GWSetGenerator)
