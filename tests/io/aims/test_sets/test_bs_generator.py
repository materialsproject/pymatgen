"""Tests the band structure input set generator"""

from __future__ import annotations

from pathlib import Path

from pymatgen.io.aims.sets.bs import BandStructureSetGenerator
from pymatgen.util.testing.aims import Si, comp_system

module_dir = Path(__file__).resolve().parents[1]
species_dir = module_dir / "species_directory"
ref_path = (module_dir / "aims_input_generator_ref").resolve()


def test_si_bs(tmp_path):
    parameters = {
        "species_dir": str(species_dir / "light"),
        "k_grid": [8, 8, 8],
    }
    comp_system(Si, parameters, "static-si-bs", tmp_path, ref_path, BandStructureSetGenerator)


def test_si_bs_output(tmp_path):
    parameters = {
        "species_dir": str(species_dir / "light"),
        "k_grid": [8, 8, 8],
        "output": ["json_log"],
    }
    comp_system(Si, parameters, "static-si-bs-output", tmp_path, ref_path, BandStructureSetGenerator)


def test_si_bs_density(tmp_path):
    parameters = {
        "species_dir": str(species_dir / "light"),
        "k_grid": [8, 8, 8],
        "k_point_density": 40,
    }
    comp_system(Si, parameters, "static-si-bs-density", tmp_path, ref_path, BandStructureSetGenerator)
