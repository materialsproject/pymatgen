"""The test of input sets generating from restart information"""

from __future__ import annotations

from pathlib import Path

from pymatgen.io.aims.sets.core import StaticSetGenerator
from pymatgen.util.testing.aims import O2, Si, comp_system

module_dir = Path(__file__).resolve().parents[1]
species_dir = module_dir / "species_directory"
ref_path = (module_dir / "aims_input_generator_ref").resolve()


def test_static_from_relax_si(tmp_path):
    user_params = {"species_dir": str(species_dir / "light")}
    comp_system(
        Si,
        user_params,
        "static-from-prev-si",
        tmp_path,
        ref_path,
        StaticSetGenerator,
        properties=["energy", "forces", "stress"],
        prev_dir=f"{ref_path}/relax-si/",
    )


def test_static_from_relax_si_no_kgrid(tmp_path):
    user_params = {"species_dir": str(species_dir / "light")}
    comp_system(
        Si,
        user_params,
        "static-from-prev-no-kgrid-si",
        tmp_path,
        ref_path,
        StaticSetGenerator,
        properties=["energy", "forces", "stress"],
        prev_dir=f"{ref_path}/relax-no-kgrid-si/",
    )


def test_static_from_relax_default_species_dir(tmp_path):
    user_params = {"species_dir": str(species_dir / "light")}
    comp_system(
        Si,
        user_params,
        "static-from-prev-si",
        tmp_path,
        ref_path,
        StaticSetGenerator,
        properties=["energy", "forces", "stress"],
        prev_dir=f"{ref_path}/relax-si/",
    )


def test_static_from_relax_o2(tmp_path):
    user_params = {"species_dir": str(species_dir / "light")}
    comp_system(
        O2,
        user_params,
        "static-from-prev-o2",
        tmp_path,
        ref_path,
        StaticSetGenerator,
        properties=["energy", "forces", "stress"],
        prev_dir=f"{ref_path}/relax-o2/",
    )


def test_static_from_relax_default_species_dir_o2(tmp_path):
    user_params = {"species_dir": str(species_dir / "light")}
    comp_system(
        O2,
        user_params,
        "static-from-prev-o2",
        tmp_path,
        ref_path,
        StaticSetGenerator,
        properties=["energy", "forces", "stress"],
        prev_dir=f"{ref_path}/relax-o2/",
    )
