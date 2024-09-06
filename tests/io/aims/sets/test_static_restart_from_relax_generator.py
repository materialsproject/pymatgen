"""The test of input sets generating from restart information"""

from __future__ import annotations

from pymatgen.io.aims.sets.core import StaticSetGenerator
from pymatgen.util.testing import TEST_FILES_DIR

from ..conftest import O2, Si, comp_system  # noqa: TID252

SPECIES_DIR = TEST_FILES_DIR / "io/aims/species_directory"
REF_PATH = (TEST_FILES_DIR / "io/aims/aims_input_generator_ref").resolve()


def test_static_from_relax_si(tmp_path):
    user_params = {"species_dir": str(SPECIES_DIR / "light")}
    comp_system(
        Si,
        user_params,
        "static-from-prev-si",
        tmp_path,
        REF_PATH,
        StaticSetGenerator,
        properties=["energy", "forces", "stress"],
        prev_dir=f"{REF_PATH}/relax-si/",
    )


def test_static_from_relax_si_no_kgrid(tmp_path):
    user_params = {"species_dir": str(SPECIES_DIR / "light")}
    comp_system(
        Si,
        user_params,
        "static-from-prev-no-kgrid-si",
        tmp_path,
        REF_PATH,
        StaticSetGenerator,
        properties=["energy", "forces", "stress"],
        prev_dir=f"{REF_PATH}/relax-no-kgrid-si/",
    )


def test_static_from_relax_default_species_dir(tmp_path):
    user_params = {"species_dir": str(SPECIES_DIR / "light")}
    comp_system(
        Si,
        user_params,
        "static-from-prev-si",
        tmp_path,
        REF_PATH,
        StaticSetGenerator,
        properties=["energy", "forces", "stress"],
        prev_dir=f"{REF_PATH}/relax-si/",
    )


def test_static_from_relax_o2(tmp_path):
    user_params = {"species_dir": str(SPECIES_DIR / "light")}
    comp_system(
        O2,
        user_params,
        "static-from-prev-o2",
        tmp_path,
        REF_PATH,
        StaticSetGenerator,
        properties=["energy", "forces", "stress"],
        prev_dir=f"{REF_PATH}/relax-o2/",
    )


def test_static_from_relax_default_species_dir_o2(tmp_path):
    user_params = {"species_dir": str(SPECIES_DIR / "light")}
    comp_system(
        O2,
        user_params,
        "static-from-prev-o2",
        tmp_path,
        REF_PATH,
        StaticSetGenerator,
        properties=["energy", "forces", "stress"],
        prev_dir=f"{REF_PATH}/relax-o2/",
    )
