from __future__ import annotations

from pymatgen.io.aims.sets.core import RelaxSetGenerator
from pymatgen.util.testing import TEST_FILES_DIR

from ..conftest import O2, Si, comp_system  # noqa: TID252

SPECIES_DIR = TEST_FILES_DIR / "io/aims/species_directory"
REF_PATH = TEST_FILES_DIR / "io/aims/aims_input_generator_ref"


def test_relax_si(tmp_path):
    params = {
        "species_dir": "light",
        "k_grid": [2, 2, 2],
    }
    comp_system(Si, params, "relax-si/", tmp_path, REF_PATH, RelaxSetGenerator)


def test_relax_si_no_kgrid(tmp_path):
    params = {"species_dir": "light"}
    comp_system(Si, params, "relax-no-kgrid-si", tmp_path, REF_PATH, RelaxSetGenerator)


def test_relax_o2(tmp_path):
    params = {"species_dir": str(SPECIES_DIR / "light")}
    comp_system(O2, params, "relax-o2", tmp_path, REF_PATH, RelaxSetGenerator)
