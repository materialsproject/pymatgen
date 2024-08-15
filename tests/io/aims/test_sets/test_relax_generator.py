from __future__ import annotations

from pathlib import Path

from pymatgen.io.aims.sets.core import RelaxSetGenerator
from pymatgen.util.testing.aims import O2, Si, comp_system

MODULE_DIR = Path(__file__).resolve().parents[1]
SPECIES_DIR = MODULE_DIR / "species_directory"
REF_PATH = (MODULE_DIR / "aims_input_generator_ref").resolve()


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
