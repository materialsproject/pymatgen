from __future__ import annotations

from pathlib import Path

from helpers.aims import O2, Si, comp_system

from pymatgen.io.aims.sets.core import RelaxSetGenerator

module_dir = Path(__file__).resolve().parents[1]
species_dir = module_dir / "species_directory"
ref_path = (module_dir / "aims_input_generator_ref").resolve()


def test_relax_si(tmp_path):
    params = {
        "species_dir": str(species_dir / "light"),
        "k_grid": [2, 2, 2],
    }
    comp_system(Si, params, "relax-si/inputs", tmp_path, ref_path, RelaxSetGenerator)


def test_relax_si_no_kgrid(tmp_path):
    params = {"species_dir": str(species_dir / "light")}
    comp_system(Si, params, "relax-no-kgrid-si", tmp_path, ref_path, RelaxSetGenerator)


def test_relax_default_species_dir(tmp_path, monkeypatch):
    monkeypatch.setenv("AIMS_SPECIES_DIR", str(species_dir / "light"))
    params = {"k_grid": [2, 2, 2]}

    comp_system(Si, params, "relax-default-sd-si", tmp_path, ref_path, RelaxSetGenerator)


def test_relax_o2(tmp_path):
    params = {"species_dir": str(species_dir / "light")}
    comp_system(O2, params, "relax-o2", tmp_path, ref_path, RelaxSetGenerator)


def test_relax_default_species_dir_o2(tmp_path, monkeypatch):
    monkeypatch.setenv("AIMS_SPECIES_DIR", str(species_dir / "light"))
    params = {"k_grid": [2, 2, 2]}

    comp_system(O2, params, "relax-default-sd-o2", tmp_path, ref_path, RelaxSetGenerator)
