"""The test of input sets generating from restart information"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

from pymatgen.io.aims.sets.core import StaticSetGenerator
from pymatgen.util.testing.aims import O2, Si, compare_files


def comp_system(atoms, prev_dir, test_name, work_path, ref_path, species_dir):
    generator = StaticSetGenerator(user_params={})
    # adjust species dir in the prev_dir
    params_file = Path(prev_dir) / "parameters.json"
    shutil.copy(params_file, Path(prev_dir) / "~parameters.json")
    with open(params_file) as f:
        params = json.load(f)
    params["species_dir"] = (species_dir / "light").as_posix()
    with open(params_file, "w") as f:
        json.dump(params, f)

    input_set = generator.get_input_set(atoms, prev_dir, properties=["energy", "forces", "stress"])
    input_set.write_input(work_path / test_name)
    compare_files(test_name, work_path, ref_path)
    shutil.move(Path(prev_dir) / "~parameters.json", params_file)


module_dir = Path(__file__).resolve().parents[1]
species_dir = module_dir / "species_directory"
ref_path = (module_dir / "aims_input_generator_ref").resolve()


def test_static_from_relax_si(tmp_path):
    comp_system(
        Si,
        f"{ref_path}/relax-si/",
        "static-from-prev-si",
        tmp_path,
        ref_path,
        species_dir,
    )


def test_static_from_relax_si_no_kgrid(tmp_path):
    comp_system(
        Si,
        f"{ref_path}/relax-no-kgrid-si/",
        "static-from-prev-no-kgrid-si",
        tmp_path,
        ref_path,
        species_dir,
    )


def test_static_from_relax_default_species_dir(tmp_path, monkeypatch):
    monkeypatch.setenv("AIMS_SPECIES_DIR", str(species_dir / "light"))

    comp_system(
        Si,
        f"{ref_path}/relax-si/",
        "static-from-prev-si",
        tmp_path,
        ref_path,
        species_dir,
    )


def test_static_from_relax_o2(tmp_path):
    comp_system(
        O2,
        f"{ref_path}/relax-o2/",
        "static-from-prev-o2",
        tmp_path,
        ref_path,
        species_dir,
    )


def test_static_from_relax_default_species_dir_o2(tmp_path, monkeypatch):
    monkeypatch.setenv("AIMS_SPECIES_DIR", str(species_dir / "light"))

    comp_system(
        O2,
        f"{ref_path}/relax-o2/",
        "static-from-prev-o2",
        tmp_path,
        ref_path,
        species_dir,
    )
