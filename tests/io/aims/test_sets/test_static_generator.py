from __future__ import annotations

import gzip
import json
import os
from glob import glob
from pathlib import Path

from pymatgen.io.aims.sets.core import StaticSetGenerator


def compare_files(test_name, work_dir, ref_dir):
    for file in glob(f"{work_dir / test_name}/*in"):
        with open(file) as test_file:
            test_lines = [line.strip() for line in test_file.readlines()[4:] if len(line.strip()) > 0]

        with gzip.open(f"{ref_dir / test_name / Path(file).name}.gz", "rt") as ref_file:
            ref_lines = [line.strip() for line in ref_file.readlines()[4:] if len(line.strip()) > 0]

        assert test_lines == ref_lines

    with open(f"{ref_dir / test_name}/parameters.json") as ref_file:
        ref = json.load(ref_file)
    ref.pop("species_dir", None)

    with open(f"{work_dir / test_name}/parameters.json") as check_file:
        check = json.load(check_file)
    check.pop("species_dir", None)

    assert ref == check


def comp_system(atoms, user_params, test_name, work_path, ref_path):
    generator = StaticSetGenerator(user_parameters=user_params)
    input_set = generator.get_input_set(atoms)
    input_set.write_input(work_path / test_name)
    compare_files(test_name, work_path, ref_path)


def test_static_si(Si, species_dir, tmp_path, ref_path):
    parameters = {
        "species_dir": str(species_dir / "light"),
        "k_grid": [2, 2, 2],
    }
    comp_system(Si, parameters, "static-si", tmp_path, ref_path)


def test_static_si_no_kgrid(Si, species_dir, tmp_path, ref_path):
    parameters = {"species_dir": str(species_dir / "light")}
    comp_system(Si, parameters, "static-no-kgrid-si", tmp_path, ref_path)


def test_static_default_species_dir(Si, species_dir, tmp_path, ref_path):
    sd_def = os.getenv("AIMS_SPECIES_DIR", None)
    os.environ["AIMS_SPECIES_DIR"] = str(species_dir / "light")
    parameters = {"k_grid": [2, 2, 2]}

    comp_system(Si, parameters, "static-default-sd-si", tmp_path, ref_path)

    if sd_def:
        os.environ["AIMS_SPECIES_DIR"] = sd_def
    else:
        os.unsetenv("AIMS_SPECIES_DIR")


def test_static_o2(O2, species_dir, tmp_path, ref_path):
    parameters = {"species_dir": str(species_dir / "light")}
    comp_system(O2, parameters, "static-o2", tmp_path, ref_path)


def test_static_default_species_dir_o2(O2, species_dir, tmp_path, ref_path):
    sd_def = os.getenv("AIMS_SPECIES_DIR", None)
    os.environ["AIMS_SPECIES_DIR"] = str(species_dir / "light")
    parameters = {"k_grid": [2, 2, 2]}

    comp_system(O2, parameters, "static-default-sd-o2", tmp_path, ref_path)

    if sd_def:
        os.environ["AIMS_SPECIES_DIR"] = sd_def
    else:
        os.unsetenv("AIMS_SPECIES_DIR")
