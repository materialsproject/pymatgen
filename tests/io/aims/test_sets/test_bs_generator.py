"""Tests the band structure input set generator"""
from __future__ import annotations

import gzip
import json
from glob import glob
from pathlib import Path

import numpy as np

from pymatgen.core import Structure
from pymatgen.io.aims.sets.bs import BandStructureSetGenerator


def check_band(test_line, ref_line):
    test_pts = [float(inp) for inp in test_line.split()[-9:-2]]
    ref_pts = [float(inp) for inp in ref_line.split()[-9:-2]]

    print(test_line)
    print(ref_line)
    return np.allclose(test_pts, ref_pts) and test_line.split()[-2:] == ref_line.split()[-2:]


def compare_files(test_name, work_dir, ref_dir):
    for file in glob(f"{work_dir / test_name}/*in"):
        with open(file) as test_file:
            test_lines = [line.strip() for line in test_file.readlines()[4:] if len(line.strip()) > 0]

        with gzip.open(f"{ref_dir / test_name / Path(file).name}.gz", "rt") as ref_file:
            ref_lines = [line.strip() for line in ref_file.readlines()[4:] if len(line.strip()) > 0]

        for test_line, ref_line in zip(test_lines, ref_lines):
            if "output" in test_line and "band" in test_line:
                assert check_band(test_line, ref_line)
            else:
                assert test_line == ref_line

    with open(f"{ref_dir / test_name}/parameters.json") as ref_file:
        ref = json.load(ref_file)
    ref.pop("species_dir", None)
    ref_output = ref.pop("output", None)

    with open(f"{work_dir / test_name}/parameters.json") as check_file:
        check = json.load(check_file)

    check.pop("species_dir", None)
    check_output = check.pop("output", None)

    assert ref == check

    if check_output:
        for ref_out, check_out in zip(ref_output, check_output):
            if "band" in check_out:
                assert check_band(check_out, ref_out)
            else:
                assert ref_out == check_out


def comp_system(atoms, user_params, test_name, work_path, ref_path):
    k_point_density = user_params.pop("k_point_density", 20)
    generator = BandStructureSetGenerator(user_parameters=user_params, k_point_density=k_point_density)
    input_set = generator.get_input_set(atoms)
    input_set.write_input(work_path / test_name)
    compare_files(test_name, work_path, ref_path)


Si = Structure(
    lattice=[[0.0, 2.715, 2.715], [2.715, 0.0, 2.715], [2.715, 2.715, 0.0]],
    species=["Si", "Si"],
    coords=[[0, 0, 0], [0.25, 0.25, 0.25]],
)
species_dir = Path(__file__).resolve().parents[1] / "species_directory"

module_dir = Path(__file__).resolve().parents[1]
ref_path = (module_dir / "aims_input_generator_ref").resolve()


def test_si_bs(tmp_path):
    parameters = {
        "species_dir": str(species_dir / "light"),
        "k_grid": [8, 8, 8],
    }
    comp_system(Si, parameters, "static-si-bs", tmp_path, ref_path)


def test_si_bs_output(tmp_path):
    parameters = {
        "species_dir": str(species_dir / "light"),
        "k_grid": [8, 8, 8],
        "output": [
            "json_log",
        ],
    }
    comp_system(Si, parameters, "static-si-bs-output", tmp_path, ref_path)


def test_si_bs_density(tmp_path):
    parameters = {
        "species_dir": str(species_dir / "light"),
        "k_grid": [8, 8, 8],
        "k_point_density": 40,
    }
    comp_system(Si, parameters, "static-si-bs-density", tmp_path, ref_path)
