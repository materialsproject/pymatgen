from __future__ import annotations

import gzip
import json
from glob import glob
from pathlib import Path

import numpy as np

from pymatgen.core import Structure


def check_band(test_line, ref_line):
    test_pts = [float(inp) for inp in test_line.split()[-9:-2]]
    ref_pts = [float(inp) for inp in ref_line.split()[-9:-2]]

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


Si = Structure(
    lattice=[[0.0, 2.715, 2.715], [2.715, 0.0, 2.715], [2.715, 2.715, 0.0]],
    species=["Si", "Si"],
    coords=[[0, 0, 0], [0.25, 0.25, 0.25]],
)

module_dir = Path(__file__).resolve().parents[1]
species_dir = module_dir / "species_directory"
ref_path = (module_dir / "aims_input_generator_ref").resolve()
