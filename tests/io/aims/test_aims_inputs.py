from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
from monty.json import MontyDecoder

from pymatgen.io.aims.inputs import AimsGeometryIn

infile_dir = Path(__file__).parent / "aims_input_files"


def test_read_write_si_in(tmpdir):
    si = AimsGeometryIn.from_file(infile_dir / "geometry.in.si")

    in_lattice = np.array([[0.0, 2.715, 2.716], [2.717, 0.0, 2.718], [2.719, 2.720, 0.0]])
    in_coords = np.array([[0.0, 0.0, 0.0], [0.25, 0.24, 0.26]])

    assert all(sp.symbol == "Si" for sp in si.structure.species)
    assert np.allclose(si.structure.lattice.matrix, in_lattice)
    assert np.allclose(si.structure.frac_coords.flatten(), in_coords.flatten())

    si_test_from_struct = AimsGeometryIn.from_structure(si.structure)
    assert si.structure == si_test_from_struct.structure

    workdir = Path.cwd()
    os.chdir(tmpdir)
    si_test_from_struct.write_file(overwrite=True)

    with open("geometry.in") as test_file:
        test_lines = test_file.readlines()[5:]

    with open(infile_dir / "geometry.in.si.ref") as ref_file:
        ref_lines = ref_file.readlines()[5:]

    for test_line, ref_line in zip(test_lines, ref_lines):
        assert test_line == ref_line

    os.chdir(workdir)

    with open(f"{infile_dir}/si_ref.json") as si_ref_json:
        si_from_dct = json.load(si_ref_json, cls=MontyDecoder)

    assert si.structure == si_from_dct.structure


def test_read_h2o_in(tmpdir):
    h2o = AimsGeometryIn.from_file(infile_dir / "geometry.in.h2o")

    in_coords = np.array([[0.0, 0.0, 0.119262], [0.0, 0.763239, -0.477047], [0.0, -0.763239, -0.477047]])

    assert all(sp.symbol == symb for sp, symb in zip(h2o.structure.species, ["O", "H", "H"]))
    assert np.allclose(h2o.structure.cart_coords.flatten(), in_coords.flatten())

    h2o_test_from_struct = AimsGeometryIn.from_structure(h2o.structure)
    assert h2o.structure == h2o_test_from_struct.structure

    workdir = Path.cwd()
    os.chdir(tmpdir)
    h2o_test_from_struct.write_file(overwrite=True)

    with open("geometry.in") as test_file:
        test_lines = test_file.readlines()[5:]

    with open(infile_dir / "geometry.in.h2o.ref") as ref_file:
        ref_lines = ref_file.readlines()[5:]

    for test_line, ref_line in zip(test_lines, ref_lines):
        assert test_line == ref_line

    os.chdir(workdir)

    with open(f"{infile_dir}/h2o_ref.json") as h2o_ref_json:
        h2o_from_dct = json.load(h2o_ref_json, cls=MontyDecoder)

    assert h2o.structure == h2o_from_dct.structure
