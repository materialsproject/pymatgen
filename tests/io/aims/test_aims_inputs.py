from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pytest
from monty.json import MontyDecoder, MontyEncoder

from pymatgen.io.aims.inputs import ALLOWED_AIMS_CUBE_TYPES, ALLOWED_AIMS_CUBE_TYPES_STATE, AimsCube, AimsGeometryIn

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

    in_coords = np.array(
        [
            [0.0, 0.0, 0.119262],
            [0.0, 0.763239, -0.477047],
            [0.0, -0.763239, -0.477047],
        ]
    )

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


def check_wrong_type_aims_cube(type, exp_err):
    with pytest.raises(ValueError, match=exp_err):
        AimsCube(type=type)


def test_aims_cube(tmpdir):
    check_wrong_type_aims_cube(type="INCORRECT_TYPE", exp_err="Cube type undefined")

    for type in ALLOWED_AIMS_CUBE_TYPES_STATE:
        check_wrong_type_aims_cube(
            type=type,
            exp_err=f"Cube of type {type} must have a state associated with it",
        )

    for type in ALLOWED_AIMS_CUBE_TYPES:
        check_wrong_type_aims_cube(
            type=f"{type} 1",
            exp_err=f"Cube of type {type} can not have a state associated with it",
        )

    with pytest.raises(
        ValueError,
        match="Cube file must have a format of cube, gOpenMol, or xsf",
    ):
        AimsCube(type=ALLOWED_AIMS_CUBE_TYPES[0], format="TEST_ERR")

    with pytest.raises(ValueError, match="Spin state must be 1 or 2"):
        AimsCube(type=ALLOWED_AIMS_CUBE_TYPES[0], spinstate=3)

    with pytest.raises(ValueError, match="The cube origin must have 3 components"):
        AimsCube(type=ALLOWED_AIMS_CUBE_TYPES[0], origin=[0])

    with pytest.raises(ValueError, match="Only three cube edges can be passed"):
        AimsCube(type=ALLOWED_AIMS_CUBE_TYPES[0], edges=[[0.0, 0.0, 0.1]])

    with pytest.raises(ValueError, match="Each cube edge must have 3 components"):
        AimsCube(
            type=ALLOWED_AIMS_CUBE_TYPES[0],
            edges=[[0.0, 0.0, 0.1], [0.1, 0.0, 0.0], [0.1, 0.0]],
        )

    with pytest.raises(ValueError, match="elf_type only used when the cube type is elf"):
        AimsCube(type=ALLOWED_AIMS_CUBE_TYPES[0], elf_type=1)

    with pytest.raises(ValueError, match="The number of points per edge must have 3 components"):
        AimsCube(type=ALLOWED_AIMS_CUBE_TYPES[0], points=[100, 100, 100, 100])

    test_cube = AimsCube(
        type="elf",
        origin=[0, 0, 0],
        edges=[[0.01, 0, 0], [0.0, 0.01, 0], [0.0, 0, 0.01]],
        points=[100, 100, 100],
        spinstate=1,
        kpoint=1,
        filename="test.cube",
        format="cube",
        elf_type=1,
    )

    test_cube_block = [
        "output cube elf",
        "    cube origin [0, 0, 0]",
        "    cube edge 100  1.000000000000e-02  0.000000000000e+00  0.000000000000e+00",
        "    cube edge 100  0.000000000000e+00  1.000000000000e-02  0.000000000000e+00",
        "    cube edge 100  0.000000000000e+00  0.000000000000e+00  1.000000000000e-02",
        "    cube format cube",
        "    cube spinstate 1",
        "    cube kpoint 1",
        "    cube filename test.cube",
        "    cube elf_type 1",
        "",
    ]
    assert test_cube.control_block == "\n".join(test_cube_block)

    test_cube_from_dict = json.loads(json.dumps(test_cube.as_dict(), cls=MontyEncoder), cls=MontyDecoder)
    assert test_cube_from_dict.control_block == test_cube.control_block
