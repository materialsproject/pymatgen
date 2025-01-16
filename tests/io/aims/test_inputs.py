from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from monty.json import MontyDecoder, MontyEncoder
from pyfhiaims.control.cube import AimsCube

from pymatgen.io.aims.inputs import AimsControlIn, AimsGeometryIn

from .conftest import compare_single_files as compare_files

TEST_DIR = Path("/home/purcellt/git/pymatgen/tests/files/io/aims/input_files")


# def test_read_write_si_in(tmp_path: Path):
#     si = AimsGeometryIn.from_file(TEST_DIR / "geometry.in.si.gz")

#     in_lattice = np.array([[0, 2.715, 2.716], [2.717, 0, 2.718], [2.719, 2.720, 0]])
#     in_coords = np.array([[0, 0, 0], [0.25, 0.24, 0.26]])

#     assert all(sp.symbol == "Si" for sp in si.structure.species)
#     assert_allclose(si.structure.lattice.matrix, in_lattice)
#     assert_allclose(si.structure.frac_coords.flatten(), in_coords.flatten())

#     si_test_from_struct = AimsGeometryIn.from_structure(si.structure)
#     assert si.structure == si_test_from_struct.structure

#     si_test_from_struct.write_file(directory=tmp_path, overwrite=True)
#     with pytest.raises(ValueError, match="geometry.in file exists in "):
#         si_test_from_struct.write_file(directory=tmp_path, overwrite=False)

#     compare_files(TEST_DIR / "geometry.in.si.ref", f"{tmp_path}/geometry.in")

#     si.structure.to(tmp_path / "si.in", fmt="aims")
#     compare_files(TEST_DIR / "geometry.in.si.ref", f"{tmp_path}/si.in")

#     si_from_file = Structure.from_file(f"{tmp_path}/geometry.in")
#     assert all(sp.symbol == "Si" for sp in si_from_file.species)

#     with gzip.open(f"{TEST_DIR}/si_ref.json.gz", mode="rt") as si_ref_json:
#         si_from_dct = json.load(si_ref_json, cls=MontyDecoder)

#     assert si.structure == si_from_dct.structure


# def test_read_h2o_in(tmp_path: Path):
#     h2o = AimsGeometryIn.from_file(TEST_DIR / "geometry.in.h2o.gz")

#     in_coords = [
#         [0, 0, 0.119262],
#         [0, 0.763239, -0.477047],
#         [0, -0.763239, -0.477047],
#     ]

#     assert all(sp.symbol == symb for sp, symb in zip(h2o.structure.species, ["O", "H", "H"], strict=True))
#     assert_allclose(h2o.structure.cart_coords, in_coords)

#     h2o_test_from_struct = AimsGeometryIn.from_structure(h2o.structure)
#     assert h2o.structure == h2o_test_from_struct.structure

#     h2o_test_from_struct.write_file(directory=tmp_path, overwrite=True)

#     with pytest.raises(ValueError, match="geometry.in file exists in "):
#         h2o_test_from_struct.write_file(directory=tmp_path, overwrite=False)

#     compare_files(TEST_DIR / "geometry.in.h2o.ref", f"{tmp_path}/geometry.in")

#     with gzip.open(f"{TEST_DIR}/h2o_ref.json.gz", mode="rt") as h2o_ref_json:
#         h2o_from_dct = json.load(h2o_ref_json, cls=MontyDecoder)

#     assert h2o.structure == h2o_from_dct.structure


# def test_write_spins(tmp_path: Path):
#     mg2mn4o8 = Structure(
#         lattice=Lattice(
#             [
#                 [5.06882343, 0.00012488, -2.66110167],
#                 [-1.39704234, 4.87249911, -2.66110203],
#                 [0.00986091, 0.01308528, 6.17649359],
#             ]
#         ),
#         species=[
#             "Mg",
#             "Mg",
#             Species("Mn", spin=5.0),
#             Species("Mn", spin=5.0),
#             Species("Mn", spin=5.0),
#             Species("Mn", spin=5.0),
#             "O",
#             "O",
#             "O",
#             "O",
#             "O",
#             "O",
#             "O",
#             "O",
#         ],
#         coords=[
#             [0.37489726, 0.62510274, 0.75000002],
#             [0.62510274, 0.37489726, 0.24999998],
#             [-0.00000000, -0.00000000, 0.50000000],
#             [-0.00000000, 0.50000000, 0.00000000],
#             [0.50000000, -0.00000000, 0.50000000],
#             [-0.00000000, -0.00000000, 0.00000000],
#             [0.75402309, 0.77826750, 0.50805882],
#             [0.77020285, 0.24594779, 0.99191316],
#             [0.22173254, 0.24597689, 0.99194116],
#             [0.24597691, 0.22173250, 0.49194118],
#             [0.24594765, 0.77020288, 0.49191313],
#             [0.22979715, 0.75405221, 0.00808684],
#             [0.75405235, 0.22979712, 0.50808687],
#             [0.77826746, 0.75402311, 0.00805884],
#         ],
#     )

#     geo_in = AimsGeometryIn.from_structure(mg2mn4o8)

#     magmom_lines = [line for line in geo_in.content.split("\n") if "initial_moment" in line]
#     assert len(magmom_lines) == 4

#     magmoms = np.array([float(line.strip().split()[-1]) for line in magmom_lines])
#     assert np.all(magmoms == 5.0)

#     mg2mn4o8 = Structure(
#         lattice=mg2mn4o8.lattice,
#         species=mg2mn4o8.species,
#         coords=mg2mn4o8.frac_coords,
#         site_properties={"magmom": np.zeros(mg2mn4o8.num_sites)},
#     )
#     with pytest.raises(
#         ValueError,
#         match="species.spin and magnetic moments don't agree. Please only define one",
#     ):
#         geo_in = AimsGeometryIn.from_structure(mg2mn4o8)


def test_aims_control_in():
    tmp_path = Path.cwd()
    parameters = {
        "cubes": [
            AimsCube(type="eigenstate 1", points=[10, 10, 10]),
            AimsCube(type="total_density", points=[10, 10, 10]),
        ],
        "xc": "LDA",
        "smearing": ["fermi-dirac", 0.01],
        "vdw_correction_hirshfeld": True,
        "compute_forces": True,
        "relax_geometry": ["trm", "1e-3"],
        "batch_size_limit": 200,
        "species_dir": f"{TEST_DIR.parent}/species_directory/light",
    }

    aims_control = AimsControlIn(parameters.copy())

    for key, val in parameters.items():
        assert aims_control[key] == val

    del aims_control["xc"]
    assert "xc" not in aims_control.parameters
    aims_control.parameters = parameters

    h2o = AimsGeometryIn.from_file(TEST_DIR / "geometry.in.h2o.gz").structure
    aims_control.write_file(h2o, directory=tmp_path, overwrite=True)
    # compare_files(TEST_DIR / "control.in.h2o", f"{tmp_path}/control.in")

    si = AimsGeometryIn.from_file(TEST_DIR / "geometry.in.si.gz").structure
    with pytest.raises(ValueError, match="k-grid must be defined for periodic systems"):
        aims_control.write_file(si, directory=tmp_path, overwrite=True)

    aims_control["k_grid"] = [1, 1, 1]
    with pytest.raises(ValueError, match="control.in file already in "):
        aims_control.write_file(si, directory=tmp_path, overwrite=False)

    aims_control["output"] = "band 0 0 0 0.5 0 0.5 10 G X"
    aims_control["output"] = "band 0 0 0 0.5 0.5 0.5 10 G L"

    aims_control_from_dict = json.loads(json.dumps(aims_control.as_dict(), cls=MontyEncoder), cls=MontyDecoder)
    for key, val in aims_control.parameters.items():
        if key in ["output", "cubes"]:
            np.all(aims_control_from_dict[key] == val)
        assert aims_control_from_dict[key] == val

    aims_control_from_dict.write_file(si, directory=tmp_path, verbose_header=True, overwrite=True)
    compare_files(TEST_DIR / "control.in.si", f"{tmp_path}/control.in")
