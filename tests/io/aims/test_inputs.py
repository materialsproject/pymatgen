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
    aims_control["xc"] = "libxc LDA_X+LDA_C_PW"

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
