from __future__ import annotations

import pytest
from numpy.testing import assert_array_equal

from pymatgen.io.shengbte import Control
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

f90nml = pytest.importorskip("f90nml")
TEST_DIR = f"{TEST_FILES_DIR}/io/shengbte"


class TestShengBTE(PymatgenTest):
    def setUp(self):
        self.filename = f"{TEST_DIR}/CONTROL-CSLD_Si"
        self.test_dict = {
            "nelements": 1,
            "natoms": 2,
            "ngrid": [25, 25, 25],
            "norientations": 0,
            "lfactor": 0.1,
            "lattvec": [
                [0.0, 2.734363999, 2.734363999],
                [2.734363999, 0.0, 2.734363999],
                [2.734363999, 2.734363999, 0.0],
            ],
            "elements": "Si",
            "types": [1, 1],
            "positions": [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]],
            "scell": [5, 5, 5],
            "t": 500,
            "scalebroad": 0.5,
            "isotopes": False,
            "onlyharmonic": False,
            "nonanalytic": False,
            "nanowires": False,
        }

    def test_from_file(self):
        io = Control.from_file(self.filename)
        assert io["nelements"] == 1
        assert io["natoms"] == 2
        assert tuple(io["ngrid"]) == (25, 25, 25)
        assert io["norientations"] == 0
        assert io["lfactor"] == 0.1
        assert io["lattvec"][0] == [0.0, 2.734363999, 2.734363999]
        assert io["lattvec"][1] == [2.734363999, 0.0, 2.734363999]
        assert io["lattvec"][2] == [2.734363999, 2.734363999, 0.0]
        assert isinstance(io["elements"], list | str)
        if isinstance(io["elements"], list):
            all_strings = all(isinstance(item, str) for item in io["elements"])
            assert all_strings
        assert isinstance(io["types"], list | int)
        if isinstance(io["types"], list):
            all_ints = all(isinstance(item, int) for item in io["types"])
            assert all_ints
        assert_array_equal(io["positions"], [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]])
        assert tuple(io["scell"]) == (5, 5, 5)
        assert io["t"] == 500
        assert io["scalebroad"] == 0.5
        assert not io["isotopes"]
        assert not io["onlyharmonic"]
        assert not io["nonanalytic"]
        assert not io["nanowires"]

        io.to_file(filename=f"{self.tmp_path}/test_control")

        with open(f"{self.tmp_path}/test_control") as file:
            test_str = file.read()
        with open(f"{TEST_DIR}/CONTROL-CSLD_Si") as reference_file:
            reference_string = reference_file.read()
        assert test_str == reference_string

    def test_from_dict(self):
        io = Control.from_dict(self.test_dict)
        io.to_file(filename=f"{self.tmp_path}/test_control")
        with open(f"{self.tmp_path}/test_control") as file:
            test_str = file.read()
        with open(f"{TEST_DIR}/CONTROL-CSLD_Si") as reference_file:
            reference_string = reference_file.read()
        assert test_str == reference_string

    def test_as_from_dict(self):
        # tests as dict and from dict methods
        ctrl_from_file = Control.from_file(self.filename)
        control_from_dict = Control.from_dict(ctrl_from_file.as_dict())
        assert control_from_dict == ctrl_from_file
        assert control_from_dict.to_json() == ctrl_from_file.to_json()
