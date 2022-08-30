# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import os

import pytest
from monty.serialization import MontyDecoder
from monty.tempfile import ScratchDir

from pymatgen.core.structure import Structure
from pymatgen.io.cif import CifParser, CifWriter
from pymatgen.io.core import InputFile, InputSet
from pymatgen.util.testing import PymatgenTest

test_dir = os.path.join(PymatgenTest.TEST_FILES_DIR)


class StructInputFile(InputFile):
    """Test implementation of an InputFile object for CIF."""

    def __init__(self, structure):
        self.structure = structure

    def get_string(self) -> str:
        cw = CifWriter(self.structure)
        return str(cw)

    @classmethod
    def from_string(cls, contents: str):
        cp = CifParser.from_string(contents)
        struct = cp.get_structures()[0]
        return cls(structure=struct)


class FakeClass:
    def __init__(self, a, b):
        self.a = a
        self.b = b

    def write_file(self):
        raise ValueError

    def __str__(self):
        return f"{self.a}\n{self.b}"


class TestInputFile:
    def test_file_io(self):
        with pytest.raises(FileNotFoundError):
            StructInputFile.from_file("fakepath.cif")

        sif = StructInputFile.from_file(os.path.join(test_dir, "Li.cif"))
        assert isinstance(sif.structure, Structure)

        with ScratchDir("."):
            sif.write_file("newLi.cif")
            assert os.path.exists("newLi.cif")

    def test_msonable(self):
        sif = StructInputFile.from_file(os.path.join(test_dir, "Li.cif"))
        sif_dict = sif.as_dict()
        decoder = MontyDecoder()
        temp_sif = decoder.process_decoded(sif_dict)
        assert isinstance(temp_sif, StructInputFile)
        assert sif.structure == temp_sif.structure


class TestInputSet:
    def test_mapping(self):
        sif1 = StructInputFile.from_file(os.path.join(test_dir, "Li.cif"))
        sif2 = StructInputFile.from_file(os.path.join(test_dir, "LiFePO4.cif"))
        sif3 = StructInputFile.from_file(os.path.join(test_dir, "Li2O.cif"))
        inp_set = InputSet(
            {
                "cif1": sif1,
                "cif2": sif2,
                "cif3": sif3,
            },
            kwarg1=1,
            kwarg2="hello",
        )

        assert len(inp_set) == 3
        assert inp_set.kwarg1 == 1
        assert inp_set.kwarg2 == "hello"
        with pytest.raises(AttributeError):
            inp_set.kwarg3
        expected = [("cif1", sif1), ("cif2", sif2), ("cif3", sif3)]

        for (fname, contents), (exp_fname, exp_contents) in zip(inp_set.items(), expected):
            assert fname == exp_fname
            assert contents is exp_contents

        assert inp_set["cif1"] is sif1
        with pytest.raises(KeyError):
            inp_set["kwarg1"]

        sif4 = StructInputFile.from_file(os.path.join(test_dir, "CuCl.cif"))
        inp_set["cif4"] = sif4
        assert inp_set.inputs["cif4"] is sif4
        assert len(inp_set) == 4

        del inp_set["cif2"]
        del inp_set["cif4"]

        assert len(inp_set) == 2
        expected = [("cif1", sif1), ("cif3", sif3)]
        for (fname, contents), (exp_fname, exp_contents) in zip(inp_set.items(), expected):
            assert fname == exp_fname
            assert contents is exp_contents

    def test_equality(self):
        sif1 = StructInputFile.from_file(os.path.join(test_dir, "Li.cif"))
        sif2 = StructInputFile.from_file(os.path.join(test_dir, "LiFePO4.cif"))
        sif3 = StructInputFile.from_file(os.path.join(test_dir, "Li2O.cif"))

        inp_set = InputSet(
            {
                "cif1": sif1,
                "cif2": sif2,
            },
            kwarg1=1,
            kwarg2="hello",
        )

        inp_set2 = InputSet(
            {
                "cif1": sif1,
                "cif2": sif2,
            },
            kwarg1=1,
            kwarg2="hello",
        )

        inp_set3 = InputSet(
            {
                "cif1": sif1,
                "cif2": sif2,
                "cif3": sif3,
            },
            kwarg1=1,
            kwarg2="hello",
        )

        inp_set4 = InputSet(
            {
                "cif1": sif1,
                "cif2": sif2,
            },
            kwarg1=1,
            kwarg2="goodbye",
        )

        inp_set5 = InputSet(
            {
                "cif1": sif1,
                "cif2": sif2,
            },
            kwarg1=1,
            kwarg2="hello",
            kwarg3="goodbye",
        )

        assert inp_set == inp_set2
        assert inp_set != inp_set3
        assert inp_set != inp_set4
        assert inp_set != inp_set5

    def test_msonable(self):
        sif1 = StructInputFile.from_file(os.path.join(test_dir, "Li.cif"))
        sif2 = StructInputFile.from_file(os.path.join(test_dir, "Li2O.cif"))
        inp_set = InputSet(
            {
                "cif1": sif1,
                "cif2": sif2,
            },
            kwarg1=1,
            kwarg2="hello",
        )

        inp_set_dict = inp_set.as_dict()
        decoder = MontyDecoder()
        temp_inp_set = decoder.process_decoded(inp_set_dict)
        assert isinstance(temp_inp_set, InputSet)
        assert temp_inp_set.kwarg1 == 1
        assert temp_inp_set.kwarg2 == "hello"
        assert temp_inp_set._kwargs == inp_set._kwargs
        for (fname, contents), (fname2, contents2) in zip(temp_inp_set.items(), inp_set.items()):
            assert fname == fname2
            assert contents.structure == contents2.structure

    def test_write(self):
        sif1 = StructInputFile.from_file(os.path.join(test_dir, "Li.cif"))
        sif2 = StructInputFile.from_file(os.path.join(test_dir, "Li2O.cif"))
        inp_set = InputSet(
            {
                "cif1": sif1,
                "cif2": sif2,
            },
            kwarg1=1,
            kwarg2="hello",
        )
        with ScratchDir("."):
            inp_set.write_input(directory="input_dir", make_dir=True, overwrite=True, zip_inputs=False)
            assert os.path.exists(os.path.join("input_dir", "cif1"))
            assert os.path.exists(os.path.join("input_dir", "cif2"))
            assert len(os.listdir("input_dir")) == 2
            with pytest.raises(FileExistsError):
                inp_set.write_input(directory="input_dir", make_dir=True, overwrite=False, zip_inputs=False)
            inp_set.write_input(directory="input_dir", make_dir=True, overwrite=True, zip_inputs=True)
            assert len(os.listdir("input_dir")) == 1
            assert os.path.exists(os.path.join("input_dir", f"{type(inp_set).__name__}.zip"))
            with pytest.raises(FileNotFoundError):
                inp_set.write_input(directory="input_dir2", make_dir=False, overwrite=True, zip_inputs=False)

        inp_set = InputSet(
            {"cif1": sif1, "file_from_str": "hello you", "file_from_strcast": FakeClass(a="Aha", b="Beh")}
        )
        with ScratchDir("."):
            inp_set.write_input(directory="input_dir", make_dir=True, overwrite=True, zip_inputs=False)
            assert os.path.exists(os.path.join("input_dir", "cif1"))
            assert os.path.exists(os.path.join("input_dir", "file_from_str"))
            assert os.path.exists(os.path.join("input_dir", "file_from_strcast"))
            assert len(os.listdir("input_dir")) == 3
            cp = CifParser(filename=os.path.join("input_dir", "cif1"))
            assert cp.get_structures()[0] == sif1.structure
            with open(os.path.join("input_dir", "file_from_str")) as f:
                file_from_str = f.read()
                assert file_from_str == "hello you"
            with open(os.path.join("input_dir", "file_from_strcast")) as f:
                file_from_strcast = f.read()
                assert file_from_strcast == "Aha\nBeh"
