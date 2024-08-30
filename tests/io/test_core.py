from __future__ import annotations

import copy
import os
from typing import TYPE_CHECKING

import pytest
from monty.serialization import MontyDecoder

from pymatgen.core.structure import Structure
from pymatgen.io.cif import CifParser, CifWriter
from pymatgen.io.core import InputFile, InputSet
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

if TYPE_CHECKING:
    from typing_extensions import Self


class StructInputFile(InputFile):
    """Test implementation of an InputFile object for CIF."""

    def __init__(self, structure: Structure):
        self.structure = structure

    def get_str(self) -> str:
        cw = CifWriter(self.structure)
        return str(cw)

    @classmethod
    def from_str(cls, contents: str) -> Self:
        struct = Structure.from_str(contents, fmt="cif")
        return cls(structure=struct)


class FakeClass:
    def __init__(self, a, b):
        self.a = a
        self.b = b

    def write_file(self):
        raise ValueError

    def __str__(self):
        return f"{self.a}\n{self.b}"


class TestInputFile(PymatgenTest):
    def test_file_io(self):
        with pytest.raises(FileNotFoundError, match="No such file or directory: 'fakepath.cif'"):
            StructInputFile.from_file("fakepath.cif")

        sif = StructInputFile.from_file(f"{TEST_FILES_DIR}/cif/Li.cif")
        assert isinstance(sif.structure, Structure)

        sif.write_file("newLi.cif")
        assert os.path.isfile("newLi.cif")

    def test_msonable(self):
        sif = StructInputFile.from_file(f"{TEST_FILES_DIR}/cif/Li.cif")
        sif_dict = sif.as_dict()
        decoder = MontyDecoder()
        temp_sif = decoder.process_decoded(sif_dict)
        assert isinstance(temp_sif, StructInputFile)
        assert sif.structure == temp_sif.structure


class TestInputSet(PymatgenTest):
    @classmethod
    def setUpClass(cls):
        cls.sif1 = StructInputFile.from_file(f"{TEST_FILES_DIR}/cif/Li.cif")
        cls.sif2 = StructInputFile.from_file(f"{TEST_FILES_DIR}/cif/LiFePO4.cif")
        cls.sif3 = StructInputFile.from_file(f"{TEST_FILES_DIR}/cif/Li2O.cif")

    def test_mapping(self):
        sif1, sif2, sif3 = self.sif1, self.sif2, self.sif3
        inp_set = InputSet({"cif1": sif1, "cif2": sif2, "cif3": sif3}, kwarg1=1, kwarg2="hello")

        assert len(inp_set) == 3
        assert inp_set.kwarg1 == 1
        assert inp_set.kwarg2 == "hello"
        with pytest.raises(AttributeError, match="has no attribute 'kwarg3'"):
            _ = inp_set.kwarg3
        expected = [("cif1", sif1), ("cif2", sif2), ("cif3", sif3)]

        for (fname, contents), (exp_fname, exp_contents) in zip(inp_set.items(), expected, strict=True):
            assert fname == exp_fname
            assert contents is exp_contents

        assert inp_set["cif1"] is sif1
        with pytest.raises(KeyError, match="'kwarg1'"):
            inp_set["kwarg1"]

        sif4 = StructInputFile.from_file(f"{TEST_FILES_DIR}/cif/CuCl.cif")
        inp_set["cif4"] = sif4
        assert inp_set.inputs["cif4"] is sif4
        assert len(inp_set) == 4

        del inp_set["cif2"]
        del inp_set["cif4"]

        assert len(inp_set) == 2
        expected = [("cif1", sif1), ("cif3", sif3)]
        for (fname, contents), (exp_fname, exp_contents) in zip(inp_set.items(), expected, strict=True):
            assert fname == exp_fname
            assert contents is exp_contents

    def test_equality(self):
        sif1, sif2, sif3 = self.sif1, self.sif2, self.sif3

        inp_set = InputSet({"cif1": sif1, "cif2": sif2}, kwarg1=1, kwarg2="hello")

        inp_set2 = InputSet({"cif1": sif1, "cif2": sif2}, kwarg1=1, kwarg2="hello")

        inp_set3 = InputSet({"cif1": sif1, "cif2": sif2, "cif3": sif3}, kwarg1=1, kwarg2="hello")

        inp_set4 = InputSet({"cif1": sif1, "cif2": sif2}, kwarg1=1, kwarg2="goodbye")

        inp_set5 = InputSet({"cif1": sif1, "cif2": sif2}, kwarg1=1, kwarg2="hello", kwarg3="goodbye")

        assert inp_set.as_dict() == inp_set2.as_dict()
        assert inp_set.as_dict() != inp_set3.as_dict()
        assert inp_set.as_dict() != inp_set4.as_dict()
        assert inp_set.as_dict() != inp_set5.as_dict()

    def test_msonable(self):
        sif1, sif2 = self.sif1, self.sif2
        inp_set = InputSet({"cif1": sif1, "cif2": sif2}, kwarg1=1, kwarg2="hello")

        inp_set_dict = inp_set.as_dict()
        decoder = MontyDecoder()
        temp_inp_set = decoder.process_decoded(inp_set_dict)
        assert isinstance(temp_inp_set, InputSet)
        assert temp_inp_set.kwarg1 == 1
        assert temp_inp_set.kwarg2 == "hello"
        assert temp_inp_set._kwargs == inp_set._kwargs
        for (fname, contents), (fname2, contents2) in zip(temp_inp_set.items(), inp_set.items(), strict=True):
            assert fname == fname2
            assert contents.structure == contents2.structure

    def test_write(self):
        inp_set = InputSet({"cif1": self.sif1, "cif2": self.sif2}, kwarg1=1, kwarg2="hello")
        inp_set.write_input(directory="input_dir", make_dir=True, overwrite=True, zip_inputs=False)
        assert os.path.isfile("input_dir/cif1")
        assert os.path.isfile("input_dir/cif2")
        assert len(os.listdir("input_dir")) == 2
        with pytest.raises(FileExistsError, match="cif1"):
            inp_set.write_input(directory="input_dir", make_dir=True, overwrite=False, zip_inputs=False)
        inp_set.write_input(directory="input_dir", make_dir=True, overwrite=True, zip_inputs=True)
        assert len(os.listdir("input_dir")) == 1
        assert os.path.isfile(f"input_dir/{type(inp_set).__name__}.zip")
        with pytest.raises(FileNotFoundError, match="input_dir2"):
            inp_set.write_input(directory="input_dir2", make_dir=False, overwrite=True, zip_inputs=False)

    def test_write_from_str(self):
        inp_set = InputSet(
            {"cif1": self.sif1, "file_from_str": "hello you", "file_from_str_cast": FakeClass(a="Aha", b="Beh")}
        )
        inp_set.write_input(directory="input_dir", make_dir=True, overwrite=True, zip_inputs=False)
        assert os.path.isfile("input_dir/cif1")
        assert os.path.isfile("input_dir/file_from_str")
        assert os.path.isfile("input_dir/file_from_str_cast")
        assert len(os.listdir("input_dir")) == 3
        parser = CifParser(filename="input_dir/cif1")
        assert parser.parse_structures()[0] == self.sif1.structure
        with open("input_dir/file_from_str") as file:
            file_from_str = file.read()
            assert file_from_str == "hello you"
        with open("input_dir/file_from_str_cast") as file:
            file_from_str_cast = file.read()
            assert file_from_str_cast == "Aha\nBeh"

    def test_copy(self):
        sif1, sif2 = self.sif1, self.sif2

        inp_set = InputSet({"cif1": sif1, "cif2": sif2}, kwarg1=1, kwarg2="hello")
        inp_set2 = copy.copy(inp_set)

        assert inp_set.as_dict() == inp_set2.as_dict()

    def test_deepcopy(self):
        sif1, sif2 = self.sif1, self.sif2

        inp_set = InputSet(
            {"cif1": sif1, "cif2": sif2},
            kwarg1=1,
            kwarg2="hello",
        )
        inp_set2 = copy.deepcopy(inp_set)

        assert inp_set.as_dict() == inp_set2.as_dict()
