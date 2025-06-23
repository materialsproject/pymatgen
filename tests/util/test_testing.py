"""Test testing utils."""

from __future__ import annotations

import os
import warnings
from pathlib import Path
from unittest.mock import patch

import orjson
import pytest
from monty.json import MontyDecoder

from pymatgen.core import Element, Structure
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.util.misc import is_np_dict_equal
from pymatgen.util.testing import (
    FAKE_POTCAR_DIR,
    STRUCTURES_DIR,
    TEST_FILES_DIR,
    VASP_IN_DIR,
    VASP_OUT_DIR,
    MatSciTest,
    PymatgenTest,
)

warnings.filterwarnings("ignore", message="PymatgenTest is deprecated", category=FutureWarning)


class TestPymatgenTest(PymatgenTest):
    """Test deprecated PymatgenTest."""

    def test_deprecated_warning(self):
        with pytest.warns(FutureWarning, match="PymatgenTest is deprecated"):
            TestPymatgenTest()


def test_paths():
    """Test paths provided in testing util."""
    assert STRUCTURES_DIR.is_dir()
    assert [f for f in os.listdir(STRUCTURES_DIR) if f.endswith(".json")]

    assert TEST_FILES_DIR.is_dir()
    assert os.path.isdir(VASP_IN_DIR)
    assert os.path.isdir(VASP_OUT_DIR)

    assert os.path.isdir(FAKE_POTCAR_DIR)
    assert any(f.startswith("POTCAR") for _root, _dir, files in os.walk(FAKE_POTCAR_DIR) for f in files)


class TestTmpDir(MatSciTest):
    def test_tmp_dir_initialization(self):
        """Test that the working directory is correctly set to a temporary directory."""
        current_dir = Path.cwd()
        assert current_dir == self.tmp_path

        assert self.tmp_path.is_dir()

    def test_tmp_dir_is_clean(self):
        """Test that the temporary directory is empty at the start of the test."""
        assert not any(self.tmp_path.iterdir())

    def test_creating_files_in_tmp_dir(self):
        """Test that files can be created in the temporary directory."""
        test_file = self.tmp_path / "test_file.txt"
        test_file.write_text("Hello, pytest!")

        assert test_file.exists()
        assert test_file.read_text() == "Hello, pytest!"


class TestAssertMSONable(MatSciTest):
    def test_valid_msonable(self):
        """Test a valid MSONable object."""
        kpts_obj = Kpoints.monkhorst_automatic((2, 2, 2), [0, 0, 0])

        result = self.assert_msonable(kpts_obj)
        serialized = orjson.loads(result)

        expected_result = {
            "@module": "pymatgen.io.vasp.inputs",
            "@class": "Kpoints",
            "comment": "Automatic kpoint scheme",
            "nkpoints": 0,
            "generation_style": "Monkhorst",
            "kpoints": [[2, 2, 2]],
            "usershift": [0, 0, 0],
            "kpts_weights": None,
            "coord_type": None,
            "labels": None,
            "tet_number": 0,
            "tet_weight": 0,
            "tet_connections": None,
        }

        assert is_np_dict_equal(serialized, expected_result)

    def test_non_msonable(self):
        non_msonable = dict(hello="world")
        # Test `test_is_subclass` is True
        with pytest.raises(TypeError, match="dict object is not MSONable"):
            self.assert_msonable(non_msonable)

        # Test `test_is_subclass` is False (dict don't have `as_dict` method)
        with pytest.raises(AttributeError, match="'dict' object has no attribute 'as_dict'"):
            self.assert_msonable(non_msonable, test_is_subclass=False)

    def test_cannot_reconstruct(self):
        """Patch the `from_dict` method of `Kpoints` to return a corrupted object"""
        kpts_obj = Kpoints.monkhorst_automatic((2, 2, 2), [0, 0, 0])

        with patch.object(Kpoints, "from_dict", side_effect=lambda d: Kpoints(comment="Corrupted Object")):
            reconstructed_obj = Kpoints.from_dict(kpts_obj.as_dict())
            assert reconstructed_obj.comment == "Corrupted Object"

            with pytest.raises(ValueError, match="Kpoints object could not be reconstructed accurately"):
                self.assert_msonable(kpts_obj)

    def test_not_round_trip(self):
        kpts_obj = Kpoints.monkhorst_automatic((2, 2, 2), [0, 0, 0])

        # Patch the MontyDecoder to return an object of a different class
        class NotAKpoints:
            pass

        with patch.object(MontyDecoder, "process_decoded", side_effect=lambda d: NotAKpoints()) as mock_decoder:
            with pytest.raises(
                TypeError,
                match="The reconstructed NotAKpoints object is not a subclass of Kpoints",
            ):
                self.assert_msonable(kpts_obj)

            mock_decoder.assert_called()


class TestMatSciTest(MatSciTest):
    def test_assert_str_content_equal(self):
        # Cases where strings are equal
        self.assert_str_content_equal("hello world", "hello world")
        self.assert_str_content_equal("  hello   world  ", "hello world")
        self.assert_str_content_equal("\nhello\tworld\n", "hello world")

        # Test whitespace handling
        self.assert_str_content_equal("", "")
        self.assert_str_content_equal("  ", "")
        self.assert_str_content_equal("hello\n", "hello")
        self.assert_str_content_equal("hello\r\n", "hello")
        self.assert_str_content_equal("hello\t", "hello")

        # Cases where strings are not equal
        with pytest.raises(AssertionError, match="Strings are not equal"):
            self.assert_str_content_equal("hello world", "hello_world")

        with pytest.raises(AssertionError, match="Strings are not equal"):
            self.assert_str_content_equal("hello", "hello world")

    def test_get_structure(self):
        # Get structure with name (string)
        structure = self.get_structure("LiFePO4")
        assert isinstance(structure, Structure)

        # Test non-existent structure
        with pytest.raises(FileNotFoundError, match="structure for non-existent doesn't exist"):
            structure = self.get_structure("non-existent")

    def test_serialize_with_pickle(self):
        # Test picklable Element
        result = self.serialize_with_pickle(Element.from_Z(1))
        assert isinstance(result, list)
        assert result[0] is Element.H
