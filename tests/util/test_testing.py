from __future__ import annotations

import os

import pytest

from pymatgen.core import Structure
from pymatgen.util.testing import (
    FAKE_POTCAR_DIR,
    STRUCTURES_DIR,
    TEST_FILES_DIR,
    VASP_IN_DIR,
    VASP_OUT_DIR,
    PymatgenTest,
)


def test_paths():
    """Test paths provided in testing util."""
    assert STRUCTURES_DIR.is_dir()
    assert [f for f in os.listdir(STRUCTURES_DIR) if f.endswith(".json")]

    assert TEST_FILES_DIR.is_dir()
    assert os.path.isdir(VASP_IN_DIR)
    assert os.path.isdir(VASP_OUT_DIR)

    assert os.path.isdir(FAKE_POTCAR_DIR)
    assert any(f.startswith("POTCAR") for _root, _dir, files in os.walk(FAKE_POTCAR_DIR) for f in files)


class TestPymatgenTest:
    def test_tmp_dir(self):
        pass

    def test_assert_msonable(self):
        pass

    def test_assert_str_content_equal(self):
        # Cases where strings are equal
        PymatgenTest.assert_str_content_equal("hello world", "hello world")
        PymatgenTest.assert_str_content_equal("  hello   world  ", "hello world")
        PymatgenTest.assert_str_content_equal("\nhello\tworld\n", "hello world")

        # Test whitespace handling
        PymatgenTest.assert_str_content_equal("", "")
        PymatgenTest.assert_str_content_equal("  ", "")
        PymatgenTest.assert_str_content_equal("hello\n", "hello")
        PymatgenTest.assert_str_content_equal("hello\r\n", "hello")
        PymatgenTest.assert_str_content_equal("hello\t", "hello")

        # Cases where strings are not equal
        with pytest.raises(AssertionError, match="Strings are not equal"):
            PymatgenTest.assert_str_content_equal("hello world", "hello_world")

        with pytest.raises(AssertionError, match="Strings are not equal"):
            PymatgenTest.assert_str_content_equal("hello", "hello world")

    def test_get_structure(self):
        # Get structure with name (string)
        structure = PymatgenTest.get_structure("LiFePO4")
        assert isinstance(structure, Structure)

        # Test non-existent structure
        with pytest.raises(FileNotFoundError, match="structure for non-existent doesn't exist"):
            structure = PymatgenTest.get_structure("non-existent")

    def test_serialize_with_pickle(self):
        pass
