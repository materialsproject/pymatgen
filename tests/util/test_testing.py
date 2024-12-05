"""Test testing utils."""

from __future__ import annotations

from pathlib import Path

from pymatgen.core import Element, Structure
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.util.testing import PymatgenTest


class TestPymatgenTest(PymatgenTest):
    def test_tmp_dir(self):
        current_dir = Path.cwd()
        assert current_dir == self.tmp_path

        assert self.tmp_path.is_dir()

        test_file = self.tmp_path / "test_file.txt"
        test_file.write_text("Hello, pytest!")

        assert test_file.exists()
        assert test_file.read_text() == "Hello, pytest!"

    def test_get_structure(self):
        assert isinstance(self.get_structure("LiFePO4"), Structure)

    def test_assert_str_content_equal(self):
        # TODO: see PR 4205, assert_str_content_equal has a minor bug
        # and would not raise AssertError but return a boolean only
        assert self.assert_str_content_equal("hi", " hi")
        assert not self.assert_str_content_equal("hi", "world")

    def test_serialize_with_pickle(self):
        result = self.serialize_with_pickle(Element.from_Z(1))
        assert isinstance(result, list)
        assert result[0] is Element.H

    def test_assert_msonable(self):
        kpoints = Kpoints.gamma_automatic((3, 3, 3), [0, 0, 0])
        kp_str = self.assert_msonable(kpoints)
        assert "Automatic kpoint scheme" in kp_str
