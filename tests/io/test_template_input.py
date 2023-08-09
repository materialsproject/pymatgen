from __future__ import annotations

import os

import pytest

from pymatgen.io.template import TemplateInputGen
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest


class TestTemplateInputGen(PymatgenTest):
    def test_write_inputs(self):
        tis = TemplateInputGen().get_input_set(
            template=f"{TEST_FILES_DIR}/template_input_file.txt",
            variables={"TEMPERATURE": 298},
            filename="hello_world.in",
        )
        tis.write_input(self.tmp_path)
        with open(os.path.join(self.tmp_path, "hello_world.in")) as f:
            assert "298" in f.read()

        with pytest.raises(FileNotFoundError, match="No such file or directory:"):
            tis.write_input(os.path.join(self.tmp_path, "temp"), make_dir=False)

        tis.write_input(os.path.join(self.tmp_path, "temp"), make_dir=True)

        tis = TemplateInputGen().get_input_set(
            template=f"{TEST_FILES_DIR}/template_input_file.txt",
            variables={"TEMPERATURE": 400},
            filename="hello_world.in",
        )

        # test len, iter, getitem
        assert len(tis.inputs) == 1
        assert len(list(tis.inputs)) == 1
        assert isinstance(tis.inputs["hello_world.in"], str)

        with pytest.raises(FileExistsError, match="hello_world.in"):
            tis.write_input(self.tmp_path, overwrite=False)

        tis.write_input(self.tmp_path, overwrite=True)

        with open(os.path.join(self.tmp_path, "hello_world.in")) as f:
            assert "400" in f.read()

        tis.write_input(self.tmp_path, zip_inputs=True)

        assert "InputSet.zip" in list(os.listdir(self.tmp_path))
