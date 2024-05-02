from __future__ import annotations

import os

import pytest

from pymatgen.io.template import TemplateInputGen
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

TEST_DIR = f"{TEST_FILES_DIR}/io"


class TestTemplateInputGen(PymatgenTest):
    def test_write_inputs(self):
        input_set = TemplateInputGen().get_input_set(
            template=f"{TEST_DIR}/template_input_file.txt",
            variables={"TEMPERATURE": 298},
            filename="hello_world.in",
        )
        input_set.write_input(self.tmp_path)
        with open(f"{self.tmp_path}/hello_world.in") as file:
            assert "298" in file.read()

        with pytest.raises(FileNotFoundError, match="No such file or directory:"):
            input_set.write_input(f"{self.tmp_path}/temp", make_dir=False)

        input_set.write_input(f"{self.tmp_path}/temp", make_dir=True)

        input_set = TemplateInputGen().get_input_set(
            template=f"{TEST_DIR}/template_input_file.txt",
            variables={"TEMPERATURE": 400},
            filename="hello_world.in",
        )

        # test len, iter, getitem
        assert len(input_set.inputs) == 1
        assert len(list(input_set.inputs)) == 1
        assert isinstance(input_set.inputs["hello_world.in"], str)

        with pytest.raises(FileExistsError, match="hello_world.in"):
            input_set.write_input(self.tmp_path, overwrite=False)

        input_set.write_input(self.tmp_path, overwrite=True)

        with open(f"{self.tmp_path}/hello_world.in") as file:
            assert "400" in file.read()

        input_set.write_input(self.tmp_path, zip_inputs=True)

        assert "InputSet.zip" in list(os.listdir(self.tmp_path))
