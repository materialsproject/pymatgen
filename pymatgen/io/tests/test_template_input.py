# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import annotations

import os
import tempfile

import pytest

from pymatgen.io.template import TemplateInputGen
from pymatgen.util.testing import PymatgenTest

test_dir = os.path.join(PymatgenTest.TEST_FILES_DIR)


class TestTemplateInputGen:
    def test_write_inputs(self):
        with tempfile.TemporaryDirectory() as scratch_dir:
            tis = TemplateInputGen().get_input_set(
                template=os.path.join(test_dir, "template_input_file.txt"),
                variables={"TEMPERATURE": 298},
                filename="hello_world.in",
            )
            tis.write_input(scratch_dir)
            with open(os.path.join(scratch_dir, "hello_world.in")) as f:
                assert "298" in f.read()

            with pytest.raises(FileNotFoundError):
                tis.write_input(os.path.join(scratch_dir, "temp"), make_dir=False)

            tis.write_input(os.path.join(scratch_dir, "temp"), make_dir=True)

            tis = TemplateInputGen().get_input_set(
                template=os.path.join(test_dir, "template_input_file.txt"),
                variables={"TEMPERATURE": 400},
                filename="hello_world.in",
            )

            # test len, iter, getitem
            assert len(tis.inputs) == 1
            assert len(list(tis.inputs)) == 1
            assert isinstance(tis.inputs["hello_world.in"], str)

            with pytest.raises(FileExistsError):
                tis.write_input(scratch_dir, overwrite=False)

            tis.write_input(scratch_dir, overwrite=True)

            with open(os.path.join(scratch_dir, "hello_world.in")) as f:
                assert "400" in f.read()

            tis.write_input(scratch_dir, zip_inputs=True)

            assert "InputSet.zip" in list(os.listdir(scratch_dir))
