# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import os
import pytest
import tempfile

from pymatgen.util.testing import PymatgenTest
from pymatgen.io.core import TemplateInputSet

test_dir = os.path.join(PymatgenTest.TEST_FILES_DIR)


class TestTemplateInputSet:
    def test_write_inputs(self):

        with tempfile.TemporaryDirectory() as scratch_dir:
            tis = TemplateInputSet(
                template=os.path.join(test_dir, "template_input_file.txt"), variables={"TEMPERATURE": 298}
            )
            tis.write_input(scratch_dir, filename="hello_world.in")
            with open(os.path.join(scratch_dir, "hello_world.in"), "r") as f:
                assert "298" in f.read()

            with pytest.raises(FileNotFoundError):
                tis.write_input(os.path.join(scratch_dir, "temp"), make_dir=False, filename="hello_world.in")

            tis.write_input(os.path.join(scratch_dir, "temp"), make_dir=True, filename="hello_world.in")

            tis = TemplateInputSet(
                template=os.path.join(test_dir, "template_input_file.txt"), variables={"TEMPERATURE": 400}
            )

            with pytest.raises(FileExistsError):
                tis.write_input(scratch_dir, overwrite=False, filename="hello_world.in")

            tis.write_input(scratch_dir, overwrite=True, filename="hello_world.in")

            with open(os.path.join(scratch_dir, "hello_world.in"), "r") as f:
                assert "400" in f.read()

            tis.write_input(scratch_dir, zip_inputs=True, filename="hello_world.in")

            assert "TemplateInputSet.zip" in [f for f in os.listdir(scratch_dir)]

    def test_from_directory(self):
        with pytest.raises(NotImplementedError):
            tis = TemplateInputSet.from_directory(test_dir)
