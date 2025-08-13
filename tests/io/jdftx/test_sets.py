"""Tests for JDFTx input sets."""

from __future__ import annotations

from pathlib import Path

from pymatgen.io.jdftx.inputs import JDFTXInfile
from pymatgen.io.jdftx.sets import JdftxInputSet
from pymatgen.util.testing import TEST_FILES_DIR

from .inputs_test_utils import assert_idential_jif, ex_infile1_fname

ex_files_dir = Path(TEST_FILES_DIR) / "io" / "jdftx" / "example_files"


def test_jdftxinputset_from_directory():
    input_set = JdftxInputSet.from_file(ex_infile1_fname)
    jdftx_inputfile = JDFTXInfile.from_file(ex_infile1_fname)
    assert_idential_jif(input_set.jdftxinput, jdftx_inputfile)


def test_jdftxinputset_write_file(tmp_path):
    jdftx_inputfile = JDFTXInfile.from_file(ex_infile1_fname)
    input_set = JdftxInputSet(jdftx_inputfile, jdftx_inputfile.structure)
    input_set.write_input(tmp_path, infile="test.in")
    written_file = tmp_path / "test.in"
    read_jdftx_inputfile = JDFTXInfile.from_file(written_file)
    assert written_file.exists()
    assert_idential_jif(read_jdftx_inputfile, jdftx_inputfile)
