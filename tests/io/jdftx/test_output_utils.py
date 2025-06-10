"""Tests for the output_utils module."""

from __future__ import annotations

import os

import pytest

from pymatgen.io.jdftx._output_utils import find_first_range_key, get_start_lines
from pymatgen.io.jdftx.joutstructures import _get_joutstructures_start_idx
from pymatgen.io.jdftx.outputs import _find_jdftx_out_file

from .outputs_test_utils import noeigstats_outfile_path
from .shared_test_utils import write_mt_file


def test_get_start_lines():
    """Test the get_start_lines function.

    This function is used to find the start of the JDFTx calculations in the outfile.
    It tests the behavior to make sure the correct errors are raised on empty files and files without
    the pattern to find the start of the JDFTx calculations.
    """
    with pytest.raises(ValueError, match="Outfile parser fed an empty file."):
        get_start_lines([])
    with pytest.raises(ValueError, match="No JDFTx calculations found in file."):
        get_start_lines(["\n", "\n"])
    with open(noeigstats_outfile_path) as f:
        outfile_text = list.copy(list(f))
    assert get_start_lines(outfile_text)[0] == 1
    assert len(get_start_lines(outfile_text)) == 1


def test_find_first_range_key():
    """Test the find_first_range_key function.

    This function is used to find the lines that begin with the key, or that start with a pound + spacing
    and then the key.
    """
    out1 = find_first_range_key("barbie", ["barbie"])
    assert len(out1) == 1
    assert out1[0] == 0
    out1 = find_first_range_key("barbie", ["# barbie"])
    assert len(out1) == 0
    out1 = find_first_range_key("barbie", ["# barbie"], skip_pound=True)
    assert len(out1) == 1
    assert out1[0] == 0
    out1 = find_first_range_key("barbie", ["barbie", "barbie", "barbie"])
    assert len(out1) == 3
    out1 = find_first_range_key("barbie", ["barbie", "ken", "barbie"])
    assert len(out1) == 2
    out1 = find_first_range_key("barbie", ["barbie", "ken", "barbie"], startline=1)
    assert len(out1) == 1


def test_get_joutstructures_start_idx():
    """Test the _get_joutstructures_start_idx function.

    This function is used to find the start of the JOutStructures in the outfile.
    This function just matches line with a matched pattern.
    """
    start_flag = "barbie"
    assert _get_joutstructures_start_idx(["ken", "barbie"], out_slice_start_flag=start_flag) == 1
    assert _get_joutstructures_start_idx(["ken", "(unrelated stuff) barbie"], out_slice_start_flag=start_flag) == 1
    assert _get_joutstructures_start_idx(["barbie", "ken"], out_slice_start_flag=start_flag) == 0
    assert _get_joutstructures_start_idx(["ken", "ken"], out_slice_start_flag=start_flag) is None


def test_find_jdftx_out_file(tmp_path):
    """Test the _find_jdftx_out_file function.

    This function is used to find the JDFTx out file in a directory.
    It tests the behavior to make sure the correct errors are raised on directories without and out file
    and directories with multiple out files. And out file must match "*.out" or "out" exactly.
    """
    with pytest.raises(FileNotFoundError, match="No JDFTx out file found in directory."):
        _find_jdftx_out_file(tmp_path)
    write_mt_file(tmp_path, "test.out")
    assert _find_jdftx_out_file(tmp_path) == tmp_path / "test.out"
    # out file has to match "*.out" or "out" exactly
    write_mt_file(tmp_path, "tinyout")
    assert _find_jdftx_out_file(tmp_path) == tmp_path / "test.out"
    os.remove(_find_jdftx_out_file(tmp_path))
    write_mt_file(tmp_path, "out")
    assert _find_jdftx_out_file(tmp_path) == tmp_path / "out"
    write_mt_file(tmp_path, "tinyout.out")
    with pytest.raises(FileNotFoundError, match="Multiple JDFTx out files found in directory."):
        _find_jdftx_out_file(tmp_path)
