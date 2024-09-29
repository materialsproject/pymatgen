from __future__ import annotations

from shutil import copy as cp
from shutil import rmtree
from typing import TYPE_CHECKING

import pytest

from pymatgen.io.jdftx.jdftxoutput import JDFTxOutput

from .conftest import (
    dump_files_dir,
    ex_files_dir,
    example_sp_outfile_known_simple,
    example_sp_outfile_path,
    jdftxoutfile_matches_known_simple,
)

if TYPE_CHECKING:
    from pathlib import Path


@pytest.fixture
def setup_tmpdir():
    calc_dir = dump_files_dir / "tmp_calc_dir"
    calc_dir.mkdir(exist_ok=True)
    yield calc_dir
    rmtree(calc_dir)


def test_jdftxoutput_from_calc_dir_initialization_exceptions(setup_tmpdir):
    calc_dir = setup_tmpdir
    # Expected error on empty directory
    with pytest.raises(FileNotFoundError):
        JDFTxOutput.from_calc_dir(calc_dir)
    outfile_src = ex_files_dir / "example_ionmin.out"
    outfile_path = calc_dir / "example_ionmin.out"
    cp(outfile_src, outfile_path)
    # Expected error on being fed a file
    with pytest.raises(
        ValueError, match=f"{outfile_path} is not a directory. To initialize from an out file, use from_out_file."
    ):
        JDFTxOutput.from_calc_dir(outfile_path)
    dup_outfile_path = calc_dir / "example_ionmin_dup.out"
    cp(outfile_src, dup_outfile_path)
    # Expected error on multiple out files
    with pytest.raises(
        ValueError,
        match=f"Multiple out files found in {calc_dir}. Please specify the out file by "
        "initializing with the from_out_file method, or by cleaning up the directory.",
    ):
        JDFTxOutput.from_calc_dir(calc_dir)


@pytest.mark.parametrize(
    ("outfile_src", "known"),
    [
        (example_sp_outfile_path, example_sp_outfile_known_simple),
    ],
)
def test_jdftxoutput_outfile_consistency(setup_tmpdir, outfile_src: Path, known: dict):
    calc_dir = setup_tmpdir
    outfile_path = calc_dir / outfile_src.name
    cp(outfile_src, outfile_path)
    for meth, var in [
        (JDFTxOutput.from_calc_dir, calc_dir),
        (JDFTxOutput.from_out_file, outfile_path),
    ]:
        jout = meth(var)
        jdftxoutfile_matches_known_simple(jout.outfile, known)
        del jout
