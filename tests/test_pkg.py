"""Check the content of source code, source distribution and binary distribution."""

from __future__ import annotations

import os
import subprocess
import tarfile
from glob import glob
from pathlib import Path
from typing import TYPE_CHECKING

import pytest
from monty.tempfile import ScratchDir

if TYPE_CHECKING:
    from pymatgen.util.typing import PathLike

SRC_DIR = Path(__file__).parent.parent


class TestCheckDistribution:
    def test_source_code(self):
        """Directly check the source code in the working directory."""
        src_txt_path = SRC_DIR / "src/pymatgen.egg-info/SOURCES.txt"

        _check_src_txt_is_complete(SRC_DIR, src_txt_path)

    def test_source_distribution(self):
        """Build the source distribution and verify its contents."""

        with ScratchDir("."):
            # Build the source distribution
            subprocess.run(["python", "-m", "build", "--sdist", SRC_DIR, "--outdir", ".", "-C--quiet"], check=True)

            # Decompress sdist
            sdist_file = next(Path(".").glob("*.tar.gz"))
            sdist_dir = sdist_file.name.removesuffix(".tar.gz")
            with tarfile.open(sdist_file, "r:gz") as tar:
                # TODO: remove attr check after only 3.12+
                if hasattr(tarfile, "data_filter"):
                    tar.extractall("", filter="data")
                else:
                    tar.extractall("")  # noqa: S202

            # Verify source distribution contents
            src_txt_path = f"{sdist_dir}/src/pymatgen.egg-info/SOURCES.txt"
            _check_src_txt_is_complete(src_dir=sdist_dir, src_txt_path=src_txt_path)

    @pytest.mark.skip(reason="WIP")
    def test_binary_distribution(self):
        """Build the binary distribution (wheels) and verify its contents."""


def _check_src_txt_is_complete(src_dir: PathLike, src_txt_path: PathLike) -> None:
    """Check that all source code and data files are listed in given SOURCES.txt.

    Args:
        src_dir (PathLike): Path to the source code directory.
        src_txt_path (PathLike): Path to the "SOURCES.txt" file.
    """
    src_dir = Path(src_dir)
    src_txt_path = Path(src_txt_path)

    assert src_dir.is_dir(), f"{src_dir} is not a directory"
    assert src_txt_path.is_file(), f"{src_txt_path} doesn't exist"

    with open(src_txt_path, encoding="utf-8") as file:
        sources = file.read()

    # Check that all files listed in "SOURCES.txt" exist
    for src_file in sources.splitlines():
        assert os.path.isfile(src_dir / src_file), f"{src_file!r} does not exist!"

    # Check that all files in src/pymatgen/ are listed in SOURCES.txt
    for ext in ("py", "json", "json.*", "yaml", "csv"):
        for filepath in glob(f"{src_dir}/src/pymatgen/**/*.{ext}", recursive=True):
            unix_path = os.path.relpath(filepath.replace("\\", "/"), start=src_dir)

            if unix_path not in sources:
                raise ValueError(f"{unix_path} not found in {src_txt_path}, check package data config")
