"""Check the content of source code, source distribution and binary distribution."""

from __future__ import annotations

import os
import subprocess
from glob import glob
from pathlib import Path

import pytest

SRC_DIR = Path(__file__).parent.parent


class TestCheckDistribution:
    def test_source_code(self):
        """Directly check the source code in the working directory."""
        src_text_path = SRC_DIR / "src/pymatgen.egg-info/SOURCES.txt"

        check_src_txt_is_complete(SRC_DIR, src_text_path)

    @pytest.mark.skip(reason="WIP")  # TODO: parameterize for bdist/sdist?
    def test_source_distribution(self):
        """Build the source distribution and verify its contents."""
        # Build the source distribution in ScratchDir
        subprocess.run(["python", "-m", "build", "--sdist"], check=True)

        # Decompress

        # Verify the distribution contents
        check_src_txt_is_complete(sdist_path, SRC_TXT_PATH)

    @pytest.mark.skip(reason="WIP")
    def test_binary_distribution(self):
        """Build the binary distribution (wheels) and verify its contents."""


def check_src_txt_is_complete(src_dir: str, src_txt_path: str) -> None:
    """Check that all source code and data files are listed in given SOURCES.txt.

    Args:
        src_dir (str): Path to the source code directory.
        src_txt_path (str): Path to the "SOURCES.txt" file.
    """
    assert os.path.isdir(src_dir)
    assert os.path.isfile(src_txt_path)

    with open(src_txt_path, encoding="utf-8") as file:
        sources = file.read()

    # Check that all files listed in "SOURCES.txt" exist
    for src_file in sources.splitlines():
        assert os.path.isfile(src_file), f"{src_file!r} does not exist!"

    # Check that all files in src/pymatgen/ are listed in SOURCES.txt
    for ext in ("py", "json", "json.*", "yaml", "csv"):
        for filepath in glob(f"{src_dir}/src/pymatgen/**/*.{ext}", recursive=True):
            unix_path = os.path.relpath(filepath.replace("\\", "/"), start=src_dir)

            if unix_path not in sources:
                raise ValueError(f"{unix_path} not found in {src_txt_path}, check package data config")
