from __future__ import annotations

import os
from glob import glob

import pytest

src_txt_path = "pymatgen.egg-info/SOURCES.txt"
src_txt_missing = not os.path.isfile(src_txt_path)


@pytest.mark.skipif(src_txt_missing, reason=f"{src_txt_path} not found. Run `pip install .` to create")
def test_egg_sources_txt_is_complete():
    """Check that all source and data files in pymatgen/ are listed in pymatgen.egg-info/SOURCES.txt."""

    with open(src_txt_path) as file:
        sources = file.read()

    # check that all files listed in SOURCES.txt exist
    for src_file in sources.splitlines():
        assert os.path.isfile(src_file), f"{src_file!r} does not exist!"

    # check that all files in pymatgen/ are listed in SOURCES.txt
    for ext in ("py", "json*", "yaml", "csv"):
        for filepath in glob(f"pymatgen/**/*.{ext}", recursive=True):
            unix_path = filepath.replace("\\", "/")
            if unix_path.endswith("dao.py"):
                continue
            if unix_path not in sources:
                raise ValueError(
                    f"{unix_path} not found in {src_txt_path}. check setup.py package_data for "
                    "outdated inclusion rules."
                )
