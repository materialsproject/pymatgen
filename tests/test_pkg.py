from __future__ import annotations

import os
from glob import glob

import pytest

SRC_TXT_PATH = "src/pymatgen.egg-info/SOURCES.txt"


@pytest.mark.skipif(
    not os.path.isfile(SRC_TXT_PATH),
    reason=f"{SRC_TXT_PATH=} not found. Run `pip install .` to create",
)
def test_egg_sources_txt_is_complete():
    """Check that all source and data files in pymatgen/ are listed in pymatgen.egg-info/SOURCES.txt."""

    with open(SRC_TXT_PATH, encoding="utf-8") as file:
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
                    f"{unix_path} not found in {SRC_TXT_PATH}. check setup.py package_data for "
                    "outdated inclusion rules."
                )
