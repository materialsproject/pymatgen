from __future__ import annotations

import os
from glob import glob

SRC_TXT_PATH = "src/pymatgen.egg-info/SOURCES.txt"


def test_egg_sources_txt_is_complete():
    """Check that all source and data files in src/pymatgen/ are listed in pymatgen.egg-info/SOURCES.txt."""

    with open(SRC_TXT_PATH, encoding="utf-8") as file:
        sources = file.read()

    # Check that all files listed in SOURCES.txt exist
    for src_file in sources.splitlines():
        assert os.path.isfile(src_file), f"{src_file!r} does not exist!"

    # Check that all files in src/pymatgen/ are listed in SOURCES.txt
    for ext in ("py", "json", "json.*", "yaml", "csv"):
        for filepath in glob(f"pymatgen/**/*.{ext}", recursive=True):
            unix_path = filepath.replace("\\", "/")
            if unix_path not in sources:
                raise ValueError(f"{unix_path} not found in {SRC_TXT_PATH}, check package data config")
