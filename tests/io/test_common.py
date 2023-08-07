from __future__ import annotations

from typing import TYPE_CHECKING

from pymatgen.io.common import VolumetricData
from pymatgen.util.testing import TEST_FILES_DIR

if TYPE_CHECKING:
    from pathlib import Path


def test_cube_io_faithful(tmp_path: Path) -> None:
    in_path = f"{TEST_FILES_DIR}/cube-gh-2817.xyz"

    cube_file = VolumetricData.from_cube(in_path)
    out_path = f"{tmp_path}/cube-gh-2817.xyz"
    cube_file.to_cube(out_path)
    out_cube = VolumetricData.from_cube(out_path)

    # structure should be preserved round-trip to/from cube file
    assert cube_file.structure.volume == out_cube.structure.volume
    assert cube_file.structure == out_cube.structure
