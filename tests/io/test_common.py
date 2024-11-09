from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from pymatgen.io.common import PMGDir, VolumetricData
from pymatgen.util.testing import TEST_FILES_DIR

if TYPE_CHECKING:
    from pathlib import Path


def test_cube_io_faithful(tmp_path: Path) -> None:
    in_path = f"{TEST_FILES_DIR}/io/cube-gh-2817.xyz"

    cube_file = VolumetricData.from_cube(in_path)
    out_path = f"{tmp_path}/cube-gh-2817.xyz"
    cube_file.to_cube(out_path)
    out_cube = VolumetricData.from_cube(out_path)

    # structure should be preserved round-trip to/from cube file
    assert cube_file.structure.volume == out_cube.structure.volume
    assert cube_file.structure == out_cube.structure


class TestPMGDir:
    def test_getitem(self):
        # Some simple testing of loading and reading since all these were tested in other classes.
        d = PMGDir(f"{TEST_FILES_DIR}/io/vasp/fixtures/relaxation")
        assert len(d) == 5
        assert d["OUTCAR"].run_stats["cores"] == 8

        d = PMGDir(f"{TEST_FILES_DIR}/io/vasp/fixtures/scan_relaxation")
        assert len(d) == 2
        assert "vasprun.xml.gz" in d
        assert "OUTCAR" in d
        assert d["vasprun.xml.gz"].incar["METAGGA"] == "R2scan"

        with pytest.raises(ValueError, match="hello not found"):
            d["hello"]

        d = PMGDir(f"{TEST_FILES_DIR}/io/pwscf")
        with pytest.warns(UserWarning, match=r"No parser defined for Si.pwscf.out"):
            assert isinstance(d["Si.pwscf.out"], str)

        # Test NEB directories.
        d = PMGDir(f"{TEST_FILES_DIR}/io/vasp/fixtures/neb_analysis/neb1/neb")

        assert len(d) == 10
        from pymatgen.io.vasp import Poscar

        assert isinstance(d["00/POSCAR"], Poscar)

        outcars = d.get_files_by_name("OUTCAR")
        assert len(outcars) == 5
        assert all("OUTCAR" for k in outcars)

        d.reset()
        for v in d._files.values():
            assert v is None
