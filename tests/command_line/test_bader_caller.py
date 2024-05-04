from __future__ import annotations

import warnings
from shutil import which
from unittest.mock import patch

import numpy as np
import pytest
from monty.shutil import copy_r
from numpy.testing import assert_allclose
from pytest import approx

from pymatgen.command_line.bader_caller import BaderAnalysis, bader_analysis_from_path
from pymatgen.util.testing import TEST_FILES_DIR, VASP_IN_DIR, VASP_OUT_DIR, PymatgenTest

TEST_DIR = f"{TEST_FILES_DIR}/command_line/bader"


@pytest.mark.skipif(not which("bader"), reason="bader executable not present")
class TestBaderAnalysis(PymatgenTest):
    def setUp(self):
        warnings.catch_warnings()

    def test_init(self):
        # test with reference file
        analysis = BaderAnalysis(
            chgcar_filename=f"{VASP_OUT_DIR}/CHGCAR.Fe3O4.gz",
            potcar_filename=f"{VASP_IN_DIR}/POTCAR_Fe3O4.gz",
            chgref_filename=f"{VASP_OUT_DIR}/CHGCAR.Fe3O4_ref.gz",
        )
        assert len(analysis.data) == 14
        assert analysis.data[0]["charge"] == approx(6.6136782, abs=1e-3)
        assert analysis.data[0]["charge"] == analysis.get_charge(0)
        assert analysis.nelectrons == 96
        assert analysis.vacuum_charge == approx(0)
        ans = [
            -1.3863218,
            -1.3812175,
            -1.3812175,
            -1.2615902,
            -1.3812175,
            -1.3862971,
            1.021523,
            1.024357,
            1.021523,
            1.021523,
            1.021523,
            1.021523,
            1.021523,
            1.024357,
        ]
        for idx in range(14):
            assert ans[idx] == approx(analysis.get_charge_transfer(idx), abs=1e-3)
        assert analysis.get_partial_charge(0) == -analysis.get_charge_transfer(0)
        struct = analysis.get_oxidation_state_decorated_structure()
        assert struct[0].specie.oxi_state == approx(1.3863218, abs=1e-3)

        # make sure bader still runs without reference file
        analysis = BaderAnalysis(chgcar_filename=f"{VASP_OUT_DIR}/CHGCAR.Fe3O4.gz")
        assert len(analysis.data) == 14

        # Test Cube file format parsing

        copy_r(TEST_DIR, self.tmp_path)
        analysis = BaderAnalysis(cube_filename=f"{TEST_DIR}/elec.cube.gz")
        assert len(analysis.data) == 9

    def test_from_path(self):
        # we need to create two copies of input files since monty decompressing files
        # deletes the compressed version which can't happen twice in same directory
        copy_r(TEST_DIR, direct_dir := f"{self.tmp_path}/direct")
        copy_r(TEST_DIR, from_path_dir := f"{self.tmp_path}/from_path")
        chgcar_path = f"{direct_dir}/CHGCAR.gz"
        chgref_path = f"{direct_dir}/_CHGCAR_sum.gz"

        analysis = BaderAnalysis(chgcar_filename=chgcar_path, chgref_filename=chgref_path)
        analysis_from_path = BaderAnalysis.from_path(from_path_dir)

        for key in analysis_from_path.summary:
            val, val_from_path = analysis.summary[key], analysis_from_path.summary[key]
            if isinstance(analysis_from_path.summary[key], (bool, str)):
                assert val == val_from_path, f"{key=}"
            elif key == "charge":
                assert_allclose(val, val_from_path, atol=1e-5)

    def test_bader_analysis_from_path(self):
        summary = bader_analysis_from_path(TEST_DIR)
        """
        Reference summary dict (with bader 1.0)
        summary_ref = {
            "magmom": [4.298761, 4.221997, 4.221997, 3.816685, 4.221997, 4.298763, 0.36292, 0.370516, 0.36292,
                0.36292, 0.36292, 0.36292, 0.36292, 0.370516],
            "min_dist": [0.835789, 0.92947, 0.92947, 0.973007, 0.92947, 0.835789, 0.94067, 0.817381, 0.94067,
                0.94067, 0.94067, 0.94067, 0.94067, 0.817381],
            "vacuum_charge": 0.0,
            "vacuum_volume": 0.0,
            "atomic_volume": [9.922887, 8.175158, 8.175158, 9.265802, 8.175158, 9.923233, 12.382546, 12.566972,
                12.382546, 12.382546, 12.382546, 12.382546, 12.382546, 12.566972],
            "charge": [12.248132, 12.26177, 12.26177, 12.600596, 12.26177, 12.248143, 7.267303, 7.256998,
                7.267303, 7.267303, 7.267303, 7.267303, 7.267303, 7.256998],
            "bader_version": 1.0,
            "reference_used": True,
        }
        """
        assert set(summary) == {
            "magmom",
            "min_dist",
            "vacuum_charge",
            "vacuum_volume",
            "atomic_volume",
            "charge",
            "bader_version",
            "reference_used",
        }
        assert summary["reference_used"]
        assert sum(summary["magmom"]) == approx(28, abs=1e-1)

    def test_atom_parsing(self):
        # test with reference file
        analysis = BaderAnalysis(
            chgcar_filename=f"{VASP_OUT_DIR}/CHGCAR.Fe3O4.gz",
            potcar_filename=f"{VASP_IN_DIR}/POTCAR_Fe3O4.gz",
            chgref_filename=f"{VASP_OUT_DIR}/CHGCAR.Fe3O4_ref.gz",
            parse_atomic_densities=True,
        )

        assert len(analysis.atomic_densities) == len(analysis.chgcar.structure)

        assert np.sum(analysis.chgcar.data["total"]) == approx(
            np.sum([np.sum(dct["data"]) for dct in analysis.atomic_densities])
        )

    def test_missing_file_bader_exe_path(self):
        pytest.skip("doesn't reliably raise RuntimeError")
        # mock which("bader") to return None so we always fall back to use bader_exe_path
        with (
            patch("shutil.which", return_value=None),
            pytest.raises(
                RuntimeError, match="BaderAnalysis requires the executable bader be in the PATH or the full path "
            ),
        ):
            BaderAnalysis(chgcar_filename=f"{VASP_OUT_DIR}/CHGCAR.Fe3O4.gz", bader_exe_path="")
