from __future__ import annotations

import os
import unittest
import warnings
from shutil import which
from unittest.mock import patch

import numpy as np
import pytest
from monty.shutil import copy_r
from pytest import approx

from pymatgen.command_line.bader_caller import BaderAnalysis, bader_analysis_from_path
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest


@unittest.skipIf(not which("bader"), "bader executable not present")
class TestBaderAnalysis(PymatgenTest):
    _multiprocess_shared_ = True

    def setUp(self):
        warnings.catch_warnings()
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_init(self):
        # test with reference file
        analysis = BaderAnalysis(
            chgcar_filename=f"{TEST_FILES_DIR}/CHGCAR.Fe3O4",
            potcar_filename=f"{TEST_FILES_DIR}/POTCAR.Fe3O4",
            chgref_filename=f"{TEST_FILES_DIR}/CHGCAR.Fe3O4_ref",
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
        for i in range(14):
            assert ans[i] == approx(analysis.get_charge_transfer(i), abs=1e-3)
        assert analysis.get_partial_charge(0) == -analysis.get_charge_transfer(0)
        struct = analysis.get_oxidation_state_decorated_structure()
        assert struct[0].specie.oxi_state == approx(1.3863218, abs=1e-3)

        # make sure bader still runs without reference file
        analysis = BaderAnalysis(chgcar_filename=f"{TEST_FILES_DIR}/CHGCAR.Fe3O4")
        assert len(analysis.data) == 14

        # Test Cube file format parsing
        test_dir = f"{TEST_FILES_DIR}/bader"
        copy_r(test_dir, self.tmp_path)
        analysis = BaderAnalysis(cube_filename=os.path.join(test_dir, "elec.cube.gz"))
        assert len(analysis.data) == 9

    def test_from_path(self):
        test_dir = f"{TEST_FILES_DIR}/bader"
        # we need to create two copies of input files since monty decompressing files
        # deletes the compressed version which can't happen twice in same directory
        copy_r(test_dir, direct_dir := f"{self.tmp_path}/direct")
        copy_r(test_dir, from_path_dir := f"{self.tmp_path}/from_path")
        chgcar_path = f"{direct_dir}/CHGCAR.gz"
        chgref_path = f"{direct_dir}/_CHGCAR_sum.gz"
        analysis = BaderAnalysis.from_path(from_path_dir)
        analysis0 = BaderAnalysis(chgcar_filename=chgcar_path, chgref_filename=chgref_path)
        charge = np.array(analysis.summary["charge"])
        charge0 = np.array(analysis0.summary["charge"])
        assert np.allclose(charge, charge0)

    def test_automatic_runner(self):
        pytest.skip("raises RuntimeError: bader exits with return code 24")
        summary = bader_analysis_from_path(f"{TEST_FILES_DIR}/bader")
        """
        Reference summary dict (with bader 1.0)
        summary_ref = {
            'magmom': [4.298761, 4.221997, 4.221997, 3.816685, 4.221997, 4.298763, 0.36292,
                       0.370516, 0.36292, 0.36292, 0.36292, 0.36292, 0.36292, 0.370516],
            'min_dist': [0.835789, 0.92947, 0.92947, 0.973007, 0.92947, 0.835789, 0.94067,
                         0.817381, 0.94067, 0.94067, 0.94067, 0.94067, 0.94067, 0.817381],
            'vacuum_charge': 0.0,
            'vacuum_volume': 0.0,
            'atomic_volume': [9.922887, 8.175158, 8.175158, 9.265802, 8.175158, 9.923233, 12.382546,
                              12.566972, 12.382546, 12.382546, 12.382546, 12.382546, 12.382546, 12.566972],
            'charge': [12.248132, 12.26177, 12.26177, 12.600596, 12.26177, 12.248143, 7.267303,
                       7.256998, 7.267303, 7.267303, 7.267303, 7.267303, 7.267303, 7.256998],
            'bader_version': 1.0,
            'reference_used': True
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
            chgcar_filename=f"{TEST_FILES_DIR}/CHGCAR.Fe3O4",
            potcar_filename=f"{TEST_FILES_DIR}/POTCAR.Fe3O4",
            chgref_filename=f"{TEST_FILES_DIR}/CHGCAR.Fe3O4_ref",
            parse_atomic_densities=True,
        )

        assert len(analysis.atomic_densities) == len(analysis.chgcar.structure)

        assert np.sum(analysis.chgcar.data["total"]) == approx(
            np.sum([np.sum(d["data"]) for d in analysis.atomic_densities])
        )

    def test_missing_file_bader_exe_path(self):
        pytest.skip("doesn't reliably raise RuntimeError")
        # mock which("bader") to return None so we always fall back to use bader_exe_path
        with patch("shutil.which", return_value=None), pytest.raises(
            RuntimeError, match="BaderAnalysis requires the executable bader be in the PATH or the full path "
        ):
            BaderAnalysis(chgcar_filename=f"{TEST_FILES_DIR}/CHGCAR.Fe3O4", bader_exe_path="")
