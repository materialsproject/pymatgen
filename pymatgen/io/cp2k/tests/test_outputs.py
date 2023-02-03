# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import annotations

import unittest
from pathlib import Path

from pytest import approx

from pymatgen.io.cp2k.outputs import Cp2kOutput
from pymatgen.util.testing import PymatgenTest


class SetTest(PymatgenTest):
    def setUp(self):
        self.TEST_FILES_DIR = Path.joinpath(self.TEST_FILES_DIR, "cp2k")
        self.out = Cp2kOutput(Path.joinpath(self.TEST_FILES_DIR, "cp2k.out"), auto_load=True)

    def test_files(self):
        """Can find files successfully"""
        self.out.parse_files()
        assert len(self.out.filenames["PDOS"]) == 1
        assert len(self.out.filenames["PDOS"]) == 1
        assert len(self.out.filenames["band_structure"]) == 1
        assert len(self.out.filenames["hyperfine_tensor"]) == 1
        assert len(self.out.filenames["chi_tensor"]) == 1
        assert len(self.out.filenames["g_tensor"]) == 1

    def test_run_info(self):
        """Can extract run info from out file"""
        assert self.out.spin_polarized is True
        assert self.out.completed is True
        assert self.out.num_warnings == [[2]]
        assert self.out.charge == 0
        assert self.out.cp2k_version == "2022.1"
        assert self.out.run_type.upper() == "ENERGY_FORCE"

    def energy_force(self):
        """Can get energy and forces"""
        assert self.out.final_energy == -197.40000341992783
        self.assertArrayAlmostEqual(
            self.out.data["forces"][0], [[-0.00000001, -0.00000001, -0.00000001], [0.00000002, 0.00000002, 0.00000002]]
        )

    def test_band(self):
        """Can parse bandstructure files"""
        assert self.out.band_structure
        assert self.out.band_structure.get_band_gap().get("energy") == 0.27940141999999923

    def test_dos(self):
        """Can parse dos files"""
        assert self.out.data["pdos"]["Si_1"]["s"]["efermi"] == approx(-6.7370756409404455)
        assert self.out.data["tdos"].energies[0] == approx(-6.781065751604123)

    def test_chi(self):
        self.out.parse_chi_tensor()
        assert len(self.out.data["chi_total"]) == 1
        assert self.out.data["PV1"][0] == approx(0.1587)
        assert self.out.data["PV2"][0] == approx(0.4582)
        assert self.out.data["PV3"][0] == approx(0.4582)
        assert self.out.data["ISO"][0] == approx(0.3584)
        assert self.out.data["ANISO"][0] == approx(0.1498)
        self.assertArrayAlmostEqual(
            self.out.data["chi_soft"][0],
            [[5.9508, -1.6579, -1.6579], [-1.6579, 5.9508, -1.6579], [-1.6579, -1.6579, 5.9508]],
        )
        self.assertArrayAlmostEqual(self.out.data["chi_local"][0], [[0, 0, 0], [0, 0, 0], [0, 0, 0]])
        self.assertArrayAlmostEqual(
            self.out.data["chi_total"][0],
            [[5.9508, -1.6579, -1.6579], [-1.6579, 5.9508, -1.6579], [-1.6579, -1.6579, 5.9508]],
        )
        self.assertArrayAlmostEqual(
            self.out.data["chi_total_ppm_cgs"][0],
            [[0.3584, -0.0998, -0.0998], [-0.0998, 0.3584, -0.0998], [-0.0998, -0.0998, 0.3584]],
        )

    def test_gtensor(self):
        self.out.parse_gtensor()
        assert len(self.out.data["gtensor_total"]) == 1
        self.assertArrayAlmostEqual(self.out.data["gmatrix_zke"][0], [[0, 0, 0], [0, 0, 0], [0, 0, 0]])
        self.assertArrayAlmostEqual(self.out.data["gmatrix_so"][0], [[0, 0, 0], [0, 0, 0], [0, 0, 0]])
        self.assertArrayAlmostEqual(self.out.data["gmatrix_soo"][0], [[0, 0, 0], [0, 0, 0], [0, 0, 0]])
        self.assertArrayAlmostEqual(
            self.out.data["gmatrix_total"][0],
            [[2.0023193044, 0.0, 0.0], [0.0, 2.0023193044, 0.0], [0.0, 0.0, 2.0023193044]],
        )
        self.assertArrayAlmostEqual(
            self.out.data["gtensor_total"][0],
            [[2.0023193044, 0.0, 0.0], [0.0, 2.0023193044, 0.0], [0.0, 0.0, 2.0023193044]],
        )
        self.assertArrayAlmostEqual(
            self.out.data["delta_g"][0],
            [
                [0.7158445077, -0.6982592888, 0.0007786468],
                [0.4071365712, 0.4164835801, -0.8128845182],
                [0.5672798720, 0.5822159333, 0.5824243761],
            ],
        )

    def test_hyperfine(self):
        self.out.parse_hyperfine()
        dat = self.out.data["hyperfine_tensor"]
        ref = [
            [
                [0.0000000000, 0.0000001288, 0.0000001288],
                [0.0000001288, -0.0000000000, 0.0000001288],
                [0.0000001288, 0.0000001288, 0.0000000000],
            ],
            [
                [0.0000000000, -0.0000001288, -0.0000001288],
                [-0.0000001288, 0.0000000000, -0.0000001288],
                [-0.0000001288, -0.0000001288, 0.0000000000],
            ],
        ]
        self.assertArrayAlmostEqual(dat[0], ref)


if __name__ == "__main__":
    unittest.main()
