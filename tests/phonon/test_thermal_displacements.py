from __future__ import annotations

import numpy as np
from numpy.testing import assert_allclose
from pytest import approx

from pymatgen.core.structure import Structure
from pymatgen.phonon.thermal_displacements import ThermalDisplacementMatrices
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest


class TestThermalDisplacement(PymatgenTest):
    """Test data from J. George's matlab code https://github.com/JaGeo/MolecularToolbox."""

    def setUp(self) -> None:
        self.thermal = ThermalDisplacementMatrices(
            thermal_displacement_matrix_cart=[
                [5.16e-03, 6.13e-03, 4.15e-03, -1.10e-04, -1.58e-03, -8.10e-04],
                [6.12e-03, 5.48e-03, 3.95e-03, 1.57e-03, -1.30e-04, -7.90e-04],
                [4.20e-03, 4.25e-03, 5.33e-03, 0.00e00, -3.00e-05, -1.39e-03],
                [5.18e-03, 6.33e-03, 3.94e-03, 1.40e-04, -1.52e-03, -6.80e-04],
                [6.16e-03, 5.22e-03, 4.14e-03, 1.53e-03, 1.00e-04, -8.90e-04],
                [5.16e-03, 6.13e-03, 4.15e-03, 1.10e-04, -1.58e-03, 8.10e-04],
                [6.12e-03, 5.48e-03, 3.95e-03, -1.57e-03, -1.30e-04, 7.90e-04],
                [4.20e-03, 4.25e-03, 5.33e-03, -0.00e00, -3.00e-05, 1.39e-03],
                [5.18e-03, 6.33e-03, 3.94e-03, -1.40e-04, -1.52e-03, 6.80e-04],
                [6.16e-03, 5.22e-03, 4.14e-03, -1.53e-03, 1.00e-04, 8.90e-04],
                [4.27e-03, 4.51e-03, 3.70e-03, 3.20e-04, -5.90e-04, -8.30e-04],
                [4.38e-03, 4.44e-03, 3.65e-03, 5.90e-04, -4.20e-04, -8.60e-04],
                [4.46e-03, 4.33e-03, 3.70e-03, 5.70e-04, -3.30e-04, -8.30e-04],
                [4.13e-03, 4.19e-03, 3.81e-03, 3.70e-04, -3.90e-04, -9.00e-04],
                [4.33e-03, 4.44e-03, 3.65e-03, 4.10e-04, -5.90e-04, -8.20e-04],
                [4.27e-03, 4.51e-03, 3.70e-03, -3.20e-04, -5.90e-04, 8.30e-04],
                [4.38e-03, 4.44e-03, 3.65e-03, -5.90e-04, -4.20e-04, 8.60e-04],
                [4.46e-03, 4.33e-03, 3.70e-03, -5.70e-04, -3.30e-04, 8.30e-04],
                [4.13e-03, 4.19e-03, 3.81e-03, -3.70e-04, -3.90e-04, 9.00e-04],
                [4.33e-03, 4.44e-03, 3.65e-03, -4.10e-04, -5.90e-04, 8.20e-04],
                [4.88e-03, 4.97e-03, 3.97e-03, 7.00e-04, -7.00e-04, -1.44e-03],
                [4.88e-03, 4.97e-03, 3.97e-03, -7.00e-04, -7.00e-04, 1.44e-03],
            ],
            structure=Structure.from_file(f"{TEST_FILES_DIR}/thermal_displacement_matrices/POSCAR"),
            temperature=0.0,
        )

        self.thermal_with_cif = ThermalDisplacementMatrices(
            thermal_displacement_matrix_cart=[
                [5.16e-03, 6.13e-03, 4.15e-03, -1.10e-04, -1.58e-03, -8.10e-04],
                [6.12e-03, 5.48e-03, 3.95e-03, 1.57e-03, -1.30e-04, -7.90e-04],
                [4.20e-03, 4.25e-03, 5.33e-03, 0.00e00, -3.00e-05, -1.39e-03],
                [5.18e-03, 6.33e-03, 3.94e-03, 1.40e-04, -1.52e-03, -6.80e-04],
                [6.16e-03, 5.22e-03, 4.14e-03, 1.53e-03, 1.00e-04, -8.90e-04],
                [5.16e-03, 6.13e-03, 4.15e-03, 1.10e-04, -1.58e-03, 8.10e-04],
                [6.12e-03, 5.48e-03, 3.95e-03, -1.57e-03, -1.30e-04, 7.90e-04],
                [4.20e-03, 4.25e-03, 5.33e-03, -0.00e00, -3.00e-05, 1.39e-03],
                [5.18e-03, 6.33e-03, 3.94e-03, -1.40e-04, -1.52e-03, 6.80e-04],
                [6.16e-03, 5.22e-03, 4.14e-03, -1.53e-03, 1.00e-04, 8.90e-04],
                [4.27e-03, 4.51e-03, 3.70e-03, 3.20e-04, -5.90e-04, -8.30e-04],
                [4.38e-03, 4.44e-03, 3.65e-03, 5.90e-04, -4.20e-04, -8.60e-04],
                [4.46e-03, 4.33e-03, 3.70e-03, 5.70e-04, -3.30e-04, -8.30e-04],
                [4.13e-03, 4.19e-03, 3.81e-03, 3.70e-04, -3.90e-04, -9.00e-04],
                [4.33e-03, 4.44e-03, 3.65e-03, 4.10e-04, -5.90e-04, -8.20e-04],
                [4.27e-03, 4.51e-03, 3.70e-03, -3.20e-04, -5.90e-04, 8.30e-04],
                [4.38e-03, 4.44e-03, 3.65e-03, -5.90e-04, -4.20e-04, 8.60e-04],
                [4.46e-03, 4.33e-03, 3.70e-03, -5.70e-04, -3.30e-04, 8.30e-04],
                [4.13e-03, 4.19e-03, 3.81e-03, -3.70e-04, -3.90e-04, 9.00e-04],
                [4.33e-03, 4.44e-03, 3.65e-03, -4.10e-04, -5.90e-04, 8.20e-04],
                [4.88e-03, 4.97e-03, 3.97e-03, 7.00e-04, -7.00e-04, -1.44e-03],
                [4.88e-03, 4.97e-03, 3.97e-03, -7.00e-04, -7.00e-04, 1.44e-03],
            ],
            structure=Structure.from_file(f"{TEST_FILES_DIR}/thermal_displacement_matrices/POSCAR"),
            temperature=0.0,
            thermal_displacement_matrix_cif=[
                [0.00457, 0.00613, 0.00415, -0.00011, -0.00081, -0.00082],
                [0.00601, 0.00548, 0.00395, 0.00157, 0.00058, -0.00049],
                [0.00423, 0.00425, 0.00533, 0.00000, 0.00092, -0.00136],
                [0.00460, 0.00633, 0.00394, 0.00014, -0.00080, -0.00065],
                [0.00613, 0.00522, 0.00414, 0.00153, 0.00083, -0.00061],
                [0.00457, 0.00613, 0.00415, 0.00011, -0.00081, 0.00082],
                [0.00601, 0.00548, 0.00395, -0.00157, 0.00058, 0.00049],
                [0.00423, 0.00425, 0.00533, -0.00000, 0.00092, 0.00136],
                [0.00460, 0.00633, 0.00394, -0.00014, -0.00080, 0.00065],
                [0.00613, 0.00522, 0.00414, -0.00153, 0.00083, 0.00061],
                [0.00405, 0.00451, 0.00370, 0.00032, 0.00008, -0.00076],
                [0.00421, 0.00444, 0.00365, 0.00059, 0.00024, -0.00074],
                [0.00432, 0.00433, 0.00370, 0.00057, 0.00033, -0.00071],
                [0.00399, 0.00419, 0.00381, 0.00037, 0.00030, -0.00082],
                [0.00410, 0.00444, 0.00365, 0.00041, 0.00007, -0.00074],
                [0.00405, 0.00451, 0.00370, -0.00032, 0.00008, 0.00076],
                [0.00421, 0.00444, 0.00365, -0.00059, 0.00024, 0.00074],
                [0.00432, 0.00433, 0.00370, -0.00057, 0.00033, 0.00071],
                [0.00399, 0.00419, 0.00381, -0.00037, 0.00030, 0.00082],
                [0.00410, 0.00444, 0.00365, -0.00041, 0.00007, 0.00074],
                [0.00461, 0.00497, 0.00397, 0.00070, 0.00002, -0.00129],
                [0.00461, 0.00497, 0.00397, -0.00070, 0.00002, 0.00129],
            ],
        )

    def test_ucart(self):
        assert self.thermal.thermal_displacement_matrix_cart[0][0] == approx(0.00516)
        # U11, U22, U33, U23, U13, U12
        assert_allclose(
            self.thermal.thermal_displacement_matrix_cart_matrixform[0],
            [[5.16e-03, -8.10e-04, -1.58e-03], [-8.10e-04, 6.13e-03, -1.10e-04], [-1.58e-03, -1.10e-04, 4.15e-03]],
            5,
        )
        assert_allclose(
            self.thermal_with_cif.thermal_displacement_matrix_cart_matrixform[0],
            [[5.16e-03, -8.10e-04, -1.58e-03], [-8.10e-04, 6.13e-03, -1.10e-04], [-1.58e-03, -1.10e-04, 4.15e-03]],
            5,
        )

    def test_u1_u2_u3(self):
        assert self.thermal.U1U2U3[0].sort() == approx(np.array([2.893872e-03, 5.691239e-03, 6.854889e-03]).sort())

    def test_ustar(self):
        Ustar = self.thermal.Ustar
        assert_allclose(
            Ustar[0],
            ThermalDisplacementMatrices.get_full_matrix(
                [[1.664527e-04, 2.287923e-04, 1.858146e-05, -1.421950e-06, -1.040138e-05, -3.009800e-05]]
            )[0],
            5,
        )

    def test_ucif(self):
        Ucif = self.thermal.Ucif
        assert_allclose(
            Ucif[0],
            ThermalDisplacementMatrices.get_full_matrix(
                [[0.004574, 0.006130, 0.004150, -0.000110, -0.000815, -0.000817]]
            )[0],
            5,
        )

    def test_b(self):
        B = self.thermal.B
        assert_allclose(
            B[0],
            ThermalDisplacementMatrices.get_full_matrix(
                [[0.361112, 0.484005, 0.327672, -0.008685, -0.064335, -0.064479]]
            )[0],
            5,
        )

    def test_beta(self):
        beta = self.thermal.beta
        assert_allclose(
            beta[0],
            ThermalDisplacementMatrices.get_full_matrix(
                [[3.285645e-03, 4.516179e-03, 3.667833e-04, -2.806818e-05, -2.053151e-04, -5.941107e-04]]
            )[0],
            5,
        )
        assert_allclose(
            beta[-1],
            ThermalDisplacementMatrices.get_full_matrix(
                [[3.308590e-03, 3.661568e-03, 3.508740e-04, -1.786229e-04, 4.787484e-06, 9.400372e-04]]
            )[0],
            5,
        )

    def test_write_file(self):
        printed = False
        self.thermal.write_cif(f"{self.tmp_path}/U.cif")
        with open(f"{self.tmp_path}/U.cif") as file:
            file.seek(0)  # set position to start of file
            lines = file.read().splitlines()  # now we won't have those newlines
            if "_atom_site_aniso_U_12" in lines:
                printed = True
        assert printed

    def test_from_ucif(self):
        thermal = ThermalDisplacementMatrices.from_Ucif(
            thermal_displacement_matrix_cif=[
                [0.00457, 0.00613, 0.00415, -0.00011, -0.00081, -0.00082],
                [0.00601, 0.00548, 0.00395, 0.00157, 0.00058, -0.00049],
                [0.00423, 0.00425, 0.00533, 0.00000, 0.00092, -0.00136],
                [0.00460, 0.00633, 0.00394, 0.00014, -0.00080, -0.00065],
                [0.00613, 0.00522, 0.00414, 0.00153, 0.00083, -0.00061],
                [0.00457, 0.00613, 0.00415, 0.00011, -0.00081, 0.00082],
                [0.00601, 0.00548, 0.00395, -0.00157, 0.00058, 0.00049],
                [0.00423, 0.00425, 0.00533, -0.00000, 0.00092, 0.00136],
                [0.00460, 0.00633, 0.00394, -0.00014, -0.00080, 0.00065],
                [0.00613, 0.00522, 0.00414, -0.00153, 0.00083, 0.00061],
                [0.00405, 0.00451, 0.00370, 0.00032, 0.00008, -0.00076],
                [0.00421, 0.00444, 0.00365, 0.00059, 0.00024, -0.00074],
                [0.00432, 0.00433, 0.00370, 0.00057, 0.00033, -0.00071],
                [0.00399, 0.00419, 0.00381, 0.00037, 0.00030, -0.00082],
                [0.00410, 0.00444, 0.00365, 0.00041, 0.00007, -0.00074],
                [0.00405, 0.00451, 0.00370, -0.00032, 0.00008, 0.00076],
                [0.00421, 0.00444, 0.00365, -0.00059, 0.00024, 0.00074],
                [0.00432, 0.00433, 0.00370, -0.00057, 0.00033, 0.00071],
                [0.00399, 0.00419, 0.00381, -0.00037, 0.00030, 0.00082],
                [0.00410, 0.00444, 0.00365, -0.00041, 0.00007, 0.00074],
                [0.00461, 0.00497, 0.00397, 0.00070, 0.00002, -0.00129],
                [0.00461, 0.00497, 0.00397, -0.00070, 0.00002, 0.00129],
            ],
            structure=Structure.from_file(f"{TEST_FILES_DIR}/thermal_displacement_matrices/POSCAR"),
            temperature=0.0,
        )
        assert_allclose(
            thermal.thermal_displacement_matrix_cart,
            [
                [5.16e-03, 6.13e-03, 4.15e-03, -1.10e-04, -1.58e-03, -8.10e-04],
                [6.12e-03, 5.48e-03, 3.95e-03, 1.57e-03, -1.30e-04, -7.90e-04],
                [4.20e-03, 4.25e-03, 5.33e-03, 0.00e00, -3.00e-05, -1.39e-03],
                [5.18e-03, 6.33e-03, 3.94e-03, 1.40e-04, -1.52e-03, -6.80e-04],
                [6.16e-03, 5.22e-03, 4.14e-03, 1.53e-03, 1.00e-04, -8.90e-04],
                [5.16e-03, 6.13e-03, 4.15e-03, 1.10e-04, -1.58e-03, 8.10e-04],
                [6.12e-03, 5.48e-03, 3.95e-03, -1.57e-03, -1.30e-04, 7.90e-04],
                [4.20e-03, 4.25e-03, 5.33e-03, -0.00e00, -3.00e-05, 1.39e-03],
                [5.18e-03, 6.33e-03, 3.94e-03, -1.40e-04, -1.52e-03, 6.80e-04],
                [6.16e-03, 5.22e-03, 4.14e-03, -1.53e-03, 1.00e-04, 8.90e-04],
                [4.27e-03, 4.51e-03, 3.70e-03, 3.20e-04, -5.90e-04, -8.30e-04],
                [4.38e-03, 4.44e-03, 3.65e-03, 5.90e-04, -4.20e-04, -8.60e-04],
                [4.46e-03, 4.33e-03, 3.70e-03, 5.70e-04, -3.30e-04, -8.30e-04],
                [4.13e-03, 4.19e-03, 3.81e-03, 3.70e-04, -3.90e-04, -9.00e-04],
                [4.33e-03, 4.44e-03, 3.65e-03, 4.10e-04, -5.90e-04, -8.20e-04],
                [4.27e-03, 4.51e-03, 3.70e-03, -3.20e-04, -5.90e-04, 8.30e-04],
                [4.38e-03, 4.44e-03, 3.65e-03, -5.90e-04, -4.20e-04, 8.60e-04],
                [4.46e-03, 4.33e-03, 3.70e-03, -5.70e-04, -3.30e-04, 8.30e-04],
                [4.13e-03, 4.19e-03, 3.81e-03, -3.70e-04, -3.90e-04, 9.00e-04],
                [4.33e-03, 4.44e-03, 3.65e-03, -4.10e-04, -5.90e-04, 8.20e-04],
                [4.88e-03, 4.97e-03, 3.97e-03, 7.00e-04, -7.00e-04, -1.44e-03],
                [4.88e-03, 4.97e-03, 3.97e-03, -7.00e-04, -7.00e-04, 1.44e-03],
            ],
            atol=1e-5,
        )

    def test_compute_directionality_quality_criterion(self):
        assert_allclose(
            self.thermal.compute_directionality_quality_criterion(self.thermal)[0]["vector0"],
            [-0.6502072, 0.67306922, 0.35243215],
        )

        assert_allclose(
            self.thermal.compute_directionality_quality_criterion(self.thermal)[0]["vector1"],
            [-0.6502072, 0.67306922, 0.35243215],
        )

        thermal = ThermalDisplacementMatrices(
            thermal_displacement_matrix_cart=[
                [6.12e-03, 5.48e-03, 3.95e-03, 1.57e-03, -1.30e-04, -7.90e-04],
                [5.16e-03, 6.13e-03, 4.15e-03, -1.10e-04, -1.58e-03, -8.10e-04],
                [4.20e-03, 4.25e-03, 5.33e-03, 0.00e00, -3.00e-05, -1.39e-03],
                [5.18e-03, 6.33e-03, 3.94e-03, 1.40e-04, -1.52e-03, -6.80e-04],
                [6.16e-03, 5.22e-03, 4.14e-03, 1.53e-03, 1.00e-04, -8.90e-04],
                [5.16e-03, 6.13e-03, 4.15e-03, 1.10e-04, -1.58e-03, 8.10e-04],
                [6.12e-03, 5.48e-03, 3.95e-03, -1.57e-03, -1.30e-04, 7.90e-04],
                [4.20e-03, 4.25e-03, 5.33e-03, -0.00e00, -3.00e-05, 1.39e-03],
                [5.18e-03, 6.33e-03, 3.94e-03, -1.40e-04, -1.52e-03, 6.80e-04],
                [6.16e-03, 5.22e-03, 4.14e-03, -1.53e-03, 1.00e-04, 8.90e-04],
                [4.27e-03, 4.51e-03, 3.70e-03, 3.20e-04, -5.90e-04, -8.30e-04],
                [4.38e-03, 4.44e-03, 3.65e-03, 5.90e-04, -4.20e-04, -8.60e-04],
                [4.46e-03, 4.33e-03, 3.70e-03, 5.70e-04, -3.30e-04, -8.30e-04],
                [4.13e-03, 4.19e-03, 3.81e-03, 3.70e-04, -3.90e-04, -9.00e-04],
                [4.33e-03, 4.44e-03, 3.65e-03, 4.10e-04, -5.90e-04, -8.20e-04],
                [4.27e-03, 4.51e-03, 3.70e-03, -3.20e-04, -5.90e-04, 8.30e-04],
                [4.38e-03, 4.44e-03, 3.65e-03, -5.90e-04, -4.20e-04, 8.60e-04],
                [4.46e-03, 4.33e-03, 3.70e-03, -5.70e-04, -3.30e-04, 8.30e-04],
                [4.13e-03, 4.19e-03, 3.81e-03, -3.70e-04, -3.90e-04, 9.00e-04],
                [4.33e-03, 4.44e-03, 3.65e-03, -4.10e-04, -5.90e-04, 8.20e-04],
                [4.88e-03, 4.97e-03, 3.97e-03, 7.00e-04, -7.00e-04, -1.44e-03],
                [4.88e-03, 4.97e-03, 3.97e-03, -7.00e-04, -7.00e-04, 1.44e-03],
            ],
            structure=Structure.from_file(f"{TEST_FILES_DIR}/thermal_displacement_matrices/POSCAR"),
            temperature=0.0,
        )
        assert self.thermal.compute_directionality_quality_criterion(self.thermal)[0]["angle"] == approx(0.0)
        assert_allclose(
            self.thermal.compute_directionality_quality_criterion(thermal)[0]["vector0"],
            self.thermal.compute_directionality_quality_criterion(thermal)[1]["vector1"],
        )

    def test_angle(self):
        assert self.thermal._angle_dot([-1, -1, -1], [1, 1, 1]) == approx(180.0)
        assert self.thermal._angle_dot([1, 1, 1], [1, 1, 1]) == approx(0.0)

    def test_ratio_prolate(self):
        assert self.thermal.ratio_prolate[0] == approx(6.854889e-03 / 2.893872e-03)

    def test_to_structure_with_site_properties(self):
        # test creation of structure with site properties
        structure = self.thermal.to_structure_with_site_properties_Ucif()
        # test reading of structure with site properties
        new_thermals = ThermalDisplacementMatrices.from_structure_with_site_properties_Ucif(structure)
        assert_allclose(
            self.thermal.thermal_displacement_matrix_cart,
            new_thermals.thermal_displacement_matrix_cart,
            atol=1e-9,
        )
        assert_allclose(self.thermal.structure.frac_coords, new_thermals.structure.frac_coords)
        assert_allclose(self.thermal.structure.volume, new_thermals.structure.volume)

    def test_visualization_directionality_criterion(self):
        # test file creation for VESTA
        printed = False
        self.thermal.visualize_directionality_quality_criterion(
            filename=f"{self.tmp_path}/U.vesta", other=self.thermal, which_structure=0
        )
        with open(f"{self.tmp_path}/U.vesta") as file:
            file.seek(0)  # set position to start of file
            lines = file.read().splitlines()  # now we won't have those newlines
            if "VECTR" in lines:
                printed = True
        assert printed

    def test_from_cif_p1(self):
        self.thermal.write_cif(f"{self.tmp_path}/U.cif")
        new_thermals = ThermalDisplacementMatrices.from_cif_P1(f"{self.tmp_path}/U.cif")
        assert_allclose(new_thermals[0].thermal_displacement_matrix_cif_matrixform, self.thermal.Ucif)
        assert_allclose(new_thermals[0].structure.frac_coords, self.thermal.structure.frac_coords)
        assert_allclose(new_thermals[0].structure.volume, self.thermal.structure.volume)
