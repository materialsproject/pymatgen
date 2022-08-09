import os

from pymatgen.core.structure import Structure
from pymatgen.phonon.thermal_displacements import ThermalDisplacementMatrices
from pymatgen.util.testing import PymatgenTest


class ThermalDisplacementTest(PymatgenTest):
    """
    Test data from J. George's matlab code https://github.com/JaGeo/MolecularToolbox
    """

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
            structure=Structure.from_file(
                os.path.join(PymatgenTest.TEST_FILES_DIR, "thermal_displacement_matrices", "POSCAR")
            ),
            temperature=0.0,
        )

    def test_Ucart(self):
        self.assertAlmostEqual(self.thermal.thermal_displacement_matrix_cart[0][0], 0.00516)
        # U11, U22, U33, U23, U13, U12
        self.assertArrayAlmostEqual(
            self.thermal.thermal_displacement_matrix_cart_matrixform[0],
            [[5.16e-03, -8.10e-04, -1.58e-03], [-8.10e-04, 6.13e-03, -1.10e-04], [-1.58e-03, -1.10e-04, 4.15e-03]],
            5,
        )

    def test_Ustar(self):
        Ustar = self.thermal.Ustar
        self.assertArrayAlmostEqual(
            Ustar[0],
            ThermalDisplacementMatrices.get_full_matrix(
                [[1.664527e-04, 2.287923e-04, 1.858146e-05, -1.421950e-06, -1.040138e-05, -3.009800e-05]]
            )[0],
            5,
        )

    def test_Ucif(self):
        Ucif = self.thermal.Ucif
        self.assertArrayAlmostEqual(
            Ucif[0],
            ThermalDisplacementMatrices.get_full_matrix(
                [[0.004574, 0.006130, 0.004150, -0.000110, -0.000815, -0.000817]]
            )[0],
            5,
        )

    def test_B(self):
        B = self.thermal.B
        self.assertArrayAlmostEqual(
            B[0],
            ThermalDisplacementMatrices.get_full_matrix(
                [[0.361112, 0.484005, 0.327672, -0.008685, -0.064335, -0.064479]]
            )[0],
            5,
        )

    def test_beta(self):
        beta = self.thermal.beta
        self.assertArrayAlmostEqual(
            beta[0],
            ThermalDisplacementMatrices.get_full_matrix(
                [[3.285645e-03, 4.516179e-03, 3.667833e-04, -2.806818e-05, -2.053151e-04, -5.941107e-04]]
            )[0],
            5,
        )
        self.assertArrayAlmostEqual(
            beta[-1],
            ThermalDisplacementMatrices.get_full_matrix(
                [[3.308590e-03, 3.661568e-03, 3.508740e-04, -1.786229e-04, 4.787484e-06, 9.400372e-04]]
            )[0],
            5,
        )
