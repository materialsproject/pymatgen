from __future__ import annotations

from unittest import TestCase

from numpy.testing import assert_allclose
from pytest import approx

from pymatgen.analysis.quasirrho import QuasiRRHO, get_avg_mom_inertia
from pymatgen.io.gaussian import GaussianOutput
from pymatgen.io.qchem.outputs import QCOutput
from pymatgen.util.testing import TEST_FILES_DIR

TEST_DIR = f"{TEST_FILES_DIR}/analysis/quasirrho"


class TestQuasiRRHO(TestCase):
    """Test class for QuasiRRHO"""

    def setUp(self):
        self.gout = GaussianOutput(f"{TEST_DIR}/quasirrho_gaufreq.log")
        self.linear_gout = GaussianOutput(f"{TEST_DIR}/co2.log.gz")
        self.qout = QCOutput(f"{TEST_DIR}/Frequency_no_equal.qout")

    def test_qrrho_gaussian(self):
        """
        Testing from a Gaussian Output file. Correct values are taken from
        Trevor Seguin's original bash script.
        """
        correct_g = -884.776886
        correct_stot = 141.584080
        qrrho = QuasiRRHO.from_gaussian_output(self.gout)
        assert correct_stot == approx(qrrho.entropy_quasiRRHO, rel=0.1), "Incorrect total entropy"
        assert correct_g == approx(qrrho.free_energy_quasiRRHO), "Incorrect Quasi-RRHO free energy"

    def test_qrrho_qchem(self):
        """
        Testing from a QChem output file. "Correct" values are from
        the version of QuasiRRHO that worked for Gaussian.
        initial run.
        """
        correct_g = -428.60768184222934
        correct_stot = 103.41012732045324
        # HO total entropy from QChem = 106.521

        qrrho = QuasiRRHO.from_qc_output(self.qout)
        assert correct_stot == approx(qrrho.entropy_quasiRRHO, rel=0.1), "Incorrect total entropy"
        assert correct_g == approx(qrrho.free_energy_quasiRRHO), "Incorrect Quasi-RRHO free energy"

    def test_rrho_manual(self):
        """
        Test manual input creation. Values from GaussianOutput
        """
        e_final = self.gout.final_energy
        mol = self.gout.final_structure
        vib_freqs = [freq["frequency"] for freq in self.gout.frequencies[-1]]

        correct_g = -884.776886
        correct_stot = 141.584080
        qrrho = QuasiRRHO(mol=mol, energy=e_final, frequencies=vib_freqs, mult=1)
        assert correct_stot == approx(qrrho.entropy_quasiRRHO, rel=0.1), "Incorrect total entropy"
        assert correct_g == approx(qrrho.free_energy_quasiRRHO), "Incorrect Quasi-RRHO free energy"

    def test_rrho_linear(self):
        """Test on a linear CO2 molecule from Gaussian Output file.
        Correct free_energy_ho is checked with Gaussian's internal calculation.
        Correct free_energy_quasirrho is compared internally in the hope of
        preventing future errors.
        """
        correct_g_ho = -187.642070
        correct_g_qrrho = -187.642725
        qrrho = QuasiRRHO.from_gaussian_output(self.linear_gout)
        assert correct_g_ho == approx(qrrho.free_energy_ho, rel=1e-5), (
            f"Incorrect harmonic oscillator free energy, {correct_g_ho} != {qrrho.free_energy_ho}"
        )
        assert correct_g_qrrho == approx(qrrho.free_energy_quasiRRHO), "Incorrect  Quasi-RRHO free energy"

    def test_extreme_temperature_and_pressure(self):
        qrrho = QuasiRRHO.from_gaussian_output(self.gout, temp=0.1, press=1e9)
        assert qrrho.temp == approx(0.1)
        assert qrrho.press == approx(1e9)

    def test_get_avg_mom_inertia(self):
        mol = self.gout.final_structure
        avg_mom_inertia, inertia_eigen_vals = get_avg_mom_inertia(mol)
        assert avg_mom_inertia == approx(0)
        assert_allclose(inertia_eigen_vals, [0, 0, 0], rtol=1e-6, atol=1e-12)
