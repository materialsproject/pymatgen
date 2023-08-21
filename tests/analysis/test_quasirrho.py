from __future__ import annotations

import unittest

import pytest

from pymatgen.analysis.quasirrho import QuasiRRHO, get_avg_mom_inertia
from pymatgen.io.gaussian import GaussianOutput
from pymatgen.io.qchem.outputs import QCOutput
from pymatgen.util.testing import TEST_FILES_DIR


class TestQuasiRRHO(unittest.TestCase):
    """Test class for QuasiRRHO"""

    def setUp(self):
        test_dir = TEST_FILES_DIR
        self.gout = GaussianOutput(f"{test_dir}/molecules/quasirrho_gaufreq.log")
        self.linear_gout = GaussianOutput(f"{test_dir}/molecules/co2.log.gz")
        self.qout = QCOutput(f"{test_dir}/molecules/new_qchem_files/Frequency_no_equal.qout")

    def test_qrrho_gaussian(self):
        """
        Testing from a Gaussian Output file. Correct values are taken from
        Trevor Seguin's original bash script.
        """
        correct_g = -884.776886
        correct_stot = 141.584080
        qrrho = QuasiRRHO.from_gaussian_output(self.gout)
        assert correct_stot == pytest.approx(qrrho.entropy_quasiRRHO, 0.1), "Incorrect total entropy"
        assert correct_g == pytest.approx(qrrho.free_energy_quasiRRHO), "Incorrect Quasi-RRHO free energy"

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
        assert correct_stot == pytest.approx(qrrho.entropy_quasiRRHO, 0.1), "Incorrect total entropy"
        assert correct_g == pytest.approx(qrrho.free_energy_quasiRRHO), "Incorrect Quasi-RRHO free energy"

    def test_rrho_manual(self):
        """
        Test manual input creation. Values from GaussianOutput
        """
        e = self.gout.final_energy
        mol = self.gout.final_structure
        vib_freqs = [f["frequency"] for f in self.gout.frequencies[-1]]

        correct_g = -884.776886
        correct_stot = 141.584080
        qrrho = QuasiRRHO(mol=mol, energy=e, frequencies=vib_freqs, mult=1)
        assert correct_stot == pytest.approx(qrrho.entropy_quasiRRHO, 0.1), "Incorrect total entropy"
        assert correct_g == pytest.approx(qrrho.free_energy_quasiRRHO), "Incorrect Quasi-RRHO free energy"

    def test_rrho_linear(self):
        """Test on a linear CO2 molecule from Gaussian Output file.
        Correct free_energy_ho is checked with Gaussian's internal calculation.
        Correct free_energy_quasirrho is compared internally in the hope of
        preventing future errors.
        """
        correct_g_ho = -187.642070
        correct_g_qrrho = -187.642725
        qrrho = QuasiRRHO.from_gaussian_output(self.linear_gout)
        assert correct_g_ho == pytest.approx(
            qrrho.free_energy_ho, rel=1e-5
        ), f"Incorrect harmonic oscillator free energy, {correct_g_ho} != {qrrho.free_energy_ho}"
        assert correct_g_qrrho == pytest.approx(qrrho.free_energy_quasiRRHO), "Incorrect  Quasi-RRHO free energy"

    def test_extreme_temperature_and_pressure(self):
        qrrho = QuasiRRHO.from_gaussian_output(self.gout, temp=0.1, press=1e9)
        assert qrrho.temp == 0.1
        assert qrrho.press == 1e9

    def test_get_avg_mom_inertia(self):
        mol = self.gout.final_structure
        avg_mom_inertia, inertia_eigen_vals = get_avg_mom_inertia(mol)
        assert avg_mom_inertia == pytest.approx(0)
        assert inertia_eigen_vals == pytest.approx([0, 0, 0])
