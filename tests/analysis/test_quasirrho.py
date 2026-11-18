from __future__ import annotations

from numpy.testing import assert_allclose
from pytest import approx

from pymatgen.analysis.quasirrho import QuasiRRHO, get_avg_mom_inertia
from pymatgen.io.gaussian import GaussianOutput
from pymatgen.io.qchem.outputs import QCOutput
from tests.testing import TEST_FILES_DIR

TEST_DIR = f"{TEST_FILES_DIR}/analysis/quasirrho"


class TestQuasiRRHO:
    """Test class for QuasiRRHO"""

    def setup_method(self):
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
        # tight tolerances: mixing up vibrational temperatures [K] and wavenumbers [cm^-1] in the
        # free-rotor interpolation gave S = 141.297 cal/mol/K and G = -884.776745 Ha
        assert qrrho.entropy_quasiRRHO == approx(correct_stot, rel=0, abs=0.01), "Incorrect total entropy"
        assert qrrho.free_energy_quasiRRHO == approx(correct_g, rel=0, abs=1e-5), "Incorrect Quasi-RRHO free energy"

    def test_qrrho_qchem(self):
        """
        Test a QChem output against independently calculated partition functions.
        """
        # Independent SI calculation from the fixture's 17 coordinates, 45 frequencies and doublet
        # multiplicity, using natural-abundance masses, T=298.15 K, p=101317 Pa and sigma_r=1:
        # sum ideal-gas translational/rotational/electronic and vibrational partition functions;
        # Grimme (2012), doi:10.1002/chem.201200497, interpolates S_HO with S_rotor using
        # mu=h/(8*pi**2*c*freq), mean principal inertia, and w=1/(1+(100/freq)**4), freq in cm^-1.
        # G=E_SCF+H_thermal-T*S, including zero-point energy and pV in H_thermal.
        correct_g = -428.6076869155551
        correct_stot = 103.42339168569977
        # HO total entropy from QChem = 106.521

        qrrho = QuasiRRHO.from_qc_output(self.qout)
        assert qrrho.entropy_quasiRRHO == approx(correct_stot, rel=0, abs=0.01), "Incorrect total entropy"
        assert qrrho.free_energy_quasiRRHO == approx(correct_g, rel=0, abs=1e-5), "Incorrect Quasi-RRHO free energy"

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
