# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
Testing for quasirrho.py
"""

import os

from pymatgen.analysis.quasirrho import QuasiRRHO
from pymatgen.io.gaussian import GaussianOutput
from pymatgen.io.qchem.outputs import QCOutput
from pymatgen.util.testing import PymatgenTest
import pytest


class TestQuasiRRHO(PymatgenTest):
    """
    Test class for QuasiRRHO
    """

    def setUp(self):
        test_dir = self.TEST_FILES_DIR
        self.gout = GaussianOutput(os.path.join(test_dir, "molecules", "quasirrho_gaufreq.log"))
        self.linear_gout = GaussianOutput(os.path.join(test_dir, "molecules", "co2.log"))
        self.qout = QCOutput(os.path.join(test_dir, "molecules", "new_qchem_files", "Frequency_no_equal.qout"))

    def test_qrrho_gaussian(self):
        """
        Testing from a Gaussian Output file. Correct values are taken from
        Trevor Seguin's original bash script.
        """
        m = 55
        correct_g_conc = -884.770084
        correct_g = -884.776886
        correct_stot = 141.584080
        qrrho = QuasiRRHO.from_GaussianOutput(self.gout, conc=m)
        assert correct_stot == pytest.approx(qrrho.entropy_quasiRRHO,0.1), \
            'Incorrect total entropy'
        assert correct_g == pytest.approx(qrrho.free_energy_quasiRRHO), \
            'Incorrect Quasi-RRHO free energy'
        assert correct_g_conc == \
               pytest.approx(qrrho.concentration_corrected_g_quasiRRHO), \
            'Incorrect concentration corrected Quasi-RRHO free energy, ' \
            '{} != {}'.format(correct_g_conc, qrrho.concentration_corrected_g_quasiRRHO)
    def test_qrrho_qchem(self):
        """
        Testing from a QChem output file. "Correct" values are from
        the version of QuasiRRHO that worked for Gaussian.
        initial run.
        """
        m = 55
        correct_g_conc = -428.60088006873156
        correct_g = -428.60768184222934
        correct_stot = 103.41012732045324
        # HO total entropy from QChem = 106.521

        qrrho = QuasiRRHO.from_QCOutput(self.qout, conc=m)
        assert correct_stot == pytest.approx(qrrho.entropy_quasiRRHO, 0.1), \
            'Incorrect total entropy'
        assert correct_g == pytest.approx(qrrho.free_energy_quasiRRHO), \
            'Incorrect Quasi-RRHO free energy'
        assert correct_g_conc == \
               pytest.approx(qrrho.concentration_corrected_g_quasiRRHO), \
            'Incorrect concentration corrected Quasi-RRHO free energy'

    def test_rrho_manual(self):
        """
        Test manual input creation. Values from GaussianOutput
        """
        rrho_in = {}
        e = self.gout.final_energy
        mol = self.gout.final_structure
        vib_freqs = [f["frequency"] for f in self.gout.frequencies[-1]]

        m = 55
        correct_g_conc = -884.770084
        correct_g = -884.776886
        correct_stot = 141.584080
        qrrho = QuasiRRHO(mol=mol, energy=e, frequencies=vib_freqs,mult=1,
                          conc=m)
        # self.assertAlmostEqual(correct_stot, qrrho.entropy_quasiRRHO, 0)
        # self.assertAlmostEqual(correct_g, qrrho.free_energy_quasiRRHO, 3)
        # self.assertAlmostEqual(correct_g_conc, qrrho.concentration_corrected_g_quasiRRHO, 3)
        assert correct_stot == pytest.approx(qrrho.entropy_quasiRRHO, 0.1), \
            'Incorrect total entropy'
        assert correct_g == pytest.approx(qrrho.free_energy_quasiRRHO), \
            'Incorrect Quasi-RRHO free energy'
        assert correct_g_conc == \
               pytest.approx(qrrho.concentration_corrected_g_quasiRRHO), \
            'Incorrect concentration corrected Quasi-RRHO free energy'

    def test_rrho_linear(self):
        """
        Testing on a linear CO2 molecule from Gaussian Output file.
        Correct free_energy_ho is checked with Gaussian's internal calculation.
        Correct free_energy_quasirrho is compared internally in the hope of
        preventing future errors.
        .
        """
        correct_g_ho = -187.642070
        correct_g_qrrho = -187.642725
        qrrho = QuasiRRHO.from_GaussianOutput(self.linear_gout)
        assert correct_g_ho == pytest.approx(qrrho.free_energy_ho, 0.0001), \
            'Incorrect harmonic oscillator free energy, {} != {}'.format(
                correct_g_ho, qrrho.free_energy_ho)
        assert correct_g_qrrho == \
               pytest.approx(qrrho.free_energy_quasiRRHO), \
            'Incorrect  Quasi-RRHO free energy'
