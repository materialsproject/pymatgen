# coding: utf-8
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
        qrrho = QuasiRRHO(self.gout, conc=m)
        self.assertAlmostEqual(correct_stot, qrrho.entropy_quasiRRHO, 0)
        self.assertAlmostEqual(correct_g, qrrho.free_energy_quasiRRHO, 3)
        self.assertAlmostEqual(correct_g_conc, qrrho.concentration_corrected_g_quasiRRHO, 3)

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

        qrrho = QuasiRRHO(self.qout, conc=m)
        self.assertAlmostEqual(correct_stot, qrrho.entropy_quasiRRHO, 0)
        self.assertAlmostEqual(correct_g, qrrho.free_energy_quasiRRHO, 3)
        self.assertAlmostEqual(correct_g_conc, qrrho.concentration_corrected_g_quasiRRHO, 3)

    def test_rrho_manual(self):
        """
        Test manual input creation. Values from GaussianOutput
        """
        rrho_in = {}
        rrho_in["mult"] = self.gout.spin_multiplicity
        rrho_in["elec_energy"] = self.gout.final_energy
        rrho_in["mol"] = self.gout.final_structure
        vib_freqs = [f["frequency"] for f in self.gout.frequencies[-1]]
        rrho_in["frequencies"] = vib_freqs

        m = 55
        correct_g_conc = -884.770084
        correct_g = -884.776886
        correct_stot = 141.584080
        qrrho = QuasiRRHO(rrho_in, conc=m)
        self.assertAlmostEqual(correct_stot, qrrho.entropy_quasiRRHO, 0)
        self.assertAlmostEqual(correct_g, qrrho.free_energy_quasiRRHO, 3)
        self.assertAlmostEqual(correct_g_conc, qrrho.concentration_corrected_g_quasiRRHO, 3)

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
        qrrho = QuasiRRHO(self.linear_gout)
        self.assertAlmostEqual(correct_g_ho, qrrho.free_energy_ho, 2)
        self.assertAlmostEqual(correct_g_qrrho, qrrho.free_energy_quasiRRHO, 2)
