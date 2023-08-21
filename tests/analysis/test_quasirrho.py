# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
Testing for quasirrho.py
"""
from __future__ import annotations

import os
import unittest

import pytest

from pymatgen.analysis.quasirrho import QuasiRRHO
from pymatgen.io.gaussian import GaussianOutput
from pymatgen.io.qchem.outputs import QCOutput
from pymatgen.util.testing import TEST_FILES_DIR


class TestQuasiRRHO(unittest.TestCase):
    """
    Test class for QuasiRRHO
    """
    def setUp(self):
        test_dir = TEST_FILES_DIR
        self.gout = GaussianOutput(os.path.join(test_dir, "molecules", "quasirrho_gaufreq.log"))
        self.linear_gout = GaussianOutput(os.path.join(test_dir, "molecules","co2.log.gz"))
        self.qout = QCOutput(os.path.join(test_dir, "molecules", "new_qchem_files", "Frequency_no_equal.qout"))

    def test_qrrho_gaussian(self):
        """
        Testing from a Gaussian Output file. Correct values are taken from
        Trevor Seguin's original bash script.
        """
        correct_g = -884.776886
        correct_stot = 141.584080
        qrrho = QuasiRRHO.from_GaussianOutput(self.gout)
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

        qrrho = QuasiRRHO.from_QCOutput(self.qout)
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
        assert correct_g_ho == pytest.approx(
            qrrho.free_energy_ho, 0.0001
        ), f"Incorrect harmonic oscillator free energy, {correct_g_ho} != {qrrho.free_energy_ho}"
        assert correct_g_qrrho == pytest.approx(qrrho.free_energy_quasiRRHO), "Incorrect  Quasi-RRHO free energy"
