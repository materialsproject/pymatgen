import unittest
import os
from pymatgen.io.gaussian import GaussianOutput
from pymatgen.io.qchem.outputs import QCOutput

from pymatgen.analysis.quasirrho import QuasiRRHO

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files', "molecules")


class TestQuasiRRHO(unittest.TestCase):

    def setUp(self):
        self.gout = GaussianOutput(
            os.path.join(test_dir, "quasirrho_gaufreq.log"))
        self.qout = QCOutput(os.path.join(test_dir,
                                          "new_qchem_files",
                                          "Frequency_no_equal.qout"))

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
        self.assertAlmostEqual(correct_g_conc,
                               qrrho.concentration_corrected_g_quasiRRHO, 3)

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
        self.assertAlmostEqual(correct_g_conc,
                               qrrho.concentration_corrected_g_quasiRRHO, 3)

    def test_rrho_manual(self):
        """
        Test manual input creation. Values from GaussianOutput
        """
        rrho_in = {}
        rrho_in["mult"] = self.gout.spin_multiplicity
        rrho_in["elec_energy"] = self.gout.final_energy
        rrho_in["mol"] = self.gout.final_structure
        vib_freqs = [f['frequency'] for f in self.gout.frequencies[-1]]
        rrho_in['frequencies'] = vib_freqs

        m = 55
        correct_g_conc = -884.770084
        correct_g = -884.776886
        correct_stot = 141.584080
        qrrho = QuasiRRHO(rrho_in, conc=m)
        self.assertAlmostEqual(correct_stot, qrrho.entropy_quasiRRHO, 0)
        self.assertAlmostEqual(correct_g, qrrho.free_energy_quasiRRHO, 3)
        self.assertAlmostEqual(correct_g_conc,
                               qrrho.concentration_corrected_g_quasiRRHO, 3)
