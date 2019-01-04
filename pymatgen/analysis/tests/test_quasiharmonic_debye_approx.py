
import unittest
import numpy as np
from pymatgen import Structure
from pymatgen.analysis.eos import EOS

from pymatgen.analysis.quasiharmonic import QuasiharmonicDebyeApprox

__author__ = 'Kiran Mathew'


class TestQuasiharmociDebyeApprox(unittest.TestCase):

    def setUp(self):
        struct = Structure.from_dict({
            'lattice': {'a': 2.5630200477817295,
                        'alpha': 59.999993839702206,
                        'b': 2.563020442699644,
                        'beta': 59.999988742674944,
                        'c': 2.56301993,
                        'gamma': 60.00000504373715,
                        'matrix': [[2.21964022, 0.0, 1.28151046],
                                   [0.73987974, 2.09269747, 1.28151046],
                                   [-0.0, -0.0, 2.56301993]],
                        'volume': 11.905318492097948},
            'sites': [{'abc': [-0.0, -0.0, -0.0],
                       'label': 'B',
                       'species': [{'element': 'B', 'occu': 1}],
                       'xyz': [0.0, 0.0, 0.0]},
                      {'abc': [0.25, 0.25, 0.25],
                       'label': 'N',
                       'species': [{'element': 'N', 'occu': 1}],
                       'xyz': [0.73987999, 0.5231743675, 1.2815102125]}]})

        self.energies = [-15.76315281, -16.11541813, -16.41784171, -16.47471523,
                         -16.63624155, -16.6741551, -16.78661144, -16.88768073,
                         -16.92450672, -17.04863261, -17.06126553, -17.15786866,
                         -17.19784976, -17.25078749, -17.30017149, -17.32578594,
                         -17.3708922, -17.38125127, -17.41231934, -17.41534352,
                         -17.42636644]
        self.volumes = [8.678977833994137, 8.971505437031707, 9.27052889309282, 15.845976281427582,
                        15.417733609491387, 9.576127994353376, 14.997270631725604, 9.888370962140854,
                        14.584523227465766, 14.179424329180256, 10.20732378093211, 13.78189117535765,
                        10.533067462993838, 13.391864274742145, 10.865663655755416, 13.009260480347871,
                        11.205193091129587, 12.634015019827533, 11.551718049704352, 12.26606042141808,
                        11.90531496343142]
        self.eos = "vinet"
        self.T = 300
        self.qhda = QuasiharmonicDebyeApprox(self.energies, self.volumes, struct, t_min=self.T,
                                             t_max=self.T, eos=self.eos)
        self.opt_vol = 11.957803302392925

    def test_bulk_modulus(self):
        eos = EOS(self.eos)
        eos_fit = eos.fit(self.volumes, self.energies)
        bulk_modulus = float(str(eos_fit.b0_GPa).split()[0])
        bulk_modulus_ans =  float(str(self.qhda.bulk_modulus).split()[0])
        np.testing.assert_almost_equal(bulk_modulus, bulk_modulus_ans, 3)

    def test_optimum_volume(self):
        opt_vol = self.qhda.optimum_volumes[0]
        np.testing.assert_almost_equal(opt_vol, self.opt_vol, 3)

    def test_debye_temperature(self):
        theta = self.qhda.debye_temperature(self.opt_vol)
        np.testing.assert_almost_equal(theta, 2559.675227, 3)

    def test_gruneisen_paramter(self):
        gamma = self.qhda.gruneisen_parameter(self.T, self.opt_vol)
        np.testing.assert_almost_equal(gamma, 1.670486, 3)

    def test_thermal_conductivity(self):
        kappa = self.qhda.thermal_conductivity(self.T, self.opt_vol)
        np.testing.assert_almost_equal(kappa, 131.736242, 1)

    def test_vibrational_internal_energy(self):
        u = self.qhda.vibrational_internal_energy(self.T, self.opt_vol)
        np.testing.assert_almost_equal(u, 0.50102, 3)

    def test_vibrational_free_energy(self):
        A = self.qhda.vibrational_free_energy(self.T, self.opt_vol)
        np.testing.assert_almost_equal(A, 0.494687, 3)

class TestAnharmonicQuasiharmociDebyeApprox(unittest.TestCase):

    def setUp(self):
        struct = Structure.from_str("""FCC Al
1.0
2.473329 0.000000 1.427977
0.824443 2.331877 1.427977
0.000000 0.000000 2.855955
Al
1
direct
0.000000 0.000000 0.000000 Al""", fmt='POSCAR')

        self.energies = [-3.69150886, -3.70788383, -3.71997361, -3.72522301,
                         -3.73569569, -3.73649743, -3.74054982]
        self.volumes = [14.824542034870653, 18.118887714656875, 15.373596786943025,
                        17.569833126580278, 15.92265868064787, 17.02077912220064,
                        16.471717630914863]
        self.eos = "vinet"
        self.T = 500
        self.qhda = QuasiharmonicDebyeApprox(self.energies, self.volumes, struct, t_min=self.T,
                                             t_max=self.T, eos=self.eos, anharmonic_contribution=True)
        self.opt_vol = 17.216094889116807

    def test_optimum_volume(self):
        opt_vol = self.qhda.optimum_volumes[0]
        np.testing.assert_almost_equal(opt_vol, self.opt_vol, 3)

    def test_debye_temperature(self):
        theta = self.qhda.debye_temperature(self.opt_vol)
        np.testing.assert_approx_equal(theta, 601.239096, 4 )

    def test_gruneisen_paramter(self):
        gamma = self.qhda.gruneisen_parameter(0, self.qhda.ev_eos_fit.v0)
        np.testing.assert_almost_equal(gamma, 2.188302, 3)

    def test_thermal_conductivity(self):
        kappa = self.qhda.thermal_conductivity(self.T, self.opt_vol)
        np.testing.assert_almost_equal(kappa, 21.810997, 1)

    def test_vibrational_internal_energy(self):
        u = self.qhda.vibrational_internal_energy(self.T, self.opt_vol)
        np.testing.assert_almost_equal(u, 0.13845, 3)

    def test_vibrational_free_energy(self):
        A = self.qhda.vibrational_free_energy(self.T, self.opt_vol)
        np.testing.assert_almost_equal(A, -0.014620, 3)
