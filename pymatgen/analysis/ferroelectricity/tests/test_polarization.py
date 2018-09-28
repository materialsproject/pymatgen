from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

import unittest

import os
from pymatgen.analysis.ferroelectricity.polarization import *
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.outputs import Outcar
from pymatgen.io.vasp.inputs import Potcar
from pymatgen.util.testing import PymatgenTest

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        'test_files/BTO_221_99_polarization')
bto_folders = ['nonpolar_polarization']
bto_folders += ['interpolation_{}_polarization'.format(str(i)) for i in range(1,9)][::-1]
bto_folders += ['polar_polarization']

structures = [Structure.from_file(test_dir+"/"+folder+"/POSCAR") for folder in bto_folders]

ions = np.matrix([[-44.363484, -44.363484, -44.79259988],
                  [-44.324764, -44.324764, -69.43452043],
                  [-44.286055, -44.286055, -69.8792077 ],
                  [-44.247335, -44.247335, -70.32475473],
                  [-44.208626, -44.208626, -70.77139856],
                  [-44.169906, -44.169906, -71.21889307],
                  [-44.131197, -44.131197, -71.66746921],
                  [-44.092477, -44.092477, -72.1168782 ],
                  [-44.053768, -44.053768, -72.56736141],
                  [-44.015048, -44.015048, -73.01874336]])


class UtilsTest(PymatgenTest):
    def setUp(self):
        self.potcar = Potcar.from_file(test_dir+"/POTCAR")
        self.zval_dict = {'Ba': 10.0, 'Ti': 10.0, 'O': 6.0}
        self.ions = ions
        self.structures = structures

    def test_zval_dict_from_potcar(self):
        zval_dict = zval_dict_from_potcar(self.potcar)
        self.assertDictEqual(self.zval_dict, zval_dict)

    def test_get_total_ionic_dipole(self):
        p_ion = get_total_ionic_dipole(self.structures[-1],self.zval_dict)
        self.assertArrayAlmostEqual(p_ion, self.ions[-1].A1.tolist())


class PolarizationTest(PymatgenTest):
    def setUp(self):
        self.p_ions = ions
        self.p_ions_outcar = np.matrix([[0.0, 0.0, 43.93437],
                                        [0.0, 0.0, 19.81697],
                                        [0.0, 0.0, 19.76076],
                                        [0.0, 0.0, 19.70306],
                                        [0.0, 0.0, 19.64372],
                                        [0.0, 0.0, -5.06619],
                                        [0.0, 0.0, -5.18997],
                                        [0.0, 0.0, -5.31457],
                                        [0.0, 0.0, -5.44026],
                                        [0.0, 0.0, -5.56684]])
        self.p_elecs = np.matrix([[4.03304, -4.03304, -3.60393],
                                  [4.02958, 4.02958, -3.77177],
                                  [4.02611, 4.02611, -3.93397],
                                  [4.02264, 4.02263, -4.08851],
                                  [4.01916, 4.01914, 3.99662],
                                  [4.01567, 4.01565, 3.90327],
                                  [4.01217, 4.01214, 3.81998],
                                  [4.00867, 4.00863, 3.74561],
                                  [4.00517, 4.00512, 3.67949],
                                  [0.00024, 0.00019, 3.61674]])
        self.same_branch = np.matrix([[  9.76948106e-05,  -9.76948108e-05,   4.59556390e-05],
                                      [ -1.36325612e-03,  -1.36325612e-03,   5.99098550e+00],
                                      [ -2.54781559e-03,  -2.54781559e-03,   1.18312234e+01],
                                      [ -3.74896442e-03,  -3.50709575e-03,   1.74695147e+01],
                                      [ -4.67728039e-03,  -4.19508654e-03,   2.28288548e+01],
                                      [ -5.38348125e-03,  -4.90281328e-03,   2.79488973e+01],
                                      [ -5.82178137e-03,  -5.10304293e-03,   3.28220345e+01],
                                      [ -6.28132190e-03,  -5.32598777e-03,   3.74721262e+01],
                                      [ -6.71430111e-03,  -5.52382219e-03,   4.19231297e+01],
                                      [ -5.69679257e-03,  -4.50996078e-03,   4.62887982e+01]])
        self.quanta = np.matrix([[  98.50186747,   98.50186747,   98.50186747],
                                 [  98.09416498,   98.09416498,   98.67403571],
                                 [  97.69065056,   97.69065056,   98.84660662],
                                 [  97.29131054,   97.29131054,   99.01967988],
                                 [  96.89603543,   96.89603543,   99.19315873],
                                 [  96.50481368,   96.50481368,   99.36714337],
                                 [  96.11753848,   96.11753848,   99.54153654],
                                 [  95.7342003 ,   95.7342003 ,   99.71643897],
                                 [  95.35469487,   95.35469487,   99.89175289],
                                 [  94.97901455,   94.97901455,  100.06757957]])
        self.structures = structures
        self.polarization = Polarization(self.p_elecs, self.p_ions, self.structures)
        self.outcars = [Outcar(test_dir+"/"+folder+"/OUTCAR") for folder in bto_folders]
        self.change = np.matrix([[ -5.79448738e-03,  -4.41226597e-03,   4.62887522e+01]])
        self.change_norm = 46.288752795325244
        self.max_jumps = [0.00021336004941047062, 0.00016254800426403291, 0.038269946959965086]
        self.smoothness = [0.00017013512377086267, 0.00013467465540412905, 0.034856268571937743]

    def test_from_outcars_and_structures(self):
        polarization = Polarization.from_outcars_and_structures(self.outcars, self.structures)
        p_elecs, p_ions = polarization.get_pelecs_and_pions(convert_to_muC_per_cm2=False)
        self.assertArrayAlmostEqual(p_elecs[0].A1.tolist(), self.p_elecs[0].A1.tolist())
        self.assertArrayAlmostEqual(p_elecs[-1].A1.tolist(), self.p_elecs[-1].A1.tolist())
        self.assertArrayAlmostEqual(p_ions[0].A1.tolist(), self.p_ions_outcar[0].A1.tolist())
        self.assertArrayAlmostEqual(p_ions[-1].A1.tolist(), self.p_ions_outcar[-1].A1.tolist())
        # Test for calc_ionic_from_zval=True
        polarization = Polarization.from_outcars_and_structures(self.outcars, self.structures,
                                                                calc_ionic_from_zval=True)
        p_elecs, p_ions = polarization.get_pelecs_and_pions(convert_to_muC_per_cm2=False)
        self.assertArrayAlmostEqual(p_elecs[0].A1.tolist(), self.p_elecs[0].A1.tolist())
        self.assertArrayAlmostEqual(p_elecs[-1].A1.tolist(), self.p_elecs[-1].A1.tolist())
        self.assertArrayAlmostEqual(p_ions[0].A1.tolist(), self.p_ions[0].A1.tolist())
        self.assertArrayAlmostEqual(p_ions[-1].A1.tolist(), self.p_ions[-1].A1.tolist())

    def test_get_same_branch_polarization_data(self):
        same_branch = self.polarization.get_same_branch_polarization_data(convert_to_muC_per_cm2=True)
        self.assertArrayAlmostEqual(same_branch[0].A1.tolist(), self.same_branch[0].A1.tolist())
        self.assertArrayAlmostEqual(same_branch[-1].A1.tolist(), self.same_branch[-1].A1.tolist())

    def test_get_lattice_quanta(self):
        quanta = self.polarization.get_lattice_quanta(convert_to_muC_per_cm2=True)
        self.assertArrayAlmostEqual(quanta[0].A1.tolist(), self.quanta[0].A1.tolist())
        self.assertArrayAlmostEqual(quanta[-1].A1.tolist(), self.quanta[-1].A1.tolist())

    def test_get_polarization_change(self):
        change = self.polarization.get_polarization_change()
        self.assertArrayAlmostEqual(change, self.change)

    def test_get_polarization_change_norm(self):
        change_norm = self.polarization.get_polarization_change_norm()
        self.assertAlmostEqual(change_norm, self.change_norm)

    def test_max_spline_jumps(self):
        max_jumps = self.polarization.max_spline_jumps()
        self.assertArrayAlmostEqual(self.max_jumps, max_jumps)

    def test_smoothness(self):
        smoothness = self.polarization.smoothness()
        self.assertArrayAlmostEqual(self.smoothness, smoothness)


class EnergyTrendTest(PymatgenTest):
    def setUp(self):
        self.energies = [-7.97738049,
                         -7.988621176,
                         -7.9793246479999995,
                         -7.987973192,
                         -7.984676138,
                         -7.982700144000001,
                         -7.986539788,
                         -7.980859048000001,
                         -7.978240114,
                         -7.977637218]
        self.energy_trend = EnergyTrend(self.energies)
        self.smoothness = 0.0029874731764648306
        self.max_jump = 0.0058893082867133018

    def test_max_spline_jump(self):
        max_jump = self.energy_trend.max_spline_jump()
        self.assertAlmostEqual(max_jump, self.max_jump)

    def test_smoothness(self):
        smoothness = self.energy_trend.smoothness()
        self.assertAlmostEqual(smoothness, self.smoothness)

    def test_endpoints_minima(self):
        endpoints = self.energy_trend.endpoints_minima(slope_cutoff=1e-2)
        self.assertDictEqual({'polar': True, 'nonpolar': True},
                             endpoints)

if __name__ == '__main__':
    unittest.main()
