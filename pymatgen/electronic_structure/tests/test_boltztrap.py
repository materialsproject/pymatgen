# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest
import os
from pymatgen.electronic_structure.boltztrap import BoltztrapAnalyzer
from pymatgen.electronic_structure.core import Spin

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class BoltztrapAnalyzerTest(unittest.TestCase):

    def setUp(self):
        self.bz = BoltztrapAnalyzer.from_files(os.path.join(test_dir, "boltztrap"))

    def test_properties(self):
        self.assertAlmostEqual(self.bz.gap, 1.6644932121620404)
        array = self.bz._cond[300][102]
        self.assertAlmostEqual(array[0][0], 7.5756518e+19)
        self.assertAlmostEqual(array[0][2], -11.14679)
        self.assertAlmostEqual(array[1][0], -88.203286)
        self.assertAlmostEqual(array[2][2], 1.7133249e+19)
        array = self.bz._seebeck[300][22]
        self.assertAlmostEqual(array[0][1], 6.4546074e-22)
        self.assertAlmostEqual(array[1][1], -0.00032073711)
        self.assertAlmostEqual(array[1][2], -2.9868424e-24)
        self.assertAlmostEqual(array[2][2], -0.0003126543)
        array = self.bz._kappa[500][300]
        self.assertAlmostEqual(array[0][1], 0.00014524309)
        self.assertAlmostEqual(array[1][1], 328834400000000.0)
        self.assertAlmostEqual(array[1][2], 3.7758069e-05)
        self.assertAlmostEqual(array[2][2], 193943750000000.0)
        self.assertAlmostEqual(self.bz._hall[400][800][1][0][0], 9.5623749e-28)
        self.assertAlmostEqual(self.bz._hall[400][68][1][2][2], 6.5106975e-10)
        self.assertAlmostEqual(self.bz.doping['p'][3], 1e18)
        self.assertAlmostEqual(self.bz.mu_doping['p'][300][2], 0.1553770018406)
        self.assertAlmostEqual(self.bz.mu_doping['n'][300][-1], 1.6486017632924719)
        self.assertAlmostEqual(self.bz._cond_doping['n'][800][3][1][1], 1.5564085e+16)
        self.assertAlmostEqual(self.bz._seebeck_doping['p'][600][2][0][1], 3.2860613e-23)
        self.assertAlmostEqual(self.bz._carrier_conc[500][67], 38.22832002)
        self.assertAlmostEqual(self.bz.vol, 612.97557323964838)
        self.assertAlmostEqual(self.bz._hall_doping['n'][700][-1][2][2][2], 5.0136483e-26)
        self.assertAlmostEqual(self.bz.dos.efermi, -0.0300005507057)
        self.assertAlmostEqual(self.bz.dos.energies[0], -2.4497049391830448)
        self.assertAlmostEqual(self.bz.dos.energies[345], -0.727088176739)
        self.assertAlmostEqual(self.bz.dos.energies[-1], 3.7569398770153524)
        self.assertAlmostEqual(self.bz.dos.densities[Spin.up][400], 118.70171)
        self.assertAlmostEqual(self.bz.dos.densities[Spin.up][200], 179.58562)
        self.assertAlmostEqual(self.bz.dos.densities[Spin.up][300], 289.43945)

    def test_get_seebeck(self):
        ref = [-768.99078999999995, -724.43919999999991, -686.84682999999973]
        for i in range(0, 3):
            self.assertAlmostEqual(self.bz.get_seebeck()['n'][800][3][i], ref[i])
        self.assertAlmostEqual(self.bz.get_seebeck(output='average')['p'][800][3], 697.608936667)
        self.assertAlmostEqual(self.bz.get_seebeck(output='average', doping_levels=False)[500][520], 1266.7056)
        self.assertAlmostEqual(self.bz.get_seebeck(output='eigs', doping_levels=False)[300][65], -36.2459389333)

    def test_get_conductivity(self):
        ref = [5.9043185000000022, 17.855599000000002, 26.462935000000002]
        for i in range(0, 3):
            self.assertAlmostEqual(self.bz.get_conductivity()['p'][600][2][i], ref[i])
        self.assertAlmostEqual(self.bz.get_conductivity(output='average')['n'][700][1], 1.58736609667)
        self.assertAlmostEqual(self.bz.get_conductivity(output='average', doping_levels=False)[300][457], 2.87163566667)
        self.assertAlmostEqual(self.bz.get_conductivity(output='eigs', doping_levels=False,
                                                        relaxation_time=1e-15)[200][63], 16573.0536667)

    def test_get_power_factor(self):
        ref = [6.2736602345523362, 17.900184232304138, 26.158282220458144]
        for i in range(0, 3):
            self.assertAlmostEqual(self.bz.get_power_factor()['p'][200][2][i], ref[i])
        self.assertAlmostEqual(self.bz.get_power_factor(output='average')['n'][600][4], 411.230962976)
        self.assertAlmostEqual(self.bz.get_power_factor(output='average', doping_levels=False,
                                                        relaxation_time=1e-15)[500][459], 6.59277148467)
        self.assertAlmostEqual(self.bz.get_power_factor(output='eigs', doping_levels=False)[800][61], 2022.67064134)

    def test_get_thermal_conductivity(self):
        ref = [2.7719565628862623e-05, 0.00010048046886793946, 0.00015874549392499391]
        for i in range(0, 3):
            self.assertAlmostEqual(self.bz.get_thermal_conductivity()['p'][300][2][i], ref[i])
        self.assertAlmostEqual(self.bz.get_thermal_conductivity(output='average', relaxation_time=1e-15)['n'][500][0],
                               1.74466575612e-07)
        self.assertAlmostEqual(self.bz.get_thermal_conductivity(output='average', doping_levels=False)[800][874],
                               8.08066254813)
        self.assertAlmostEqual(self.bz.get_thermal_conductivity(output='eigs', doping_levels=False)[200][32],
                               0.0738961845832)

    def test_get_zt(self):
        ref = [0.0002228294548133532, 0.00081441896388844142, 0.00085232847622913053]
        for i in range(0, 3):
            self.assertAlmostEqual(self.bz.get_zt()['n'][400][0][i], ref[i])
        self.assertAlmostEqual(self.bz.get_zt(output='average', kl=0.5)['p'][700][2], 0.0170001879916)
        self.assertAlmostEqual(self.bz.get_zt(output='average', doping_levels=False, relaxation_time=1e-15)[300][240],
                               0.00953842615332)
        self.assertAlmostEqual(self.bz.get_zt(output='eigs', doping_levels=False)[700][65], 0.335990406091)

    def test_get_average_eff_mass(self):
        ref = [0.76045816788363574, 0.96181142990667101, 2.9428428773308628]
        for i in range(0, 3):
            self.assertAlmostEqual(self.bz.get_average_eff_mass()['p'][300][2][i], ref[i])
        ref = [1.1295783824744523, 1.3898454041924351, 5.2459984671977935]
        for i in range(0, 3):
            self.assertAlmostEqual(self.bz.get_average_eff_mass()['n'][600][1][i], ref[i])
        ref = [[9.61811430e-01, -8.25159596e-19, -4.70319444e-19],
               [-8.25159596e-19, 2.94284288e+00, 3.00368916e-18],
               [-4.70319444e-19, 3.00368916e-18, 7.60458168e-01]]
        for i in range(0, 3):
            for j in range(0, 3):
                self.assertAlmostEqual(self.bz.get_average_eff_mass(output='tensor')['p'][300][2][i][j], ref[i][j])
        self.assertAlmostEqual(self.bz.get_average_eff_mass(output='average')['n'][300][2], 1.53769093989)

    def test_get_carrier_concentration(self):
        self.assertAlmostEqual(self.bz.get_carrier_concentration()[300][39], 6.480515652533115e+22)
        self.assertAlmostEqual(self.bz.get_carrier_concentration()[300][693], -6590800965604750.0)

    def test_get_hall_carrier_concentration(self):
        self.assertAlmostEqual(self.bz.get_hall_carrier_concentration()[600][120], 6.773394626767555e+21)
        self.assertAlmostEqual(self.bz.get_hall_carrier_concentration()[500][892], -9.136803845741777e+21)

if __name__ == '__main__':
    unittest.main()
