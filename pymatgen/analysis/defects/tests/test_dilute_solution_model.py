# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest2 as unittest
import json
import os

from monty.json import MontyDecoder
from pymatgen.analysis.defects.dilute_solution_model import *
import random

try:
    import sympy
except ImportError:
    sympy = None

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        'test_files')
with open(
        os.path.join(test_dir, 'mp1048_defect_formation_energies.json')) as fp:
    formation_energy_dict = json.load(fp, cls=MontyDecoder)
with open(os.path.join(test_dir, 'mp1048_raw_defect_energies.json')) as fp:
    raw_energy_dict = json.load(fp, cls=MontyDecoder)
with open(os.path.join(test_dir, 'mp1487_raw_defect_energies.json')) as fp:
    mp1487_raw_energy_dict = json.load(fp, cls=MontyDecoder)


# TODO (from SP): You MUST redo this entire test. The whole tset is
# monstrously slow. It takes more than 10 mins to get through this test alone.


@unittest.skipIf((not sympy) or random.randint(0, 10) % 10 != 0,
                 "sympy not present or random skip.")
class DiluteSolutionModelTest(unittest.TestCase):
    def setUp(self):
        """
        Setup mandatory inputs for dilute_solution_model
        """
        self.e0 = raw_energy_dict['bulk_energy']
        self.asites = raw_energy_dict['antisites']
        self.vac = raw_energy_dict['vacancies']
        self.struct = raw_energy_dict['structure']
        self.T = 600
        self.trial_mu = formation_energy_dict[str(self.T)]['chemical_potential']

    def test_formation_energies_without_chem_pot(self):
        """
        Should generate formation energies without input chempot
        """
        energies, chem_pot = dilute_solution_model(
            self.struct, self.e0, self.vac, self.asites, self.T,
            generate='energy')
        self.assertIsNotNone(energies)
        self.assertIsNotNone(chem_pot)

    def test_formation_energies_with_chem_pot(self):
        energies, chem_pot = dilute_solution_model(
            self.struct, self.e0, self.vac, self.asites, self.T,
            trial_chem_pot=self.trial_mu, generate='energy')
        self.assertIsNotNone(energies)
        self.assertIsNotNone(chem_pot)

    def test_plot_data_without_chem_pot(self):
        conc_data, en_data, mu_data = dilute_solution_model(
            self.struct, self.e0, self.vac, self.asites, self.T,
            generate='plot')
        self.assertIsNotNone(conc_data)
        self.assertIsNotNone(en_data)
        self.assertIsNotNone(mu_data)
        for key, value in conc_data.items():
            self.assertIsNotNone(value)
        for key, value in mu_data.items():
            self.assertIsNotNone(value)
        for key, value in en_data.items():
            self.assertIsNotNone(value)

    def test_plot_data_with_chem_pot(self):
        conc_data, en_data, mu_data = dilute_solution_model(
            self.struct, self.e0, self.vac, self.asites, self.T,
            trial_chem_pot=self.trial_mu, generate='plot')
        self.assertIsNotNone(conc_data)
        self.assertIsNotNone(en_data)
        self.assertIsNotNone(mu_data)
        for key, value in conc_data.items():
            self.assertIsNotNone(value)
        for key, value in mu_data.items():
            self.assertIsNotNone(value)
        for key, value in en_data.items():
            self.assertIsNotNone(value)
            # print(plot_data['y'])


@unittest.skipIf((not sympy) or random.randint(0, 10) % 10 != 0,
                 "sympy not present or random skip.")
class SoluteSiteFinderTest(unittest.TestCase):
    def setUp(self):
        """
        Setup mandatory inputs for dilute_solution_model
        """
        self.e0 = mp1487_raw_energy_dict['bulk_energy']
        self.asites = mp1487_raw_energy_dict['antisites']
        self.vac = mp1487_raw_energy_dict['vacancies']
        self.solutes = mp1487_raw_energy_dict['solutes']
        self.struct = mp1487_raw_energy_dict['structure']
        self.T = 1000

    def test_plot_data_without_chem_pot(self):
        plot_data = solute_site_preference_finder(
            self.struct, self.e0, self.T, self.vac, self.asites, self.solutes,
            solute_concen=0.01)
        self.assertIsNotNone(plot_data)

    def still_wait_plot_data_with_chem_pot(self):
        plot_data = dilute_solution_model(
            self.struct, self.e0, self.vac, self.asites, self.T,
            trial_chem_pot=self.trial_mu, generate='plot')
        self.assertIsNotNone(plot_data)
        for key, value in plot_data.items():
            self.assertIsNotNone(value)


if __name__ == "__main__":
    unittest.main()
