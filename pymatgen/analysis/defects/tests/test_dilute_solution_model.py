#!/usr/bin/python

import unittest
import json
import os

from pymatgen.serializers.json_coders import PMGJSONDecoder
from pymatgen.analysis.defects.dilute_solution_model import dilute_solution_model


try:
    import sympy
except ImportError:
    sympy = None

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        'test_files')
with open(os.path.join(test_dir,'mp1048_defect_formation_energies.json')) as fp:
    formation_energy_dict = json.load(fp,cls=PMGJSONDecoder)
with open(os.path.join(test_dir,'mp1048_raw_defect_energies.json')) as fp:
    raw_energy_dict = json.load(fp,cls=PMGJSONDecoder)


@unittest.skipIf(not sympy, "sympy not present.")
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
            self.struct,self.e0,self.vac,self.asites,self.T,generate='energy')
        self.assertIsNotNone(energies)
        self.assertIsNotNone(chem_pot)

    def test_formation_energies_with_chem_pot(self):
        energies, chem_pot = dilute_solution_model(
            self.struct,self.e0,self.vac,self.asites,self.T,
            trial_chem_pot=self.trial_mu,generate='energy')
        self.assertIsNotNone(energies)
        self.assertIsNotNone(chem_pot)

    def test_plot_data_without_chem_pot(self):
        plot_data = dilute_solution_model(
            self.struct,self.e0,self.vac,self.asites,self.T,generate='plot')
        print plot_data.keys()
        self.assertIsNotNone(plot_data)

    def test_plot_data_with_chem_pot(self):
        plot_data = dilute_solution_model(
            self.struct,self.e0,self.vac,self.asites,self.T,
            trial_chem_pot=self.trial_mu,generate='plot')
        self.assertIsNotNone(plot_data)
        for key,value in plot_data.items():
            self.assertIsNotNone(value)
        print plot_data['y']



if __name__ == "__main__":
    unittest.main()
