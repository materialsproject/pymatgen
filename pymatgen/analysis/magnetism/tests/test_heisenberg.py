# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import warnings

import os
import unittest
import pandas as pd

from monty.serialization import loadfn

import pymatgen.analysis.magnetism.heisenberg as heisenberg
from pymatgen import Structure

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        'test_files', 'magnetic_orderings')

formulas = ['CoS2', 'MnNi2Sn', 'MnSi', 'LiMnPO4', 'Ta2CoO6', 'Ni(SbO3)2']

class HeisenbergMapperTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.df = pd.read_json(os.path.join(test_dir, 'mag_orderings_test_cases.json'))

        df = cls.df

        cls.CoS2 = df[df['formula_pretty']=='CoS2']
        cls.MnNi2Sn = df[df['formula_pretty']=='MnNi2Sn']
        cls.MnSi = df[df['formula_pretty']=='MnSi']
        cls.LiMnPO4 = df[df['formula_pretty']=='LiMnPO4']
        cls.Ta2CoO6 = df[df['formula_pretty']=='Ta2CoO6']
        cls.Ni_SbO3_2 = df[df['formula_pretty']=='Ni(SbO3)2']

        # cls.compounds = [cls.CoS2, cls.MnNi2Sn, cls.MnSi,
        # cls.LiMnPO4, cls.Ta2CoO6, cls.Ni_SbO3_2]
        cls.compounds = [cls.CoS2]

        cls.hms = []
        for c in cls.compounds:
            ordered_structures = list(c['structure'])
            ordered_structures = [Structure.from_dict(d) for d in ordered_structures]
            epa = list(c['energy_per_atom'])
            energies = [e*len(s) for (e, s) in zip(epa, ordered_structures)]

            hm = heisenberg.HeisenbergMapper(ordered_structures, energies)
            cls.hms.append(hm)

    def setUp(self):
        pass

    def tearDown(self):
        warnings.simplefilter('default')

    def test_graphs(self):
        for hm in self.hms:
            sgraphs = hm.sgraphs
        pass

    def test_sites(self):
        pass

    def test_nn_interactions(self):
        pass

    def test_exchange_matrix(self):
        pass

    def test_exchange_params(self):
        pass

    def test_mean_field(self):
        pass

if __name__ == '__main__':
    unittest.main()