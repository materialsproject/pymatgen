# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest
import json
import os

from pymatgen.core.periodic_table import Specie
from pymatgen import Composition
from pymatgen.structure_prediction.substitution_probability \
    import SubstitutionProbability, SubstitutionPredictor


def get_table():
    """
    Loads a lightweight lambda table for use in unit tests to reduce
    initialization time, and make unit tests insensitive to changes in the
    default lambda table.
    """
    data_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                            'test_files', "struct_predictor")

    json_file = os.path.join(data_dir, 'test_lambda.json')
    with open(json_file) as f:
        lambda_table = json.load(f)
    return lambda_table


class SubstitutionProbabilityTest(unittest.TestCase):

    def test_full_lambda_table(self):
        """
        This test tests specific values in the data folder. If the
        json is updated, these tests will have to be as well
        """
        sp = SubstitutionProbability(alpha= -5.)
        sp1 = Specie('Fe', 4)
        sp3 = Specie('Mn', 3)
        prob1 = sp.prob(sp1, sp3)
        self.assertAlmostEqual(prob1, 1.69243954552e-05, 5
                               , "probability isn't correct")
        sp2 = Specie('Pt', 4)
        sp4 = Specie('Pd', 4)
        prob2 = sp.prob(sp2, sp4)
        self.assertAlmostEqual(prob2, 4.7174906021e-05, 5
                               , "probability isn't correct")
        corr = sp.pair_corr(Specie("Cu", 2), Specie("Fe", 2))
        self.assertAlmostEqual(corr, 6.82496631637, 5
                               , "probability isn't correct")
        prob3 = sp.cond_prob_list([sp1, sp2], [sp3, sp4])
        self.assertAlmostEqual(prob3, 0.000300298841302, 6
                               , "probability isn't correct")

    def test_mini_lambda_table(self):
        sp = SubstitutionProbability(lambda_table=get_table(), alpha= -5.)
        o2 = Specie('O', -2)
        s2 = Specie('S', -2)
        li1 = Specie('Li', 1)
        na1 = Specie('Na', 1)
        self.assertAlmostEqual(sp.prob(s2, o2), 0.124342317272, 5
                               , "probability isn't correct")
        self.assertAlmostEqual(sp.pair_corr(li1, na1), 1.65425296864, 5
                               , "correlation isn't correct")
        prob = sp.cond_prob_list([o2, li1], [na1, li1])
        self.assertAlmostEqual(prob, 0.00102673915742, 5
                               , "probability isn't correct")


class SubstitutionPredictorTest(unittest.TestCase):

    def test_prediction(self):
        sp = SubstitutionPredictor(threshold = 8e-3)
        result = sp.list_prediction(['Na+', 'Cl-'], to_this_composition = True)[5]
        cprob = sp.p.cond_prob_list(result['substitutions'].keys(),
                                    result['substitutions'].values())
        self.assertAlmostEqual(result['probability'], cprob)
        self.assertEqual(set(result['substitutions'].values()), set(['Na+', 'Cl-']))

        result = sp.list_prediction(['Na+', 'Cl-'], to_this_composition = False)[5]
        cprob = sp.p.cond_prob_list(result['substitutions'].keys(),
                                    result['substitutions'].values())
        self.assertAlmostEqual(result['probability'], cprob)
        self.assertNotEqual(set(result['substitutions'].values()),
                            set(['Na+', 'Cl-']))

        c = Composition({'Ag2+' : 1, 'Cl-' : 2})
        result = sp.composition_prediction(c, to_this_composition = True)[2]
        self.assertEqual(set(result['substitutions'].values()), set(c.elements))
        result = sp.composition_prediction(c, to_this_composition = False)[2]
        self.assertEqual(set(result['substitutions'].keys()), set(c.elements))


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
