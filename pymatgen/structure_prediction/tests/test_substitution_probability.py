#!/usr/bin/env python

import unittest

from pymatgen.core.periodic_table import Specie
from pymatgen.structure_prediction.substitution_probability \
    import SubstitutionProbability, test_table


class SubstitutionProbabilityTest(unittest.TestCase):
    def test_full_lambda_table(self):
        """
        This test tests specific values in the data folder. If the 
        json is updated, these tests will have to be as well
        """
        sp = SubstitutionProbability()
        sp1 = Specie('Fe', 4)
        sp3 = Specie('Mn', 3)
        prob1 = sp.prob(sp1, sp3)
        self.assertAlmostEqual(prob1, 1.09474411449e-05, 5
                               , "probability isn't correct")
        sp2 = Specie('Pt', 4)
        sp4 = Specie('Pd', 4)
        prob2 = sp.prob(sp2, sp4)
        self.assertAlmostEqual(prob2, 3.0514797917e-05, 5
                               , "probability isn't correct")
        corr = sp.pair_corr(Specie("Cu", 2), Specie("Fe", 2))
        self.assertAlmostEqual(corr, 5.43050207593, 5
                               , "probability isn't correct")                        
        prob3 = sp.cond_prob_list([sp1, sp2], [sp3, sp4])
        self.assertAlmostEqual(prob3, 4.39320136113e-05, 5
                               , "probability isn't correct")
        
    def test_mini_lambda_table(self):
        sp = SubstitutionProbability(lambda_table = test_table())
        o2 = Specie('O', -2)
        s2 = Specie('S', -2)
        li1 = Specie('Li', 1)
        na1 = Specie('Na', 1)
        
        self.assertAlmostEqual(sp.prob(s2, o2), 0.0887053176007, 5
                               , "probability isn't correct")
        self.assertAlmostEqual(sp.pair_corr(li1, na1), 1.07252419202, 5
                               , "correlation isn't correct")
        prob = sp.cond_prob_list([o2, li1], [na1, li1])
        self.assertAlmostEqual(prob, 0.0704874560872, 5
                               , "probability isn't correct")


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
