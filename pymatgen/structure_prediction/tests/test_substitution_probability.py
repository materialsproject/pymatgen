#!/usr/bin/env python

import unittest

from pymatgen.core.periodic_table import Specie
from pymatgen.structure_prediction.substitution_probability \
    import SubstitutionProbability



class SubstitutionProbabilityTest(unittest.TestCase):
    def test_substitution_probability(self):
        """
        This test tests specific values in the data folder. If the 
        json is updated, these tests will have to be as well
        """
        sp = SubstitutionProbability.from_defaults()
        sp1 = Specie('Fe', 4)
        sp3 = Specie('Mn', 3)
        prob1 = sp.probability(sp1, sp3)
        self.assertAlmostEqual(prob1, 1.09474411449e-05, 5
                               , "probability isn't correct")
        sp2 = Specie('Pt', 4)
        sp4 = Specie('Pd', 4)
        prob2 = sp.probability(sp2, sp4)
        self.assertAlmostEqual(prob2, 3.0514797917e-05, 5
                               , "probability isn't correct")
        corr = sp.pair_correlation(Specie("Cu", 2), Specie("Fe", 2))
        self.assertAlmostEqual(corr, 5.43050207593, 5
                               , "probability isn't correct")                        
        prob3 = sp.conditional_probability_list([sp1, sp2], [sp3, sp4])
        self.assertAlmostEqual(prob3, 4.39320136113e-05, 5
                               , "probability isn't correct")
        
        print sp.pair_correlation(sp2, sp2)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
