#!/usr/bin/env python

import unittest

from pymatgen.core.periodic_table import Specie
from pymatgen.structure_prediction.substitution_probability import SubstitutionProbability


class SubstitutionProbabilityTest(unittest.TestCase):
    def test_substitution_probability(self):
        """
        This test tests specific values in the data folder. If the json is updated,
        these tests will have to be as well
        """
        pc = SubstitutionProbability.from_defaults()
        sp1 = Specie('Ag', 2)
        sp3 = Specie('Cl', -1)
        prob1 = pc.list_probability([sp1, sp3], [sp1, sp3])
        print prob1
        sp1 = Specie('Fe', 4)
        sp3 = Specie('Mn', 3)
        prob1 = pc.probability(sp1, sp3)
        self.assertAlmostEqual(prob1, 1.09474411449e-05, 5, "probability isn't correct")
        sp2 = Specie('Pt', 4)
        sp4 = Specie('Pd', 4)
        prob2 = pc.probability(sp2, sp4)
        self.assertAlmostEqual(prob2, 3.0514797917e-05, 5, "probability isn't correct")
        prob3 = pc.list_probability([sp1, sp2], [sp3, sp4])
        self.assertAlmostEqual(prob3, prob1 * prob2, 5, "probability isn't correct")

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
