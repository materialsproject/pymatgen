#!/usr/bin/env python

import unittest

from pymatgen.core.periodic_table import Specie
from pymatgen.core.structure import Composition
from pymatgen.structure_prediction.substitution_probability import test_table
from pymatgen.structure_prediction.substitutor import Substitutor


class SubstitutorTest(unittest.TestCase):

    def setUp(self):
        self.s = Substitutor(threshold=1e-3, lambda_table=test_table(), alpha= -5.)

    def test_substitutor(self):
        s_list = [Specie('O', -2), Specie('Li', 1)]
        subs = self.s.pred_from_list(s_list)
        self.assertEqual(len(subs), 4
                         , 'incorrect number of substitutions')
        c = Composition({'O2-': 1, 'Li1+': 2})
        subs = self.s.pred_from_comp(c)
        self.assertEqual(len(subs), 4
                         , 'incorrect number of substitutions')

    def test_to_dict(self):
        Substitutor.from_dict(self.s.to_dict)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
