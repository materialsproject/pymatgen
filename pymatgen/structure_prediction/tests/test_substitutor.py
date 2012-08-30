#!/usr/bin/env python

import unittest

from pymatgen.core.periodic_table import Specie
from pymatgen.core.structure import Composition
from pymatgen.structure_prediction.substitutor import Substitutor

class SubstitutorTest(unittest.TestCase):
    def test_substitutor(self):
        s = Substitutor()
        s_list = [Specie('Ag', 2), Specie('Cl', -1)]
        subs =  s.pred_from_list(s_list)
        self.assertEqual(len(subs), 197
                         , 'incorrect number of substitutions')
        c = Composition({'Ag2+': 1, 'Cl1-':2})
        subs = s.pred_from_comp(c)
        self.assertEqual(len(subs), 43
                         , 'incorrect number of substitutions')


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
