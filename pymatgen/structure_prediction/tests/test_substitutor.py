#!/usr/bin/env python

import unittest

from pymatgen.core.periodic_table import Specie
from pymatgen.core.structure import Composition
from pymatgen.structure_prediction.substitutor import Substitutor

class SubstitutorTest(unittest.TestCase):
    def setUp(self):
        self.s = Substitutor()
    
    def test_substitutor(self):
        s_list = [Specie('Ag', 2), Specie('Cl', -1)]
        subs =  self.s.pred_from_list(s_list)
        self.assertEqual(len(subs), 197
                         , 'incorrect number of substitutions')
        c = Composition({'Ag2+': 1, 'Cl1-':2})
        subs = self.s.pred_from_comp(c)
        self.assertEqual(len(subs), 43
                         , 'incorrect number of substitutions')
        
    def test_to_dict(self):
        Substitutor.from_dict(self.s.to_dict)
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
