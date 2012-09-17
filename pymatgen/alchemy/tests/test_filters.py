#!/usr/bin/env python
from pymatgen.alchemy.filters import ContainsSpecieFilter
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.periodic_table import Element, Specie

import unittest

class ContainsSpecieFilterTest(unittest.TestCase):
    def test_filtering(self):
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.75, 0.75])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.25, 0.25, 0.25])
        lattice = Lattice([[ 3.0, 0.0, 0.0]
                           , [1.0, 3.0, 0.00]
                           , [0.00, -2.0, 3.0]])
        s = Structure(lattice
                      , [{"Si4+": 0.5, "O2-": 0.25, "P5+": 0.25}
                         , {"Si4+": 0.5, "O2-": 0.25, "P5+": 0.25}
                         , {"Si4+": 0.5, "O2-": 0.25, "P5+": 0.25}
                         , {"Si4+": 0.5, "O2-": 0.25, "P5+": 0.25}] 
                      , coords)
        
        species1 = [Specie('Si', 5), Specie('Mg', 2)]
        f1 = ContainsSpecieFilter(species1, strict_compare = True
                                  , AND = False)
        self.assertFalse(f1.test(s), 'Incorrect filter')
        f2 = ContainsSpecieFilter(species1, strict_compare = False
                                  , AND = False)
        self.assertTrue(f2.test(s), 'Incorrect filter')
        species2 = [Specie('Si', 4), Specie('Mg', 2)]
        f3 = ContainsSpecieFilter(species2, strict_compare = True
                                  , AND = False)
        self.assertTrue(f3.test(s), 'Incorrect filter')
        f4 = ContainsSpecieFilter(species2, strict_compare = False
                                  , AND = False)
        self.assertTrue(f4.test(s), 'Incorrect filter')
        
        species3 = [Specie('Si', 5), Specie('O', -2)]
        f5 = ContainsSpecieFilter(species3, strict_compare = True
                                  , AND = True)
        self.assertFalse(f5.test(s), 'Incorrect filter')
        f6 = ContainsSpecieFilter(species3, strict_compare = False
                                  , AND = True)
        self.assertTrue(f6.test(s), 'Incorrect filter')
        species4 = [Specie('Si', 4), Specie('Mg', 2)]
        f7 = ContainsSpecieFilter(species4, strict_compare = True
                                  , AND = True)
        self.assertFalse(f7.test(s), 'Incorrect filter')
        f8 = ContainsSpecieFilter(species4, strict_compare = False
                                  , AND = True)
        self.assertFalse(f8.test(s), 'Incorrect filter')
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
