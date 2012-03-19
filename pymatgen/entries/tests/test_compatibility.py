#!/usr/bin/env python

'''
Created on Mar 19, 2012
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Mar 19, 2012"

import unittest

from pymatgen.entries.compatibility import MaterialsProjectCompatibility, MITCompatibility
from pymatgen.core.structure import Composition
from pymatgen.entries.computed_entries import ComputedEntry

class MaterialsProjectCompatibilityTest(unittest.TestCase):

    def test_process_entry(self):
        compat = MaterialsProjectCompatibility()

        #Correct parameters
        entry = ComputedEntry(Composition.from_formula('Fe2O3'), -1, 0.0, parameters = {'is_hubbard':True, 'hubbards':{'Fe':5.3, 'O':0}, 'run_type':'GGA+U', 'potcar_symbols':['PAW_PBE Fe_pv 06Sep2000', 'PAW_PBE O 08Apr2002']})
        self.assertIsNotNone(compat.process_entry(entry))

        #Wrong U value
        entry = ComputedEntry(Composition.from_formula('Fe2O3'), -1, 0.0, parameters = {'is_hubbard':True, 'hubbards':{'Fe':5.2, 'O':0}, 'run_type':'GGA+U', 'potcar_symbols':['PAW_PBE Fe_pv 06Sep2000', 'PAW_PBE O 08Apr2002']})
        self.assertIsNone(compat.process_entry(entry))

        #GGA run
        entry = ComputedEntry(Composition.from_formula('Fe2O3'), -1, 0.0, parameters = {'is_hubbard':False, 'hubbards':None, 'run_type':'GGA', 'potcar_symbols':['PAW_PBE Fe_pv 06Sep2000', 'PAW_PBE O 08Apr2002']})
        self.assertIsNone(compat.process_entry(entry))

        #Wrong psp
        entry = ComputedEntry(Composition.from_formula('Fe2O3'), -1, 0.0, parameters = {'is_hubbard':True, 'hubbards':{'Fe':5.3, 'O':0}, 'run_type':'GGA+U', 'potcar_symbols':['PAW_PBE Fe 06Sep2000', 'PAW_PBE O 08Apr2002']})
        self.assertIsNone(compat.process_entry(entry))

        #Testing processing of elements.
        entry = ComputedEntry(Composition.from_formula('O'), -1, 0.0, parameters = {'is_hubbard':False, 'hubbards':{}, 'potcars':['O'], 'run_type':'GGA'})
        entry = compat.process_entry(entry)
        self.assertEqual(entry.structureid, -8)
        self.assertAlmostEqual(entry.energy, -4.22986844926)

class MITCompatibilityTest(unittest.TestCase):

    def test_process_entry(self):
        compat = MITCompatibility()

        #Correct parameters
        entry = ComputedEntry(Composition.from_formula('Fe2O3'), -1, 0.0, parameters = {'is_hubbard':True, 'hubbards':{'Fe':4.0, 'O':0}, 'run_type':'GGA+U', 'potcar_symbols':['PAW_PBE Fe 06Sep2000', 'PAW_PBE O 08Apr2002']})
        self.assertIsNotNone(compat.process_entry(entry))

        #Wrong U value
        entry = ComputedEntry(Composition.from_formula('Fe2O3'), -1, 0.0, parameters = {'is_hubbard':True, 'hubbards':{'Fe':5.2, 'O':0}, 'run_type':'GGA+U', 'potcar_symbols':['PAW_PBE Fe 06Sep2000', 'PAW_PBE O 08Apr2002']})
        self.assertIsNone(compat.process_entry(entry))

        #GGA run
        entry = ComputedEntry(Composition.from_formula('Fe2O3'), -1, 0.0, parameters = {'is_hubbard':False, 'hubbards':None, 'run_type':'GGA', 'potcar_symbols':['PAW_PBE Fe 06Sep2000', 'PAW_PBE O 08Apr2002']})
        self.assertIsNone(compat.process_entry(entry))

        #Wrong psp
        entry = ComputedEntry(Composition.from_formula('Fe2O3'), -1, 0.0, parameters = {'is_hubbard':True, 'hubbards':{'Fe':4.0, 'O':0}, 'run_type':'GGA+U', 'potcar_symbols':['PAW_PBE Fe_pv 06Sep2000', 'PAW_PBE O 08Apr2002']})
        self.assertIsNone(compat.process_entry(entry))

        #Testing processing of elements.
        entry = ComputedEntry(Composition.from_formula('O'), -1, 0.0, parameters = {'is_hubbard':False, 'hubbards':{}, 'potcars':['O'], 'run_type':'GGA'})
        entry = compat.process_entry(entry)
        self.assertEqual(entry.structureid, -8)
        self.assertAlmostEqual(entry.energy, -4.25915626315)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
