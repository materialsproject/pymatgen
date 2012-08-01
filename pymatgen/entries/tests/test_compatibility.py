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

from pymatgen.entries.compatibility import MaterialsProjectCompatibility, \
    MITCompatibility
from pymatgen.entries.computed_entries import ComputedEntry


class MaterialsProjectCompatibilityTest(unittest.TestCase):

    def test_process_entry(self):
        compat = MaterialsProjectCompatibility()

        #Correct parameters
        entry = ComputedEntry('Fe2O3', -1, 0.0,
                              parameters={'is_hubbard': True,
                                          'hubbards': {'Fe': 5.3, 'O': 0},
                                          'run_type': 'GGA+U',
                                          'potcar_symbols':
        ['PAW_PBE Fe_pv 06Sep2000', 'PAW_PBE O 08Apr2002']})
        self.assertIsNotNone(compat.process_entry(entry))

        #Check actual correction
        self.assertAlmostEqual(compat.process_entry(entry).correction,
                               - 2.733 * 2)

        entry = ComputedEntry('FeF3', -2, 0.0,
                              parameters={'is_hubbard': True,
                                          'hubbards': {'Fe': 5.3, 'F': 0},
                                          'run_type': 'GGA+U',
                                          'potcar_symbols':
        ['PAW_PBE Fe_pv 06Sep2000', 'PAW_PBE F 08Apr2002']})
        self.assertIsNotNone(compat.process_entry(entry))

        #Check actual correction
        self.assertAlmostEqual(compat.process_entry(entry).correction, -2.733)

        #Wrong U value
        entry = ComputedEntry('Fe2O3', -1, 0.0,
                              parameters={'is_hubbard': True,
                                          'hubbards': {'Fe': 5.2, 'O': 0},
                                          'run_type': 'GGA+U',
                                          'potcar_symbols':
        ['PAW_PBE Fe_pv 06Sep2000', 'PAW_PBE O 08Apr2002']})
        self.assertIsNone(compat.process_entry(entry))

        #GGA run of U
        entry = ComputedEntry('Fe2O3', -1, 0.0,
                              parameters={'is_hubbard': False,
                                          'hubbards': None,
                                          'run_type': 'GGA',
                                          'potcar_symbols':
        ['PAW_PBE Fe_pv 06Sep2000', 'PAW_PBE O 08Apr2002']})
        self.assertIsNone(compat.process_entry(entry))

        #GGA+U run of non-U
        entry = ComputedEntry('Al2O3', -1, 0.0,
                              parameters={'is_hubbard': True,
                                          'hubbards': {'Al': 5.3, 'O': 0},
                                          'run_type': 'GGA+U',
                                          'potcar_symbols':
        ['PAW_PBE Al 06Sep2000', 'PAW_PBE O 08Apr2002']})
        self.assertIsNone(compat.process_entry(entry))

        #Materials project should not have a U for sulfides
        entry = ComputedEntry('FeS2', -2, 0.0,
                              parameters={'is_hubbard': True,
                                          'hubbards': {'Fe': 5.3, 'S': 0},
                                          'run_type': 'GGA+U',
                                          'potcar_symbols':
        ['PAW_PBE Fe_pv 06Sep2000', 'PAW_PBE S 08Apr2002']})
        self.assertIsNone(compat.process_entry(entry))

        #Wrong psp
        entry = ComputedEntry('Fe2O3', -1, 0.0,
                              parameters={'is_hubbard': True,
                                          'hubbards': {'Fe': 5.3, 'O': 0},
                                          'run_type': 'GGA+U',
                                          'potcar_symbols':
        ['PAW_PBE Fe 06Sep2000', 'PAW_PBE O 08Apr2002']})
        self.assertIsNone(compat.process_entry(entry))

        #Testing processing of elements.
        entry = ComputedEntry('O', -1, 0.0,
                              parameters={'is_hubbard': False, 'hubbards': {},
                                          'potcar_symbols':
        ['PAW_PBE O 08Apr2002'], 'run_type': 'GGA'})
        entry = compat.process_entry(entry)
        self.assertEqual(entry.structureid, -8)
        self.assertAlmostEqual(entry.energy, -4.22986844926)

    def test_requires_hubbard(self):
        compat = MaterialsProjectCompatibility()
        self.assertTrue(compat.requires_hubbard("Fe2O3"))
        self.assertTrue(compat.requires_hubbard("FeSO4"))
        self.assertFalse(compat.requires_hubbard("FeS2"))
        self.assertFalse(compat.requires_hubbard("Li2O"))
        self.assertTrue(compat.requires_hubbard("FeOF"))


class MITCompatibilityTest(unittest.TestCase):

    def test_process_entry(self):
        compat = MITCompatibility()

        #Correct parameters
        entry = ComputedEntry('Fe2O3', -1, 0.0,
                              parameters={'is_hubbard': True,
                                          'hubbards': {'Fe': 4.0, 'O': 0},
                                          'run_type': 'GGA+U',
                                          'potcar_symbols':
        ['PAW_PBE Fe 06Sep2000', 'PAW_PBE O 08Apr2002']})
        self.assertIsNotNone(compat.process_entry(entry))
        self.assertAlmostEqual(compat.process_entry(entry).correction,
                               - 1.723 * 2)

        entry = ComputedEntry('FeF3', -2, 0.0,
                              parameters={'is_hubbard': True,
                                          'hubbards': {'Fe': 4.0, 'F': 0},
                                          'run_type': 'GGA+U',
                                          'potcar_symbols':
        ['PAW_PBE Fe 06Sep2000', 'PAW_PBE F 08Apr2002']})
        self.assertIsNotNone(compat.process_entry(entry))

        #Check actual correction
        self.assertAlmostEqual(compat.process_entry(entry).correction, -1.723)

        #MIT should not have a U for sulfides
        entry = ComputedEntry('FeS2', -2, 0.0,
                              parameters={'is_hubbard': True,
                                          'hubbards': {'Fe': 1.9, 'S': 0},
                                          'run_type': 'GGA+U',
                                          'potcar_symbols':
        ['PAW_PBE Fe 06Sep2000', 'PAW_PBE S 08Apr2002']})
        self.assertIsNotNone(compat.process_entry(entry))

        self.assertAlmostEqual(compat.process_entry(entry).correction, -1.113)

        #Wrong U value
        entry = ComputedEntry('Fe2O3', -1, 0.0,
                              parameters={'is_hubbard': True,
                                          'hubbards': {'Fe': 5.2, 'O': 0},
                                          'run_type': 'GGA+U',
                                          'potcar_symbols':
        ['PAW_PBE Fe 06Sep2000', 'PAW_PBE O 08Apr2002']})
        self.assertIsNone(compat.process_entry(entry))

        #GGA run
        entry = ComputedEntry('Fe2O3', -1, 0.0,
                              parameters={'is_hubbard': False,
                                          'hubbards': None,
                                          'run_type': 'GGA',
                                          'potcar_symbols':
        ['PAW_PBE Fe 06Sep2000', 'PAW_PBE O 08Apr2002']})
        self.assertIsNone(compat.process_entry(entry))

        #Wrong psp
        entry = ComputedEntry('Fe2O3', -1, 0.0,
                              parameters={'is_hubbard': True,
                                          'hubbards': {'Fe': 4.0, 'O': 0},
                                          'run_type': 'GGA+U',
                                          'potcar_symbols':
        ['PAW_PBE Fe_pv 06Sep2000', 'PAW_PBE O 08Apr2002']})
        self.assertIsNone(compat.process_entry(entry))

        #Testing processing of elements.
        entry = ComputedEntry('O', -1, 0.0,
                              parameters={'is_hubbard': False, 'hubbards': {},
                                          'potcar_symbols':
        ['PAW_PBE O 08Apr2002'], 'run_type': 'GGA'})
        entry = compat.process_entry(entry)
        self.assertEqual(entry.structureid, -8)
        self.assertAlmostEqual(entry.energy, -4.25915626315)

    def test_requires_hubbard(self):
        compat = MITCompatibility()
        self.assertTrue(compat.requires_hubbard("Fe2O3"))
        self.assertTrue(compat.requires_hubbard("FeSO4"))
        self.assertTrue(compat.requires_hubbard("FeS2"))
        self.assertFalse(compat.requires_hubbard("Li2O"))
        self.assertTrue(compat.requires_hubbard("FeOF"))

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
