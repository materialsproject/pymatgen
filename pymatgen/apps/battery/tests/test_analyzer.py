import unittest2 as unittest
import os
from pymatgen.apps.battery.analyzer import BatteryAnalyzer

from pymatgen.util.testing import PymatgenTest
from pymatgen.core.structure import Structure

import pymatgen

module_dir = os.path.join(os.path.dirname(pymatgen.__file__), '..', 'test_files')


class BatteryAnalyzerTest(PymatgenTest):

    def load_from_cif(self, filename, oxidations, cation="Li"):
        s = Structure.from_file(os.path.join(module_dir, filename))
        s.add_oxidation_state_by_element(oxidations)
        return BatteryAnalyzer(s, cation)

    def load_from_internal(self, name, oxidations, cation="Li"):
        s = self.get_structure(name).copy()
        s.add_oxidation_state_by_element(oxidations)
        return BatteryAnalyzer(s, cation)

    def setUp(self):
        self.lifepo4 = self.load_from_internal("LiFePO4",{'Li': 1, 'Fe': 2, 'P': 5, 'O': -2})
        self.nafepo4 = self.load_from_internal("NaFePO4", {'Na': 1, 'Fe': 2, 'P': 5, 'O': -2}, cation='Na')
        self.fepo4 = self.load_from_cif("FePO4a.cif", {'Fe': 3, 'P': 5, 'O': -2})
        self.lifemnpo4 = self.load_from_cif("Li4Fe3Mn1(PO4)4.cif", {'Li': 1, 'Fe': 2, 'Mn': 2, 'P': 5, 'O': -2})
        self.li8nicofe208 = self.load_from_cif("Li8Fe2NiCoO8.cif", {'Li': 1, 'Fe': 2, 'Mn': 2, 'Co': 2, 'Ni': 2, 'O': -2})
        self.li3v2p3o12 = self.load_from_internal("Li3V2(PO4)3", {'Li': 1, 'V': 3, 'O': -2, 'P': 5})

    def test_oxid_check(self):
        s = self.get_structure("LiFePO4")
        self.assertRaises(ValueError, BatteryAnalyzer, s, 'Li')

    def test_capacitygrav_calculations(self):
        lifepo4_cap = 169.89053 # same as fepo4 cap
        nafepo4_cap = 154.20331
        li3v2p3o12_cap_remove = 197.25339
        li3v2p3o12_cap_insert = 127.17129

        self.assertAlmostEqual(self.lifepo4.get_max_capgrav(), lifepo4_cap, 3)
        self.assertEqual(self.lifepo4.get_max_capgrav(remove = False), 0)
        self.assertAlmostEqual(self.lifepo4.get_max_capgrav(insert = False), lifepo4_cap, 3)

        self.assertAlmostEqual(self.nafepo4.get_max_capgrav(), nafepo4_cap, 3)
        self.assertEqual(self.nafepo4.get_max_capgrav(remove = False), 0)

        self.assertAlmostEqual(self.fepo4.get_max_capgrav(), lifepo4_cap, 3)
        self.assertEqual(self.fepo4.get_max_capgrav(insert = False), 0)

        self.assertAlmostEqual(self.li3v2p3o12.get_max_capgrav(insert=False), li3v2p3o12_cap_remove, 3)
        self.assertAlmostEqual(self.li3v2p3o12.get_max_capgrav(remove=False), li3v2p3o12_cap_insert, 3)

    def test_capacityvol_calculations(self):

        lifepo4_cap = 594.17518
        nafepo4_cap = 542.86104
        fepo4_cap = 624.82289  #this is different than lifepo4 cap if lifepo4 volume not known

        self.assertAlmostEqual(self.lifepo4.get_max_capvol(), lifepo4_cap, 3)
        self.assertEqual(self.lifepo4.get_max_capvol(remove = False), 0)
        self.assertAlmostEqual(self.lifepo4.get_max_capvol(insert = False), lifepo4_cap, 3)

        self.assertAlmostEqual(self.nafepo4.get_max_capvol(), nafepo4_cap, 3)
        self.assertEqual(self.nafepo4.get_max_capvol(remove = False), 0)
        self.assertAlmostEqual(self.nafepo4.get_max_capvol(insert = False), nafepo4_cap, 3)

        self.assertAlmostEqual(self.fepo4.get_max_capvol(), fepo4_cap, 3)
        self.assertAlmostEqual(self.fepo4.get_max_capvol(remove = False), fepo4_cap, 3)
        self.assertEqual(self.fepo4.get_max_capvol(insert = False), 0)

        #give the lifepo4 volume, should get lifepo4 capacity
        self.assertAlmostEqual(self.fepo4.get_max_capvol(volume = self.lifepo4.struc_oxid.volume), lifepo4_cap, 3)

    def test_delithiation(self):

        self.assertEqual(self.lifemnpo4.get_removals_int_oxid(), set([1.0, 2.0, 3.0, 4.0]))

        self.assertEqual(self.li8nicofe208.get_removals_int_oxid(), set([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]))

        self.assertEqual(self.li3v2p3o12.get_removals_int_oxid(), set([4.0, 6.0]))


if __name__ == '__main__':
    unittest.main()