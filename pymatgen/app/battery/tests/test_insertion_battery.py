#!/usr/bin/env python

'''
Created on Jan 25, 2012
'''

__author__ = "Anubhav Jain"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Anubhav Jain"
__email__ = "ajain@lbl.gov"
__date__ = "Jan 25, 2012"

import unittest
import os

from pymatgen.entries.computed_entries import computed_entries_from_json
from pymatgen.app.battery.insertion_battery import InsertionElectrode

module_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class InsertionElectrodeTest(unittest.TestCase):

    def setUp(self):
        data_Li = '[{"correction": 0.0, "composition": {"Li":1}, "formulasum": "Li1", "energy": -1.90753119}]'
        self.entry_Li = computed_entries_from_json(data_Li)[0]

        data_LTO = open(os.path.join(module_dir, "tests", "LiTiO2.json")).readline()
        self.entries_LTO = computed_entries_from_json(data_LTO)

        self.ie_LTO = InsertionElectrode(self.entries_LTO, self.entry_Li)

    def test_voltage(self):
        #test basic voltage
        self.assertAlmostEqual(self.ie_LTO.max_voltage, 2.78583901, 3)
        self.assertAlmostEqual(self.ie_LTO.min_voltage, 0.89702381, 3)
        self.assertAlmostEqual(self.ie_LTO.get_average_voltage(), 1.84143141, 3)
        #test voltage range selectors
        self.assertAlmostEqual(self.ie_LTO.get_average_voltage(0, 1), 0.89702381, 3)
        self.assertAlmostEqual(self.ie_LTO.get_average_voltage(2, 3), 2.78583901, 3)
        #test non-existing voltage range
        self.assertAlmostEqual(self.ie_LTO.get_average_voltage(0, 0.1), 0, 3)
        self.assertAlmostEqual(self.ie_LTO.get_average_voltage(4, 5), 0, 3)

    def test_capacities(self):
        #test basic capacity
        self.assertAlmostEqual(self.ie_LTO.get_capacity_grav(), 308.74865045, 3)
        self.assertAlmostEqual(self.ie_LTO.get_capacity_vol(), 1205.99391136, 3)

        #test capacity selector
        self.assertAlmostEqual(self.ie_LTO.get_capacity_grav(1, 3), 154.374325225, 3)

        #test alternate normalization option
        self.assertAlmostEqual(self.ie_LTO.get_capacity_grav(1, 3, False), 160.803169506, 3)
        self.assertIsNotNone(self.ie_LTO.to_dict_summary(True))

    def test_entries(self):
        #test that the proper number of sub-electrodes are returned
        self.assertEqual(len(self.ie_LTO.get_sub_electrodes(False, True)), 3)
        self.assertEqual(len(self.ie_LTO.get_sub_electrodes(True, True)), 2)


if __name__ == '__main__':
    unittest.main()
