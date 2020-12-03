# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import json
import os
import unittest

from pymatgen import MontyDecoder, MontyEncoder
from pymatgen.apps.battery.insertion_battery import InsertionElectrode
from pymatgen.entries.computed_entries import ComputedEntry

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "test_files")


class InsertionElectrodeTest(unittest.TestCase):
    def setUp(self):
        self.entry_Li = ComputedEntry("Li", -1.90753119)
        self.entry_Ca = ComputedEntry("Ca", -1.99689568)

        with open(os.path.join(test_dir, "LiTiO2_batt.json"), "r") as f:
            self.entries_LTO = json.load(f, cls=MontyDecoder)

        with open(os.path.join(test_dir, "MgVO_batt.json"), "r") as file:
            self.entries_MVO = json.load(file, cls=MontyDecoder)

        with open(os.path.join(test_dir, "Mg_batt.json"), "r") as file:
            self.entry_Mg = json.load(file, cls=MontyDecoder)

        with open(os.path.join(test_dir, "CaMoO2_batt.json"), "r") as f:
            self.entries_CMO = json.load(f, cls=MontyDecoder)

        self.ie_LTO = InsertionElectrode(self.entries_LTO, self.entry_Li)
        self.ie_MVO = InsertionElectrode(self.entries_MVO, self.entry_Mg)
        self.ie_CMO = InsertionElectrode(self.entries_CMO, self.entry_Ca)

    def test_voltage(self):
        # test basic voltage
        self.assertAlmostEqual(self.ie_LTO.max_voltage, 2.78583901, 3)
        self.assertAlmostEqual(self.ie_LTO.min_voltage, 0.89702381, 3)
        self.assertAlmostEqual(self.ie_LTO.get_average_voltage(), 1.84143141, 3)
        # test voltage range selectors
        self.assertAlmostEqual(self.ie_LTO.get_average_voltage(0, 1), 0.89702381, 3)
        self.assertAlmostEqual(self.ie_LTO.get_average_voltage(2, 3), 2.78583901, 3)
        # test non-existing voltage range
        self.assertAlmostEqual(self.ie_LTO.get_average_voltage(0, 0.1), 0, 3)
        self.assertAlmostEqual(self.ie_LTO.get_average_voltage(4, 5), 0, 3)

        self.assertAlmostEqual(self.ie_MVO.get_average_voltage(), 2.513767, 3)

    def test_capacities(self):
        # test basic capacity
        self.assertAlmostEqual(self.ie_LTO.get_capacity_grav(), 308.74865045, 3)
        self.assertAlmostEqual(self.ie_LTO.get_capacity_vol(), 1205.99391136, 3)

        # test capacity selector
        self.assertAlmostEqual(self.ie_LTO.get_capacity_grav(1, 3), 154.374325225, 3)

        # test alternate normalization option
        self.assertAlmostEqual(
            self.ie_LTO.get_capacity_grav(1, 3, False), 160.803169506, 3
        )
        self.assertIsNotNone(self.ie_LTO.as_dict_summary(True))

        self.assertAlmostEqual(self.ie_MVO.get_capacity_grav(), 281.845548242, 3)
        self.assertAlmostEqual(self.ie_MVO.get_capacity_vol(), 1145.80087994, 3)

    def test_get_instability(self):
        self.assertIsNone(self.ie_LTO.get_max_instability())
        self.assertAlmostEqual(self.ie_MVO.get_max_instability(), 0.7233711650000014)
        self.assertAlmostEqual(self.ie_MVO.get_min_instability(), 0.4913575099999994)

    def test_get_muO2(self):
        self.assertIsNone(self.ie_LTO.get_max_muO2())
        self.assertAlmostEqual(self.ie_MVO.get_max_muO2(), -4.93552791875)
        self.assertAlmostEqual(self.ie_MVO.get_min_muO2(), -11.06599657)

    def test_entries(self):
        # test that the proper number of sub-electrodes are returned
        self.assertEqual(len(self.ie_LTO.get_sub_electrodes(False, True)), 3)
        self.assertEqual(len(self.ie_LTO.get_sub_electrodes(True, True)), 2)

    def test_get_all_entries(self):
        self.ie_LTO.get_all_entries()

    def test_to_from_dict(self):
        d = self.ie_LTO.as_dict()
        ie = InsertionElectrode.from_dict(d)
        self.assertAlmostEqual(ie.max_voltage, 2.78583901, 3)
        self.assertAlmostEqual(ie.min_voltage, 0.89702381, 3)
        self.assertAlmostEqual(ie.get_average_voltage(), 1.84143141, 3)

        # Just to make sure json string works.
        json_str = json.dumps(self.ie_LTO, cls=MontyEncoder)
        ie = json.loads(json_str, cls=MontyDecoder)
        self.assertAlmostEqual(ie.max_voltage, 2.78583901, 3)
        self.assertAlmostEqual(ie.min_voltage, 0.89702381, 3)
        self.assertAlmostEqual(ie.get_average_voltage(), 1.84143141, 3)

    def test_voltage_pair(self):
        vpair = self.ie_LTO[0]
        self.assertAlmostEqual(vpair.voltage, 2.78583901)
        self.assertAlmostEqual(vpair.mAh, 13400.7411749, 2)
        self.assertAlmostEqual(vpair.mass_charge, 79.8658)
        self.assertAlmostEqual(vpair.mass_discharge, 83.3363)
        self.assertAlmostEqual(vpair.vol_charge, 37.553684467)
        self.assertAlmostEqual(vpair.vol_discharge, 37.917719932)
        self.assertAlmostEqual(vpair.frac_charge, 0.0)
        self.assertAlmostEqual(vpair.frac_discharge, 0.14285714285714285)

    def test_as_dict_summary(self):
        d = self.ie_CMO.as_dict_summary()
        self.assertAlmostEqual(d["stability_charge"], 0.2346574583333325)
        self.assertAlmostEqual(d["stability_discharge"], 0.33379544031249786)
        self.assertAlmostEqual(
            d["muO2_data"]["mp-714969"][0]["chempot"], -4.93552791875
        )


if __name__ == "__main__":
    unittest.main()
