# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import json
import os
import unittest

from monty.json import MontyDecoder, MontyEncoder

from pymatgen.apps.battery.insertion_battery import InsertionElectrode, InsertionVoltagePair
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.util.testing import PymatgenTest


class InsertionElectrodeTest(unittest.TestCase):
    def setUp(self):
        self.entry_Li = ComputedEntry("Li", -1.90753119)
        self.entry_Ca = ComputedEntry("Ca", -1.99689568)

        with open(os.path.join(PymatgenTest.TEST_FILES_DIR, "LiTiO2_batt.json")) as f:
            self.entries_LTO = json.load(f, cls=MontyDecoder)

        with open(os.path.join(PymatgenTest.TEST_FILES_DIR, "MgVO_batt.json")) as file:
            self.entries_MVO = json.load(file, cls=MontyDecoder)

        with open(os.path.join(PymatgenTest.TEST_FILES_DIR, "Mg_batt.json")) as file:
            self.entry_Mg = json.load(file, cls=MontyDecoder)

        with open(os.path.join(PymatgenTest.TEST_FILES_DIR, "CaMoO2_batt.json")) as f:
            self.entries_CMO = json.load(f, cls=MontyDecoder)

        self.ie_LTO = InsertionElectrode.from_entries(self.entries_LTO, self.entry_Li)
        self.ie_MVO = InsertionElectrode.from_entries(self.entries_MVO, self.entry_Mg)
        self.ie_CMO = InsertionElectrode.from_entries(self.entries_CMO, self.entry_Ca)

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
        self.assertAlmostEqual(self.ie_LTO.get_capacity_grav(1, 3, False), 160.803169506, 3)
        self.assertIsNotNone(self.ie_LTO.get_summary_dict(True))

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
        self.assertAlmostEqual(vpair.x_charge, 0.0)
        self.assertAlmostEqual(vpair.x_discharge, 0.5)
        # reconstruct the voltage pair
        dd = vpair.as_dict()
        vv = InsertionVoltagePair.from_dict(dd)
        self.assertAlmostEqual(vv.entry_charge.energy, -105.53608265)
        self.assertAlmostEqual(vv.voltage, 2.78583901)

    def test_get_summary_dict(self):
        d = self.ie_CMO.get_summary_dict()
        self.assertAlmostEqual(d["stability_charge"], 0.2346574583333325)
        self.assertAlmostEqual(d["stability_discharge"], 0.33379544031249786)
        self.assertAlmostEqual(d["muO2_data"]["mp-714969"][0]["chempot"], -4.93552791875)

        self.assertAlmostEqual(d["adj_pairs"][0]["muO2_data"]["mp-714969"][0]["chempot"], -4.93552791875)
        self.assertAlmostEqual(d["framework_formula"], "MoO2")
        self.assertAlmostEqual(d["adj_pairs"][1]["framework_formula"], "MoO2")

    def test_init_no_structure(self):
        def remove_structure(entries):
            ents = []
            for ient in entries:
                dd = ient.as_dict()
                ent = ComputedEntry.from_dict(dd)
                ent.data["volume"] = ient.structure.volume
                ents.append(ent)
            return ents

        ie_CMO_no_struct = InsertionElectrode.from_entries(remove_structure(self.entries_CMO), self.entry_Ca)
        d = ie_CMO_no_struct.get_summary_dict()
        self.assertAlmostEqual(d["stability_charge"], 0.2346574583333325)
        self.assertAlmostEqual(d["stability_discharge"], 0.33379544031249786)
        self.assertAlmostEqual(d["muO2_data"]["mp-714969"][0]["chempot"], -4.93552791875)

        ie_LTO_no_struct = InsertionElectrode.from_entries(self.entries_LTO, self.entry_Li, strip_structures=True)
        vols_no_struct = [ient.data["volume"] for ient in ie_LTO_no_struct.get_all_entries()]
        vols_struct = [ient.structure.volume for ient in self.ie_LTO.get_all_entries()]
        self.assertAlmostEqual(vols_no_struct, vols_struct)


if __name__ == "__main__":
    unittest.main()
