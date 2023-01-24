# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


from __future__ import annotations

import json
import os
import unittest

import pytest
from monty.json import MontyDecoder, MontyEncoder

from pymatgen.apps.battery.insertion_battery import (
    InsertionElectrode,
    InsertionVoltagePair,
)
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
        assert self.ie_LTO.max_voltage == pytest.approx(2.78583901)
        assert self.ie_LTO.min_voltage == pytest.approx(0.89702381)
        assert self.ie_LTO.get_average_voltage() == pytest.approx(1.84143141)
        # test voltage range selectors
        assert self.ie_LTO.get_average_voltage(0, 1) == pytest.approx(0.89702381)
        assert self.ie_LTO.get_average_voltage(2, 3) == pytest.approx(2.78583901)
        # test non-existing voltage range
        assert self.ie_LTO.get_average_voltage(0, 0.1) == pytest.approx(0)
        assert self.ie_LTO.get_average_voltage(4, 5) == pytest.approx(0)

        assert self.ie_MVO.get_average_voltage() == pytest.approx(2.513767)

    def test_capacities(self):
        # test basic capacity
        assert self.ie_LTO.get_capacity_grav() == pytest.approx(308.74865045)
        assert self.ie_LTO.get_capacity_vol() == pytest.approx(1205.99391136)

        # test capacity selector
        assert self.ie_LTO.get_capacity_grav(1, 3) == pytest.approx(154.374325225)

        # test alternate normalization option
        assert self.ie_LTO.get_capacity_grav(1, 3, False) == pytest.approx(160.803169506)
        assert self.ie_LTO.get_summary_dict(True) is not None

        assert self.ie_MVO.get_capacity_grav() == pytest.approx(281.845548242)
        assert self.ie_MVO.get_capacity_vol() == pytest.approx(1145.80087994)

    def test_get_instability(self):
        assert self.ie_LTO.get_max_instability() is None
        assert self.ie_MVO.get_max_instability() == pytest.approx(0.7233711650000014)
        assert self.ie_MVO.get_min_instability() == pytest.approx(0.4913575099999994)

    def test_get_muO2(self):
        assert self.ie_LTO.get_max_muO2() is None
        assert self.ie_MVO.get_max_muO2() == pytest.approx(-4.93552791875)
        assert self.ie_MVO.get_min_muO2() == pytest.approx(-11.06599657)

    def test_entries(self):
        # test that the proper number of sub-electrodes are returned
        assert len(self.ie_LTO.get_sub_electrodes(False, True)) == 3
        assert len(self.ie_LTO.get_sub_electrodes(True, True)) == 2

    def test_get_all_entries(self):
        self.ie_LTO.get_all_entries()

    def test_to_from_dict(self):
        d = self.ie_LTO.as_dict()
        ie = InsertionElectrode.from_dict(d)
        assert ie.max_voltage == pytest.approx(2.78583901)
        assert ie.min_voltage == pytest.approx(0.89702381)
        assert ie.get_average_voltage() == pytest.approx(1.84143141)

        # Just to make sure json string works.
        json_str = json.dumps(self.ie_LTO, cls=MontyEncoder)
        ie = json.loads(json_str, cls=MontyDecoder)
        assert ie.max_voltage == pytest.approx(2.78583901)
        assert ie.min_voltage == pytest.approx(0.89702381)
        assert ie.get_average_voltage() == pytest.approx(1.84143141)

    def test_voltage_pair(self):
        vpair = self.ie_LTO[0]
        assert vpair.voltage == pytest.approx(2.78583901)
        assert vpair.mAh == pytest.approx(13400.7411749)
        assert vpair.mass_charge == pytest.approx(79.8658)
        assert vpair.mass_discharge == pytest.approx(83.3363)
        assert vpair.vol_charge == pytest.approx(37.553684467)
        assert vpair.vol_discharge == pytest.approx(37.917719932)
        assert vpair.frac_charge == pytest.approx(0.0)
        assert vpair.frac_discharge == pytest.approx(0.14285714285714285)
        assert vpair.x_charge == pytest.approx(0.0)
        assert vpair.x_discharge == pytest.approx(0.5)
        # reconstruct the voltage pair
        dd = vpair.as_dict()
        vv = InsertionVoltagePair.from_dict(dd)
        assert vv.entry_charge.energy == pytest.approx(-105.53608265)
        assert vv.voltage == pytest.approx(2.78583901)

    def test_get_summary_dict(self):
        d = self.ie_CMO.get_summary_dict()
        assert d["stability_charge"] == pytest.approx(0.2346574583333325)
        assert d["stability_discharge"] == pytest.approx(0.33379544031249786)
        assert d["muO2_data"]["mp-714969"][0]["chempot"] == pytest.approx(-4.93552791875)

        assert d["adj_pairs"][0]["muO2_data"]["mp-714969"][0]["chempot"] == pytest.approx(-4.93552791875)
        assert d["framework_formula"] == "MoO2"
        assert d["adj_pairs"][1]["framework_formula"] == "MoO2"

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
        assert d["stability_charge"] == pytest.approx(0.2346574583333325)
        assert d["stability_discharge"] == pytest.approx(0.33379544031249786)
        assert d["muO2_data"]["mp-714969"][0]["chempot"] == pytest.approx(-4.93552791875)

        ie_LTO_no_struct = InsertionElectrode.from_entries(self.entries_LTO, self.entry_Li, strip_structures=True)
        vols_no_struct = [ient.data["volume"] for ient in ie_LTO_no_struct.get_all_entries()]
        vols_struct = [ient.structure.volume for ient in self.ie_LTO.get_all_entries()]
        assert vols_no_struct == pytest.approx(vols_struct, rel=1e-3)


if __name__ == "__main__":
    unittest.main()
