from __future__ import annotations

import json

import pytest
from monty.json import MontyDecoder, MontyEncoder
from pytest import approx

from pymatgen.apps.battery.insertion_battery import InsertionElectrode, InsertionVoltagePair
from pymatgen.core import Composition
from pymatgen.entries.computed_entries import ComputedEntry
from tests.testing import TEST_FILES_DIR

TEST_DIR = f"{TEST_FILES_DIR}/apps/battery"


class TestInsertionElectrode:
    def setup_method(self):
        self.entry_Li = ComputedEntry("Li", -1.90753119)
        self.entry_Ca = ComputedEntry("Ca", -1.99689568)

        with open(f"{TEST_DIR}/LiTiO2_batt.json", encoding="utf-8") as file:
            self.entries_LTO = json.load(file, cls=MontyDecoder)

        with open(f"{TEST_DIR}/MgVO_batt.json", encoding="utf-8") as file:
            self.entries_MVO = json.load(file, cls=MontyDecoder)

        with open(f"{TEST_DIR}/Mg_batt.json", encoding="utf-8") as file:
            self.entry_Mg = json.load(file, cls=MontyDecoder)

        with open(f"{TEST_DIR}/CaMoO2_batt.json", encoding="utf-8") as file:
            self.entries_CMO = json.load(file, cls=MontyDecoder)

        self.ie_LTO = InsertionElectrode.from_entries(self.entries_LTO, self.entry_Li)
        self.ie_MVO = InsertionElectrode.from_entries(self.entries_MVO, self.entry_Mg)
        self.ie_CMO = InsertionElectrode.from_entries(self.entries_CMO, self.entry_Ca)

    def test_voltage(self):
        # test basic voltage
        assert self.ie_LTO.max_voltage == approx(2.78583901)
        assert self.ie_LTO.min_voltage == approx(0.89702381)
        assert self.ie_LTO.get_average_voltage() == approx(1.84143141)
        # test voltage range selectors
        assert self.ie_LTO.get_average_voltage(0, 1) == approx(0.89702381)
        assert self.ie_LTO.get_average_voltage(2, 3) == approx(2.78583901)
        # test non-existing voltage range
        assert self.ie_LTO.get_average_voltage(0, 0.1) == approx(0)
        assert self.ie_LTO.get_average_voltage(4, 5) == approx(0)

        assert self.ie_MVO.get_average_voltage() == approx(2.513767)

    def test_capacities(self):
        # test basic capacity
        assert self.ie_LTO.get_capacity_grav() == approx(308.74865045)
        assert self.ie_LTO.get_capacity_vol() == approx(1205.99391136)

        # test capacity selector
        assert self.ie_LTO.get_capacity_grav(1, 3) == approx(154.374325225)

        # test alternate normalization option
        assert self.ie_LTO.get_capacity_grav(1, 3, use_overall_normalization=False) == approx(160.803169506)
        assert self.ie_LTO.get_summary_dict() is not None

        assert self.ie_MVO.get_capacity_grav() == approx(281.845548242)
        assert self.ie_MVO.get_capacity_vol() == approx(1145.80087994)

    def test_get_instability(self):
        assert self.ie_LTO.get_max_instability() is None
        assert self.ie_MVO.get_max_instability() == approx(0.7233711650000014)
        assert self.ie_MVO.get_min_instability() == approx(0.4913575099999994)

    def test_get_muO2(self):
        assert self.ie_LTO.get_max_muO2() is None
        assert self.ie_MVO.get_max_muO2() == approx(-4.93552791875)
        assert self.ie_MVO.get_min_muO2() == approx(-11.06599657)

    def test_entries(self):
        # test that the proper number of sub-electrodes are returned
        assert len(self.ie_LTO.get_sub_electrodes(adjacent_only=False, include_myself=True)) == 3
        assert len(self.ie_LTO.get_sub_electrodes(adjacent_only=True, include_myself=True)) == 2

    def test_get_all_entries(self):
        self.ie_LTO.get_all_entries()

    def test_as_from_dict(self):
        dct = self.ie_LTO.as_dict()
        ie = InsertionElectrode.from_dict(dct)
        assert ie.max_voltage == approx(2.78583901)
        assert ie.min_voltage == approx(0.89702381)
        assert ie.get_average_voltage() == approx(1.84143141)

        # Just to make sure JSON string works.
        json_str = json.dumps(self.ie_LTO, cls=MontyEncoder)
        ie = json.loads(json_str, cls=MontyDecoder)
        assert ie.max_voltage == approx(2.78583901)
        assert ie.min_voltage == approx(0.89702381)
        assert ie.get_average_voltage() == approx(1.84143141)

    def test_voltage_pair(self):
        volt_pair = self.ie_LTO[0]
        assert volt_pair.voltage == approx(2.78583901)
        assert volt_pair.mAh == approx(13400.7411749)
        assert volt_pair.mass_charge == approx(79.8658)
        assert volt_pair.mass_discharge == approx(83.3363)
        assert volt_pair.vol_charge == approx(37.553684467)
        assert volt_pair.vol_discharge == approx(37.917719932)
        assert volt_pair.frac_charge == approx(0)
        assert volt_pair.frac_discharge == approx(0.14285714285714285)
        assert volt_pair.x_charge == approx(0)
        assert volt_pair.x_discharge == approx(0.5)
        # reconstruct the voltage pair
        vv = InsertionVoltagePair.from_dict(volt_pair.as_dict())
        assert vv.entry_charge.energy == approx(-105.53608265)
        assert vv.voltage == approx(2.78583901)

    def test_get_summary_dict(self):
        dct = self.ie_CMO.get_summary_dict()
        assert dct["stability_charge"] == approx(0.2346574583333325)
        assert dct["stability_discharge"] == approx(0.33379544031249786)
        assert dct["muO2_data"]["mp-714969"][0]["chempot"] == approx(-4.93552791875)

        assert dct["adj_pairs"][0]["muO2_data"]["mp-714969"][0]["chempot"] == approx(-4.93552791875)
        assert dct["framework_formula"] == "MoO2"
        assert dct["adj_pairs"][1]["framework_formula"] == "MoO2"

    def test_init_no_structure(self):
        def remove_structure(entries):
            ents = []
            for entry in entries:
                dd = entry.as_dict()
                ent = ComputedEntry.from_dict(dd)
                ent.data["volume"] = entry.structure.volume
                ents.append(ent)
            return ents

        ie_CMO_no_struct = InsertionElectrode.from_entries(remove_structure(self.entries_CMO), self.entry_Ca)
        dct = ie_CMO_no_struct.get_summary_dict()
        assert dct["stability_charge"] == approx(0.2346574583333325)
        assert dct["stability_discharge"] == approx(0.33379544031249786)
        assert dct["muO2_data"]["mp-714969"][0]["chempot"] == approx(-4.93552791875)

        ie_LTO_no_struct = InsertionElectrode.from_entries(self.entries_LTO, self.entry_Li, strip_structures=True)
        vols_no_struct = [ient.data["volume"] for ient in ie_LTO_no_struct.get_all_entries()]
        vols_struct = [ient.structure.volume for ient in self.ie_LTO.get_all_entries()]
        assert vols_no_struct == approx(vols_struct, rel=1e-3)

    def test_entry_ordering(self):
        def formulas(entries):
            return [entry.reduced_formula for entry in entries]

        assert formulas(self.ie_LTO.get_stable_entries()) == ["TiO2", "LiTi2O4", "LiTiO2"]
        assert formulas(self.ie_LTO.get_stable_entries(charge_to_discharge=False)) == ["LiTiO2", "LiTi2O4", "TiO2"]
        assert formulas(self.ie_LTO.get_unstable_entries()) == ["TiO2"] * 3
        assert formulas(self.ie_LTO.get_unstable_entries(charge_to_discharge=False)) == ["TiO2"] * 3

        all_entries = self.ie_LTO.get_all_entries()
        assert len(all_entries) == 6
        assert self.ie_LTO.get_all_entries(charge_to_discharge=False) == all_entries[::-1]
        assert self.ie_LTO.fully_charged_entry.reduced_formula == "TiO2"
        assert self.ie_LTO.fully_discharged_entry.reduced_formula == "LiTiO2"

    def test_repr(self):
        lines = repr(self.ie_LTO).splitlines()
        assert lines[0] == "InsertionElectrode with endpoints at TiO2 and LiTiO2"
        assert lines[1].startswith("Avg. volt. = 1.84143")
        assert lines[2].startswith("Grav. cap. = 308.748")
        assert lines[3].startswith("Vol. cap. = 1205.99")

        lines = repr(self.ie_LTO[0]).splitlines()
        assert lines[0] == "Insertion voltage pair with working ion Li"
        assert lines[1].startswith("V = 2.78583")
        assert len(lines) == 5

    def test_as_from_dict_legacy(self):
        dct = self.ie_LTO.as_dict_legacy()
        assert dct["@class"] == "InsertionElectrode"
        assert len(dct["entries"]) == 6
        assert dct["working_ion_entry"]["energy"] == approx(-1.90753119)

        ie = InsertionElectrode.from_dict_legacy(dct)
        assert isinstance(ie, InsertionElectrode)
        assert ie.num_steps == 2
        assert ie.get_average_voltage() == approx(1.84143141)
        assert ie.get_capacity_grav() == approx(308.74865045)
        assert ie.framework_formula == "TiO2"

    def test_sub_electrodes_exclude_myself(self):
        subs = self.ie_LTO.get_sub_electrodes(adjacent_only=False, include_myself=False)
        assert [sub.num_steps for sub in subs] == [1, 1]
        assert sorted(sub.get_average_voltage() for sub in subs) == approx([0.89702381, 2.78583901])

    def test_voltage_range_selection(self):
        assert self.ie_MVO.get_max_instability(min_voltage=10) is None
        assert self.ie_MVO.get_min_muO2(max_voltage=-1) is None
        assert self.ie_LTO.get_capacity_grav(min_voltage=1, max_voltage=3) == approx(154.374325225)
        assert self.ie_LTO.get_specific_energy(min_voltage=1, max_voltage=3) == approx(154.374325225 * 2.78583901)


class TestInsertionVoltagePair:
    def setup_method(self):
        self.entry_Li = ComputedEntry("Li", -1.90753119)
        self.entry_TiO2 = ComputedEntry("TiO2", -26.38402066)
        self.entry_LiTiO2 = ComputedEntry("LiTiO2", -30.13297575)

    def test_from_entries_order_independent(self):
        pair = InsertionVoltagePair.from_entries(self.entry_TiO2, self.entry_LiTiO2, self.entry_Li)
        swapped = InsertionVoltagePair.from_entries(self.entry_LiTiO2, self.entry_TiO2, self.entry_Li)
        for vpair in (pair, swapped):
            assert vpair.entry_charge is self.entry_TiO2
            assert vpair.entry_discharge is self.entry_LiTiO2
            # V = [E(TiO2) + E(Li) - E(LiTiO2)] / e
            assert vpair.voltage == approx(-26.38402066 - 1.90753119 + 30.13297575)
            assert vpair.frac_charge == approx(0)
            assert vpair.frac_discharge == approx(0.25)
            assert vpair.x_discharge == approx(1)
            assert vpair.framework == Composition("TiO2")
            assert vpair.mass_charge == approx(Composition("TiO2").weight)
            assert vpair.mass_discharge == approx(Composition("LiTiO2").weight)

    @pytest.mark.parametrize(
        ("entry1", "entry2", "ion", "match"),
        [
            ("TiO2", "LiTiO2", ("LiO", -5.0), "working ion specified must be an element"),
            ("TiO2", "Ti2O4", ("Li", -1.9), "working ion must be present in one of the entries"),
            ("LiTiO2", "Li2Ti2O4", ("Li", -1.9), "atomic percentage cannot be the same"),
            ("TiO2", "LiVO2", ("Li", -1.9), "must have the same compositional framework"),
        ],
    )
    def test_from_entries_invalid(self, entry1, entry2, ion, match):
        with pytest.raises(ValueError, match=match):
            InsertionVoltagePair.from_entries(
                ComputedEntry(entry1, -20.0), ComputedEntry(entry2, -25.0), ComputedEntry(*ion)
            )
