from __future__ import annotations

import json

import pytest
from monty.json import MontyDecoder
from pytest import approx

from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.apps.battery.conversion_battery import ConversionElectrode, ConversionVoltagePair
from pymatgen.core import Element
from pymatgen.core.composition import Composition
from tests.testing import TEST_FILES_DIR

TEST_DIR = f"{TEST_FILES_DIR}/apps/battery"


class TestConversionElectrode:
    def setup_method(self):
        self.formulas = ["LiCoO2", "FeF3", "MnO2"]
        self.conversion_electrodes = {}
        for formula in self.formulas:
            with open(f"{TEST_DIR}/{formula}_batt.json", encoding="utf-8") as fid:
                entries = json.load(fid, cls=MontyDecoder)
            if formula in ["LiCoO2", "FeF3"]:
                working_ion = "Li"
            elif formula == "MnO2":
                working_ion = "Mg"
            conv_elec = ConversionElectrode.from_composition_and_entries(
                Composition(formula), entries, working_ion_symbol=working_ion
            )
            self.conversion_electrodes[formula] = {
                "working_ion": working_ion,
                "CE": conv_elec,
            }

        self.expected_properties = {
            "LiCoO2": {
                "average_voltage": 2.26940307125,
                "capacity_grav": 903.19752911225669,
                "capacity_vol": 2903.35804724,
                "specific_energy": 2049.7192465127678,
                "energy_density": 6588.8896693479574,
            },
            "FeF3": {
                "average_voltage": 3.06179925889,
                "capacity_grav": 601.54508701578118,
                "capacity_vol": 2132.2069115142394,
                "specific_energy": 1841.8103016131706,
                "energy_density": 6528.38954147,
            },
            "MnO2": {
                "average_voltage": 1.7127027687901726,
                "capacity_grav": 790.9142070034802,
                "capacity_vol": 3543.202003526853,
                "specific_energy": 1354.6009522103434,
                "energy_density": 6068.451881823329,
            },
        }

        # expected composite upon discharge process, of which entry object has been simplified to reduced formula
        self.expected_composite = {
            "LiCoO2": {
                "entries_charge": [
                    ["CoO2"],
                    ["Li(CoO2)2"],
                    ["LiCoO2"],
                    ["Li6CoO4", "CoO"],
                    ["Li6CoO4", "Co"],
                ],
                "entries_discharge": [
                    ["Li(CoO2)2"],
                    ["LiCoO2"],
                    ["Li6CoO4", "CoO"],
                    ["Li6CoO4", "Co"],
                    ["Co", "Li2O"],
                ],
            },
            "FeF3": {
                "entries_charge": [["FeF3"], ["FeF2", "LiF"]],
                "entries_discharge": [["FeF2", "LiF"], ["Fe", "LiF"]],
            },
            "MnO2": {
                "entries_charge": [["MnO2"]],
                "entries_discharge": [["Mn", "MgO"]],
            },
        }

    def test_init(self):
        # both 'LiCoO2' and "FeF3" are using Li+ as working ion; MnO2 is for the multivalent Mg2+ ion
        for formula in self.formulas:
            conv_electrode = self.conversion_electrodes[formula]["CE"]

            assert len(conv_electrode.get_sub_electrodes(adjacent_only=True)) == conv_electrode.num_steps
            assert len(conv_electrode.get_sub_electrodes(adjacent_only=False)) == sum(
                range(1, conv_electrode.num_steps + 1)
            )
            props = self.expected_properties[formula]

            for key, val in props.items():
                assert getattr(conv_electrode, f"get_{key}")() == approx(val, abs=1e-2)

            assert {*conv_electrode.get_summary_dict(print_subelectrodes=True)} == {
                "adj_pairs",
                "reactions",
                "energy_vol",
                "max_voltage_step",
                "fracA_discharge",
                "energy_grav",
                "nsteps",
                "average_voltage",
                "fracA_charge",
                "working_ion",
                "reactant_compositions",
                "max_delta_volume",
                "framework_formula",
                "all_pairs",
                "min_voltage",
                "max_voltage",
                "capacity_grav",
                "capacity_vol",
            }

            # try to export/import a voltage pair via a dict
            pair = conv_electrode.voltage_pairs[0]
            dct = pair.as_dict()
            pair2 = ConversionVoltagePair.from_dict(dct)
            for prop in ["voltage", "mass_charge", "mass_discharge"]:
                assert getattr(pair, prop) == getattr(pair2, prop)  # 2

            # try to create an electrode from a dict and test methods
            dct = conv_electrode.as_dict()
            electrode = ConversionElectrode.from_dict(dct)
            for key, val in props.items():
                assert getattr(electrode, f"get_{key}")() == approx(val, abs=1e-2)

    def test_repr(self):
        conv_electrode = self.conversion_electrodes[self.formulas[0]]["CE"]
        assert (
            repr(conv_electrode)
            == "ConversionElectrode with formula='LiCoO2' and n_steps=5, avg_voltage=2.269 V, min_voltage=1.656 V, "
            "max_voltage=3.492 V\nCapacity (grav.) 903.197 mAh/g, capacity (vol.) 2903.358 Ah/l\n"
            "Specific energy 2049.719 Wh/kg, energy density 6588.890 Wh/l"
        )

    def test_summary(self):
        key_map = {"specific_energy": "energy_grav", "energy_density": "energy_vol"}
        for formula in self.formulas:
            conv_elec = self.conversion_electrodes[formula]["CE"]
            dct = conv_elec.get_summary_dict()
            props = self.expected_properties[formula]
            for k, v in props.items():
                summary_key = key_map.get(k, k)
                assert dct[summary_key] == approx(v, abs=1e-2)

    def test_composite(self):
        # check entries in charged/discharged state
        for formula in self.formulas:
            CE = self.conversion_electrodes[formula]["CE"]
            for step, volt_pair in enumerate(CE.voltage_pairs):
                # entries_charge/entries_discharge attributes should return entries equal with the expected
                composite_dict = self.expected_composite[formula]
                for attr in ["entries_charge", "entries_discharge"]:
                    # composite at each discharge step, of which entry object is simplified to reduced formula
                    entries_formula_list = [entry.reduced_formula for entry in getattr(volt_pair, attr)]
                    assert entries_formula_list == composite_dict[attr][step]


@pytest.fixture(scope="module")
def fef3_pd() -> PhaseDiagram:
    with open(f"{TEST_DIR}/FeF3_batt.json", encoding="utf-8") as file:
        return PhaseDiagram(json.load(file, cls=MontyDecoder))


@pytest.fixture(scope="module")
def fef3_electrode(fef3_pd) -> ConversionElectrode:
    return ConversionElectrode.from_composition_and_pd(Composition("FeF3"), fef3_pd)


@pytest.fixture(scope="module")
def licoo2_electrode() -> ConversionElectrode:
    with open(f"{TEST_DIR}/LiCoO2_batt.json", encoding="utf-8") as file:
        entries = json.load(file, cls=MontyDecoder)
    return ConversionElectrode.from_composition_and_entries(Composition("LiCoO2"), entries)


class TestConversionElectrodeFromPD:
    def test_from_composition_and_pd(self, fef3_electrode):
        assert fef3_electrode.initial_comp == Composition("FeF3")
        assert fef3_electrode.framework_formula == "FeF3"
        assert fef3_electrode.working_ion.symbol == "Li"
        assert fef3_electrode.num_steps == 2
        assert [pair.voltage for pair in fef3_electrode] == approx([3.07574726, 3.05482526])
        assert fef3_electrode.get_average_voltage() == approx(3.06179926)
        assert fef3_electrode.get_capacity_grav() == approx(601.545087, abs=1e-3)
        assert fef3_electrode.get_specific_energy() == approx(1841.810302, abs=1e-3)

    @pytest.mark.parametrize(
        ("formula", "n_steps", "avg_voltage"),
        [("FeF2", 1, 3.05482526), ("LiF", 1, 6.36220186), ("Li3FeF6", 3, 4.71200056)],
    )
    def test_allow_unstable(self, fef3_pd, formula, n_steps, avg_voltage):
        electrode = ConversionElectrode.from_composition_and_pd(Composition(formula), fef3_pd, allow_unstable=True)
        assert electrode.num_steps == n_steps
        assert electrode.get_average_voltage() == approx(avg_voltage)

    def test_unstable_composition_raises(self, fef3_pd):
        with pytest.raises(ValueError, match="Not stable compound found at composition"):
            ConversionElectrode.from_composition_and_pd(Composition("Li3FeF6"), fef3_pd)

    def test_single_point_profile_returns_none(self, fef3_pd):
        # Fe cannot be lithiated, so the element profile has a single step
        assert ConversionElectrode.from_composition_and_pd(Composition("Fe"), fef3_pd) is None

    def test_sub_electrodes(self, fef3_electrode):
        adjacent = fef3_electrode.get_sub_electrodes(adjacent_only=True)
        assert [sub.get_average_voltage() for sub in adjacent] == approx([3.07574726, 3.05482526])
        assert all(sub.initial_comp_formula == "FeF3" for sub in adjacent)

        all_subs = fef3_electrode.get_sub_electrodes(adjacent_only=False)
        assert [sub.num_steps for sub in all_subs] == [1, 2, 1]
        assert all_subs[1].get_capacity_grav() == approx(fef3_electrode.get_capacity_grav())

    def test_is_super_electrode(self, licoo2_electrode):
        subs = licoo2_electrode.get_sub_electrodes(adjacent_only=False)
        assert len(subs) == 15
        assert all(licoo2_electrode.is_super_electrode(sub) for sub in subs)
        assert not subs[0].is_super_electrode(licoo2_electrode)

    def test_eq_and_hash(self, fef3_electrode):
        first, second = fef3_electrode.get_sub_electrodes(adjacent_only=True)
        assert fef3_electrode == fef3_electrode.get_sub_electrodes(adjacent_only=False)[1]
        assert fef3_electrode != first  # different number of steps
        assert first != second  # same number of steps, different reactions
        assert hash(first) == hash(second)

    def test_get_summary_dict(self, fef3_electrode):
        dct = fef3_electrode.get_summary_dict(print_subelectrodes=False)
        assert "adj_pairs" not in dct
        assert "all_pairs" not in dct
        assert dct["working_ion"] == "Li"
        assert dct["nsteps"] == 2
        assert dct["reactions"] == ["FeF3 + Li -> LiF + FeF2", "FeF2 + 2 Li -> 2 LiF + Fe"]
        assert dct["reactant_compositions"] == [
            {"Fe": 1.0, "F": 3.0},
            {"Li": 1.0, "F": 1.0},
            {"Fe": 1.0, "F": 2.0},
            {"Fe": 1.0, "F": 2.0},
            {"Li": 1.0, "F": 1.0},
            {"Fe": 1.0},
        ]
        assert dct["fracA_charge"] == approx(0)
        assert dct["fracA_discharge"] == approx(3 / 7)
        assert dct["max_voltage_step"] == approx(3.07574726 - 3.05482526)

        dct = fef3_electrode.get_summary_dict(print_subelectrodes=True)
        assert len(dct["adj_pairs"]) == 2
        assert len(dct["all_pairs"]) == 3
        assert dct["adj_pairs"][1]["reactions"] == ["FeF2 + 2 Li -> 2 LiF + Fe"]


class TestConversionVoltagePair:
    def test_from_steps(self, fef3_pd, fef3_electrode):
        profile = fef3_pd.get_element_profile(Element("Li"), Composition("FeF3"))[::-1]
        pair = ConversionVoltagePair.from_steps(profile[1], profile[2], {Element("Fe"): 1, Element("F"): 3}, "FeF3")
        assert str(pair.rxn) == "FeF2 + 2 Li -> 2 LiF + Fe"
        assert pair.voltage == approx(fef3_electrode[1].voltage)
        assert pair.voltage == approx(3.05482526)
        assert pair.mAh == approx(53602.962291, abs=1e-3)
        assert pair.mass_charge == approx(119.7812096)
        assert pair.mass_discharge == approx(133.6632096)
        assert pair.vol_charge == approx(54.81127194)
        assert pair.vol_discharge == approx(62.61808506)
        assert pair.frac_charge == approx(0.2)
        assert pair.frac_discharge == approx(3 / 7)
        assert pair.x_charge == approx(1)
        assert pair.x_discharge == approx(3)
        assert pair.framework == Composition("FeF3")
        assert pair.working_ion == Element("Li")
        assert [entry.reduced_formula for entry in pair.entries_discharge] == ["Fe", "LiF"]

    def test_repr(self, fef3_electrode):
        lines = repr(fef3_electrode[0]).splitlines()
        assert lines[0] == "Conversion voltage pair with working ion Li"
        assert lines[1] == "Reaction : FeF3 + Li -> LiF + FeF2"
        assert lines[2].startswith("V = 3.07574")
        assert len(lines) == 6

    def test_electrode_sequence_protocol(self, fef3_electrode):
        first = fef3_electrode[0]
        assert first in fef3_electrode
        assert list(fef3_electrode) == list(fef3_electrode.voltage_pairs)
        assert len(fef3_electrode) == 2
