from __future__ import annotations

import json
import os
import unittest

from monty.json import MontyDecoder
from pytest import approx

from pymatgen.apps.battery.conversion_battery import ConversionElectrode, ConversionVoltagePair
from pymatgen.core.composition import Composition
from pymatgen.util.testing import TEST_FILES_DIR


class TestConversionElectrode(unittest.TestCase):
    def setUp(self):
        self.formulas = ["LiCoO2", "FeF3", "MnO2"]
        self.conversion_eletrodes = {}
        for f in self.formulas:
            with open(os.path.join(TEST_FILES_DIR, f + "_batt.json")) as fid:
                entries = json.load(fid, cls=MontyDecoder)
            if f in ["LiCoO2", "FeF3"]:
                working_ion = "Li"
            elif f in ["MnO2"]:
                working_ion = "Mg"
            c = ConversionElectrode.from_composition_and_entries(
                Composition(f), entries, working_ion_symbol=working_ion
            )
            self.conversion_eletrodes[f] = {"working_ion": working_ion, "CE": c}

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
                "entries_charge": [["CoO2"], ["Li(CoO2)2"], ["LiCoO2"], ["Li6CoO4", "CoO"], ["Li6CoO4", "Co"]],
                "entries_discharge": [["Li(CoO2)2"], ["LiCoO2"], ["Li6CoO4", "CoO"], ["Li6CoO4", "Co"], ["Co", "Li2O"]],
            },
            "FeF3": {
                "entries_charge": [["FeF3"], ["FeF2", "LiF"]],
                "entries_discharge": [["FeF2", "LiF"], ["Fe", "LiF"]],
            },
            "MnO2": {"entries_charge": [["MnO2"]], "entries_discharge": [["Mn", "MgO"]]},
        }

    def test_init(self):
        # both 'LiCoO2' and "FeF3" are using Li+ as working ion; MnO2 is for the multivalent Mg2+ ion
        for f in self.formulas:
            c = self.conversion_eletrodes[f]["CE"]

            assert len(c.get_sub_electrodes(True)) == c.num_steps
            assert len(c.get_sub_electrodes(False)) == sum(range(1, c.num_steps + 1))
            assert str(c) is not None
            p = self.expected_properties[f]

            for k, v in p.items():
                assert getattr(c, "get_" + k)() == approx(v, abs=1e-2)

            assert c.get_summary_dict(True) is not None

            # try to export/import a voltage pair via a dict
            pair = c.voltage_pairs[0]
            d = pair.as_dict()
            pair2 = ConversionVoltagePair.from_dict(d)
            for prop in ["voltage", "mass_charge", "mass_discharge"]:
                assert getattr(pair, prop) == getattr(pair2, prop), 2

            # try to create an electrode from a dict and test methods
            d = c.as_dict()
            electrode = ConversionElectrode.from_dict(d)
            for k, v in p.items():
                assert getattr(electrode, "get_" + k)() == approx(v, abs=1e-2)

    def test_summary(self):
        kmap = {"specific_energy": "energy_grav", "energy_density": "energy_vol"}
        for f in self.formulas:
            c = self.conversion_eletrodes[f]["CE"]
            d = c.get_summary_dict()
            p = self.expected_properties[f]
            for k, v in p.items():
                summary_key = kmap.get(k, k)
                assert d[summary_key] == approx(v, abs=1e-2)

    def test_composite(self):
        # check entries in charged/discharged state
        for formula in self.formulas:
            CE = self.conversion_eletrodes[formula]["CE"]
            for step, vpair in enumerate(CE.voltage_pairs):
                # entries_charge/entries_discharge attributes should return entries equal with the expected
                composite_dict = self.expected_composite[formula]
                for attri in ["entries_charge", "entries_discharge"]:
                    # composite at each discharge step, of which entry object is simplified to reduced formula
                    entries_formula_list = [entry.composition.reduced_formula for entry in getattr(vpair, attri)]
                    assert entries_formula_list == composite_dict[attri][step]
