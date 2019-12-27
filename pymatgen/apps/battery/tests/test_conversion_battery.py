# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import unittest
import os
import json

from monty.json import MontyDecoder

from pymatgen import Composition
from pymatgen.apps.battery.conversion_battery import ConversionElectrode, \
    ConversionVoltagePair

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        'test_files')


class ConversionElectrodeTest(unittest.TestCase):

    def setUp(self):
        pass

    def test_init(self):
        # both 'LiCoO2' and "FeF3" are using Li+ as working ion; MnO2 is for the multivalent Mg2+ ion
        formulas = ['LiCoO2', "FeF3", "MnO2"]
        expected_properties = {}
        expected_properties['LiCoO2'] = {'average_voltage': 2.26940307125,
                                         'capacity_grav': 903.19752911225669,
                                         'capacity_vol': 2903.35804724,
                                         'specific_energy': 2049.7192465127678,
                                         'energy_density': 6588.8896693479574}
        expected_properties['FeF3'] = {'average_voltage': 3.06179925889,
                                       'capacity_grav': 601.54508701578118,
                                       'capacity_vol': 2132.2069115142394,
                                       'specific_energy': 1841.8103016131706,
                                       'energy_density': 6528.38954147}
        expected_properties['MnO2'] = {'average_voltage': 1.7127027687901726,
                                       'capacity_grav': 790.9142070034802,
                                       'capacity_vol': 3543.202003526853,
                                       'specific_energy': 1354.6009522103434,
                                       'energy_density': 6068.451881823329}
        for f in formulas:

            with open(os.path.join(test_dir, f + "_batt.json"), 'r') as fid:
                entries = json.load(fid, cls=MontyDecoder)

                # entries = computed_entries_from_json(fid.read())

            # with open(os.path.join(test_dir, f + "_batt.json"), 'w') as fid:
            # json.dump(entries, fid, cls=MontyEncoder)
            if f in ['LiCoO2', "FeF3"]:
                working_ion = "Li"
            elif f in ["MnO2"]:
                working_ion = "Mg"
            c = ConversionElectrode.from_composition_and_entries(
                Composition(f), entries, working_ion_symbol=working_ion)
            self.assertEqual(len(c.get_sub_electrodes(True)), c.num_steps)
            self.assertEqual(len(c.get_sub_electrodes(False)),
                             sum(range(1, c.num_steps + 1)))
            self.assertIsNotNone(str(c))
            p = expected_properties[f]

            for k, v in p.items():
                self.assertAlmostEqual(getattr(c, "get_" + k).__call__(), v, 2)

            self.assertIsNotNone(c.get_summary_dict(True))

            # Test pair to dict

            pair = c.voltage_pairs[0]
            d = pair.as_dict()
            pair2 = ConversionVoltagePair.from_dict(d)
            for prop in ['voltage', 'mass_charge', 'mass_discharge']:
                self.assertEqual(getattr(pair, prop), getattr(pair2, prop), 2)

            # Test
            d = c.as_dict()
            electrode = ConversionElectrode.from_dict(d)
            for k, v in p.items():
                self.assertAlmostEqual(getattr(electrode,
                                               "get_" + k).__call__(), v, 2)


if __name__ == "__main__":
    unittest.main()
