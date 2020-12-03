# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import json
import os
import unittest

from pymatgen import Composition, MontyDecoder
from pymatgen.apps.battery.conversion_battery import ConversionElectrode
from pymatgen.apps.battery.insertion_battery import InsertionElectrode
from pymatgen.apps.battery.plotter import VoltageProfilePlotter
from pymatgen.entries.computed_entries import ComputedEntry

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "test_files")


class VoltageProfilePlotterTest(unittest.TestCase):
    def testName(self):
        entry_Li = ComputedEntry("Li", -1.90753119)

        with open(os.path.join(test_dir, "LiTiO2_batt.json"), "r") as f:
            entries_LTO = json.load(f, cls=MontyDecoder)
            ie_LTO = InsertionElectrode(entries_LTO, entry_Li)

        with open(os.path.join(test_dir, "FeF3_batt.json"), "r") as fid:
            entries = json.load(fid, cls=MontyDecoder)
            ce_FF = ConversionElectrode.from_composition_and_entries(
                Composition("FeF3"), entries
            )

        plotter = VoltageProfilePlotter(xaxis="frac_x")
        plotter.add_electrode(ie_LTO, "LTO insertion")
        plotter.add_electrode(ce_FF, "FeF3 conversion")
        self.assertIsNotNone(plotter.get_plot_data(ie_LTO))
        self.assertIsNotNone(plotter.get_plot_data(ce_FF))


if __name__ == "__main__":
    unittest.main()
