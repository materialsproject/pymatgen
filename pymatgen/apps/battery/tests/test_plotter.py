# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


from __future__ import annotations

import json
import os
import unittest

from monty.json import MontyDecoder

from pymatgen.apps.battery.conversion_battery import ConversionElectrode
from pymatgen.apps.battery.insertion_battery import InsertionElectrode
from pymatgen.apps.battery.plotter import VoltageProfilePlotter
from pymatgen.core.composition import Composition
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.util.testing import PymatgenTest


class VoltageProfilePlotterTest(unittest.TestCase):
    def setUp(self):
        entry_Li = ComputedEntry("Li", -1.90753119)

        with open(os.path.join(PymatgenTest.TEST_FILES_DIR, "LiTiO2_batt.json")) as f:
            entries_LTO = json.load(f, cls=MontyDecoder)
            self.ie_LTO = InsertionElectrode.from_entries(entries_LTO, entry_Li)

        with open(os.path.join(PymatgenTest.TEST_FILES_DIR, "FeF3_batt.json")) as fid:
            entries = json.load(fid, cls=MontyDecoder)
            self.ce_FF = ConversionElectrode.from_composition_and_entries(Composition("FeF3"), entries)

    def testName(self):
        plotter = VoltageProfilePlotter(xaxis="frac_x")
        plotter.add_electrode(self.ie_LTO, "LTO insertion")
        plotter.add_electrode(self.ce_FF, "FeF3 conversion")
        self.assertIsNotNone(plotter.get_plot_data(self.ie_LTO))
        self.assertIsNotNone(plotter.get_plot_data(self.ce_FF))

    def testPlotly(self):
        plotter = VoltageProfilePlotter(xaxis="frac_x")
        plotter.add_electrode(self.ie_LTO, "LTO insertion")
        plotter.add_electrode(self.ce_FF, "FeF3 conversion")
        fig = plotter.get_plotly_figure()
        self.assertEqual(fig.layout.xaxis.title.text, "Atomic Fraction of Li")
        plotter = VoltageProfilePlotter(xaxis="x_form")
        plotter.add_electrode(self.ce_FF, "FeF3 conversion")
        fig = plotter.get_plotly_figure()
        self.assertEqual(fig.layout.xaxis.title.text, "x in Li<sub>x</sub>FeF3")

        plotter.add_electrode(self.ie_LTO, "LTO insertion")
        fig = plotter.get_plotly_figure()
        self.assertEqual(fig.layout.xaxis.title.text, "x Workion Ion per Host F.U.")


if __name__ == "__main__":
    unittest.main()
