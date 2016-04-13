# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

'''
Created on Jul 15, 2012
'''


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Jul 15, 2012"

import unittest2 as unittest
import json
import os

from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.apps.battery.insertion_battery import InsertionElectrode
from pymatgen.apps.battery.conversion_battery import ConversionElectrode
from pymatgen import MontyDecoder, Composition
from pymatgen.apps.battery.plotter import VoltageProfilePlotter

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        'test_files')


class VoltageProfilePlotterTest(unittest.TestCase):

    def testName(self):
        entry_Li = ComputedEntry("Li", -1.90753119)

        with open(os.path.join(test_dir, "LiTiO2_batt.json"), "r") as f:
            entries_LTO = json.load(f, cls=MontyDecoder)
            ie_LTO = InsertionElectrode(entries_LTO, entry_Li)

        with open(os.path.join(test_dir, "FeF3_batt.json"), 'r') as fid:
            entries = json.load(fid, cls=MontyDecoder)
            ce_FF = ConversionElectrode.from_composition_and_entries(
                Composition("FeF3"),
                entries)

        plotter = VoltageProfilePlotter(xaxis="frac_x")
        plotter.add_electrode(ie_LTO, "LTO insertion")
        plotter.add_electrode(ce_FF, "FeF3 conversion")
        self.assertIsNotNone(plotter.get_plot_data(ie_LTO))
        self.assertIsNotNone(plotter.get_plot_data(ce_FF))


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
