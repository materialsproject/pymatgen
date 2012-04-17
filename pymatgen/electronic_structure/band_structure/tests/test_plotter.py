#!/usr/bin/python

import unittest
import os
import pymatgen.electronic_structure.band_structure
from pymatgen.electronic_structure.band_structure.plotter import BSPlotter

module_dir = os.path.dirname(os.path.abspath(__file__))

class BSPlotterTest(unittest.TestCase):

    def setUp(self):
        import json
        with open(os.path.join(module_dir, "Cao_2605.json"), "rb") as f:
            d = json.loads(f.read())
            self.bs = pymatgen.electronic_structure.band_structure.band_structure.BandStructureSymmLine.from_dict(d)
            self.plotter = BSPlotter(self.bs)

    def test_bs_plot_data(self):
        self.assertEqual(len(self.plotter.bs_plot_data['distances']), 160, "wrong number of distances")
        self.assertEqual(self.plotter.bs_plot_data['ticks']['label'][5], "K", "wrong tick label")
        self.assertEqual(len(self.plotter.bs_plot_data['ticks']['label']), 19, "wrong number of tick labels")

if __name__ == '__main__':
    unittest.main()
