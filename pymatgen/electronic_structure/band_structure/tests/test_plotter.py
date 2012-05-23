#!/usr/bin/python

import unittest
import os
import json

from pymatgen.electronic_structure.band_structure.band_structure import BandStructureSymmLine
from pymatgen.electronic_structure.band_structure.plotter import BSPlotter

import pymatgen

test_dir = os.path.join(os.path.dirname(os.path.abspath(pymatgen.__file__)), '..', 'test_files')

class BSPlotterTest(unittest.TestCase):

    def setUp(self):
        with open(os.path.join(test_dir, "CaO_2605_bandstructure.json"), "rb") as f:
            d = json.loads(f.read())
            self.bs = BandStructureSymmLine.from_dict(d)
            self.plotter = BSPlotter(self.bs)


    def test_bs_plot_data(self):

        self.assertEqual(len(self.plotter.bs_plot_data()['distances']), 160, "wrong number of distances")
        self.assertEqual(self.plotter.bs_plot_data()['ticks']['label'][5], "K", "wrong tick label")
        self.assertEqual(len(self.plotter.bs_plot_data()['ticks']['label']), 19, "wrong number of tick labels")

if __name__ == '__main__':
    unittest.main()
