# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
Created on May 1, 2012
"""


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "May 1, 2012"

import unittest
import os
import json

from io import open

from pymatgen.electronic_structure.dos import CompleteDos
from pymatgen.electronic_structure.plotter import DosPlotter, BSPlotter, \
    _qvertex_target
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')

import scipy


class DosPlotterTest(unittest.TestCase):

    def setUp(self):
        with open(os.path.join(test_dir, "complete_dos.json"), "r",
                  encoding='utf-8') as f:
            self.dos = CompleteDos.from_dict(json.load(f))
            self.plotter = DosPlotter(sigma=0.2, stack=True)

    def test_add_dos_dict(self):
        d = self.plotter.get_dos_dict()
        self.assertEqual(len(d), 0)
        self.plotter.add_dos_dict(self.dos.get_element_dos(),
                                  key_sort_func=lambda x: x.X)
        d = self.plotter.get_dos_dict()
        self.assertEqual(len(d), 4)

    def test_get_dos_dict(self):
        self.plotter.add_dos_dict(self.dos.get_element_dos(),
                                  key_sort_func=lambda x: x.X)
        d = self.plotter.get_dos_dict()
        for el in ["Li", "Fe", "P", "O"]:
            self.assertIn(el, d)


class BSPlotterTest(unittest.TestCase):

    def setUp(self):
        with open(os.path.join(test_dir, "CaO_2605_bandstructure.json"),
                  "r", encoding='utf-8') as f:
            d = json.loads(f.read())
            self.bs = BandStructureSymmLine.from_dict(d)
            self.plotter = BSPlotter(self.bs)

    def test_bs_plot_data(self):
        self.assertEqual(len(self.plotter.bs_plot_data()['distances'][0]), 16,
                         "wrong number of distances in the first branch")
        self.assertEqual(len(self.plotter.bs_plot_data()['distances']), 10,
                         "wrong number of branches")
        self.assertEqual(
            sum([len(e) for e in self.plotter.bs_plot_data()['distances']]),
            160, "wrong number of distances")
        self.assertEqual(self.plotter.bs_plot_data()['ticks']['label'][5], "K",
                         "wrong tick label")
        self.assertEqual(len(self.plotter.bs_plot_data()['ticks']['label']),
                         19, "wrong number of tick labels")

    def test_qvertex_target(self):
        results = _qvertex_target([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
                                   [1.0, 1.0, 0.0], [0.0, 1.0, 0.0],
                                   [0.0, 0.0, 1.0], [1.0, 0.0, 1.0],
                                   [1.0, 1.0, 1.0], [0.0, 1.0, 1.0],
                                   [0.5, 0.5, 0.5]], 8)
        self.assertEqual(len(results), 6)
        self.assertEqual(results[3][1], 0.5)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
