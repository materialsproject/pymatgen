from __future__ import division, unicode_literals

import unittest
import os
import json
import scipy
from io import open

from pymatgen.phonon.dos import CompletePhononDos
from pymatgen.phonon.plotter import PhononDosPlotter, PhononBSPlotter
from pymatgen.phonon.bandstructure import PhononBandStructureSymmLine


test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class PhononDosPlotterTest(unittest.TestCase):

    def setUp(self):
        with open(os.path.join(test_dir, "NaCl_complete_ph_dos.json"), "r") as f:
            self.dos = CompletePhononDos.from_dict(json.load(f))
            self.plotter = PhononDosPlotter(sigma=0.2, stack=True)

    def test_add_dos_dict(self):
        d = self.plotter.get_dos_dict()
        self.assertEqual(len(d), 0)
        self.plotter.add_dos_dict(self.dos.get_element_dos(),
                                  key_sort_func=lambda x: x.X)
        d = self.plotter.get_dos_dict()
        self.assertEqual(len(d), 2)

    def test_get_dos_dict(self):
        self.plotter.add_dos_dict(self.dos.get_element_dos(),
                                  key_sort_func=lambda x: x.X)
        d = self.plotter.get_dos_dict()
        for el in ["Na", "Cl"]:
            self.assertIn(el, d)


class PhononBSPlotterTest(unittest.TestCase):

    def setUp(self):
        with open(os.path.join(test_dir, "NaCl_phonon_bandstructure.json"), "r") as f:
            d = json.loads(f.read())
            self.bs = PhononBandStructureSymmLine.from_dict(d)
            self.plotter = PhononBSPlotter(self.bs)

    def test_bs_plot_data(self):
        self.assertEqual(len(self.plotter.bs_plot_data()['distances'][0]), 51,
                         "wrong number of distances in the first branch")
        self.assertEqual(len(self.plotter.bs_plot_data()['distances']), 4,
                         "wrong number of branches")
        self.assertEqual(
            sum([len(e) for e in self.plotter.bs_plot_data()['distances']]),
            204, "wrong number of distances")
        self.assertEqual(self.plotter.bs_plot_data()['ticks']['label'][4], "Y",
                         "wrong tick label")
        self.assertEqual(len(self.plotter.bs_plot_data()['ticks']['label']),
                         8, "wrong number of tick labels")


if __name__ == "__main__":
    unittest.main()
