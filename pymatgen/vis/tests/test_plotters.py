# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import os
import unittest
import json

from monty.json import MontyDecoder
import numpy as np
import matplotlib
matplotlib.use("pdf")

from pymatgen.util.testing import PymatgenTest
from pymatgen.analysis.xas.spectrum import XANES
from pymatgen.vis.plotters import SpectrumPlotter

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files/spectrum_test")

with open(os.path.join(test_dir, 'Pd2O.json')) as fp:
    spect_data_dict = json.load(fp, cls=MontyDecoder)


class SpectrumPlotterTest(PymatgenTest):
    def setUp(self):
        self.xanes = XANES.from_dict(spect_data_dict)

    def test_get_plot(self):
        self.plotter = SpectrumPlotter(yshift=0.2)
        self.plotter.add_spectrum("Pd2O", self.xanes)
        xanes = self.xanes.copy()
        xanes.y += np.random.randn(len(xanes.y)) * 0.005
        self.plotter.add_spectrum("Pd2O + noise", xanes)
        self.plotter.add_spectrum("Pd2O - replot", xanes, "k")
        plt = self.plotter.get_plot()
        self.plotter.save_plot("spectrum_plotter_test.eps")
        os.remove("spectrum_plotter_test.eps")

    def test_get_stacked_plot(self):
        self.plotter = SpectrumPlotter(yshift=0.2, stack=True)
        self.plotter.add_spectrum("Pd2O", self.xanes, "b")
        xanes = self.xanes.copy()
        xanes.y += np.random.randn(len(xanes.y)) * 0.005
        self.plotter.add_spectrum("Pd2O + noise", xanes, "r")
        plt = self.plotter.get_plot()

if __name__ == '__main__':
    unittest.main()