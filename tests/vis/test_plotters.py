from __future__ import annotations

import json
import os

import matplotlib.pyplot as plt
import numpy as np
from monty.json import MontyDecoder

from pymatgen.analysis.xas.spectrum import XAS
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest
from pymatgen.vis.plotters import SpectrumPlotter

with open(f"{TEST_FILES_DIR}/analysis/spectrum_test/LiCoO2_k_xanes.json") as file:
    spect_data_dict = json.load(file, cls=MontyDecoder)


class TestSpectrumPlotter(PymatgenTest):
    def setUp(self):
        self.xanes = XAS.from_dict(spect_data_dict)

    def test_get_plot(self):
        self.plotter = SpectrumPlotter(yshift=0.2)
        self.plotter.add_spectrum("LiCoO2", self.xanes)
        xanes = self.xanes.copy()
        xanes.y += np.random.default_rng().standard_normal(len(xanes.y)) * 0.005
        self.plotter.add_spectrum("LiCoO2 + noise", xanes)
        self.plotter.add_spectrum("LiCoO2 - replot", xanes, "k")
        ax = self.plotter.get_plot()
        assert len(ax.lines) == 3
        img_path = f"{self.tmp_path}/spectrum_plotter_test.png"
        self.plotter.save_plot(img_path)
        assert os.path.isfile(img_path)

    def test_get_stacked_plot(self):
        self.plotter = SpectrumPlotter(yshift=0.2, stack=True)
        self.plotter.add_spectrum("LiCoO2", self.xanes, "b")
        xanes = self.xanes.copy()
        xanes.y += np.random.default_rng().standard_normal(len(xanes.y)) * 0.005
        self.plotter.add_spectrum("LiCoO2 + noise", xanes, "r")
        ax = self.plotter.get_plot()
        assert isinstance(ax, plt.Axes)
        assert len(ax.lines) == 0

    def test_get_plot_with_add_spectrum(self):
        # create spectra_dict
        spectra_dict = {"LiCoO2": self.xanes}
        xanes = self.xanes.copy()
        xanes.y += np.random.default_rng().standard_normal(len(xanes.y)) * 0.005
        spectra_dict["LiCoO2 + noise"] = spectra_dict["LiCoO2 - replot"] = xanes

        self.plotter = SpectrumPlotter(yshift=0.2)
        self.plotter.add_spectra(spectra_dict)
        ax = self.plotter.get_plot()
        assert isinstance(ax, plt.Axes)
        assert len(ax.lines) == 3
        img_path = f"{self.tmp_path}/spectrum_plotter_test2.pdf"
        self.plotter.save_plot(img_path)
        assert os.path.isfile(img_path)
