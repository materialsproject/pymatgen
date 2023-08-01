from __future__ import annotations

import json
import os
import unittest

from pymatgen.phonon.bandstructure import PhononBandStructureSymmLine
from pymatgen.phonon.dos import CompletePhononDos
from pymatgen.phonon.plotter import PhononBSPlotter, PhononDosPlotter, ThermoPlotter
from pymatgen.util.testing import PymatgenTest


class PhononDosPlotterTest(unittest.TestCase):
    def setUp(self):
        with open(os.path.join(PymatgenTest.TEST_FILES_DIR, "NaCl_complete_ph_dos.json")) as f:
            self.dos = CompletePhononDos.from_dict(json.load(f))
            self.plotter = PhononDosPlotter(sigma=0.2, stack=True)
            self.plotter_nostack = PhononDosPlotter(sigma=0.2, stack=False)

    def test_add_dos_dict(self):
        d = self.plotter.get_dos_dict()
        assert len(d) == 0
        self.plotter.add_dos_dict(self.dos.get_element_dos(), key_sort_func=lambda x: x.X)
        d = self.plotter.get_dos_dict()
        assert len(d) == 2

    def test_get_dos_dict(self):
        self.plotter.add_dos_dict(self.dos.get_element_dos(), key_sort_func=lambda x: x.X)
        d = self.plotter.get_dos_dict()
        for el in ["Na", "Cl"]:
            assert el in d

    def test_plot(self):
        # Disabling latex for testing.
        from matplotlib import axes, rc

        rc("text", usetex=False)
        self.plotter.add_dos("Total", self.dos)
        self.plotter.get_plot(units="mev")
        self.plotter_nostack.add_dos("Total", self.dos)
        plt = self.plotter_nostack.get_plot(units="mev")

        ax = plt.gca()
        assert isinstance(ax, axes.Axes)
        assert ax.get_ylabel() == "$\\mathrm{Density\\ of\\ states}$"
        assert ax.get_xlabel() == "$\\mathrm{Frequencies\\ (meV)}$"


class PhononBSPlotterTest(unittest.TestCase):
    def setUp(self):
        with open(os.path.join(PymatgenTest.TEST_FILES_DIR, "NaCl_phonon_bandstructure.json")) as f:
            d = json.loads(f.read())
            self.bs = PhononBandStructureSymmLine.from_dict(d)
            self.plotter = PhononBSPlotter(self.bs)
        with open(os.path.join(PymatgenTest.TEST_FILES_DIR, "SrTiO3_phonon_bandstructure.json")) as f:
            d = json.loads(f.read())
            self.bs_sto = PhononBandStructureSymmLine.from_dict(d)
            self.plotter_sto = PhononBSPlotter(self.bs_sto)

    def test_bs_plot_data(self):
        assert len(self.plotter.bs_plot_data()["distances"][0]) == 51, "wrong number of distances in the first branch"
        assert len(self.plotter.bs_plot_data()["distances"]) == 4, "wrong number of branches"
        assert sum(len(e) for e in self.plotter.bs_plot_data()["distances"]) == 204, "wrong number of distances"
        assert self.plotter.bs_plot_data()["ticks"]["label"][4] == "Y", "wrong tick label"
        assert len(self.plotter.bs_plot_data()["ticks"]["label"]) == 8, "wrong number of tick labels"

    def test_plot(self):
        # Disabling latex for testing.
        from matplotlib import rc

        rc("text", usetex=False)
        self.plotter.get_plot(units="mev")

    def test_proj_plot(self):
        # Disabling latex for testing.
        from matplotlib import rc

        rc("text", usetex=False)
        self.plotter.get_proj_plot(units="mev")
        self.plotter.get_proj_plot(units="mev", ylim=(15, 30), rgb_labels=("NA", "CL"))
        self.plotter.get_proj_plot(units="mev", site_comb=[[0], [1]])
        self.plotter.get_proj_plot(units="mev", site_comb=[[0], [1]])

        self.plotter_sto.get_proj_plot()
        self.plotter_sto.get_proj_plot(ylim=(-2.5, 5), site_comb=[[0], [1], [2, 3, 4]])
        self.plotter_sto.get_proj_plot(site_comb=[[0], [1], [2, 3, 4]], rgb_labels=("SR", "TI", "O"))
        self.plotter_sto.get_proj_plot(site_comb=[[0], [1], [2], [3, 4]])

    def test_plot_compare(self):
        # Disabling latex for testing.
        from matplotlib import rc

        rc("text", usetex=False)
        self.plotter.plot_compare(self.plotter, units="mev")


class ThermoPlotterTest(unittest.TestCase):
    def setUp(self):
        with open(os.path.join(PymatgenTest.TEST_FILES_DIR, "NaCl_complete_ph_dos.json")) as f:
            self.dos = CompletePhononDos.from_dict(json.load(f))
            self.plotter = ThermoPlotter(self.dos, self.dos.structure)

    def test_plot_functions(self):
        # Disabling latex for testing.
        from matplotlib import rc

        rc("text", usetex=False)
        self.plotter.plot_cv(5, 100, 5, show=False)
        self.plotter.plot_entropy(5, 100, 5, show=False)
        self.plotter.plot_internal_energy(5, 100, 5, show=False)
        self.plotter.plot_helmholtz_free_energy(5, 100, 5, show=False)
        self.plotter.plot_thermodynamic_properties(5, 100, 5, show=False, fig_close=True)


# Gruneisen plotter is already tested in test_gruneisen


if __name__ == "__main__":
    unittest.main()
