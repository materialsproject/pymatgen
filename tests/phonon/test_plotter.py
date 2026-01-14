from __future__ import annotations

import matplotlib.pyplot as plt
import orjson
import pytest
from numpy.testing import assert_allclose

from pymatgen.phonon import CompletePhononDos, PhononBandStructureSymmLine
from pymatgen.phonon.plotter import PhononBSPlotter, PhononDosPlotter, ThermoPlotter
from pymatgen.util.testing import TEST_FILES_DIR

TEST_DIR = f"{TEST_FILES_DIR}/phonon/dos"

plt.rc("text", usetex=False)  # Disabling latex for testing


class TestPhononDosPlotter:
    def setup_method(self):
        with open(f"{TEST_DIR}/NaCl_complete_ph_dos.json", "rb") as file:
            self.dos = CompletePhononDos.from_dict(orjson.loads(file.read()))
        self.plotter = PhononDosPlotter(sigma=0.2, stack=True)
        self.plotter_no_stack = PhononDosPlotter(sigma=0.2, stack=False)
        self.plotter_no_sigma = PhononDosPlotter(sigma=None, stack=False)

    def test_add_dos_dict(self):
        dct = self.plotter.get_dos_dict()
        assert len(dct) == 0
        self.plotter.add_dos_dict(self.dos.get_element_dos(), key_sort_func=lambda x: x.X)
        dct = self.plotter.get_dos_dict()
        assert len(dct) == 2

    def test_get_dos_dict(self):
        self.plotter.add_dos_dict(self.dos.get_element_dos(), key_sort_func=lambda x: x.X)
        dct = self.plotter.get_dos_dict()
        assert {*dct} >= {"Na", "Cl"}

    def test_plot(self):
        self.plotter.add_dos("Total", self.dos)
        self.plotter.get_plot(units="mev")
        self.plotter_no_stack.add_dos("Total", self.dos)
        ax = self.plotter_no_stack.get_plot(units="mev")
        assert isinstance(ax, plt.Axes)
        assert ax.get_ylabel() == "$\\mathrm{Density\\ of\\ states}$"
        assert ax.get_xlabel() == "$\\mathrm{Frequencies\\ (meV)}$"
        self.plotter_no_sigma.add_dos("Total", self.dos)
        ax2 = self.plotter_no_sigma.get_plot(units="mev")
        assert_allclose(ax2.get_ylim(), (min(self.dos.densities), max(self.dos.densities)))
        ax3 = self.plotter_no_sigma.get_plot(units="mev", invert_axes=True)
        assert ax3.get_ylabel() == "$\\mathrm{Frequencies\\ (meV)}$"
        assert ax3.get_xlabel() == "$\\mathrm{Density\\ of\\ states}$"
        assert_allclose(ax3.get_xlim(), (min(self.dos.densities), max(self.dos.densities)))
        assert ax3.get_ylim() == ax.get_xlim()


class TestPhononBSPlotter:
    def setup_method(self):
        with open(
            f"{TEST_FILES_DIR}/electronic_structure/bandstructure/NaCl_phonon_bandstructure.json", encoding="utf-8"
        ) as file:
            dct = orjson.loads(file.read())
        self.bs = PhononBandStructureSymmLine.from_dict(dct)
        self.plotter = PhononBSPlotter(self.bs, label="NaCl")
        with open(
            f"{TEST_FILES_DIR}/electronic_structure/bandstructure/SrTiO3_phonon_bandstructure.json", encoding="utf-8"
        ) as file:
            dct = orjson.loads(file.read())
        self.bs_sto = PhononBandStructureSymmLine.from_dict(dct)
        self.plotter_sto = PhononBSPlotter(self.bs_sto)

    def test_bs_plot_data(self):
        assert len(self.plotter.bs_plot_data()["distances"][0]) == 51, "wrong number of distances in the first branch"
        assert len(self.plotter.bs_plot_data()["distances"]) == 4, "wrong number of branches"
        assert sum(len(dist) for dist in self.plotter.bs_plot_data()["distances"]) == 204, "wrong number of distances"
        assert self.plotter.bs_plot_data()["ticks"]["label"][4] == "Y", "wrong tick label"
        assert len(self.plotter.bs_plot_data()["ticks"]["label"]) == 8, "wrong number of tick labels"

    def test_plot(self):
        ax = self.plotter.get_plot(units="mev")
        assert isinstance(ax, plt.Axes)
        assert ax.get_ylabel() == "$\\mathrm{Frequencies\\ (meV)}$"
        assert ax.get_xlabel() == "$\\mathrm{Wave\\ Vector}$"

    def test_proj_plot(self):
        self.plotter.get_proj_plot(units="mev")
        self.plotter.get_proj_plot(units="mev", ylim=(15, 30), rgb_labels=("NA", "CL"))
        self.plotter.get_proj_plot(units="mev", site_comb=[[0], [1]])
        self.plotter.get_proj_plot(units="mev", site_comb=[[0], [1]])

        self.plotter_sto.get_proj_plot()
        self.plotter_sto.get_proj_plot(ylim=(-2.5, 5), site_comb=[[0], [1], [2, 3, 4]])
        self.plotter_sto.get_proj_plot(site_comb=[[0], [1], [2, 3, 4]], rgb_labels=("SR", "TI", "O"))
        self.plotter_sto.get_proj_plot(site_comb=[[0], [1], [2], [3, 4]])

    def test_plot_compare(self):
        labels = ("NaCl", "NaCl 2")
        ax = self.plotter.plot_compare({labels[1]: self.plotter}, units="mev")
        assert isinstance(ax, plt.Axes)
        assert ax.get_ylabel() == "$\\mathrm{Frequencies\\ (meV)}$"
        assert ax.get_xlabel() == "$\\mathrm{Wave\\ Vector}$"
        assert ax.get_title() == ""
        assert [itm.get_text() for itm in ax.get_legend().get_texts()] == list(labels)
        ax = self.plotter.plot_compare(self.plotter, units="mev")
        assert [itm.get_text() for itm in ax.get_legend().get_texts()] == [
            "NaCl",
            "NaCl",
        ]
        labels = ("NaCl", "NaCl 2", "NaCl 3")
        ax = self.plotter.plot_compare({labels[1]: self.plotter, labels[2]: self.plotter}, units="mev")
        assert [itm.get_text() for itm in ax.get_legend().get_texts()] == list(labels)
        colors = tuple(itm.get_color() for itm in ax.get_legend().get_lines())
        assert colors == ("blue", "red", "green")
        with pytest.raises(ValueError, match="The two band structures are not compatible."):
            self.plotter.plot_compare(self.plotter_sto)
        ax = self.plotter.plot_compare(self.plotter_sto, on_incompatible="ignore")
        assert ax is None


class TestThermoPlotter:
    def setup_method(self):
        with open(f"{TEST_DIR}/NaCl_complete_ph_dos.json", "rb") as file:
            self.dos = CompletePhononDos.from_dict(orjson.loads(file.read()))
        self.plotter = ThermoPlotter(self.dos, self.dos.structure)

    @pytest.mark.filterwarnings("ignore:More than 20 figures have been opened")
    def test_plot_functions(self):
        fig = self.plotter.plot_cv(5, 100, 5, show=False)
        assert isinstance(fig, plt.Figure)
        fig = self.plotter.plot_entropy(5, 100, 5, show=False)
        assert isinstance(fig, plt.Figure)
        fig = self.plotter.plot_internal_energy(5, 100, 5, show=False)
        assert isinstance(fig, plt.Figure)
        fig = self.plotter.plot_helmholtz_free_energy(5, 100, 5, show=False)
        assert isinstance(fig, plt.Figure)
        fig = self.plotter.plot_thermodynamic_properties(5, 100, 5, show=False, fig_close=True)
        assert isinstance(fig, plt.Figure)


# Gruneisen plotter is already tested in test_gruneisen
