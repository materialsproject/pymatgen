# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import annotations

import json
import os
import unittest
import warnings
from shutil import which

import numpy as np
import scipy
from pytest import approx

from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from pymatgen.electronic_structure.boltztrap import BoltztrapAnalyzer
from pymatgen.electronic_structure.cohp import CompleteCohp
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.dos import CompleteDos
from pymatgen.electronic_structure.plotter import (
    BoltztrapPlotter,
    BSDOSPlotter,
    BSPlotter,
    BSPlotterProjected,
    CohpPlotter,
    DosPlotter,
    fold_point,
    plot_brillouin_zone,
    plot_ellipsoid,
)
from pymatgen.io.vasp import Vasprun
from pymatgen.util.testing import PymatgenTest


class DosPlotterTest(unittest.TestCase):
    def setUp(self):
        with open(os.path.join(PymatgenTest.TEST_FILES_DIR, "complete_dos.json")) as f:
            self.dos = CompleteDos.from_dict(json.load(f))
            self.plotter = DosPlotter(sigma=0.2, stack=True)
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_add_dos_dict(self):
        d = self.plotter.get_dos_dict()
        assert len(d) == 0
        self.plotter.add_dos_dict(self.dos.get_element_dos(), key_sort_func=lambda x: x.X)
        d = self.plotter.get_dos_dict()
        assert len(d) == 4

    def test_get_dos_dict(self):
        self.plotter.add_dos_dict(self.dos.get_element_dos(), key_sort_func=lambda x: x.X)
        d = self.plotter.get_dos_dict()
        for el in ["Li", "Fe", "P", "O"]:
            assert el in d

    # Minimal baseline testing for get_plot. not a true test. Just checks that
    # it can actually execute.
    def test_get_plot(self):
        # Disabling latex is needed for this test to work.
        from matplotlib import rc

        rc("text", usetex=False)
        self.plotter.add_dos_dict(self.dos.get_element_dos(), key_sort_func=lambda x: x.X)
        plt = self.plotter.get_plot()
        self.plotter.save_plot("dosplot.png")
        assert os.path.isfile("dosplot.png")
        os.remove("dosplot.png")
        plt.close("all")


class BSPlotterTest(unittest.TestCase):
    def setUp(self):
        with open(os.path.join(PymatgenTest.TEST_FILES_DIR, "CaO_2605_bandstructure.json")) as f:
            d = json.loads(f.read())
            self.bs = BandStructureSymmLine.from_dict(d)
            self.plotter = BSPlotter(self.bs)

        assert len(self.plotter._bs) == 1, "wrong number of band objects"

        with open(os.path.join(PymatgenTest.TEST_FILES_DIR, "N2_12103_bandstructure.json")) as f:
            d = json.loads(f.read())
            self.sbs_sc = BandStructureSymmLine.from_dict(d)

        with open(os.path.join(PymatgenTest.TEST_FILES_DIR, "C_48_bandstructure.json")) as f:
            d = json.loads(f.read())
            self.sbs_met = BandStructureSymmLine.from_dict(d)

        self.plotter_multi = BSPlotter([self.sbs_sc, self.sbs_met])
        assert len(self.plotter_multi._bs) == 2, "wrong number of band objects"
        assert self.plotter_multi._nb_bands == [96, 96], "wrong number of bands"
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_add_bs(self):
        self.plotter_multi.add_bs(self.sbs_sc)
        assert len(self.plotter_multi._bs) == 3, "wrong number of band objects"
        assert self.plotter_multi._nb_bands == [96, 96, 96], "wrong number of bands"

    def test_get_branch_steps(self):
        steps_idx = BSPlotter._get_branch_steps(self.sbs_sc.branches)
        assert steps_idx == [0, 121, 132, 143], "wrong list of steps idx"

    def test_rescale_distances(self):
        rescaled_distances = self.plotter_multi._rescale_distances(self.sbs_sc, self.sbs_met)
        assert len(rescaled_distances) == len(self.sbs_met.distance), "wrong length of distances list"
        assert rescaled_distances[-1] == 6.5191398067252875, "wrong last distance value"
        assert rescaled_distances[148] == self.sbs_sc.distance[19], "wrong distance at high symm k-point"

    def test_interpolate_bands(self):
        data = self.plotter.bs_plot_data()
        d = data["distances"]
        en = data["energy"]["1"]
        int_distances, int_energies = self.plotter._interpolate_bands(d, en)

        assert len(int_distances) == 10, "wrong length of distances list"
        assert len(int_distances[0]) == 100, "wrong length of distances in a branch"
        assert len(int_energies) == 10, "wrong length of distances list"
        assert int_energies[0].shape == (16, 100), "wrong length of distances list"

    def test_bs_plot_data(self):
        assert len(self.plotter.bs_plot_data()["distances"]) == 10, "wrong number of sequences of branches"
        assert (
            len(self.plotter.bs_plot_data()["distances"][0]) == 16
        ), "wrong number of distances in the first sequence of branches"
        assert sum(len(e) for e in self.plotter.bs_plot_data()["distances"]) == 160, "wrong number of distances"

        length = len(self.plotter.bs_plot_data(split_branches=False)["distances"][0])
        assert length == 144, "wrong number of distances in the first sequence of branches"

        length = len(self.plotter.bs_plot_data(split_branches=False)["distances"])
        assert length == 2, "wrong number of distances in the first sequence of branches"

        assert self.plotter.bs_plot_data()["ticks"]["label"][5] == "K", "wrong tick label"
        assert len(self.plotter.bs_plot_data()["ticks"]["label"]) == 19, "wrong number of tick labels"

    def test_get_ticks(self):
        assert self.plotter.get_ticks()["label"][5] == "K", "wrong tick label"
        assert self.plotter.get_ticks()["distance"][5] == 2.406607625322699, "wrong tick distance"

    # Minimal baseline testing for get_plot. not a true test. Just checks that
    # it can actually execute.
    def test_get_plot(self):
        # zero_to_efermi = True, ylim = None, smooth = False,
        # vbm_cbm_marker = False, smooth_tol = None

        # Disabling latex is needed for this test to work.
        from matplotlib import rc

        rc("text", usetex=False)

        plt = self.plotter.get_plot()
        assert plt.ylim() == (-4.0, 7.6348), "wrong ylim"
        plt = self.plotter.get_plot(smooth=True)
        plt = self.plotter.get_plot(vbm_cbm_marker=True)
        self.plotter.save_plot("bsplot.png")
        assert os.path.isfile("bsplot.png")
        os.remove("bsplot.png")
        plt.close("all")

        # test plotter with 2 bandstructures
        plt = self.plotter_multi.get_plot()
        assert len(plt.gca().get_lines()) == 874, "wrong number of lines"
        assert plt.ylim() == (-10.0, 10.0), "wrong ylim"
        plt = self.plotter_multi.get_plot(zero_to_efermi=False)
        assert plt.ylim() == (-15.2379, 12.67141266), "wrong ylim"
        plt = self.plotter_multi.get_plot(smooth=True)
        self.plotter_multi.save_plot("bsplot.png")
        assert os.path.isfile("bsplot.png")
        os.remove("bsplot.png")
        plt.close("all")


class BSPlotterProjectedTest(unittest.TestCase):
    def setUp(self):
        with open(os.path.join(PymatgenTest.TEST_FILES_DIR, "Cu2O_361_bandstructure.json")) as f:
            d = json.load(f)
            self.bs = BandStructureSymmLine.from_dict(d)
            self.plotter = BSPlotterProjected(self.bs)
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    # Minimal baseline testing for get_plot. not a true test. Just checks that
    # it can actually execute.
    def test_methods(self):
        self.plotter.get_elt_projected_plots().close()
        self.plotter.get_elt_projected_plots_color().close()
        self.plotter.get_projected_plots_dots({"Cu": ["d", "s"], "O": ["p"]}).close()
        self.plotter.get_projected_plots_dots_patom_pmorb(
            {"Cu": ["dxy", "s", "px"], "O": ["px", "py", "pz"]},
            {"Cu": [3, 5], "O": [1]},
        ).close()


class BSDOSPlotterTest(unittest.TestCase):
    def setUp(self):
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    # Minimal baseline testing for get_plot. not a true test. Just checks that
    # it can actually execute.
    def test_methods(self):
        v = Vasprun(os.path.join(PymatgenTest.TEST_FILES_DIR, "vasprun_Si_bands.xml"))
        p = BSDOSPlotter()
        plt = p.get_plot(
            v.get_band_structure(kpoints_filename=os.path.join(PymatgenTest.TEST_FILES_DIR, "KPOINTS_Si_bands"))
        )
        plt.close()
        plt = p.get_plot(
            v.get_band_structure(kpoints_filename=os.path.join(PymatgenTest.TEST_FILES_DIR, "KPOINTS_Si_bands")),
            v.complete_dos,
        )
        plt.close("all")

        with open(os.path.join(PymatgenTest.TEST_FILES_DIR, "SrBa2Sn2O7.json")) as f:
            bandstr_dict = json.load(f)
        # generate random projections
        data_structure = [[[[0 for _ in range(12)] for _ in range(9)] for _ in range(70)] for _ in range(90)]
        bandstr_dict["projections"]["1"] = data_structure
        d = bandstr_dict["projections"]["1"]
        for i in range(len(d)):
            for j in range(len(d[i])):
                for k in range(len(d[i][j])):
                    for m in range(len(d[i][j][k])):
                        d[i][j][k][m] = 0
                        # d[i][j][k][m] = np.random.rand()
                    # generate random number for two atoms
                    a = np.random.randint(0, 7)
                    b = np.random.randint(0, 7)
                    # c = np.random.randint(0,7)
                    d[i][j][k][a] = np.random.rand()
                    d[i][j][k][b] = np.random.rand()
                    # d[i][j][k][c] = np.random.rand()
        bandstr = BandStructureSymmLine.from_dict(bandstr_dict)
        plt = p.get_plot(bandstr)
        plt.show()


class PlotBZTest(unittest.TestCase):
    def setUp(self):
        self.rec_latt = Structure.from_file(
            os.path.join(PymatgenTest.TEST_FILES_DIR, "Si.cssr")
        ).lattice.reciprocal_lattice
        self.kpath = [[[0.0, 0.0, 0.0], [0.5, 0.0, 0.5], [0.5, 0.25, 0.75], [0.375, 0.375, 0.75]]]
        self.labels = {
            "\\Gamma": [0.0, 0.0, 0.0],
            "K": [0.375, 0.375, 0.75],
            "L": [0.5, 0.5, 0.5],
            "U": [0.625, 0.25, 0.625],
            "W": [0.5, 0.25, 0.75],
            "X": [0.5, 0.0, 0.5],
        }
        self.hessian = [
            [17.64757034, 3.90159625, -4.77845607],
            [3.90159625, 14.88874142, 6.75776076],
            [-4.77845607, 6.75776076, 12.12987493],
        ]
        self.center = [0.41, 0.0, 0.41]
        self.points = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_bz_plot(self):
        _, ax = plot_ellipsoid(self.hessian, self.center, lattice=self.rec_latt)
        plot_brillouin_zone(
            self.rec_latt,
            lines=self.kpath,
            labels=self.labels,
            kpoints=self.points,
            ax=ax,
            show=False,
        )

    def test_fold_point(self):
        assert scipy.allclose(
            fold_point([0.0, -0.5, 0.5], lattice=self.rec_latt),
            self.rec_latt.get_cartesian_coords([0.0, 0.5, 0.5]),
        )
        assert scipy.allclose(
            fold_point([0.1, -0.6, 0.2], lattice=self.rec_latt),
            self.rec_latt.get_cartesian_coords([0.1, 0.4, 0.2]),
        )


x_trans = which("x_trans")


@unittest.skipIf(not x_trans, "No x_trans.")
class BoltztrapPlotterTest(unittest.TestCase):
    def setUp(self):
        bz = BoltztrapAnalyzer.from_files(os.path.join(PymatgenTest.TEST_FILES_DIR, "boltztrap/transp/"))
        self.plotter = BoltztrapPlotter(bz)
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_plot_carriers(self):
        plt = self.plotter.plot_carriers()
        assert len(plt.gca().get_lines()) == 7, "wrong number of lines"
        assert plt.gca().get_lines()[0].get_data()[0][0] == -2.0702422655947665, "wrong 0 data in line 0"
        assert plt.gca().get_lines()[0].get_data()[1][0] == 6.525490122298364e22, "wrong 1 data in line 0"
        plt.close()

    def test_plot_complexity_factor_mu(self):
        plt = self.plotter.plot_complexity_factor_mu()
        assert len(plt.gca().get_lines()) == 2, "wrong number of lines"
        assert plt.gca().get_lines()[0].get_data()[0][0] == -2.0702422655947665, "wrong 0 data in line 0"
        assert plt.gca().get_lines()[0].get_data()[1][0] == 0.004708835456903449, "wrong 1 data in line 0"
        plt.close()

    def test_plot_conductivity_dop(self):
        plt = self.plotter.plot_conductivity_dop()
        assert len(plt.gca().get_lines()) == 8, "wrong number of lines"
        assert plt.gca().get_lines()[0].get_data()[0][0] == 1000000000000000.0, "wrong 0 data in line 0"
        assert plt.gca().get_lines()[0].get_data()[1][0] == 0.3801957596666667, "wrong 1 data in line 0"
        plt.close()

    def test_plot_conductivity_mu(self):
        plt = self.plotter.plot_conductivity_mu()
        assert len(plt.gca().get_lines()) == 9, "wrong number of lines"
        assert plt.gca().get_lines()[0].get_data()[0][0] == -2.0702422655947665, "wrong 0 data in line 0"
        assert plt.gca().get_lines()[0].get_data()[1][0] == 1965.1306, "wrong 1 data in line 0"
        plt.close()

    def test_plot_conductivity_temp(self):
        plt = self.plotter.plot_conductivity_temp()
        assert len(plt.gca().get_lines()) == 6, "wrong number of lines"
        assert plt.gca().get_lines()[0].get_data()[0][0] == 100, "wrong 0 data in line 0"
        assert plt.gca().get_lines()[0].get_data()[1][0] == 0.3801957596666667, "wrong 1 data in line 0"
        plt.close()

    def test_plot_dos(self):
        plt = self.plotter.plot_dos()
        assert len(plt.gca().get_lines()) == 3, "wrong number of lines"
        assert plt.gca().get_lines()[0].get_data()[0][0] == -2.4197044934588674, "wrong 0 data in line 0"
        assert plt.gca().get_lines()[0].get_data()[1][0] == 0.0, "wrong 1 data in line 0"
        plt.close()

    def test_plot_eff_mass_dop(self):
        plt = self.plotter.plot_eff_mass_dop()
        assert len(plt.gca().get_lines()) == 8, "wrong number of lines"
        assert plt.gca().get_lines()[0].get_data()[0][0] == 1000000000000000.0, "wrong 0 data in line 0"
        assert plt.gca().get_lines()[0].get_data()[1][0] == 1.4231240011719886, "wrong 1 data in line 0"
        plt.close()

    def test_plot_eff_mass_temp(self):
        plt = self.plotter.plot_eff_mass_temp()
        assert len(plt.gca().get_lines()) == 6, "wrong number of lines"
        assert plt.gca().get_lines()[0].get_data()[0][0] == 100, "wrong 0 data in line 0"
        assert plt.gca().get_lines()[0].get_data()[1][0] == 1.4231240011719886, "wrong 1 data in line 0"
        plt.close()

    def test_plot_hall_carriers(self):
        plt = self.plotter.plot_hall_carriers()
        assert len(plt.gca().get_lines()) == 7, "wrong number of lines"
        assert plt.gca().get_lines()[0].get_data()[0][0] == -2.0702422655947665, "wrong 0 data in line 0"
        assert plt.gca().get_lines()[0].get_data()[1][0] == 9.538187273102463e17, "wrong 1 data in line 0"
        plt.close()

    def test_plot_power_factor_dop(self):
        plt = self.plotter.plot_power_factor_dop()
        assert len(plt.gca().get_lines()) == 8, "wrong number of lines"
        assert plt.gca().get_lines()[0].get_data()[0][0] == 1000000000000000.0, "wrong 0 data in line 0"
        assert plt.gca().get_lines()[0].get_data()[1][0] == 0.40606868935796925, "wrong 1 data in line 0"
        plt.close()

    def test_plot_power_factor_mu(self):
        plt = self.plotter.plot_power_factor_mu()
        assert len(plt.gca().get_lines()) == 9, "wrong number of lines"
        assert plt.gca().get_lines()[0].get_data()[0][0] == -2.0702422655947665, "wrong 0 data in line 0"
        assert plt.gca().get_lines()[0].get_data()[1][0] == 365.5514594136157, "wrong 1 data in line 0"
        plt.close()

    def test_plot_power_factor_temp(self):
        plt = self.plotter.plot_power_factor_temp()
        assert len(plt.gca().get_lines()) == 6, "wrong number of lines"
        assert plt.gca().get_lines()[0].get_data()[0][0] == 100, "wrong 0 data in line 0"
        assert plt.gca().get_lines()[0].get_data()[1][0] == 0.40606868935796925, "wrong 1 data in line 0"
        plt.close()

    def test_plot_seebeck_dop(self):
        plt = self.plotter.plot_seebeck_dop()
        assert len(plt.gca().get_lines()) == 8, "wrong number of lines"
        assert plt.gca().get_lines()[0].get_data()[0][0] == 1000000000000000.0, "wrong 0 data in line 0"
        assert plt.gca().get_lines()[0].get_data()[1][0] == 1050.8197666666667, "wrong 1 data in line 0"
        plt.close()

    def test_plot_seebeck_eff_mass_mu(self):
        plt = self.plotter.plot_seebeck_eff_mass_mu()
        assert len(plt.gca().get_lines()) == 2, "wrong number of lines"
        assert plt.gca().get_lines()[0].get_data()[0][0] == -2.0702422655947665, "wrong 0 data in line 0"
        assert plt.gca().get_lines()[0].get_data()[1][0] == 6412.881888198197, "wrong 1 data in line 0"
        plt.close()

    def test_plot_seebeck_mu(self):
        plt = self.plotter.plot_seebeck_mu()
        assert len(plt.gca().get_lines()) == 9, "wrong number of lines"
        assert plt.gca().get_lines()[0].get_data()[0][0] == -2.0702422655947665, "wrong 0 data in line 0"
        assert plt.gca().get_lines()[0].get_data()[1][0] == -433.11096000000003, "wrong 1 data in line 0"
        plt.close()

    def test_plot_seebeck_temp(self):
        plt = self.plotter.plot_seebeck_temp()
        assert len(plt.gca().get_lines()) == 6, "wrong number of lines"
        assert plt.gca().get_lines()[0].get_data()[0][0] == 100, "wrong 0 data in line 0"
        assert plt.gca().get_lines()[0].get_data()[1][0] == 1050.8197666666667, "wrong 1 data in line 0"
        plt.close()

    def test_plot_zt_dop(self):
        plt = self.plotter.plot_zt_dop()
        assert len(plt.gca().get_lines()) == 8, "wrong number of lines"
        assert plt.gca().get_lines()[0].get_data()[0][0] == 1000000000000000.0, "wrong 0 data in line 0"
        assert plt.gca().get_lines()[0].get_data()[1][0] == 4.060682863129955e-05, "wrong 1 data in line 0"
        plt.close()

    def test_plot_zt_mu(self):
        plt = self.plotter.plot_zt_mu()
        assert len(plt.gca().get_lines()) == 9, "wrong number of lines"
        assert plt.gca().get_lines()[0].get_data()[0][0] == -2.0702422655947665, "wrong 0 data in line 0"
        assert plt.gca().get_lines()[0].get_data()[1][0] == 0.2153839699235254, "wrong 1 data in line 0"
        plt.close()

    def test_plot_zt_temp(self):
        plt = self.plotter.plot_zt_temp()
        assert len(plt.gca().get_lines()) == 6, "wrong number of lines"
        assert plt.gca().get_lines()[0].get_data()[0][0] == 100, "wrong 0 data in line 0"
        assert plt.gca().get_lines()[0].get_data()[1][0] == 4.060682863129955e-05, "wrong 1 data in line 0"
        plt.close()


class CohpPlotterTest(PymatgenTest):
    def setUp(self):
        path = os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "complete_cohp_lobster.json")
        with open(os.path.join(path)) as f:
            self.cohp = CompleteCohp.from_dict(json.load(f))
        path = os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "complete_coop_lobster.json")
        with open(os.path.join(path)) as f:
            self.coop = CompleteCohp.from_dict(json.load(f))
        self.cohp_plot = CohpPlotter(zero_at_efermi=False)
        self.coop_plot = CohpPlotter(are_coops=True)
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_attributes(self):
        assert not self.cohp_plot.are_coops
        assert self.coop_plot.are_coops
        assert not self.cohp_plot.zero_at_efermi
        assert self.coop_plot.zero_at_efermi
        self.cohp_plot.add_cohp_dict(self.cohp.all_cohps)
        cohp_energies = self.cohp_plot._cohps["1"]["energies"]
        assert len(cohp_energies) == 301
        assert cohp_energies[0] == approx(-0.27768)
        assert cohp_energies[-1] == approx(14.77248)
        self.coop_plot.add_cohp_dict(self.coop.all_cohps)
        coop_energies = self.coop_plot._cohps["10"]["energies"]
        assert len(coop_energies) == 241
        assert coop_energies[0] == approx(-6.02510)
        assert coop_energies[-1] == approx(6.02510)

    def test_add_cohp_dict(self):
        # Sorts the populations by z-coordinates of the sites
        def sortkeys(sites):
            return sites[0].z, sites[1].z

        sorted_keys = ["3", "4", "7", "8", "9", "10", "11", "6", "5", "2", "1"]

        d_coop = self.coop_plot.get_cohp_dict()
        assert len(d_coop) == 0
        bonds = self.coop.bonds
        self.coop_plot.add_cohp_dict(self.coop.all_cohps, key_sort_func=lambda x: sortkeys(bonds[x]["sites"]))
        d_coop = self.coop_plot.get_cohp_dict()
        assert len(d_coop) == 11
        assert list(self.coop_plot._cohps) == sorted_keys

    def test_get_cohp_dict(self):
        self.cohp_plot.add_cohp_dict(self.cohp.all_cohps)
        d_cohp = self.cohp_plot.get_cohp_dict()
        for bond in ["1", "2"]:
            assert bond in d_cohp

    def test_get_plot(self):
        self.cohp_plot.add_cohp_dict(self.cohp.all_cohps)
        plt_cohp = self.cohp_plot.get_plot()
        ax_cohp = plt_cohp.gca()
        assert ax_cohp.get_xlabel() == "-COHP"
        assert ax_cohp.get_ylabel() == "$E$ (eV)"
        legend_labels = ax_cohp.get_legend_handles_labels()[1]
        assert len(self.cohp_plot._cohps) == len(legend_labels)
        assert ax_cohp.lines[0].get_linestyle() == "-"
        assert ax_cohp.lines[1].get_linestyle() == "--"
        for label in legend_labels:
            assert label in self.cohp_plot._cohps
        linesindex = legend_labels.index("1")
        linestyles = {Spin.up: "-", Spin.down: "--"}
        cohp_fe_fe = self.cohp.all_cohps["1"]
        for s, spin in enumerate([Spin.up, Spin.down]):
            lines = ax_cohp.lines[2 * linesindex + s]
            self.assertArrayAlmostEqual(lines.get_xdata(), -cohp_fe_fe.cohp[spin])
            self.assertArrayAlmostEqual(lines.get_ydata(), self.cohp.energies)
            assert lines.get_linestyle() == linestyles[spin]
        plt_cohp.close()

        plt_cohp = self.cohp_plot.get_plot(invert_axes=False, plot_negative=False)
        ax_cohp = plt_cohp.gca()
        assert ax_cohp.get_xlabel() == "$E$ (eV)"
        assert ax_cohp.get_ylabel() == "COHP"
        for s, spin in enumerate([Spin.up, Spin.down]):
            lines = ax_cohp.lines[2 * linesindex + s]
            self.assertArrayAlmostEqual(lines.get_xdata(), self.cohp.energies)
            self.assertArrayAlmostEqual(lines.get_ydata(), cohp_fe_fe.cohp[spin])
        plt_cohp.close()

        plt_cohp = self.cohp_plot.get_plot(integrated=True)
        ax_cohp = plt_cohp.gca()
        assert ax_cohp.get_xlabel() == "-ICOHP (eV)"
        for s, spin in enumerate([Spin.up, Spin.down]):
            lines = ax_cohp.lines[2 * linesindex + s]
            self.assertArrayAlmostEqual(lines.get_xdata(), -cohp_fe_fe.icohp[spin])

        coop_dict = {"Bi5-Bi6": self.coop.all_cohps["10"]}
        self.coop_plot.add_cohp_dict(coop_dict)
        plt_coop = self.coop_plot.get_plot()
        ax_coop = plt_coop.gca()
        assert ax_coop.get_xlabel() == "COOP"
        assert ax_coop.get_ylabel() == "$E - E_f$ (eV)"
        lines_coop = ax_coop.get_lines()[0]
        self.assertArrayAlmostEqual(lines_coop.get_ydata(), self.coop.energies - self.coop.efermi)
        coop_bi_bi = self.coop.all_cohps["10"].cohp[Spin.up]
        self.assertArrayAlmostEqual(lines_coop.get_xdata(), coop_bi_bi)

        # Cleanup.
        plt_cohp.close()
        plt_coop.close("all")

    def test_save_plot(self):
        self.cohp_plot.add_cohp_dict(self.cohp.all_cohps)
        plt_cohp = self.cohp_plot.get_plot()
        self.cohp_plot.save_plot("cohpplot.png")
        assert os.path.isfile("cohpplot.png")
        os.remove("cohpplot.png")
        plt_cohp.close("all")


if __name__ == "__main__":
    unittest.main()
