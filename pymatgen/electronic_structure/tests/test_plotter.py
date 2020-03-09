# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import unittest
import os
import json
import warnings
from io import open

import scipy

from monty.os.path import which
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.cohp import CompleteCohp
from pymatgen.electronic_structure.dos import CompleteDos
from pymatgen.electronic_structure.boltztrap import BoltztrapAnalyzer
from pymatgen.electronic_structure.plotter import DosPlotter, BSPlotter, \
    plot_ellipsoid, fold_point, plot_brillouin_zone, BSPlotterProjected, \
    BSDOSPlotter, CohpPlotter, BoltztrapPlotter
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from pymatgen.core.structure import Structure
from pymatgen.io.vasp import Vasprun
from pymatgen.util.testing import PymatgenTest

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class DosPlotterTest(unittest.TestCase):
    def setUp(self):
        with open(os.path.join(test_dir, "complete_dos.json"), "r",
                  encoding='utf-8') as f:
            self.dos = CompleteDos.from_dict(json.load(f))
            self.plotter = DosPlotter(sigma=0.2, stack=True)
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

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

    # Minimal baseline testing for get_plot. not a true test. Just checks that
    # it can actually execute.
    def test_get_plot(self):
        # Disabling latex is needed for this test to work.
        from matplotlib import rc
        rc('text', usetex=False)
        self.plotter.add_dos_dict(self.dos.get_element_dos(),
                                  key_sort_func=lambda x: x.X)
        plt = self.plotter.get_plot()
        self.plotter.save_plot("dosplot.png")
        self.assertTrue(os.path.isfile("dosplot.png"))
        os.remove("dosplot.png")
        plt.close("all")


class BSPlotterTest(unittest.TestCase):
    def setUp(self):
        with open(os.path.join(test_dir, "CaO_2605_bandstructure.json"),
                  "r", encoding='utf-8') as f:
            d = json.loads(f.read())
            self.bs = BandStructureSymmLine.from_dict(d)
            self.plotter = BSPlotter(self.bs)
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

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

    # Minimal baseline testing for get_plot. not a true test. Just checks that
    # it can actually execute.
    def test_get_plot(self):
        # zero_to_efermi = True, ylim = None, smooth = False,
        # vbm_cbm_marker = False, smooth_tol = None

        # Disabling latex is needed for this test to work.
        from matplotlib import rc
        rc('text', usetex=False)

        plt = self.plotter.get_plot()
        plt = self.plotter.get_plot(smooth=True)
        plt = self.plotter.get_plot(vbm_cbm_marker=True)
        self.plotter.save_plot("bsplot.png")
        self.assertTrue(os.path.isfile("bsplot.png"))
        os.remove("bsplot.png")
        plt.close("all")


class BSPlotterProjectedTest(unittest.TestCase):
    def setUp(self):
        with open(os.path.join(test_dir, "Cu2O_361_bandstructure.json"),
                  "r", encoding='utf-8') as f:
            d = json.load(f)
            self.bs = BandStructureSymmLine.from_dict(d)
            self.plotter = BSPlotterProjected(self.bs)
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    # Minimal baseline testing for get_plot. not a true test. Just checks that
    # it can actually execute.
    def test_methods(self):
        pass
        # self.plotter.get_elt_projected_plots().close()
        # self.plotter.get_elt_projected_plots_color().close()
        # self.plotter.get_projected_plots_dots({'Cu': ['d', 's'], 'O': ['p']}).close()
        # self.plotter.get_projected_plots_dots_patom_pmorb(
        #     {'Cu': ['dxy', 's', 'px'], 'O': ['px', 'py', 'pz']},
        #     {'Cu': [3, 5], 'O': [1]}
        # ).close()


class BSDOSPlotterTest(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    # Minimal baseline testing for get_plot. not a true test. Just checks that
    # it can actually execute.
    def test_methods(self):
        v = Vasprun(os.path.join(test_dir, "vasprun_Si_bands.xml"))
        p = BSDOSPlotter()
        plt = p.get_plot(v.get_band_structure(
            kpoints_filename=os.path.join(test_dir, "KPOINTS_Si_bands")))
        plt.close()
        plt = p.get_plot(v.get_band_structure(
            kpoints_filename=os.path.join(test_dir, "KPOINTS_Si_bands")),
            v.complete_dos)
        plt.close("all")


class PlotBZTest(unittest.TestCase):
    def setUp(self):
        self.rec_latt = Structure.from_file(
            os.path.join(test_dir, "Si.cssr")).lattice.reciprocal_lattice
        self.kpath = [[[0., 0., 0.], [0.5, 0., 0.5], [0.5, 0.25, 0.75],
                       [0.375, 0.375, 0.75]]]
        self.labels = {'\\Gamma': [0., 0., 0.], 'K': [0.375, 0.375, 0.75],
                       u'L': [0.5, 0.5, 0.5],
                       'U': [0.625, 0.25, 0.625], 'W': [0.5, 0.25, 0.75],
                       'X': [0.5, 0., 0.5]}
        self.hessian = [[17.64757034, 3.90159625, -4.77845607],
                        [3.90159625, 14.88874142, 6.75776076],
                        [-4.77845607, 6.75776076, 12.12987493]]
        self.center = [0.41, 0., 0.41]
        self.points = [[0., 0., 0.], [0.5, 0.5, 0.5]]
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_bz_plot(self):
        fig, ax = plot_ellipsoid(self.hessian, self.center,
                                 lattice=self.rec_latt)
        fig = plot_brillouin_zone(self.rec_latt, lines=self.kpath, labels=self.labels,
                                  kpoints=self.points, ax=ax, show=False)

    def test_fold_point(self):
        self.assertTrue(
            scipy.allclose(fold_point([0., -0.5, 0.5], lattice=self.rec_latt),
                           self.rec_latt.get_cartesian_coords([0., 0.5, 0.5])))
        self.assertTrue(
            scipy.allclose(fold_point([0.1, -0.6, 0.2], lattice=self.rec_latt),
                           self.rec_latt.get_cartesian_coords([0.1, 0.4, 0.2])))


x_trans = which("x_trans")


@unittest.skipIf(not x_trans, "No x_trans.")
class BoltztrapPlotterTest(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_plots(self):
        bz = BoltztrapAnalyzer.from_files(
            os.path.join(test_dir, "boltztrap/transp/"))
        plotter = BoltztrapPlotter(bz)
        plotter.plot_seebeck_eff_mass_mu().close()
        plotter.plot_complexity_factor_mu().close()
        plotter.plot_conductivity_mu().close()
        plotter.plot_power_factor_mu().close()
        plotter.plot_zt_mu().close()
        plotter.plot_dos().close()

        # TODO: These tests fail. Whoever is responsible for the
        # BoltztrapPlotter needs to fix these. The fact that there are not tests
        # for the plotter is atrocious. I will reject all future additions to
        # the plotter until these are fixed.
        # plotter.plot_seebeck_temp()
        # plotter.plot_seebeck_dop()
        # plotter.plot_carriers()
        # plotter.plot_conductivity_dop()
        # plotter.plot_conductivity_temp()
        # plotter.plot_power_factor_dop()
        # plotter.plot_power_factor_temp()
        # plotter.plot_eff_mass_dop()
        # plotter.plot_zt_dop()
        # plotter.plot_zt_temp()


class CohpPlotterTest(PymatgenTest):
    def setUp(self):
        path = os.path.join(test_dir, "cohp", "complete_cohp_lobster.json")
        with open(os.path.join(path), "r") as f:
            self.cohp = CompleteCohp.from_dict(json.load(f))
        path = os.path.join(test_dir, "cohp", "complete_coop_lobster.json")
        with open(os.path.join(path), "r") as f:
            self.coop = CompleteCohp.from_dict(json.load(f))
        self.cohp_plot = CohpPlotter(zero_at_efermi=False)
        self.coop_plot = CohpPlotter(are_coops=True)
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_attributes(self):
        self.assertFalse(self.cohp_plot.are_coops)
        self.assertTrue(self.coop_plot.are_coops)
        self.assertFalse(self.cohp_plot.zero_at_efermi)
        self.assertTrue(self.coop_plot.zero_at_efermi)
        self.cohp_plot.add_cohp_dict(self.cohp.all_cohps)
        cohp_energies = self.cohp_plot._cohps["1"]["energies"]
        self.assertEqual(len(cohp_energies), 301)
        self.assertAlmostEqual(cohp_energies[0], -0.27768)
        self.assertAlmostEqual(cohp_energies[-1], 14.77248)
        self.coop_plot.add_cohp_dict(self.coop.all_cohps)
        coop_energies = self.coop_plot._cohps["10"]["energies"]
        self.assertEqual(len(coop_energies), 241)
        self.assertAlmostEqual(coop_energies[0], -6.02510)
        self.assertAlmostEqual(coop_energies[-1], 6.02510)

    def test_add_cohp_dict(self):
        # Sorts the populations by z-coordinates of the sites
        def sortkeys(sites):
            return sites[0].z, sites[1].z

        sorted_keys = ["3", "4", "7", "8",
                       "9", "10", "11", "6",
                       "5", "2", "1"]

        d_coop = self.coop_plot.get_cohp_dict()
        self.assertEqual(len(d_coop), 0)
        bonds = self.coop.bonds
        self.coop_plot.add_cohp_dict(self.coop.all_cohps,
                                     key_sort_func=lambda x:
                                     sortkeys(bonds[x]["sites"]))
        d_coop = self.coop_plot.get_cohp_dict()
        self.assertEqual(len(d_coop), 11)
        self.assertEqual(list(self.coop_plot._cohps.keys()), sorted_keys)

    def test_get_cohp_dict(self):
        self.cohp_plot.add_cohp_dict(self.cohp.all_cohps)
        d_cohp = self.cohp_plot.get_cohp_dict()
        for bond in ["1", "2"]:
            self.assertIn(bond, d_cohp)

    def test_get_plot(self):
        self.cohp_plot.add_cohp_dict(self.cohp.all_cohps)
        plt_cohp = self.cohp_plot.get_plot()
        ax_cohp = plt_cohp.gca()
        self.assertEqual(ax_cohp.get_xlabel(), "-COHP")
        self.assertEqual(ax_cohp.get_ylabel(), "$E$ (eV)")
        legend_labels = ax_cohp.get_legend_handles_labels()[1]
        self.assertEqual(len(self.cohp_plot._cohps), len(legend_labels))
        self.assertEqual(ax_cohp.lines[0].get_linestyle(), "-")
        self.assertEqual(ax_cohp.lines[1].get_linestyle(), "--")
        for label in legend_labels:
            self.assertIn(label, self.cohp_plot._cohps)
        linesindex = legend_labels.index("1")
        linestyles = {Spin.up: '-', Spin.down: '--'}
        cohp_fe_fe = self.cohp.all_cohps["1"]
        for s, spin in enumerate([Spin.up, Spin.down]):
            lines = ax_cohp.lines[2 * linesindex + s]
            self.assertArrayAlmostEqual(lines.get_xdata(),
                                        -cohp_fe_fe.cohp[spin])
            self.assertArrayAlmostEqual(lines.get_ydata(), self.cohp.energies)
            self.assertEqual(lines.get_linestyle(), linestyles[spin])
        plt_cohp.close()

        plt_cohp = self.cohp_plot.get_plot(invert_axes=False,
                                           plot_negative=False)
        ax_cohp = plt_cohp.gca()
        self.assertEqual(ax_cohp.get_xlabel(), "$E$ (eV)")
        self.assertEqual(ax_cohp.get_ylabel(), "COHP")
        for s, spin in enumerate([Spin.up, Spin.down]):
            lines = ax_cohp.lines[2 * linesindex + s]
            self.assertArrayAlmostEqual(lines.get_xdata(), self.cohp.energies)
            self.assertArrayAlmostEqual(lines.get_ydata(),
                                        cohp_fe_fe.cohp[spin])
        plt_cohp.close()

        plt_cohp = self.cohp_plot.get_plot(integrated=True)
        ax_cohp = plt_cohp.gca()
        self.assertEqual(ax_cohp.get_xlabel(), "-ICOHP (eV)")
        for s, spin in enumerate([Spin.up, Spin.down]):
            lines = ax_cohp.lines[2 * linesindex + s]
            self.assertArrayAlmostEqual(lines.get_xdata(),
                                        -cohp_fe_fe.icohp[spin])

        coop_dict = {"Bi5-Bi6": self.coop.all_cohps["10"]}
        self.coop_plot.add_cohp_dict(coop_dict)
        plt_coop = self.coop_plot.get_plot()
        ax_coop = plt_coop.gca()
        self.assertEqual(ax_coop.get_xlabel(), "COOP")
        self.assertEqual(ax_coop.get_ylabel(), "$E - E_f$ (eV)")
        lines_coop = ax_coop.get_lines()[0]
        self.assertArrayAlmostEqual(lines_coop.get_ydata(),
                                    self.coop.energies - self.coop.efermi)
        coop_bi_bi = self.coop.all_cohps["10"].cohp[Spin.up]
        self.assertArrayAlmostEqual(lines_coop.get_xdata(), coop_bi_bi)

        # Cleanup.
        plt_cohp.close()
        plt_coop.close("all")

    def test_save_plot(self):
        self.cohp_plot.add_cohp_dict(self.cohp.all_cohps)
        plt_cohp = self.cohp_plot.get_plot()
        self.cohp_plot.save_plot("cohpplot.png")
        self.assertTrue(os.path.isfile("cohpplot.png"))
        os.remove("cohpplot.png")
        plt_cohp.close("all")


if __name__ == "__main__":
    unittest.main()
