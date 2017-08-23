# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals
import unittest
import os
import json

from io import open

import matplotlib
matplotlib.use("pdf")  # Use non-graphical display backend during test.

from pymatgen.electronic_structure.dos import CompleteDos
from pymatgen.electronic_structure.plotter import DosPlotter, BSPlotter, \
    plot_ellipsoid, fold_point, plot_brillouin_zone, BSPlotterProjected, \
    BSDOSPlotter
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from pymatgen.core.structure import Structure
from pymatgen.io.vasp import Vasprun

"""
Created on May 1, 2012
"""

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "May 1, 2012"

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


class BSPlotterProjectedTest(unittest.TestCase):
    def setUp(self):
        with open(os.path.join(test_dir, "Cu2O_361_bandstructure.json"),
                  "r", encoding='utf-8') as f:
            d = json.load(f)
            self.bs = BandStructureSymmLine.from_dict(d)
            self.plotter = BSPlotterProjected(self.bs)

    # Minimal baseline testing for get_plot. not a true test. Just checks that
    # it can actually execute.
    def test_methods(self):
        self.plotter.get_elt_projected_plots()
        self.plotter.get_elt_projected_plots_color()
        self.plotter.get_projected_plots_dots({'Cu': ['d', 's'], 'O': ['p']})
        # self.plotter.get_projected_plots_dots_patom_pmorb(
        #     {'Cu': ['dxy', 's', 'px'], 'O': ['px', 'py', 'pz']},
        #     {'Cu': [3, 5], 'O': [1]}
        # )


class BSDOSPlotterTest(unittest.TestCase):
    # Minimal baseline testing for get_plot. not a true test. Just checks that
    # it can actually execute.
    def test_methods(self):
        v = Vasprun(os.path.join(test_dir, "vasprun_Si_bands.xml"))
        p = BSDOSPlotter()
        plt = p.get_plot(v.get_band_structure(
            kpoints_filename=os.path.join(test_dir, "KPOINTS_Si_bands")))
        plt = p.get_plot(v.get_band_structure(
            kpoints_filename=os.path.join(test_dir, "KPOINTS_Si_bands")),
            v.complete_dos)


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

    def test_bz_plot(self):
        fig, ax = plot_ellipsoid(self.hessian, self.center,
                                 lattice=self.rec_latt)
        plot_brillouin_zone(self.rec_latt, lines=self.kpath, labels=self.labels,
                            kpoints=self.points, ax=ax, show=False)

    def test_fold_point(self):
        self.assertTrue(
            scipy.allclose(fold_point([0., -0.5, 0.5], lattice=self.rec_latt),
                           self.rec_latt.get_cartesian_coords([0., 0.5, 0.5])))
        self.assertTrue(
            scipy.allclose(fold_point([0.1, -0.6, 0.2], lattice=self.rec_latt),
                           self.rec_latt.get_cartesian_coords([0.1, 0.4, 0.2])))


if __name__ == "__main__":
    unittest.main()
