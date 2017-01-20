# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest2 as unittest
import os
import numpy as np

from pymatgen import Element
from pymatgen.phasediagram.entries import PDEntryIO, PDEntry
from pymatgen.phasediagram.maker import PhaseDiagram
from pymatgen.phasediagram.plotter import PDPlotter, uniquelines, \
     triangular_coord, tet_coord


module_dir = os.path.dirname(os.path.abspath(__file__))


class PDPlotterTest(unittest.TestCase):

    def setUp(self):
        (elements, entries) = PDEntryIO.from_csv(os.path.join(module_dir, "pdentries_test.csv"))
        self.pd = PhaseDiagram(entries)
        self.plotter = PDPlotter(self.pd, show_unstable=True)
        entrieslio = [e for e in entries
                   if len(e.composition) < 3 and ("Fe" not in e.composition)]

        self.pd_formation = PhaseDiagram(entrieslio)
        self.plotter_formation = PDPlotter(self.pd_formation, show_unstable=True)
        entries.append(PDEntry("C", 0))
        self.pd3d = PhaseDiagram(entries)
        self.plotter3d = PDPlotter(self.pd3d, show_unstable=True)


    def test_pd_plot_data(self):
        (lines, labels, unstable_entries) = self.plotter.pd_plot_data
        self.assertEqual(len(lines), 22)
        self.assertEqual(len(labels), len(self.pd.stable_entries), "Incorrect number of lines generated!")
        self.assertEqual(len(unstable_entries), len(self.pd.all_entries) - len(self.pd.stable_entries), "Incorrect number of lines generated!")
        (lines, labels, unstable_entries) = self.plotter3d.pd_plot_data
        self.assertEqual(len(lines), 33)
        self.assertEqual(len(labels), len(self.pd3d.stable_entries))
        self.assertEqual(len(unstable_entries),
                         len(self.pd3d.all_entries) - len(self.pd3d.stable_entries))
        (lines, labels, unstable_entries) = self.plotter_formation.pd_plot_data
        self.assertEqual(len(lines), 3)
        self.assertEqual(len(labels), len(self.pd_formation.stable_entries))

    def test_get_plot(self):
        # Some very basic non-tests. Just to make sure the methods are callable.
        import matplotlib
        matplotlib.use("pdf")
        self.plotter.get_plot()
        self.plotter3d.get_plot()
        self.plotter.get_plot(energy_colormap="Reds", process_attributes=True)
        plt = self.plotter3d.get_plot(energy_colormap="Reds", process_attributes=True)
        self.plotter.get_plot(energy_colormap="Reds", process_attributes=False)
        plt = self.plotter3d.get_plot(energy_colormap="Reds",
                                      process_attributes=False)
        self.plotter.get_chempot_range_map_plot([Element("Li"), Element("O")])
        # self.plotter.get_contour_pd_plot()


class UtilityFunctionTest(unittest.TestCase):

    def test_unique_lines(self):
        testdata = [[5, 53, 353], [399, 20, 52], [399, 400, 20], [13, 399, 52],
                    [21, 400, 353], [393, 5, 353], [400, 393, 353],
                    [393, 400, 399], [393, 13, 5], [13, 393, 399],
                    [400, 17, 20], [21, 17, 400]]
        expected_ans = set([(5, 393), (21, 353), (353, 400), (5, 13), (17, 20),
                            (21, 400), (17, 400), (52, 399), (393, 399),
                            (20, 52), (353, 393), (5, 353), (5, 53), (13, 399),
                            (393, 400), (13, 52), (53, 353), (17, 21),
                            (13, 393), (20, 399), (399, 400), (20, 400)])
        self.assertEqual(uniquelines(testdata), expected_ans)

    def test_triangular_coord(self):
        coord = [0.5, 0.5]
        coord = triangular_coord(coord)
        self.assertTrue(np.allclose(coord, [ 0.75, 0.4330127]))

    def test_tet_coord(self):
        coord = [0.5, 0.5, 0.5]
        coord = tet_coord(coord)
        self.assertTrue(np.allclose(coord, [ 1., 0.57735027, 0.40824829]))

if __name__ == '__main__':
    unittest.main()
