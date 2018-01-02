# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest
import os
from monty.serialization import loadfn
import warnings
from pymatgen.analysis.pourbaix.maker import PourbaixDiagram
from pymatgen.analysis.pourbaix.entry import PourbaixEntryIO
from pymatgen.analysis.pourbaix.plotter import PourbaixPlotter

test_dir = os.path.join(os.path.dirname(__file__), '..', '..', '..', '..',
                        'test_files')


class TestPourbaixPlotter(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore")

        module_dir = os.path.dirname(os.path.abspath(__file__))
        (elements, entries) = PourbaixEntryIO.from_csv(os.path.join(module_dir,
                                                    "test_entries.csv"))
        self.num_simplices = {"Zn(s)": 7, "ZnO2(s)": 7, "Zn[2+]": 4, "ZnO2[2-]": 4, "ZnHO2[-]": 4}
        self.e_above_hull_test = {"ZnHO[+]": 0.0693, "ZnO(aq)": 0.0624}
        self.decomp_test = {"ZnHO[+]": {"ZnO(s)": 0.5, "Zn[2+]": 0.5}, "ZnO(aq)": {"ZnO(s)": 1.0}}
        self.pd = PourbaixDiagram(entries)
        self.multi_data = loadfn(os.path.join(test_dir, 'multicomp_pbx.json'))
        self.plotter = PourbaixPlotter(self.pd)

    def tearDown(self):
        warnings.resetwarnings()

    def test_plot_pourbaix(self):
        # Default limits
        plt = self.plotter.get_pourbaix_plot()
        # Non-standard limits
        plt = self.plotter.get_pourbaix_plot(limits=[[-5, 4], [-2, 2]])
        
        # Try 3-D plot
        plot_3d = self.plotter._get_plot()
        plot_3d_unstable = self.plotter._get_plot(label_unstable=True)

        plt.close()
        plot_3d.close()
        plot_3d_unstable.close()

    def test_plot_entry_stability(self):
        entry = self.pd.all_entries[0]
        plt = self.plotter.plot_entry_stability(entry, limits=[[-2, 14], [-3, 3]])

        # binary system
        pd_binary = PourbaixDiagram(self.multi_data['binary'],
                                    comp_dict = {"Ag": 0.5, "Te": 0.5})
        binary_plotter = PourbaixPlotter(pd_binary)
        test_entry = pd_binary._unprocessed_entries[0]
        plt = binary_plotter.plot_entry_stability(test_entry)
        plt.close()

if __name__ == '__main__':
    unittest.main()
