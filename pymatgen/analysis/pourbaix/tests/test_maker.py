# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals


import unittest
import os
from monty.serialization import loadfn
import warnings

from pymatgen.analysis.pourbaix.maker import PourbaixDiagram, PourbaixEntry,\
    MultiEntry, PourbaixPlotter, IonEntry
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.core.ion import Ion


class TestPourbaixEntry(unittest.TestCase):
    """
    Test all functions using a fictitious entry
    """
    def setUp(self):
        # comp = Composition("Mn2O3")
        self.solentry = ComputedEntry("Mn2O3", 49)
        ion = Ion.from_formula("MnO4-")
        self.ionentry = IonEntry(ion, 25)
        self.PxIon = PourbaixEntry(self.ionentry)
        self.PxSol = PourbaixEntry(self.solentry)
        self.PxIon.concentration = 1e-4

    def test_pourbaix_entry(self):
        self.assertEqual(self.PxIon.entry.energy, 25, "Wrong Energy!")
        self.assertEqual(self.PxIon.entry.name,\
                          "MnO4[-]", "Wrong Entry!")
        self.assertEqual(self.PxSol.entry.energy, 49, "Wrong Energy!")
        self.assertEqual(self.PxSol.entry.name,\
                           "Mn2O3", "Wrong Entry!")
        # self.assertEqual(self.PxIon.energy, 25, "Wrong Energy!")
        # self.assertEqual(self.PxSol.energy, 49, "Wrong Energy!")
        self.assertEqual(self.PxIon.concentration, 1e-4, "Wrong concentration!")

    def test_calc_coeff_terms(self):
        self.assertEqual(self.PxIon.npH, -8, "Wrong npH!")
        self.assertEqual(self.PxIon.nPhi, -7, "Wrong nPhi!")
        self.assertEqual(self.PxIon.nH2O, 4, "Wrong nH2O!")

        self.assertEqual(self.PxSol.npH, -6, "Wrong npH!")
        self.assertEqual(self.PxSol.nPhi, -6, "Wrong nPhi!")
        self.assertEqual(self.PxSol.nH2O, 3, "Wrong nH2O!")

    def test_to_from_dict(self):
        d = self.PxIon.as_dict()
        ion_entry = self.PxIon.from_dict(d)
        self.assertEqual(ion_entry.entry.name, "MnO4[-]", "Wrong Entry!")

        d = self.PxSol.as_dict()
        sol_entry = self.PxSol.from_dict(d)
        self.assertEqual(sol_entry.name, "Mn2O3(s)", "Wrong Entry!")
        self.assertEqual(sol_entry.energy, self.PxSol.energy,
                         "as_dict and from_dict energies unequal")

class TestPourbaixDiagram(unittest.TestCase):
    def setUp(self):
        module_dir = os.path.dirname(os.path.abspath(__file__))
        self.test_data = loadfn(os.path.join(module_dir, 'pourbaix_test_data.json'))

    def test_pourbaix_diagram(self):
        # Single
        pbx = PourbaixDiagram(self.test_data['Zn'])
        self.assertEqual(set([e.name for e in pbx.stable_entries]),
                         {"ZnO(s)", "Zn[2+]", "ZnHO2[-]", "ZnO2[2-]", "Zn(s)"},
                         "List of stable entries does not match")

        pbx_nofilter = PourbaixDiagram(self.test_data['Zn'], filter_solids=False)
        self.assertEqual(set([e.name for e in pbx_nofilter.stable_entries]),
                         {"ZnO(s)", "Zn[2+]", "ZnHO2[-]", "ZnO2[2-]", "Zn(s)",
                          "ZnO2(s)", "ZnH(s)"},
                         "List of stable entries for unfiltered pbx does not match")

        pbx_lowconc = PourbaixDiagram(self.test_data['Zn'], conc_dict={"Zn": 1e-8})
        self.assertEqual(set([e.name for e in pbx_lowconc.stable_entries]),
                         {"Zn(HO)2(aq)", "Zn[2+]", "ZnHO2[-]", "ZnO2[2-]", "Zn(s)"})

    def test_get_pourbaix_domains(self):
        domains = PourbaixDiagram.get_pourbaix_domains(self.test_data['Zn'])

    def test_get_decomposition(self):
        for entry in [entry for entry in self.pd.all_entries
                      if entry not in self.pd.stable_entries]:
            decomp_entries = self.analyzer.get_decomposition(entry)
            for entr in decomp_entries:
                self.assertEqual(decomp_entries[entr], self.decomp_test[entry.name][entr.name])
            e_above_hull = self.analyzer.get_e_above_hull(entry)
            self.assertAlmostEqual(e_above_hull, self.e_above_hull_test[entry.name], 3)

    def test_mpr_pipeline(self):
        from pymatgen import MPRester
        mpr = MPRester()
        data = mpr.get_pourbaix_entries(["Zn"])
        pbx = PourbaixDiagram(data, filter_solids=True, conc_dict={"Zn": 1e-8})
        pbx.find_stable_entry(10, 0)


class TestPourbaixPlotter(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore")

        module_dir = os.path.dirname(os.path.abspath(__file__))
        self.test_data = loadfn(os.path.join(module_dir, "pourbaix_test_data.json"))
        self.pd = PourbaixDiagram(self.test_data["Zn"])
        self.plotter = PourbaixPlotter(self.pd)

    def tearDown(self):
        warnings.resetwarnings()

    def test_plot_pourbaix(self):
        plotter = PourbaixPlotter(self.pd)
        # Default limits
        plt = plotter.get_pourbaix_plot()
        # Non-standard limits
        plt = plotter.get_pourbaix_plot(limits=[[-5, 4], [-2, 2]])

    def test_plot_entry_stability(self):
        entry = self.pd.all_entries[0]
        plt = self.plotter.plot_entry_stability(entry, limits=[[-2, 14], [-3, 3]])

        # binary system
        pd_binary = PourbaixDiagram(self.test_data['Ag-Te'],
                                    comp_dict = {"Ag": 0.5, "Te": 0.5})
        binary_plotter = PourbaixPlotter(pd_binary)
        test_entry = pd_binary._unprocessed_entries[0]
        plt = binary_plotter.plot_entry_stability(test_entry)
        plt.close()

if __name__ == '__main__':
    unittest.main()
