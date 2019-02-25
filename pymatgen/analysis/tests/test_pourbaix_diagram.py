# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.



import unittest
import os
from monty.serialization import loadfn
import warnings
import numpy as np
import multiprocessing

from pymatgen.analysis.pourbaix_diagram import PourbaixDiagram, PourbaixEntry,\
    PourbaixPlotter, IonEntry, MultiEntry
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.core.ion import Ion
from pymatgen import SETTINGS

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class PourbaixEntryTest(unittest.TestCase):
    _multiprocess_shared_ = True
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
        self.assertEqual(self.PxIon.entry.name,
                         "MnO4[-]", "Wrong Entry!")
        self.assertEqual(self.PxSol.entry.energy, 49, "Wrong Energy!")
        self.assertEqual(self.PxSol.entry.name,
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

    def test_energy_functions(self):
        # TODO: test these for values
        self.PxSol.energy_at_conditions(10, 0)
        self.PxSol.energy_at_conditions(np.array([1, 2, 3]), 0)
        self.PxSol.energy_at_conditions(10, np.array([1, 2, 3]))
        self.PxSol.energy_at_conditions(np.array([1, 2, 3]),
                                        np.array([1, 2, 3]))

    def test_multi_entry(self):
        # TODO: More robust multientry test
        m_entry = MultiEntry([self.PxSol, self.PxIon])
        for attr in ['energy', 'composition', 'nPhi']:
            self.assertEqual(getattr(m_entry, attr),
                             getattr(self.PxSol, attr) + getattr(self.PxIon, attr))

        # As dict, from dict
        m_entry_dict = m_entry.as_dict()
        m_entry_new = MultiEntry.from_dict(m_entry_dict)
        self.assertEqual(m_entry_new.energy, m_entry.energy)


class PourbaixDiagramTest(unittest.TestCase):
    _multiprocess_shared_ = True

    @classmethod
    def setUpClass(cls):
        cls.test_data = loadfn(os.path.join(test_dir, 'pourbaix_test_data.json'))
        cls.pbx = PourbaixDiagram(cls.test_data['Zn'], filter_solids=True)
        cls.pbx_nofilter = PourbaixDiagram(cls.test_data['Zn'],
                                           filter_solids=False)

    def test_pourbaix_diagram(self):
        self.assertEqual(set([e.name for e in self.pbx.stable_entries]),
                         {"ZnO(s)", "Zn[2+]", "ZnHO2[-]", "ZnO2[2-]", "Zn(s)"},
                         "List of stable entries does not match")

        self.assertEqual(set([e.name for e in self.pbx_nofilter.stable_entries]),
                         {"ZnO(s)", "Zn[2+]", "ZnHO2[-]", "ZnO2[2-]", "Zn(s)",
                          "ZnO2(s)", "ZnH(s)"},
                         "List of stable entries for unfiltered pbx does not match")

        pbx_lowconc = PourbaixDiagram(self.test_data['Zn'], conc_dict={"Zn": 1e-8},
                                      filter_solids=True)
        self.assertEqual(set([e.name for e in pbx_lowconc.stable_entries]),
                         {"Zn(HO)2(aq)", "Zn[2+]", "ZnHO2[-]", "ZnO2[2-]", "Zn(s)"})

    def test_properties(self):
        self.assertEqual(len(self.pbx.unstable_entries), 2)

    def test_multicomponent(self):
        # Assure no ions get filtered at high concentration
        ag_n = [e for e in self.test_data['Ag-Te-N']
                if not "Te" in e.composition]
        highconc = PourbaixDiagram(ag_n, filter_solids=True,
                                   conc_dict={"Ag": 1e-5, "N": 1})
        entry_sets = [set(e.entry_id) for e in highconc.stable_entries]
        self.assertIn({"mp-124", "ion-17"}, entry_sets)

        # Binary system
        pd_binary = PourbaixDiagram(self.test_data['Ag-Te'], filter_solids=True,
                                    comp_dict={"Ag": 0.5, "Te": 0.5},
                                    conc_dict={"Ag": 1e-8, "Te": 1e-8})
        self.assertEqual(len(pd_binary.stable_entries), 30)
        test_entry = pd_binary.find_stable_entry(8, 2)
        self.assertTrue("mp-499" in test_entry.entry_id)

        # Find a specific multientry to test
        self.assertEqual(pd_binary.get_decomposition_energy(test_entry, 8, 2), 0)
        self.assertEqual(pd_binary.get_decomposition_energy(
            test_entry.entry_list[0], 8, 2), 0)

        pd_ternary = PourbaixDiagram(self.test_data['Ag-Te-N'], filter_solids=True)
        self.assertEqual(len(pd_ternary.stable_entries), 49)

        ag = self.test_data['Ag-Te-N'][30]
        self.assertAlmostEqual(pd_ternary.get_decomposition_energy(ag, 2, -1), 0)
        self.assertAlmostEqual(pd_ternary.get_decomposition_energy(ag, 10, -2), 0)

        # Test invocation of pourbaix diagram from ternary data
        new_ternary = PourbaixDiagram(pd_ternary.all_entries)
        self.assertEqual(len(new_ternary.stable_entries), 49)
        self.assertAlmostEqual(new_ternary.get_decomposition_energy(ag, 2, -1), 0)
        self.assertAlmostEqual(new_ternary.get_decomposition_energy(ag, 10, -2), 0)

    def test_get_pourbaix_domains(self):
        domains = PourbaixDiagram.get_pourbaix_domains(self.test_data['Zn'])
        self.assertEqual(len(domains[0]), 7)

    def test_get_decomposition(self):
        # Test a stable entry to ensure that it's zero in the stable region
        entry = self.test_data['Zn'][12] # Should correspond to mp-2133
        self.assertAlmostEqual(self.pbx.get_decomposition_energy(entry, 10, 1),
                               0.0, 5, "Decomposition energy of ZnO is not 0.")

        # Test an unstable entry to ensure that it's never zero
        entry = self.test_data['Zn'][11]
        ph, v = np.meshgrid(np.linspace(0, 14), np.linspace(-2, 4))
        result = self.pbx_nofilter.get_decomposition_energy(entry, ph, v)
        self.assertTrue((result >= 0).all(),
                        "Unstable energy has hull energy of 0 or less")

        # Test an unstable hydride to ensure HER correction works
        self.assertAlmostEqual(self.pbx.get_decomposition_energy(entry, -3, -2),
                               11.093744395)
        # Test a list of pHs
        self.pbx.get_decomposition_energy(entry, np.linspace(0, 2, 5), 2)

        # Test a list of Vs
        self.pbx.get_decomposition_energy(entry, 4, np.linspace(-3, 3, 10))

        # Test a set of matching arrays
        ph, v = np.meshgrid(np.linspace(0, 14), np.linspace(-3, 3))
        self.pbx.get_decomposition_energy(entry, ph, v)

    def test_multielement_parallel(self):
        # Simple test to ensure that multiprocessing is working
        test_entries = self.test_data["Ag-Te-N"]
        nproc = multiprocessing.cpu_count()
        pbx = PourbaixDiagram(test_entries, filter_solids=True, nproc=nproc)
        self.assertEqual(len(pbx.stable_entries), 49)

    @unittest.skipIf(not SETTINGS.get("PMG_MAPI_KEY"),
                     "PMG_MAPI_KEY environment variable not set.")
    def test_mpr_pipeline(self):
        from pymatgen import MPRester
        mpr = MPRester()
        data = mpr.get_pourbaix_entries(["Zn"])
        pbx = PourbaixDiagram(data, filter_solids=True, conc_dict={"Zn": 1e-8})
        pbx.find_stable_entry(10, 0)

        data = mpr.get_pourbaix_entries(["Ag", "Te"])
        pbx = PourbaixDiagram(data, filter_solids=True,
                              conc_dict={"Ag": 1e-8, "Te": 1e-8})
        self.assertEqual(len(pbx.stable_entries), 30)
        test_entry = pbx.find_stable_entry(8, 2)
        self.assertAlmostEqual(test_entry.energy, 2.3894017960000009, 3)

        # Test custom ions
        entries = mpr.get_pourbaix_entries(["Sn", "C", "Na"])
        ion = IonEntry(Ion.from_formula("NaO28H80Sn12C24+"), -161.676)
        custom_ion_entry = PourbaixEntry(ion, entry_id='my_ion')
        pbx = PourbaixDiagram(entries + [custom_ion_entry], filter_solids=True,
                              comp_dict={"Na": 1, "Sn": 12, "C": 24})
        self.assertAlmostEqual(pbx.get_decomposition_energy(custom_ion_entry, 5, 2),
                               8.31202738629504, 2)

    def test_nofilter(self):
        entries = self.test_data['Ag-Te']
        pbx = PourbaixDiagram(entries)
        pbx.get_decomposition_energy(entries[0], 0, 0)

    def test_solid_filter(self):
        entries = self.test_data['Ag-Te-N']
        pbx = PourbaixDiagram(entries, filter_solids=True)
        pbx.get_decomposition_energy(entries[0], 0, 0)

    def test_serialization(self):
        d = self.pbx.as_dict()
        new = PourbaixDiagram.from_dict(d)
        self.assertEqual(set([e.name for e in new.stable_entries]),
                         {"ZnO(s)", "Zn[2+]", "ZnHO2[-]", "ZnO2[2-]", "Zn(s)"},
                         "List of stable entries does not match")

        # Test with unprocessed entries included, this should result in the
        # previously filtered entries being included
        d = self.pbx.as_dict(include_unprocessed_entries=True)
        new = PourbaixDiagram.from_dict(d)
        self.assertEqual(
            set([e.name for e in new.stable_entries]),
            {"ZnO(s)", "Zn[2+]", "ZnHO2[-]", "ZnO2[2-]", "Zn(s)", "ZnO2(s)", "ZnH(s)"},
            "List of stable entries for unfiltered pbx does not match")

        pd_binary = PourbaixDiagram(self.test_data['Ag-Te'], filter_solids=True,
                                    comp_dict={"Ag": 0.5, "Te": 0.5},
                                    conc_dict={"Ag": 1e-8, "Te": 1e-8})
        new_binary = PourbaixDiagram.from_dict(pd_binary.as_dict())
        self.assertEqual(len(pd_binary.stable_entries),
                         len(new_binary.stable_entries))


class PourbaixPlotterTest(unittest.TestCase):
    def setUp(self):
        warnings.simplefilter("ignore")
        self.test_data = loadfn(os.path.join(test_dir, "pourbaix_test_data.json"))
        self.pd = PourbaixDiagram(self.test_data["Zn"])
        self.plotter = PourbaixPlotter(self.pd)

    def tearDown(self):
        warnings.simplefilter("default")

    def test_plot_pourbaix(self):
        plotter = PourbaixPlotter(self.pd)
        # Default limits
        plotter.get_pourbaix_plot()
        # Non-standard limits
        plotter.get_pourbaix_plot(limits=[[-5, 4], [-2, 2]])

    def test_plot_entry_stability(self):
        entry = self.pd.all_entries[0]
        self.plotter.plot_entry_stability(entry, limits=[[-2, 14], [-3, 3]])

        # binary system
        pd_binary = PourbaixDiagram(self.test_data['Ag-Te'],
                                    comp_dict = {"Ag": 0.5, "Te": 0.5})
        binary_plotter = PourbaixPlotter(pd_binary)
        test_entry = pd_binary._unprocessed_entries[0]
        plt = binary_plotter.plot_entry_stability(test_entry)
        plt.close()


if __name__ == '__main__':
    unittest.main()
