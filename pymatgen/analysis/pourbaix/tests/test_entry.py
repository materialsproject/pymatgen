# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest
import os

from pymatgen.analysis.pourbaix.entry import PourbaixEntry, IonEntry, MultiEntry
from pymatgen.analysis.pourbaix.entry import PourbaixEntryIO
from pymatgen.phasediagram.entries import PDEntry
from pymatgen.core.ion import Ion

from pymatgen.core.structure import Composition
from pymatgen.core.periodic_table import Element


class TestPourbaixEntry(unittest.TestCase):
    """
    Test all functions using a fictitious entry
    """
    def setUp(self):
        comp = Composition("Mn2O3")
        self.solentry = PDEntry(comp, 49)
        ion = Ion.from_formula("MnO4-")
        self.ionentry = IonEntry(ion, 25)
        self.PxIon = PourbaixEntry(self.ionentry)
        self.PxSol = PourbaixEntry(self.solentry)
        self.PxIon.conc = 1e-4

    def test_pourbaix_entry(self):
        self.assertEqual(self.PxIon.entry.energy, 25, "Wrong Energy!")
        self.assertEqual(self.PxIon.entry.name,\
                          "MnO4[-]", "Wrong Entry!")
        self.assertEqual(self.PxSol.entry.energy, 49, "Wrong Energy!")
        self.assertEqual(self.PxSol.entry.name,\
                           "Mn2O3", "Wrong Entry!")
        self.assertEqual(self.PxIon.g0, 25, "Wrong Energy!")
        self.assertEqual(self.PxSol.g0, 49, "Wrong Energy!")
        self.assertEqual(self.PxIon.conc, 1e-4, "Wrong concentration!")

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


class MultiEntryTest(unittest.TestCase):
    """
    Test MultiEntry using fictitious entries
    """

    def setUp(self):
        entrylist = list()
        weights = list()
        comp = Composition("Mn2O3")
        entry = PDEntry(comp, 49)
        entrylist.append(PourbaixEntry(entry))
        weights.append(1.0)
        comp = Ion.from_formula("MnO4[-]")
        entry = IonEntry(comp, 25)
        entrylist.append(PourbaixEntry(entry))
        weights.append(0.25)
        comp = Composition("Fe2O3")
        entry = PDEntry(comp, 50)
        entrylist.append(PourbaixEntry(entry))
        weights.append(0.5)
        comp = Ion.from_formula("Fe[2+]")
        entry = IonEntry(comp, 15)
        entrylist.append(PourbaixEntry(entry))
        weights.append(2.5)
        comp = Ion.from_formula("Fe[3+]")
        entry = IonEntry(comp, 20)
        entrylist.append(PourbaixEntry(entry))
        weights.append(1.5)
        self.weights = weights
        self.entrylist = entrylist
        self.multientry = MultiEntry(entrylist, weights)

    def test_multi_entry(self):
        sum_g0 = 0.0
        sum_npH = 0.0
        sum_nPhi = 0.0
        sum_nH2O = 0.0
        for w, e in zip(self.weights, self.entrylist):
            sum_g0 += w * e.g0
            sum_npH += w * e.npH
            sum_nPhi += w * e.nPhi
            sum_nH2O += w * e.nH2O
        self.assertAlmostEqual(sum_g0, self.multientry.g0, "g0 doesn't match")
        self.assertAlmostEqual(sum_npH, self.multientry.npH, "npH doesn't match")
        self.assertAlmostEqual(sum_nPhi, self.multientry.nPhi, "nPhi doesn't match")
        self.assertAlmostEqual(sum_nH2O, self.multientry.nH2O, "nH2O doesn't match")


class IonEntryTest(unittest.TestCase):
    """
    Test IonEntry using fictitious entry
    """
    def setUp(self):
        ion = Ion.from_formula("MnO4[-]")
        self.entry = IonEntry(ion, 49)

    def test_get_energy(self):
        self.assertEqual(self.entry.energy, 49, "Wrong energy!")

    def test_get_name(self):
        self.assertEqual(self.entry.name, 'MnO4[-]', "Wrong name!")

    def test_get_composition(self):
        comp = self.entry.composition
        expected_comp = Ion.from_formula('MnO4[-]')
        self.assertEqual(comp, expected_comp, "Wrong composition!")

    def test_to_from_dict(self):
        d = self.entry.as_dict()
        entry = IonEntry.from_dict(d)

        self.assertEqual(entry.name, 'MnO4[-]', "Wrong name!")
        self.assertEqual(entry.energy_per_atom, 49.0 / 5)


class TestPourbaixEntryIO(unittest.TestCase):
    """
    Test Pourbaix Entry IO class
    """

    def test_read_write_csv(self):
        Zn_solids = ["Zn", "ZnO", "ZnO2"]
        sol_g = [0.0, -3.338, -1.315]
        Zn_ions = ["Zn[2+]", "ZnOH[+]", "HZnO2[-]", "ZnO2[2-]", "ZnO"]
        liq_g = [-1.527, -3.415, -4.812, -4.036, -2.921]
        liq_conc = [1e-6, 1e-6, 1e-6, 1e-6, 1e-6]
        solid_entry = list()
        for sol in Zn_solids:
            comp = Composition(sol)
            delg = sol_g[Zn_solids.index(sol)]
            solid_entry.append(PourbaixEntry(PDEntry(comp, delg)))
        ion_entry = list()
        for ion in Zn_ions:
            comp_ion = Ion.from_formula(ion)
            delg = liq_g[Zn_ions.index(ion)]
            conc = liq_conc[Zn_ions.index(ion)]
            PoE = PourbaixEntry(IonEntry(comp_ion, delg))
            PoE.conc = conc
            ion_entry.append(PoE)
        entries = solid_entry + ion_entry
        PourbaixEntryIO.to_csv("pourbaix_test_entries.csv", entries)

        (elements, entries) = PourbaixEntryIO.from_csv(
            "pourbaix_test_entries.csv")
        self.assertEqual(elements,
                         [Element('Zn'), Element('H'), Element('O')],
                         "Wrong elements!")
        self.assertEqual(len(entries), 8, "Wrong number of entries!")
        os.remove("pourbaix_test_entries.csv")


if __name__ == '__main__':
    unittest.main()
