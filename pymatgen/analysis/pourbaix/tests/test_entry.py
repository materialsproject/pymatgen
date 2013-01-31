#!/usr/bin/env python

import unittest
import os

from pymatgen.analysis.pourbaix.entry import PourbaixEntry, IonEntry
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
        self.PxIon.set_conc(1e-4)

    def test_pourbaix_entry(self):
        self.assertEquals(self.PxIon.entry.energy, 25, "Wrong Energy!")
        self.assertEquals(self.PxIon.entry.name,\
                          "MnO4[-]", "Wrong Entry!")
        self.assertEquals(self.PxSol.entry.energy, 49, "Wrong Energy!")
        self.assertEquals(self.PxSol.entry.name,\
                           "Mn2O3", "Wrong Entry!")
        self.assertEquals(self.PxIon.g0, 25, "Wrong Energy!")
        self.assertEquals(self.PxSol.g0, 49, "Wrong Energy!")
        self.assertEquals(self.PxIon.conc, 1e-4, "Wrong concentration!")

    def test_calc_coeff_terms(self):
        self.assertEquals(self.PxIon.npH, -8, "Wrong npH!")
        self.assertEquals(self.PxIon.nPhi, -7, "Wrong nPhi!")
        self.assertEquals(self.PxIon.nH2O, 4, "Wrong nH2O!")

        self.assertEquals(self.PxSol.npH, -6, "Wrong npH!")
        self.assertEquals(self.PxSol.nPhi, -6, "Wrong nPhi!")
        self.assertEquals(self.PxSol.nH2O, 3, "Wrong nH2O!")

    def test_to_from_dict(self):
        d = self.PxIon.to_dict
        ion_entry = self.PxIon.from_dict(d)
        self.assertEquals(ion_entry.entry.name, "MnO4[-]", "Wrong Entry!")


class IonEntryTest(unittest.TestCase):
    """
    Test IonEntry using fictitious entry
    """
    def setUp(self):
        ion = Ion.from_formula("MnO4[-]")
        self.entry = IonEntry(ion, 49)

    def test_get_energy(self):
        self.assertEquals(self.entry.energy, 49, "Wrong energy!")

    def test_get_name(self):
        self.assertEquals(self.entry.name, 'MnO4[-]', "Wrong name!")

    def test_get_composition(self):
        comp = self.entry.composition
        expected_comp = Ion.from_formula('MnO4[-]')
        self.assertEquals(comp, expected_comp, "Wrong composition!")

    def test_to_from_dict(self):
        d = self.entry.to_dict
        entry = IonEntry.from_dict(d)

        self.assertEquals(entry.name, 'MnO4[-]', "Wrong name!")
        self.assertEquals(entry.energy_per_atom, 49.0 / 5)


class TestPourbaixEntryIO(unittest.TestCase):
    """
    Test Pourbaix Entry IO class
    """

    def test_write_csv(self):
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
            PoE.set_conc(conc)
            ion_entry.append(PoE)
        entries = solid_entry + ion_entry
        PourbaixEntryIO.to_csv("test_entries.csv", entries)

    def test_read_csv(self):
        module_dir = os.path.dirname(os.path.abspath(__file__))
        (elements, entries) = PourbaixEntryIO.from_csv(os.path.join(module_dir,
                                                    "test_entries.csv"))
        self.assertEqual(elements,
                         [Element('Zn'), Element('H'), Element('O')],
                         "Wrong elements!")
        self.assertEqual(len(entries), 8, "Wrong number of entries!")


if __name__ == '__main__':
    unittest.main()