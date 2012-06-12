#!/usr/bin/env python

'''
Created on Jun 9, 2012
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Jun 9, 2012"

import unittest

from pymatgen.matproj.rest import MPRestAdaptor
from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Structure, Composition
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.electronic_structure.dos import CompleteDos
from pymatgen.electronic_structure.band_structure.band_structure import BandStructureSymmLine


class MPRestAdaptorTest(unittest.TestCase):

    def setUp(self):
        self.adaptor = MPRestAdaptor("test_dev")

    def test_get_data(self):
        props = ["energy", "energy_per_atom", "formation_energy_per_atom",
                 "nsites", "formula", "pretty_formula", "is_hubbard",
                 "elements", "nelements", "e_above_hull", "hubbards", "is_compatible"]
        expected_vals = [-191.33404309, -6.83335868179, -2.5574286372085706, 28, {u'P': 4, u'Fe': 4, u'O': 16, u'Li': 4},
                "LiFePO4", True, [u'Li', u'O', u'P', u'Fe'], 4, 0.0, {u'Fe': 5.3, u'Li': 0.0, u'O': 0.0, u'P': 0.0}, True]

        for (i, prop) in enumerate(props):
            if prop != 'hubbards':
                self.assertAlmostEqual(expected_vals[i], self.adaptor.get_data(19017, prop)[0][prop])
            else:
                self.assertEqual(expected_vals[i], self.adaptor.get_data(19017, prop)[0][prop])

        props = ['structure', 'initial_structure', 'final_structure', 'entry']
        for prop in props:
            obj = self.adaptor.get_data(19017, prop)[0][prop]
            if prop.endswith("structure"):
                self.assertIsInstance(obj, Structure)
            elif prop == "entry":
                obj = self.adaptor.get_data(19017, prop)[0][prop]
                self.assertIsInstance(obj, ComputedEntry)

        #Test get all data.
        self.assertEqual(len(self.adaptor.get_data(19017)[0]), 13)

        #Test formula search
        data = self.adaptor.get_data('Fe2O3', 'formula')
        self.assertTrue(len(data) > 1)
        for d in data:
            self.assertEqual(Composition(d['formula']).reduced_formula, 'Fe2O3')

        #Test chemsys search
        data = self.adaptor.get_data('Fe-Li-O', 'formula')
        self.assertTrue(len(data) > 1)
        elements = set([Element("Li"), Element("Fe"), Element("O")])
        for d in data:
            self.assertTrue(set(Composition(d['formula']).elements).issubset(elements))

    def test_get_entries_in_chemsys(self):
        syms = ["Li", "Fe", "O"]
        all_entries = self.adaptor.get_entries_in_chemsys(syms, False)
        entries = self.adaptor.get_entries_in_chemsys(syms)
        self.assertTrue(len(entries) <= len(all_entries))
        elements = set([Element(sym) for sym in syms])
        for e in entries:
            self.assertIsInstance(e, ComputedEntry)
            self.assertTrue(set(e.composition.elements).issubset(elements))

    def test_get_structure_by_material_id(self):
        s1 = self.adaptor.get_structure_by_material_id(19017)
        s2 = self.adaptor.get_structure_by_material_id(19017, False)
        self.assertEqual(s1.formula, "Li4 Fe4 P4 O16")
        self.assertEqual(s2.formula, "Li4 Fe4 P4 O16")
        #Initial and final structures have different volumes 
        self.assertNotEqual(s1.volume, s2.volume)

    def test_get_entry_by_material_id(self):
        e = self.adaptor.get_entry_by_material_id(19017)
        self.assertIsInstance(e, ComputedEntry)
        self.assertTrue(e.composition.reduced_formula, "LiFePO4")

    def test_mpquery(self):
        criteria = {'elements':{'$in':['Li', 'Na', 'K'], '$all': ['O']}}
        props = ['formula', 'energy']

        data = self.adaptor.mpquery(criteria=criteria, properties=props)
        self.assertTrue(data['response']['num_results'] > 0)
        self.assertTrue(len(data['response']['results']) > 0)
        self.assertTrue(len(data['response']['results']) <= data['response']['num_results'])

    def test_get_exp_data(self):
        data = self.adaptor.get_exp_data("Fe2O3")
        self.assertTrue(len(data) > 0)
        for d in data:
            self.assertEqual(d.formula, "Fe2O3")

    def test_get_dos_by_id(self):
        dos = self.adaptor.get_dos_by_material_id(2254)
        self.assertIsInstance(dos, CompleteDos)

    def test_get_bandstructure_by_material_id(self):
        bs = self.adaptor.get_bandstructure_by_material_id(2254)
        self.assertIsInstance(bs, BandStructureSymmLine)


if __name__ == "__main__":
    unittest.main()
