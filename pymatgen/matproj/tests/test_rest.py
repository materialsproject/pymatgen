#!/usr/bin/env python

"""
Created on Jun 9, 2012
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Jun 9, 2012"

import unittest
import os

from pymatgen.matproj.rest import MPRester
from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Structure, Composition
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.electronic_structure.dos import CompleteDos
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from nose.exc import SkipTest

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class MPResterTest(unittest.TestCase):

    def setUp(self):
        if "MAPI_KEY" not in os.environ:
            raise SkipTest("MAPI_KEY environment variable not set. Skipping...")
        self.adaptor = MPRester()

    def test_get_data(self):
        props = ["energy", "energy_per_atom", "formation_energy_per_atom",
                 "nsites", "unit_cell_formula", "pretty_formula", "is_hubbard",
                 "elements", "nelements", "e_above_hull", "hubbards",
                 "is_compatible", "task_ids",
                 "density", "icsd_id", "total_magnetization"
                 ]
        expected_vals = [-191.33812137, -6.833504334642858, -2.5532843405078913, 28,
                         {u'P': 4, u'Fe': 4, u'O': 16, u'Li': 4}, "LiFePO4",
                         True, [u'Li', u'O', u'P', u'Fe'], 4, 0.0,
                         {u'Fe': 5.3, u'Li': 0.0, u'O': 0.0, u'P': 0.0}, True,
                         [540081, 19017], 3.4660546589064176,56291, 16.0002716]

        for (i, prop) in enumerate(props):
            if prop not in ['hubbards', 'unit_cell_formula', 'elements']:
                val =  self.adaptor.get_data(540081, prop=prop)[0][prop]
                self.assertAlmostEqual(expected_vals[i], val)
            elif prop == "elements":
                self.assertEqual(set(expected_vals[i]),
                                 set(self.adaptor.get_data(540081,
                                                        prop=prop)[0][prop]))
            else:
                self.assertEqual(expected_vals[i],
                                 self.adaptor.get_data(540081, prop=prop)[0][prop])

        props = ['structure', 'initial_structure', 'final_structure', 'entry']
        for prop in props:
            obj = self.adaptor.get_data(540081, prop=prop)[0][prop]
            if prop.endswith("structure"):
                self.assertIsInstance(obj, Structure)
            elif prop == "entry":
                obj = self.adaptor.get_data(540081, prop=prop)[0][prop]
                self.assertIsInstance(obj, ComputedEntry)

        #Test chemsys search
        data = self.adaptor.get_data('Fe-Li-O', prop='unit_cell_formula')
        self.assertTrue(len(data) > 1)
        elements = set([Element("Li"), Element("Fe"), Element("O")])
        for d in data:
            self.assertTrue(set(Composition(d['unit_cell_formula']).elements).issubset(elements))

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
        s1 = self.adaptor.get_structure_by_material_id(1)
        self.assertEqual(s1.formula, "Cs1")

    def test_get_entry_by_material_id(self):
        e = self.adaptor.get_entry_by_material_id(540081)
        self.assertIsInstance(e, ComputedEntry)
        self.assertTrue(e.composition.reduced_formula, "LiFePO4")

    def test_mpquery(self):
        criteria = {'elements':{'$in':['Li', 'Na', 'K'], '$all': ['O']}}
        props = ['formula', 'energy']
        data = self.adaptor.mpquery(criteria=criteria, properties=props)
        self.assertTrue(len(data) > 6)

    def test_get_exp_thermo_data(self):
        data = self.adaptor.get_exp_thermo_data("Fe2O3")
        self.assertTrue(len(data) > 0)
        for d in data:
            self.assertEqual(d.formula, "Fe2O3")

    def test_get_dos_by_id(self):
        dos = self.adaptor.get_dos_by_material_id(2254)
        self.assertIsInstance(dos, CompleteDos)

    def test_get_bandstructure_by_material_id(self):
        bs = self.adaptor.get_bandstructure_by_material_id(2254)
        self.assertIsInstance(bs, BandStructureSymmLine)

    def test_get_structures(self):
        structs = self.adaptor.get_structures("Mn3O4")
        self.assertTrue(len(structs) > 0)
        for s in structs:
            self.assertEqual(s.composition.reduced_formula, "Mn3O4")

    def test_get_entries(self):
        entries = self.adaptor.get_entries("TiO2")
        self.assertTrue(len(entries) > 1)
        for e in entries:
            self.assertEqual(e.composition.reduced_formula, "TiO2")

    def test_get_exp_entry(self):
        entry = self.adaptor.get_exp_entry("Fe2O3")
        self.assertEqual(entry.energy, -825.5)


if __name__ == "__main__":
    unittest.main()
