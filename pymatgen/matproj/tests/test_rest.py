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
import os
import json

from pymatgen.serializers.json_coders import PMGJSONDecoder
from pymatgen.matproj.rest import MPRester
from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Structure, Composition
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.electronic_structure.dos import CompleteDos
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine

import pymatgen

test_dir = os.path.join(os.path.dirname(os.path.abspath(pymatgen.__file__)), '..', 'test_files')


class MPResterMock(MPRester):
    """
    A mock subclass for the REST adaptor to allow secure testing.
    """

    def get_data(self, chemsys_formula_id, prop=""):

        if chemsys_formula_id == "Fe-Li-O":
            filename = os.path.join(test_dir, "{}_{}.json".format("Fe-Li-O", prop))
            return json.load(open(filename, "r"), cls=PMGJSONDecoder)

        props = ["energy", "energy_per_atom", "formation_energy_per_atom",
                 "nsites", "formula", "pretty_formula", "is_hubbard",
                 "elements", "nelements", "e_above_hull", "hubbards", "is_compatible"]
        expected_vals = [-191.33404309, -6.83335868179, -2.5574286372085706, 28,
                         {u'P': 4, u'Fe': 4, u'O': 16, u'Li': 4}, "LiFePO4",
                         True, [u'Li', u'O', u'P', u'Fe'], 4, 0.0,
                         {u'Fe': 5.3, u'Li': 0.0, u'O': 0.0, u'P': 0.0}, True]
        try:
            i = props.index(prop)
            return [{prop:expected_vals[i]}]
        except ValueError:
            if prop in ("structure", "initial_structure", "final_structure"):
                obj = Structure([[1, 0, 0], [0, 1, 0], [0, 0, 1]], ["Si"], [[0, 0, 0]])
            elif prop == "bandstructure":
                filename = os.path.join(test_dir, "CaO_2605_bandstructure.json")
                d = json.load(open(filename, "r"))
                obj = BandStructureSymmLine.from_dict(d)
            elif prop == "entry":
                obj = ComputedEntry("Fe2O3", 0, 0, parameters={'hubbards':{}, "is_hubbard":False, "potcar_symbols":["PBE Fe_pv", "PBE O"]})
            elif prop == "dos":
                filename = os.path.join(test_dir, "complete_dos.json")
                d = json.load(open(filename, "r"))
                obj = CompleteDos.from_dict(d)
            return [{prop: obj}]

    def get_exp_thermo_data(self, formula):
        filename = os.path.join(test_dir, "Fe2O3_exp.json")
        return json.load(open(filename, "r"), cls=PMGJSONDecoder)

    def mpquery(self, criteria, properties):
        filename = os.path.join(test_dir, "mpquery.json")
        return json.load(open(filename, "r"), cls=PMGJSONDecoder)


class MPResterTest(unittest.TestCase):

    def setUp(self):
        self.adaptor = MPResterMock("")

    def test_get_data(self):
        props = ["energy", "energy_per_atom", "formation_energy_per_atom",
                 "nsites", "formula", "pretty_formula", "is_hubbard",
                 "elements", "nelements", "e_above_hull", "hubbards", "is_compatible"]
        expected_vals = [-191.33404309, -6.83335868179, -2.5574286372085706, 28,
                         {u'P': 4, u'Fe': 4, u'O': 16, u'Li': 4}, "LiFePO4",
                         True, [u'Li', u'O', u'P', u'Fe'], 4, 0.0,
                         {u'Fe': 5.3, u'Li': 0.0, u'O': 0.0, u'P': 0.0}, True]

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
        self.assertEqual(s1.formula, "Si1")

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


if __name__ == "__main__":
    unittest.main()
