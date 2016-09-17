# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import unittest2 as unittest
import os

from pymatgen import SETTINGS
from pymatgen.matproj.rest import MPRester, MPRestError
from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Structure, Composition
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.electronic_structure.dos import CompleteDos
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from pymatgen.entries.compatibility import MaterialsProjectCompatibility
from pymatgen.phasediagram.maker import PhaseDiagram
from pymatgen.phasediagram.analyzer import PDAnalyzer
from pymatgen.io.cif import CifParser

"""
Created on Jun 9, 2012
"""


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Jun 9, 2012"


test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


@unittest.skipIf("MAPI_KEY" not in SETTINGS,
                 "MAPI_KEY environment variable not set.")
class MPResterTest(unittest.TestCase):

    def setUp(self):
        self.rester = MPRester()

    def test_get_data(self):
        props = ["energy", "energy_per_atom", "formation_energy_per_atom",
                 "nsites", "unit_cell_formula", "pretty_formula", "is_hubbard",
                 "elements", "nelements", "e_above_hull", "hubbards",
                 "is_compatible", "task_ids",
                 "density", "icsd_ids", "total_magnetization"]
        # unicode literals have been reintroduced in py>3.2
        expected_vals = [-191.33812137, -6.833504334642858, -2.551358929370749,
                         28, {k: v for k, v in {'P': 4, 'Fe': 4, 'O': 16, 'Li': 4}.items()},
                         "LiFePO4", True, ['Li', 'O', 'P', 'Fe'], 4, 0.0,
                         {k: v for k, v in {'Fe': 5.3, 'Li': 0.0, 'O': 0.0, 'P': 0.0}.items()}, True,
                         [u'mp-601412', u'mp-19017', u'mp-796535', u'mp-797820',
                          u'mp-540081', u'mp-797269'],
                         3.4662026991351147,
                         [159107, 154117, 160776, 99860, 181272, 166815,
                          260571, 92198, 165000, 155580, 38209, 161479, 153699,
                          260569, 260570, 200155, 260572, 181341, 181342,
                          72545, 56291, 97764, 162282, 155635],
                         16.0002716]

        for (i, prop) in enumerate(props):
            if prop not in ['hubbards', 'unit_cell_formula', 'elements',
                            'icsd_ids', 'task_ids']:
                val = self.rester.get_data("mp-19017", prop=prop)[0][prop]
                self.assertAlmostEqual(expected_vals[i], val)
            elif prop in ["elements", "icsd_ids", "task_ids"]:
                self.assertEqual(set(expected_vals[i]),
                                 set(self.rester.get_data("mp-19017",
                                                          prop=prop)[0][prop]))
            else:
                self.assertEqual(expected_vals[i],
                                 self.rester.get_data("mp-19017",
                                                      prop=prop)[0][prop])

        props = ['structure', 'initial_structure', 'final_structure', 'entry']
        for prop in props:
            print(prop)
            obj = self.rester.get_data("mp-19017", prop=prop)[0][prop]
            if prop.endswith("structure"):
                self.assertIsInstance(obj, Structure)
            elif prop == "entry":
                obj = self.rester.get_data("mp-19017", prop=prop)[0][prop]
                self.assertIsInstance(obj, ComputedEntry)

        #Test chemsys search
        data = self.rester.get_data('Fe-Li-O', prop='unit_cell_formula')
        self.assertTrue(len(data) > 1)
        elements = {Element("Li"), Element("Fe"), Element("O")}
        for d in data:
            self.assertTrue(
                set(Composition(d['unit_cell_formula']).elements).issubset(
                    elements))

        self.assertRaises(MPRestError, self.rester.get_data, "Fe2O3",
                          "badmethod")

    def test_get_materials_id_from_task_id(self):
        self.assertEqual(self.rester.get_materials_id_from_task_id(
            "mp-540081"), "mp-19017")

    def test_get_materials_id_references(self):
        # nosetests pymatgen/matproj/tests/test_rest.py:MPResterTest.test_get_materials_id_references
        m = MPRester()
        data = m.get_materials_id_references('mp-123')
        self.assertTrue(len(data) > 1000)

    def test_find_structure(self):
        # nosetests pymatgen/matproj/tests/test_rest.py:MPResterTest.test_find_structure
        m = MPRester()
        ciffile = os.path.join(test_dir, 'Fe3O4.cif')
        data = m.find_structure(ciffile)
        self.assertTrue(len(data) > 1)
        s = CifParser(ciffile).get_structures()[0]
        data = m.find_structure(s)
        self.assertTrue(len(data) > 1)

    def test_get_entries_in_chemsys(self):
        syms = ["Li", "Fe", "O"]
        entries = self.rester.get_entries_in_chemsys(syms)
        elements = set([Element(sym) for sym in syms])
        for e in entries:
            self.assertIsInstance(e, ComputedEntry)
            self.assertTrue(set(e.composition.elements).issubset(elements))

    def test_get_structure_by_material_id(self):
        s1 = self.rester.get_structure_by_material_id("mp-1")
        self.assertEqual(s1.formula, "Cs1")

    def test_get_entry_by_material_id(self):
        e = self.rester.get_entry_by_material_id("mp-19017")
        self.assertIsInstance(e, ComputedEntry)
        self.assertTrue(e.composition.reduced_formula, "LiFePO4")

    def test_query(self):
        criteria = {'elements': {'$in': ['Li', 'Na', 'K'], '$all': ['O']}}
        props = ['pretty_formula', 'energy']
        data = self.rester.query(criteria=criteria, properties=props)
        self.assertTrue(len(data) > 6)
        data = self.rester.query(criteria="*2O", properties=props)
        self.assertGreaterEqual(len(data), 52)
        self.assertIn("Li2O", (d["pretty_formula"] for d in data))

    def test_get_exp_thermo_data(self):
        data = self.rester.get_exp_thermo_data("Fe2O3")
        self.assertTrue(len(data) > 0)
        for d in data:
            self.assertEqual(d.formula, "Fe2O3")

    def test_get_dos_by_id(self):
        dos = self.rester.get_dos_by_material_id("mp-2254")
        self.assertIsInstance(dos, CompleteDos)

    def test_get_bandstructure_by_material_id(self):
        bs = self.rester.get_bandstructure_by_material_id("mp-2254")
        self.assertIsInstance(bs, BandStructureSymmLine)

    def test_get_structures(self):
        structs = self.rester.get_structures("Mn3O4")
        self.assertTrue(len(structs) > 0)

    def test_get_entries(self):
        entries = self.rester.get_entries("TiO2")
        self.assertTrue(len(entries) > 1)
        for e in entries:
            self.assertEqual(e.composition.reduced_formula, "TiO2")

        entries = self.rester.get_entries("TiO2", inc_structure="final")
        self.assertTrue(len(entries) > 1)
        for e in entries:
            self.assertEqual(e.structure.composition.reduced_formula, "TiO2")

        all_entries = self.rester.get_entries("Fe", compatible_only=False)
        entries = self.rester.get_entries("Fe", compatible_only=True)
        self.assertTrue(len(entries) < len(all_entries))

        entries = self.rester.get_entries("Fe", compatible_only=True,
                                          property_data=["cif"])
        self.assertIn("cif", entries[0].data)

    def test_get_exp_entry(self):
        entry = self.rester.get_exp_entry("Fe2O3")
        self.assertEqual(entry.energy, -825.5)

    def test_submit_query_delete_snl(self):
        s = Structure([[5, 0, 0], [0, 5, 0], [0, 0, 5]], ["Fe"], [[0, 0, 0]])
        # d = self.rester.submit_snl(
        #     [s, s], remarks=["unittest"],
        #     authors="Test User <test@materialsproject.com>")
        # self.assertEqual(len(d), 2)
        # data = self.rester.query_snl({"about.remarks": "unittest"})
        # self.assertEqual(len(data), 2)
        # snlids = [d["_id"] for d in data]
        # self.rester.delete_snl(snlids)
        # data = self.rester.query_snl({"about.remarks": "unittest"})
        # self.assertEqual(len(data), 0)

    def test_get_stability(self):
        entries = self.rester.get_entries_in_chemsys(["Fe", "O"])
        modified_entries = []
        for entry in entries:
            # Create modified entries with energies that are 0.01eV higher
            # than the corresponding entries.
            if entry.composition.reduced_formula == "Fe2O3":
                modified_entries.append(
                    ComputedEntry(entry.composition,
                                  entry.uncorrected_energy + 0.01,
                                  parameters=entry.parameters,
                                  entry_id="mod_{}".format(entry.entry_id)))
        rest_ehulls = self.rester.get_stability(modified_entries)
        all_entries = entries + modified_entries
        compat = MaterialsProjectCompatibility()
        all_entries = compat.process_entries(all_entries)
        pd = PhaseDiagram(all_entries)
        a = PDAnalyzer(pd)
        for e in all_entries:
            if str(e.entry_id).startswith("mod"):
                for d in rest_ehulls:
                    if d["entry_id"] == e.entry_id:
                        data = d
                        break
                self.assertAlmostEqual(a.get_e_above_hull(e),
                                       data["e_above_hull"])

    def test_get_reaction(self):
        rxn = self.rester.get_reaction(["Li", "O"], ["Li2O"])
        self.assertIn("Li2O", rxn["Experimental_references"])

    def test_get_substrates(self):
        substrate_data = self.rester.get_substrates('mp-123', 5, [1, 0, 0])
        substrates = [sub_dict['sub_id'] for sub_dict in substrate_data]
        self.assertIn("mp-2534", substrates)

    def test_parse_criteria(self):
        crit = MPRester.parse_criteria("mp-1234 Li-*")
        self.assertIn("Li-O", crit["$or"][1]["chemsys"]["$in"])
        self.assertIn({"task_id": "mp-1234"}, crit["$or"])

        crit = MPRester.parse_criteria("Li2*")
        self.assertIn("Li2O", crit["pretty_formula"]["$in"])
        self.assertIn("Li2I", crit["pretty_formula"]["$in"])
        self.assertIn("CsLi2", crit["pretty_formula"]["$in"])

        crit = MPRester.parse_criteria("Li-*-*")
        self.assertIn("Li-Re-Ru", crit["chemsys"]["$in"])
        self.assertNotIn("Li-Li", crit["chemsys"]["$in"])

        comps = MPRester.parse_criteria("**O3")["pretty_formula"]["$in"]
        for c in comps:
            self.assertEqual(len(Composition(c)), 3, "Failed in %s" % c)

        chemsys = MPRester.parse_criteria("{Fe,Mn}-O")["chemsys"]["$in"]
        self.assertEqual(len(chemsys), 2)
        comps = MPRester.parse_criteria("{Fe,Mn,Co}O")["pretty_formula"]["$in"]
        self.assertEqual(len(comps), 3, comps)

        #Let's test some invalid symbols

        self.assertRaises(ValueError, MPRester.parse_criteria, "li-fe")
        self.assertRaises(ValueError, MPRester.parse_criteria, "LO2")

        crit = MPRester.parse_criteria("POPO2")
        self.assertIn("P2O3", crit["pretty_formula"]["$in"])


if __name__ == "__main__":
    unittest.main()
