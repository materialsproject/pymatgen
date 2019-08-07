# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
import platform
import re
import unittest
import warnings
import random
import sys
from pymatgen import SETTINGS, __version__ as pmg_version
from pymatgen.ext.matproj import MPRester, MPRestError
from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Structure, Composition
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.electronic_structure.dos import CompleteDos
from pymatgen.electronic_structure.bandstructure import (
    BandStructureSymmLine, BandStructure)
from pymatgen.entries.compatibility import MaterialsProjectCompatibility
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.analysis.pourbaix_diagram import PourbaixEntry, PourbaixDiagram
from pymatgen.analysis.wulff import WulffShape
from pymatgen.analysis.reaction_calculator import Reaction
from pymatgen.io.cif import CifParser
from pymatgen.phonon.bandstructure import PhononBandStructureSymmLine
from pymatgen.phonon.dos import CompletePhononDos
from pymatgen.util.testing import PymatgenTest

"""
Created on Jun 9, 2012
"""


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Jun 9, 2012"


@unittest.skipIf(not SETTINGS.get("PMG_MAPI_KEY"),
                 "PMG_MAPI_KEY environment variable not set.")
class MPResterTest(PymatgenTest):
    _multiprocess_shared_ = True

    def setUp(self):
        self.rester = MPRester()
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_get_all_materials_ids_doc(self):
        mids = self.rester.get_materials_ids("Al2O3")
        random.shuffle(mids)
        doc = self.rester.get_doc(mids.pop(0))
        self.assertEqual(doc["pretty_formula"], "Al2O3")

    def test_get_xas_data(self):
        # Test getting XAS data
        data = self.rester.get_xas_data("mp-19017", "Li")
        self.assertEqual("mp-19017,Li", data['mid_and_el'])
        self.assertAlmostEqual(data['spectrum']['x'][0], 55.178, places=2)
        self.assertAlmostEqual(data['spectrum']['y'][0], 0.0164634, places=2)
        
    def test_get_data(self):
        props = ["energy", "energy_per_atom", "formation_energy_per_atom",
                 "nsites", "unit_cell_formula", "pretty_formula", "is_hubbard",
                 "elements", "nelements", "e_above_hull", "hubbards",
                 "is_compatible", "task_ids",
                 "density", "icsd_ids", "total_magnetization"]
        # unicode literals have been reintroduced in py>3.2

        expected_vals = [-191.3359011, -6.833425039285714, -2.5515769497278913,
                         28, {'P': 4, 'Fe': 4, 'O': 16, 'Li': 4},
                         "LiFePO4", True, ['Li', 'O', 'P', 'Fe'], 4, 0.0,
                         {'Fe': 5.3, 'Li': 0.0, 'O': 0.0, 'P': 0.0}, True,
                         {'mp-19017', 'mp-540081', 'mp-601412'},
                         3.464840709092822,
                         [159107, 154117, 160776, 99860, 181272, 166815,
                          260571, 92198, 165000, 155580, 38209, 161479, 153699,
                          260569, 260570, 200155, 260572, 181341, 181342,
                          72545, 56291, 97764, 162282, 155635],
                         15.9996841]

        for (i, prop) in enumerate(props):
            if prop not in ['hubbards', 'unit_cell_formula', 'elements',
                            'icsd_ids', 'task_ids']:
                val = self.rester.get_data("mp-19017", prop=prop)[0][prop]
                self.assertAlmostEqual(expected_vals[i], val, places=2)
            elif prop in ["elements", "icsd_ids", "task_ids"]:
                upstream_vals = set(
                    self.rester.get_data("mp-19017", prop=prop)[0][prop])
                self.assertLessEqual(set(expected_vals[i]), upstream_vals)
            else:
                self.assertEqual(expected_vals[i],
                                 self.rester.get_data("mp-19017",
                                                      prop=prop)[0][prop])

        props = ['structure', 'initial_structure', 'final_structure', 'entry']
        for prop in props:
            obj = self.rester.get_data("mp-19017", prop=prop)[0][prop]
            if prop.endswith("structure"):
                self.assertIsInstance(obj, Structure)
            elif prop == "entry":
                obj = self.rester.get_data("mp-19017", prop=prop)[0][prop]
                self.assertIsInstance(obj, ComputedEntry)

        # Test chemsys search
        data = self.rester.get_data('Fe-Li-O', prop='unit_cell_formula')
        self.assertTrue(len(data) > 1)
        elements = {Element("Li"), Element("Fe"), Element("O")}
        for d in data:
            self.assertTrue(
                set(Composition(d['unit_cell_formula']).elements).issubset(
                    elements))

        self.assertRaises(MPRestError, self.rester.get_data, "Fe2O3",
                          "badmethod")

    def test_get_data(self):
        # Test getting supported properties
        self.assertNotEqual(self.rester.get_task_data("mp-30"), [])
        # Test aliasing
        data = self.rester.get_task_data("mp-30", "energy")
        self.assertAlmostEqual(data[0]["energy"], -4.09929227, places=2)

    def test_get_materials_id_from_task_id(self):
        self.assertEqual(self.rester.get_materials_id_from_task_id(
            "mp-540081"), "mp-19017")

    def test_get_materials_id_references(self):
        # nosetests pymatgen/matproj/tests/test_matproj.py:MPResterTest.test_get_materials_id_references
        m = MPRester()
        data = m.get_materials_id_references('mp-123')
        self.assertTrue(len(data) > 1000)

    def test_find_structure(self):
        # nosetests pymatgen/matproj/tests/test_matproj.py:MPResterTest.test_find_structure
        m = MPRester()
        ciffile = self.TEST_FILES_DIR / 'Fe3O4.cif'
        data = m.find_structure(str(ciffile))
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
        data = self.rester.query(
            criteria=criteria, properties=props, chunk_size=0)
        self.assertTrue(len(data) > 6)
        data = self.rester.query(
            criteria="*2O", properties=props, chunk_size=0)
        self.assertGreaterEqual(len(data), 52)
        self.assertIn("Li2O", (d["pretty_formula"] for d in data))

    def test_query_chunk_size(self):
        criteria = {"nelements": 2, "elements": "O"}
        props = ['pretty_formula']
        data1 = self.rester.query(
            criteria=criteria, properties=props, chunk_size=0)
        data2 = self.rester.query(
            criteria=criteria, properties=props, chunk_size=500)
        self.assertEqual({d['pretty_formula'] for d in data1},
                         {d['pretty_formula'] for d in data2})
        self.assertIn("Al2O3", {d['pretty_formula'] for d in data1})

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
        bs_unif = self.rester.get_bandstructure_by_material_id(
            "mp-2254", line_mode=False)
        self.assertIsInstance(bs_unif, BandStructure)
        self.assertNotIsInstance(bs_unif, BandStructureSymmLine)

    def test_get_phonon_data_by_material_id(self):
        bs = self.rester.get_phonon_bandstructure_by_material_id("mp-661")
        self.assertIsInstance(bs, PhononBandStructureSymmLine)
        dos = self.rester.get_phonon_dos_by_material_id("mp-661")
        self.assertIsInstance(dos, CompletePhononDos)
        ddb_str = self.rester.get_phonon_ddb_by_material_id("mp-661")
        self.assertIsInstance(ddb_str, str)

    def test_get_structures(self):
        structs = self.rester.get_structures("Mn3O4")
        self.assertTrue(len(structs) > 0)

    def test_get_entries(self):
        entries = self.rester.get_entries("TiO2")
        self.assertTrue(len(entries) > 1)
        for e in entries:
            self.assertEqual(e.composition.reduced_formula, "TiO2")

        entries = self.rester.get_entries("TiO2", inc_structure=True)
        self.assertTrue(len(entries) > 1)
        for e in entries:
            self.assertEqual(e.structure.composition.reduced_formula, "TiO2")

        # all_entries = self.rester.get_entries("Fe", compatible_only=False)
        # entries = self.rester.get_entries("Fe", compatible_only=True)
        # self.assertTrue(len(entries) < len(all_entries))

        entries = self.rester.get_entries("Fe", compatible_only=True,
                                          property_data=["cif"])
        self.assertIn("cif", entries[0].data)

        for e in self.rester.get_entries("CdO2", inc_structure=False):
            self.assertIsNotNone(e.data["oxide_type"])

        # test if it will retrieve the conventional unit cell of Ni
        entry = self.rester.get_entry_by_material_id(
            "mp-23", inc_structure=True, conventional_unit_cell=True)
        Ni = entry.structure
        self.assertEqual(Ni.lattice.a, Ni.lattice.b)
        self.assertEqual(Ni.lattice.a, Ni.lattice.c)
        self.assertEqual(Ni.lattice.alpha, 90)
        self.assertEqual(Ni.lattice.beta, 90)
        self.assertEqual(Ni.lattice.gamma, 90)

        # Ensure energy per atom is same
        primNi = self.rester.get_entry_by_material_id(
            "mp-23", inc_structure=True, conventional_unit_cell=False)
        self.assertEqual(primNi.energy_per_atom, entry.energy_per_atom)

        Ni = self.rester.get_structure_by_material_id(
            "mp-23", conventional_unit_cell=True)
        self.assertEqual(Ni.lattice.a, Ni.lattice.b)
        self.assertEqual(Ni.lattice.a, Ni.lattice.c)
        self.assertEqual(Ni.lattice.alpha, 90)
        self.assertEqual(Ni.lattice.beta, 90)
        self.assertEqual(Ni.lattice.gamma, 90)

        # Test case where convs are different from initial and final
        th = self.rester.get_structure_by_material_id(
            "mp-37", conventional_unit_cell=True)
        th_entry = self.rester.get_entry_by_material_id(
            "mp-37", inc_structure=True, conventional_unit_cell=True)
        th_entry_initial = self.rester.get_entry_by_material_id(
            "mp-37", inc_structure="initial", conventional_unit_cell=True)
        self.assertEqual(th, th_entry.structure)
        self.assertEqual(len(th_entry.structure), 4)
        self.assertEqual(len(th_entry_initial.structure), 2)

        # Test if the polymorphs of Fe are properly sorted
        # by e_above_hull when sort_by_e_above_hull=True
        Fe_entries = self.rester.get_entries("Fe", sort_by_e_above_hull=True)
        self.assertEqual(Fe_entries[0].data["e_above_hull"], 0)

    def test_get_pourbaix_entries(self):
        pbx_entries = self.rester.get_pourbaix_entries(["Fe", "Cr"])
        for pbx_entry in pbx_entries:
            self.assertTrue(isinstance(pbx_entry, PourbaixEntry))
        # Ensure entries are pourbaix compatible
        pbx = PourbaixDiagram(pbx_entries)

        # Try binary system
        #pbx_entries = self.rester.get_pourbaix_entries(["Fe", "Cr"])
        #pbx = PourbaixDiagram(pbx_entries)

        # TODO: Shyue Ping: I do not understand this test. You seem to
        # be grabbing Zn-S system, but I don't see proper test for anything,
        # including Na ref. This test also takes a long time.

        # Test Zn-S, which has Na in reference solids
        # pbx_entries = self.rester.get_pourbaix_entries(["Zn", "S"])

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
        for e in all_entries:
            if str(e.entry_id).startswith("mod"):
                for d in rest_ehulls:
                    if d["entry_id"] == e.entry_id:
                        data = d
                        break
                self.assertAlmostEqual(pd.get_e_above_hull(e),
                                       data["e_above_hull"])

    def test_get_reaction(self):
        rxn = self.rester.get_reaction(["Li", "O"], ["Li2O"])
        self.assertIn("Li2O", rxn["Experimental_references"])

    def test_get_substrates(self):
        substrate_data = self.rester.get_substrates('mp-123', 5, [1, 0, 0])
        substrates = [sub_dict['sub_id'] for sub_dict in substrate_data]
        self.assertIn("mp-2534", substrates)

    def test_get_surface_data(self):
        data = self.rester.get_surface_data("mp-126") # Pt
        one_surf = self.rester.get_surface_data('mp-129', miller_index=[-2, -3, 1])
        self.assertAlmostEqual(one_surf['surface_energy'], 2.99156963, places=2)
        self.assertArrayAlmostEqual(one_surf['miller_index'], [3, 2, 1])
        self.assertIn("surfaces", data)
        surfaces = data["surfaces"]
        self.assertTrue(len(surfaces) > 0)
        surface = surfaces.pop()
        self.assertIn("miller_index", surface)
        self.assertIn("surface_energy", surface)
        self.assertIn("is_reconstructed", surface)
        data_inc = self.rester.get_surface_data("mp-126", inc_structures=True)
        self.assertIn("structure", data_inc["surfaces"][0])

    def test_get_wulff_shape(self):
        ws = self.rester.get_wulff_shape("mp-126")
        self.assertTrue(isinstance(ws, WulffShape))

    def test_get_cohesive_energy(self):
        ecoh = self.rester.get_cohesive_energy("mp-13")
        self.assertTrue(ecoh, 5.04543279)

    def test_get_gb_data(self):
        mo_gbs = self.rester.get_gb_data(chemsys='Mo')
        self.assertEqual(len(mo_gbs), 10)
        mo_gbs_s5 = self.rester.get_gb_data(pretty_formula='Mo', sigma=5)
        self.assertEqual(len(mo_gbs_s5), 3)
        mo_s3_112 = self.rester.get_gb_data(material_id='mp-129', sigma=3,
                                            gb_plane=[1, -1, -2],
                                            include_work_of_separation=True)
        self.assertEqual(len(mo_s3_112), 1)
        gb_f = mo_s3_112[0]['final_structure']
        self.assertArrayAlmostEqual(gb_f.rotation_axis, [1, 1, 0])
        self.assertAlmostEqual(gb_f.rotation_angle, 109.47122, places=4)
        self.assertAlmostEqual(mo_s3_112[0]['gb_energy'], 0.47965, places=2)
        self.assertAlmostEqual(mo_s3_112[0]['work_of_separation'], 6.318144, places=2)
        self.assertIn("Mo24", gb_f.formula)
        hcp_s7 = self.rester.get_gb_data(material_id='mp-87', gb_plane=[0, 0, 0, 1],
                                         include_work_of_separation=True)
        self.assertAlmostEqual(hcp_s7[0]['gb_energy'], 1.12, places=2)
        self.assertAlmostEqual(hcp_s7[0]['work_of_separation'], 2.46, places=2)


    def test_get_interface_reactions(self):
        kinks = self.rester.get_interface_reactions("LiCoO2", "Li3PS4")
        self.assertTrue(len(kinks) > 0)
        kink = kinks[0]
        self.assertIn("energy", kink)
        self.assertIn("ratio_atomic", kink)
        self.assertIn("rxn", kink)
        self.assertTrue(isinstance(kink['rxn'], Reaction))
        kinks_open_O = self.rester.get_interface_reactions(
            "LiCoO2", "Li3PS4", open_el="O", relative_mu=-1)
        self.assertTrue(len(kinks_open_O) > 0)
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings("always", message="The reactant.+")
            self.rester.get_interface_reactions("LiCoO2", "MnO9")
            self.assertTrue("The reactant" in str(w[-1].message))

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

    def test_include_user_agent(self):
        headers = self.rester.session.headers
        self.assertIn("user-agent", headers, msg="Include user-agent header by default")
        m = re.match(
            r"pymatgen/(\d+)\.(\d+)\.(\d+) \(Python/(\d+)\.(\d)+\.(\d+) ([^\/]*)/([^\)]*)\)",
            headers['user-agent'])
        self.assertIsNotNone(m, msg="Unexpected user-agent value {}".format(headers['user-agent']))
        self.assertEqual(m.groups()[:3], tuple(pmg_version.split(".")))
        self.assertEqual(
            m.groups()[3:6],
            tuple(str(n) for n in (sys.version_info.major, sys.version_info.minor, sys.version_info.micro))
        )
        self.rester = MPRester(include_user_agent=False)
        self.assertNotIn("user-agent", self.rester.session.headers, msg="user-agent header unwanted")


if __name__ == "__main__":
    unittest.main()
