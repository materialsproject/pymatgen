# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
import random
import re
import unittest
import warnings

import requests

from pymatgen.core import SETTINGS, SETTINGS_FILE, yaml
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.analysis.pourbaix_diagram import PourbaixDiagram, PourbaixEntry
from pymatgen.analysis.reaction_calculator import Reaction
from pymatgen.analysis.wulff import WulffShape
from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Composition, Structure
from pymatgen.electronic_structure.bandstructure import (
    BandStructure,
    BandStructureSymmLine,
)
from pymatgen.electronic_structure.dos import CompleteDos
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.ext.matproj import MPRester, MPRestError, TaskType
from pymatgen.io.cif import CifParser
from pymatgen.phonon.bandstructure import PhononBandStructureSymmLine
from pymatgen.phonon.dos import CompletePhononDos
from pymatgen.util.testing import PymatgenTest


website_is_up = requests.get("https://www.materialsproject.org").status_code == 200


@unittest.skipIf(
    (not SETTINGS.get("PMG_MAPI_KEY")) or (not website_is_up),
    "PMG_MAPI_KEY environment variable not set or MP is down.",
)
class MPResterTest(PymatgenTest):
    _multiprocess_shared_ = True

    def setUp(self):
        self.rester = MPRester()
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")
        self.rester.session.close()

    def test_get_all_materials_ids_doc(self):
        mids = self.rester.get_materials_ids("Al2O3")
        random.shuffle(mids)
        doc = self.rester.get_doc(mids.pop(0))
        self.assertEqual(doc["pretty_formula"], "Al2O3")

    def test_get_xas_data(self):
        # Test getting XAS data
        data = self.rester.get_xas_data("mp-19017", "Li")
        self.assertEqual("mp-19017,Li", data["mid_and_el"])
        self.assertAlmostEqual(data["spectrum"]["x"][0], 55.178, places=2)
        self.assertAlmostEqual(data["spectrum"]["y"][0], 0.0164634, places=2)

    def test_get_data(self):
        props = {
            "energy",
            "energy_per_atom",
            "formation_energy_per_atom",
            "nsites",
            "unit_cell_formula",
            "pretty_formula",
            "is_hubbard",
            "elements",
            "nelements",
            "e_above_hull",
            "hubbards",
            "is_compatible",
            "task_ids",
            "density",
            "icsd_ids",
            "total_magnetization",
        }
        mpid = "mp-1143"
        vals = requests.get(f"http://www.materialsproject.org/materials/{mpid}/json/")
        expected_vals = vals.json()

        for prop in props:
            if prop not in [
                "hubbards",
                "unit_cell_formula",
                "elements",
                "icsd_ids",
                "task_ids",
            ]:
                val = self.rester.get_data(mpid, prop=prop)[0][prop]
                if prop in ["energy", "energy_per_atom"]:
                    prop = "final_" + prop
                self.assertAlmostEqual(expected_vals[prop], val, 2, "Failed with property %s" % prop)
            elif prop in ["elements", "icsd_ids", "task_ids"]:
                upstream_vals = set(self.rester.get_data(mpid, prop=prop)[0][prop])
                self.assertLessEqual(set(expected_vals[prop]), upstream_vals)
            else:
                self.assertEqual(
                    expected_vals[prop],
                    self.rester.get_data(mpid, prop=prop)[0][prop],
                )

        props = ["structure", "initial_structure", "final_structure", "entry"]
        for prop in props:
            obj = self.rester.get_data(mpid, prop=prop)[0][prop]
            if prop.endswith("structure"):
                self.assertIsInstance(obj, Structure)
            elif prop == "entry":
                obj = self.rester.get_data(mpid, prop=prop)[0][prop]
                self.assertIsInstance(obj, ComputedEntry)

        # Test chemsys search
        data = self.rester.get_data("Fe-Li-O", prop="unit_cell_formula")
        self.assertTrue(len(data) > 1)
        elements = {Element("Li"), Element("Fe"), Element("O")}
        for d in data:
            self.assertTrue(set(Composition(d["unit_cell_formula"]).elements).issubset(elements))

        self.assertRaises(MPRestError, self.rester.get_data, "Fe2O3", "badmethod")

    def test_get_materials_id_from_task_id(self):
        self.assertEqual(self.rester.get_materials_id_from_task_id("mp-540081"), "mp-19017")

    def test_get_materials_id_references(self):
        # nosetests pymatgen/matproj/tests/test_matproj.py:MPResterTest.test_get_materials_id_references
        m = MPRester()
        data = m.get_materials_id_references("mp-123")
        self.assertTrue(len(data) > 1000)

    def test_find_structure(self):
        # nosetests pymatgen/matproj/tests/test_matproj.py:MPResterTest.test_find_structure
        m = MPRester()
        ciffile = self.TEST_FILES_DIR / "Fe3O4.cif"
        data = m.find_structure(str(ciffile))
        self.assertTrue(len(data) > 1)
        s = CifParser(ciffile).get_structures()[0]
        data = m.find_structure(s)
        self.assertTrue(len(data) > 1)

    def test_get_entries_in_chemsys(self):
        syms = ["Li", "Fe", "O"]
        syms2 = "Li-Fe-O"
        entries = self.rester.get_entries_in_chemsys(syms)
        entries2 = self.rester.get_entries_in_chemsys(syms2)
        elements = set([Element(sym) for sym in syms])
        for e in entries:
            self.assertIsInstance(e, ComputedEntry)
            self.assertTrue(set(e.composition.elements).issubset(elements))

        e1 = set([i.entry_id for i in entries])
        e2 = set([i.entry_id for i in entries2])
        self.assertTrue(e1 == e2)

    def test_get_structure_by_material_id(self):
        s1 = self.rester.get_structure_by_material_id("mp-1")
        self.assertEqual(s1.formula, "Cs1")

        # requesting via task-id instead of mp-id
        self.assertWarns(Warning, self.rester.get_structure_by_material_id, "mp-698856")

        # requesting unknown mp-id
        self.assertRaises(MPRestError, self.rester.get_structure_by_material_id, "mp-does-not-exist")

    def test_get_entry_by_material_id(self):
        e = self.rester.get_entry_by_material_id("mp-19017")
        self.assertIsInstance(e, ComputedEntry)
        self.assertTrue(e.composition.reduced_formula, "LiFePO4")

    def test_query(self):
        criteria = {"elements": {"$in": ["Li", "Na", "K"], "$all": ["O"]}}
        props = ["pretty_formula", "energy"]
        data = self.rester.query(criteria=criteria, properties=props, chunk_size=0)
        self.assertTrue(len(data) > 6)
        data = self.rester.query(criteria="*2O", properties=props, chunk_size=0)
        self.assertGreaterEqual(len(data), 52)
        self.assertIn("Li2O", (d["pretty_formula"] for d in data))

    def test_query_chunk_size(self):
        criteria = {"nelements": 2, "elements": "O"}
        props = ["pretty_formula"]
        data1 = self.rester.query(criteria=criteria, properties=props, chunk_size=0)
        data2 = self.rester.query(criteria=criteria, properties=props, chunk_size=500)
        self.assertEqual({d["pretty_formula"] for d in data1}, {d["pretty_formula"] for d in data2})
        self.assertIn("Al2O3", {d["pretty_formula"] for d in data1})

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
        bs_unif = self.rester.get_bandstructure_by_material_id("mp-2254", line_mode=False)
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
        entries = self.rester.get_entries("Fe", compatible_only=True, property_data=["cif"])
        self.assertIn("cif", entries[0].data)

        for e in self.rester.get_entries("CdO2", inc_structure=False):
            self.assertIsNotNone(e.data["oxide_type"])

        # test if it will retrieve the conventional unit cell of Ni
        entry = self.rester.get_entry_by_material_id("mp-23", inc_structure=True, conventional_unit_cell=True)
        Ni = entry.structure
        self.assertEqual(Ni.lattice.a, Ni.lattice.b)
        self.assertEqual(Ni.lattice.a, Ni.lattice.c)
        self.assertEqual(Ni.lattice.alpha, 90)
        self.assertEqual(Ni.lattice.beta, 90)
        self.assertEqual(Ni.lattice.gamma, 90)

        # Ensure energy per atom is same
        primNi = self.rester.get_entry_by_material_id("mp-23", inc_structure=True, conventional_unit_cell=False)
        self.assertEqual(primNi.energy_per_atom, entry.energy_per_atom)

        Ni = self.rester.get_structure_by_material_id("mp-23", conventional_unit_cell=True)
        self.assertEqual(Ni.lattice.a, Ni.lattice.b)
        self.assertEqual(Ni.lattice.a, Ni.lattice.c)
        self.assertEqual(Ni.lattice.alpha, 90)
        self.assertEqual(Ni.lattice.beta, 90)
        self.assertEqual(Ni.lattice.gamma, 90)

        # Test case where convs are different from initial and final
        # th = self.rester.get_structure_by_material_id(
        #     "mp-37", conventional_unit_cell=True)
        # th_entry = self.rester.get_entry_by_material_id(
        #     "mp-37", inc_structure=True, conventional_unit_cell=True)
        # th_entry_initial = self.rester.get_entry_by_material_id(
        #     "mp-37", inc_structure="initial", conventional_unit_cell=True)
        # self.assertEqual(th, th_entry.structure)
        # self.assertEqual(len(th_entry.structure), 4)
        # self.assertEqual(len(th_entry_initial.structure), 2)

        # Test if the polymorphs of Fe are properly sorted
        # by e_above_hull when sort_by_e_above_hull=True
        Fe_entries = self.rester.get_entries("Fe", sort_by_e_above_hull=True)
        self.assertEqual(Fe_entries[0].data["e_above_hull"], 0)

    def test_get_pourbaix_entries(self):
        # test input chemsys as a list of elements
        pbx_entries = self.rester.get_pourbaix_entries(["Fe", "Cr"])
        for pbx_entry in pbx_entries:
            self.assertTrue(isinstance(pbx_entry, PourbaixEntry))

        # test input chemsys as a string
        pbx_entries = self.rester.get_pourbaix_entries("Fe-Cr")
        for pbx_entry in pbx_entries:
            self.assertTrue(isinstance(pbx_entry, PourbaixEntry))

        # fe_two_plus = [e for e in pbx_entries if e.entry_id == "ion-0"][0]
        # self.assertAlmostEqual(fe_two_plus.energy, -1.12369, places=3)
        #
        # feo2 = [e for e in pbx_entries if e.entry_id == "mp-25332"][0]
        # self.assertAlmostEqual(feo2.energy, 3.56356, places=3)
        #
        # # Test S, which has Na in reference solids
        # pbx_entries = self.rester.get_pourbaix_entries(["S"])
        # so4_two_minus = pbx_entries[9]
        # self.assertAlmostEqual(so4_two_minus.energy, 0.301511, places=3)

        # Ensure entries are pourbaix compatible
        PourbaixDiagram(pbx_entries)

    def test_get_exp_entry(self):
        entry = self.rester.get_exp_entry("Fe2O3")
        self.assertEqual(entry.energy, -825.5)

    # def test_submit_query_delete_snl(self):
    # s = Structure([[5, 0, 0], [0, 5, 0], [0, 0, 5]], ["Fe"], [[0, 0, 0]])
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
                    ComputedEntry(
                        entry.composition,
                        entry.uncorrected_energy + 0.01,
                        parameters=entry.parameters,
                        entry_id="mod_{}".format(entry.entry_id),
                    )
                )
        rest_ehulls = self.rester.get_stability(modified_entries)
        all_entries = entries + modified_entries
        compat = MaterialsProject2020Compatibility()
        all_entries = compat.process_entries(all_entries)
        pd = PhaseDiagram(all_entries)
        for e in all_entries:
            if str(e.entry_id).startswith("mod"):
                for d in rest_ehulls:
                    if d["entry_id"] == e.entry_id:
                        data = d
                        break
                self.assertAlmostEqual(pd.get_e_above_hull(e), data["e_above_hull"])

    def test_get_reaction(self):
        rxn = self.rester.get_reaction(["Li", "O"], ["Li2O"])
        self.assertIn("Li2O", rxn["Experimental_references"])

    def test_get_substrates(self):
        substrate_data = self.rester.get_substrates("mp-123", 5, [1, 0, 0])
        substrates = [sub_dict["sub_id"] for sub_dict in substrate_data]
        self.assertIn("mp-2534", substrates)

    def test_get_surface_data(self):
        data = self.rester.get_surface_data("mp-126")  # Pt
        one_surf = self.rester.get_surface_data("mp-129", miller_index=[-2, -3, 1])
        self.assertAlmostEqual(one_surf["surface_energy"], 2.99156963, places=2)
        self.assertArrayAlmostEqual(one_surf["miller_index"], [3, 2, 1])
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
        mo_gbs = self.rester.get_gb_data(chemsys="Mo")
        self.assertEqual(len(mo_gbs), 10)
        mo_gbs_s5 = self.rester.get_gb_data(pretty_formula="Mo", sigma=5)
        self.assertEqual(len(mo_gbs_s5), 3)
        mo_s3_112 = self.rester.get_gb_data(
            material_id="mp-129",
            sigma=3,
            gb_plane=[1, -1, -2],
            include_work_of_separation=True,
        )
        self.assertEqual(len(mo_s3_112), 1)
        gb_f = mo_s3_112[0]["final_structure"]
        self.assertArrayAlmostEqual(gb_f.rotation_axis, [1, 1, 0])
        self.assertAlmostEqual(gb_f.rotation_angle, 109.47122, places=4)
        self.assertAlmostEqual(mo_s3_112[0]["gb_energy"], 0.47965, places=2)
        self.assertAlmostEqual(mo_s3_112[0]["work_of_separation"], 6.318144, places=2)
        self.assertIn("Mo24", gb_f.formula)
        hcp_s7 = self.rester.get_gb_data(material_id="mp-87", gb_plane=[0, 0, 0, 1], include_work_of_separation=True)
        self.assertAlmostEqual(hcp_s7[0]["gb_energy"], 1.12, places=2)
        self.assertAlmostEqual(hcp_s7[0]["work_of_separation"], 2.47, places=2)

    def test_get_interface_reactions(self):
        kinks = self.rester.get_interface_reactions("LiCoO2", "Li3PS4")
        self.assertTrue(len(kinks) > 0)
        kink = kinks[0]
        self.assertIn("energy", kink)
        self.assertIn("ratio_atomic", kink)
        self.assertIn("rxn", kink)
        self.assertTrue(isinstance(kink["rxn"], Reaction))
        kinks_open_O = self.rester.get_interface_reactions("LiCoO2", "Li3PS4", open_el="O", relative_mu=-1)
        self.assertTrue(len(kinks_open_O) > 0)
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings("always", message="The reactant.+")
            self.rester.get_interface_reactions("LiCoO2", "MnO9")
            self.assertTrue("The reactant" in str(w[-1].message))

    def test_download_info(self):
        material_ids = ["mp-32800", "mp-23494"]
        task_types = [TaskType.GGA_OPT, TaskType.GGA_UNIFORM]
        file_patterns = ["vasprun*", "OUTCAR*"]
        meta, urls = self.rester.get_download_info(material_ids, task_types=task_types, file_patterns=file_patterns)
        self.assertDictEqual(
            dict(meta),
            {
                "mp-23494": [{"task_id": "mp-1752825", "task_type": "GGA NSCF Uniform"}],
                "mp-32800": [{"task_id": "mp-739635", "task_type": "GGA NSCF Uniform"}],
            },
        )
        prefix = "http://labdev-nomad.esc.rzg.mpg.de/fairdi/nomad/mp/api/raw/query?"
        # previous test
        # ids = 'mp-23494,mp-688563,mp-32800,mp-746913'
        ids = "mp-1752825,mp-739635"
        self.assertEqual(
            urls[0],
            f"{prefix}file_pattern=vasprun*&file_pattern=OUTCAR*&external_id={ids}",
        )

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

        # Let's test some invalid symbols

        self.assertRaises(ValueError, MPRester.parse_criteria, "li-fe")
        self.assertRaises(ValueError, MPRester.parse_criteria, "LO2")

        crit = MPRester.parse_criteria("POPO2")
        self.assertIn("P2O3", crit["pretty_formula"]["$in"])

    def test_include_user_agent(self):
        headers = self.rester.session.headers
        self.assertIn("user-agent", headers, msg="Include user-agent header by default")
        m = re.match(
            r"pymatgen/(\d+)\.(\d+)\.(\d+)\.?(\d+)? \(Python/(\d+)\.(\d)+\.(\d+) ([^\/]*)/([^\)]*)\)",
            headers["user-agent"],
        )
        self.assertIsNotNone(m, msg="Unexpected user-agent value {}".format(headers["user-agent"]))
        self.rester = MPRester(include_user_agent=False)
        self.assertNotIn("user-agent", self.rester.session.headers, msg="user-agent header unwanted")

    def test_database_version(self):

        with MPRester(notify_db_version=True) as mpr:
            db_version = mpr.get_database_version()

        self.assertIsInstance(db_version, str)

        with open(SETTINGS_FILE, "rt") as f:
            d = yaml.safe_load(f)

        self.assertEqual(d["MAPI_DB_VERSION"]["LAST_ACCESSED"], db_version)
        self.assertIsInstance(d["MAPI_DB_VERSION"]["LOG"][db_version], int)

    def test_pourbaix_heavy(self):

        entries = self.rester.get_pourbaix_entries(["Li", "Mg", "Sn", "Pd"])
        pbx = PourbaixDiagram(entries, nproc=4, filter_solids=False)
        entries = self.rester.get_pourbaix_entries(["Ba", "Ca", "V", "Cu", "F"])
        pbx = PourbaixDiagram(entries, nproc=4, filter_solids=False)
        entries = self.rester.get_pourbaix_entries(["Ba", "Ca", "V", "Cu", "F", "Fe"])
        pbx = PourbaixDiagram(entries, nproc=4, filter_solids=False)
        entries = self.rester.get_pourbaix_entries(["Na", "Ca", "Nd", "Y", "Ho", "F"])
        pbx = PourbaixDiagram(entries, nproc=4, filter_solids=False)

    def test_pourbaix_mpr_pipeline(self):

        data = self.rester.get_pourbaix_entries(["Zn"])
        pbx = PourbaixDiagram(data, filter_solids=True, conc_dict={"Zn": 1e-8})
        pbx.find_stable_entry(10, 0)

        data = self.rester.get_pourbaix_entries(["Ag", "Te"])
        pbx = PourbaixDiagram(data, filter_solids=True, conc_dict={"Ag": 1e-8, "Te": 1e-8})
        self.assertEqual(len(pbx.stable_entries), 29)
        test_entry = pbx.find_stable_entry(8, 2)
        self.assertEqual(sorted(test_entry.entry_id), ["ion-10", "mp-499"])

        # Test against ion sets with multiple equivalent ions (Bi-V regression)
        entries = self.rester.get_pourbaix_entries(["Bi", "V"])
        pbx = PourbaixDiagram(entries, filter_solids=True, conc_dict={"Bi": 1e-8, "V": 1e-8})
        self.assertTrue(all(["Bi" in entry.composition and "V" in entry.composition for entry in pbx.all_entries]))


if __name__ == "__main__":
    unittest.main()
