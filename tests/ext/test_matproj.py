from __future__ import annotations

import random
import re
from unittest.mock import patch

import pytest
import requests
from numpy.testing import assert_allclose
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.analysis.pourbaix_diagram import PourbaixDiagram, PourbaixEntry
from pymatgen.analysis.reaction_calculator import Reaction
from pymatgen.analysis.wulff import WulffShape
from pymatgen.core import SETTINGS, Composition, Element, Structure
from pymatgen.electronic_structure.bandstructure import BandStructure, BandStructureSymmLine
from pymatgen.electronic_structure.dos import CompleteDos
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.ext.matproj import MP_LOG_FILE, _MPResterBasic
from pymatgen.ext.matproj_legacy import MPRestError, TaskType, _MPResterLegacy
from pymatgen.phonon.bandstructure import PhononBandStructureSymmLine
from pymatgen.phonon.dos import CompletePhononDos
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest
from pytest import approx
from ruamel.yaml import YAML

PMG_MAPI_KEY = SETTINGS.get("PMG_MAPI_KEY", "")
if (10 < len(PMG_MAPI_KEY) <= 20) and "PMG_MAPI_KEY" in SETTINGS:
    MP_URL = "https://legacy.materialsproject.org"
elif len(PMG_MAPI_KEY) > 20:
    MP_URL = "https://api.materialsproject.org"
else:
    MP_URL = "https://materialsproject.org"
try:
    skip_mprester_tests = requests.get(MP_URL, timeout=600).status_code != 200

except (ModuleNotFoundError, ImportError, requests.exceptions.ConnectionError):
    # Skip all MPRester tests if some downstream problem on the website, mp-api or whatever.
    skip_mprester_tests = True


@pytest.mark.skipif(
    skip_mprester_tests or (not 10 < len(PMG_MAPI_KEY) <= 20),
    reason="Legacy PMG_MAPI_KEY environment variable not set or MP API is down.",
)
class TestMPResterOld(PymatgenTest):
    def setUp(self):
        self.rester = _MPResterLegacy()

    def test_get_all_materials_ids_doc(self):
        mids = self.rester.get_materials_ids("Al2O3")
        random.shuffle(mids)
        doc = self.rester.get_doc(mids.pop(0))
        assert doc["pretty_formula"] == "Al2O3"

    def test_get_xas_data(self):
        # Test getting XAS data
        data = self.rester.get_xas_data("mp-19017", "Li")
        assert data["mid_and_el"] == "mp-19017,Li"
        assert data["spectrum"]["x"][0] == approx(55.178)
        assert data["spectrum"]["y"][0] == approx(0.0164634)

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
        mp_id = "mp-1143"
        vals = requests.get(f"http://legacy.materialsproject.org/materials/{mp_id}/json/", timeout=600)
        expected_vals = vals.json()

        for prop in props:
            if prop not in ["hubbards", "unit_cell_formula", "elements", "icsd_ids", "task_ids"]:
                val = self.rester.get_data(mp_id, prop=prop)[0][prop]
                if prop in ["energy", "energy_per_atom"]:
                    prop = f"final_{prop}"
                assert expected_vals[prop] == approx(val), f"Failed with property {prop}"
            elif prop in ["elements", "icsd_ids", "task_ids"]:
                upstream_vals = set(self.rester.get_data(mp_id, prop=prop)[0][prop])
                assert set(expected_vals[prop]) <= upstream_vals
            else:
                assert expected_vals[prop] == self.rester.get_data(mp_id, prop=prop)[0][prop]

        props = ["structure", "initial_structure", "final_structure", "entry"]
        for prop in props:
            obj = self.rester.get_data(mp_id, prop=prop)[0][prop]
            if prop.endswith("structure"):
                assert isinstance(obj, Structure)
            elif prop == "entry":
                obj = self.rester.get_data(mp_id, prop=prop)[0][prop]
                assert isinstance(obj, ComputedEntry)

        # Test chemsys search
        data = self.rester.get_data("Fe-Li-O", prop="unit_cell_formula")
        assert len(data) > 1
        elements = {Element("Li"), Element("Fe"), Element("O")}
        for d in data:
            assert set(Composition(d["unit_cell_formula"]).elements).issubset(elements)

        with pytest.raises(MPRestError, match="REST query returned with error status code 404"):
            self.rester.get_data("Fe2O3", "bad-method")

    def test_get_materials_id_from_task_id(self):
        assert self.rester.get_materials_id_from_task_id("mp-540081") == "mp-19017"

    def test_get_materials_id_references(self):
        mpr = _MPResterLegacy()
        data = mpr.get_materials_id_references("mp-123")
        assert len(data) > 1000

    def test_find_structure(self):
        mpr = _MPResterLegacy()
        cif_file = f"{TEST_FILES_DIR}/cif/Fe3O4.cif"
        data = mpr.find_structure(str(cif_file))
        assert len(data) > 1
        struct = Structure.from_file(cif_file, primitive=True)
        data = mpr.find_structure(struct)
        assert len(data) > 1

    def test_get_entries_in_chemsys(self):
        syms = ["Li", "Fe", "O"]
        syms2 = "Li-Fe-O"
        entries = self.rester.get_entries_in_chemsys(syms)
        entries2 = self.rester.get_entries_in_chemsys(syms2)
        elements = {Element(sym) for sym in syms}
        for entry in entries:
            assert isinstance(entry, ComputedEntry)
            assert set(entry.elements).issubset(elements)

        e1 = {ent.entry_id for ent in entries}
        e2 = {ent.entry_id for ent in entries2}
        assert e1 == e2

        stable_entries = self.rester.get_entries_in_chemsys(syms, additional_criteria={"e_above_hull": {"$lte": 0.001}})
        assert len(stable_entries) < len(entries)

    def test_get_structure_by_material_id(self):
        s1 = self.rester.get_structure_by_material_id("mp-1")
        assert s1.formula == "Cs1"

        # requesting via task-id instead of mp-id
        with pytest.warns(
            Warning,
            match="calculation task mp-698856 is mapped to canonical mp-id mp-1394, so structure for mp-1394 returned",
        ):
            self.rester.get_structure_by_material_id("mp-698856")

        # requesting unknown mp-id
        # TODO (janosh) this seems like the wrong error message for this case
        with pytest.raises(MPRestError, match="'id' is not a valid Element"):
            self.rester.get_structure_by_material_id("id-does-not-exist")

    def test_get_entry_by_material_id(self):
        entry = self.rester.get_entry_by_material_id("mp-19017")
        assert isinstance(entry, ComputedEntry)
        assert entry.reduced_formula == "LiFePO4"

        with pytest.raises(MPRestError, match="material_id = 'mp-2022' does not exist"):
            self.rester.get_entry_by_material_id("mp-2022")  # "mp-2022" does not exist

    def test_query(self):
        criteria = {"elements": {"$in": ["Li", "Na", "K"], "$all": ["O"]}}
        props = ["pretty_formula", "energy"]
        data = self.rester.query(criteria=criteria, properties=props, chunk_size=0)
        assert len(data) > 6
        data = self.rester.query(criteria="*2O", properties=props, chunk_size=0)
        assert len(data) >= 52
        assert "Li2O" in (d["pretty_formula"] for d in data)

    def test_query_chunk_size(self):
        criteria = {"nelements": 2, "elements": "O"}
        props = ["pretty_formula"]
        data1 = self.rester.query(criteria=criteria, properties=props, chunk_size=0)
        data2 = self.rester.query(criteria=criteria, properties=props, chunk_size=500)
        assert {d["pretty_formula"] for d in data1} == {d["pretty_formula"] for d in data2}
        assert "Al2O3" in {d["pretty_formula"] for d in data1}

    def test_get_exp_thermo_data(self):
        data = self.rester.get_exp_thermo_data("Fe2O3")
        assert len(data) > 0
        for d in data:
            assert d.formula == "Fe2O3"

    def test_get_dos_by_id(self):
        dos = self.rester.get_dos_by_material_id("mp-2254")
        assert isinstance(dos, CompleteDos)

    def test_get_bandstructure_by_material_id(self):
        bs = self.rester.get_bandstructure_by_material_id("mp-2254")
        assert isinstance(bs, BandStructureSymmLine)
        bs_unif = self.rester.get_bandstructure_by_material_id("mp-2254", line_mode=False)
        assert isinstance(bs_unif, BandStructure)
        assert not isinstance(bs_unif, BandStructureSymmLine)

    def test_get_phonon_data_by_material_id(self):
        bs = self.rester.get_phonon_bandstructure_by_material_id("mp-661")
        assert isinstance(bs, PhononBandStructureSymmLine)
        dos = self.rester.get_phonon_dos_by_material_id("mp-661")
        assert isinstance(dos, CompletePhononDos)
        ddb_str = self.rester.get_phonon_ddb_by_material_id("mp-661")
        assert isinstance(ddb_str, str)

    def test_get_structures(self):
        structs = self.rester.get_structures("Mn3O4")
        assert len(structs) > 0

    def test_get_entries(self):
        entries = self.rester.get_entries("TiO2")
        assert len(entries) > 1
        for entry in entries:
            assert entry.reduced_formula == "TiO2"

        entries = self.rester.get_entries("TiO2", inc_structure=True)
        assert len(entries) > 1
        for entry in entries:
            assert entry.structure.reduced_formula == "TiO2"

        # all_entries = self.rester.get_entries("Fe", compatible_only=False)
        # entries = self.rester.get_entries("Fe", compatible_only=True)
        # assert len(entries) < len(all_entries)
        entries = self.rester.get_entries("Fe", compatible_only=True, property_data=["cif"])
        assert "cif" in entries[0].data

        for entry in self.rester.get_entries("CdO2", inc_structure=False):
            assert entry.data["oxide_type"] is not None

        # test if it will retrieve the conventional unit cell of Ni
        entry = self.rester.get_entry_by_material_id("mp-23", inc_structure=True, conventional_unit_cell=True)
        Ni = entry.structure
        assert Ni.lattice.a == Ni.lattice.b
        assert Ni.lattice.a == Ni.lattice.c
        assert Ni.lattice.alpha == 90
        assert Ni.lattice.beta == 90
        assert Ni.lattice.gamma == 90

        # Ensure energy per atom is same
        primNi = self.rester.get_entry_by_material_id("mp-23", inc_structure=True, conventional_unit_cell=False)
        assert primNi.energy_per_atom == entry.energy_per_atom

        Ni = self.rester.get_structure_by_material_id("mp-23", conventional_unit_cell=True)
        assert Ni.lattice.a == Ni.lattice.b
        assert Ni.lattice.a == Ni.lattice.c
        assert Ni.lattice.alpha == 90
        assert Ni.lattice.beta == 90
        assert Ni.lattice.gamma == 90

        # Test case where convs are different from initial and final
        # th = self.rester.get_structure_by_material_id("mp-37", conventional_unit_cell=True)
        # th_entry = self.rester.get_entry_by_material_id("mp-37", inc_structure=True, conventional_unit_cell=True)
        # th_entry_initial = self.rester.get_entry_by_material_id(
        #     "mp-37", inc_structure="initial", conventional_unit_cell=True
        # )
        # assert th == th_entry.structure
        # assert len(th_entry.structure) == 4
        # assert len(th_entry_initial.structure) == 2

        # Test if the polymorphs of Fe are properly sorted
        # by e_above_hull when sort_by_e_above_hull=True
        Fe_entries = self.rester.get_entries("Fe", sort_by_e_above_hull=True)
        assert Fe_entries[0].data["e_above_hull"] == 0

    def test_get_pourbaix_entries(self):
        # test input chemsys as a list of elements
        pbx_entries = self.rester.get_pourbaix_entries(["Fe", "Cr"])
        for pbx_entry in pbx_entries:
            assert isinstance(pbx_entry, PourbaixEntry)

        # test input chemsys as a string
        pbx_entries = self.rester.get_pourbaix_entries("Fe-Cr")
        for pbx_entry in pbx_entries:
            assert isinstance(pbx_entry, PourbaixEntry)

        # fe_two_plus = next(entry for entry in pbx_entries if entry.entry_id == "ion-0")
        # assert fe_two_plus.energy == approx(-1.12369, abs=1e-3)

        # feo2 = next(entry for entry in pbx_entries if entry.entry_id == "mp-25332")
        # assert feo2.energy == approx(3.56356, abs=1e-3)

        # # Test S, which has Na in reference solids
        # pbx_entries = self.rester.get_pourbaix_entries(["S"])
        # so4_two_minus = pbx_entries[9]
        # assert so4_two_minus.energy == approx(0.301511, abs=1e-3)

        # Ensure entries are Pourbaix compatible
        PourbaixDiagram(pbx_entries)

    def test_get_exp_entry(self):
        entry = self.rester.get_exp_entry("Fe2O3")
        assert entry.energy == -825.5

    # def test_submit_query_delete_snl(self):
    #     struct = Structure(np.eye(3) * 5, ["Fe"], [[0, 0, 0]])
    #     submission_ids = self.rester.submit_snl(
    #         [struct, struct], remarks=["unittest"], authors="Test User <test@materialsproject.com>"
    #     )
    #     assert len(submission_ids) == 2
    #     data = self.rester.query_snl({"about.remarks": "unittest"})
    #     assert len(data) == 2
    #     snl_ids = [d["_id"] for d in data]
    #     self.rester.delete_snl(snl_ids)
    #     data = self.rester.query_snl({"about.remarks": "unittest"})
    #     assert len(data) == 0

    def test_get_stability(self):
        entries = self.rester.get_entries_in_chemsys(["Fe", "O"])
        modified_entries = [
            ComputedEntry(
                entry.composition,
                entry.uncorrected_energy + 0.01,
                parameters=entry.parameters,
                entry_id=f"mod_{entry.entry_id}",
            )
            for entry in entries
            if entry.reduced_formula == "Fe2O3"
        ]
        rester_ehulls = self.rester.get_stability(modified_entries)
        all_entries = entries + modified_entries
        compat = MaterialsProject2020Compatibility()
        all_entries = compat.process_entries(all_entries)
        pd = PhaseDiagram(all_entries)
        for entry in all_entries:
            if str(entry.entry_id).startswith("mod"):
                for dct in rester_ehulls:
                    if dct["entry_id"] == entry.entry_id:
                        data = dct
                        break
                assert pd.get_e_above_hull(entry) == approx(data["e_above_hull"])

    def test_get_reaction(self):
        rxn = self.rester.get_reaction(["Li", "O"], ["Li2O"])
        assert "Li2O" in rxn["Experimental_references"]

    def test_get_substrates(self):
        substrate_data = self.rester.get_substrates("mp-123", 5, [1, 0, 0])
        substrates = [sub_dict["sub_id"] for sub_dict in substrate_data]
        assert "mp-2534" in substrates

    def test_get_surface_data(self):
        data = self.rester.get_surface_data("mp-126")  # Pt
        one_surf = self.rester.get_surface_data("mp-129", miller_index=[-2, -3, 1])
        assert one_surf["surface_energy"] == approx(2.99156963)
        assert_allclose(one_surf["miller_index"], [3, 2, 1])
        assert "surfaces" in data
        surfaces = data["surfaces"]
        assert len(surfaces) > 0
        surface = surfaces.pop()
        assert "miller_index" in surface
        assert "surface_energy" in surface
        assert "is_reconstructed" in surface
        data_inc = self.rester.get_surface_data("mp-126", inc_structures=True)
        assert "structure" in data_inc["surfaces"][0]

    def test_get_wulff_shape(self):
        ws = self.rester.get_wulff_shape("mp-126")
        assert isinstance(ws, WulffShape)

    def test_get_cohesive_energy(self):
        e_coh = self.rester.get_cohesive_energy("mp-13")
        assert e_coh, 5.04543279

    def test_get_gb_data(self):
        mo_gbs = self.rester.get_gb_data(chemsys="Mo")
        assert len(mo_gbs) == 10
        mo_gbs_s5 = self.rester.get_gb_data(pretty_formula="Mo", sigma=5)
        assert len(mo_gbs_s5) == 3
        mo_s3_112 = self.rester.get_gb_data(
            material_id="mp-129",
            sigma=3,
            gb_plane=[1, -1, -2],
            include_work_of_separation=True,
        )
        assert len(mo_s3_112) == 1
        gb_f = mo_s3_112[0]["final_structure"]
        assert_allclose(gb_f.rotation_axis, [1, 1, 0])
        assert gb_f.rotation_angle == approx(109.47122)
        assert mo_s3_112[0]["gb_energy"] == approx(0.47965, rel=1e-4)
        assert mo_s3_112[0]["work_of_separation"] == approx(6.318144)
        assert "Mo24" in gb_f.formula
        hcp_s7 = self.rester.get_gb_data(material_id="mp-87", gb_plane=[0, 0, 0, 1], include_work_of_separation=True)
        assert hcp_s7[0]["gb_energy"] == approx(1.1206, rel=1e-4)
        assert hcp_s7[0]["work_of_separation"] == approx(2.4706, rel=1e-4)

    def test_get_interface_reactions(self):
        kinks = self.rester.get_interface_reactions("LiCoO2", "Li3PS4")
        assert len(kinks) > 0
        kink = kinks[0]
        assert "energy" in kink
        assert "ratio_atomic" in kink
        assert "rxn" in kink
        assert isinstance(kink["rxn"], Reaction)
        kinks_open_O = self.rester.get_interface_reactions("LiCoO2", "Li3PS4", open_el="O", relative_mu=-1)
        assert len(kinks_open_O) > 0
        with pytest.warns(
            UserWarning,
            match="The reactant MnO9 has no matching entry with negative formation energy, "
            "instead convex hull energy for this composition will be used for reaction energy calculation.",
        ):
            self.rester.get_interface_reactions("LiCoO2", "MnO9")

    def test_download_info(self):
        material_ids = ["mvc-2970"]
        task_types = [TaskType.GGA_OPT, TaskType.GGAU_UNIFORM]
        file_patterns = ["vasprun*", "OUTCAR*"]
        meta, urls = self.rester.get_download_info(material_ids, task_types=task_types, file_patterns=file_patterns)
        assert dict(meta) == {
            "mvc-2970": [{"task_id": "mp-1738602", "task_type": "GGA+U NSCF Uniform"}],
        }
        assert (
            urls[0] == "https://nomad-lab.eu/prod/rae/api/raw/query?file_pattern=vasprun*"
            "&file_pattern=OUTCAR*&external_id=mp-1738602"
        )

    def test_parse_criteria(self):
        crit = _MPResterLegacy.parse_criteria("mp-1234 Li-*")
        assert "Li-O" in crit["$or"][1]["chemsys"]["$in"]
        assert {"task_id": "mp-1234"} in crit["$or"]

        crit = _MPResterLegacy.parse_criteria("Li2*")
        assert "Li2O" in crit["pretty_formula"]["$in"]
        assert "Li2I" in crit["pretty_formula"]["$in"]
        assert "CsLi2" in crit["pretty_formula"]["$in"]

        crit = _MPResterLegacy.parse_criteria("Li-*-*")
        assert "Li-Re-Ru" in crit["chemsys"]["$in"]
        assert "Li-Li" not in crit["chemsys"]["$in"]

        comps = _MPResterLegacy.parse_criteria("**O3")["pretty_formula"]["$in"]
        for comp in comps:
            assert len(Composition(comp)) == 3, f"Failed in {comp}"

        chemsys = _MPResterLegacy.parse_criteria("{Fe,Mn}-O")["chemsys"]["$in"]
        assert len(chemsys) == 2
        comps = _MPResterLegacy.parse_criteria("{Fe,Mn,Co}O")["pretty_formula"]["$in"]
        assert len(comps) == 3, comps

        # Let's test some invalid symbols

        with pytest.raises(ValueError, match="'li' is not a valid Element"):
            _MPResterLegacy.parse_criteria("li-fe")
        with pytest.raises(ValueError, match="'L' is not a valid Element"):
            _MPResterLegacy.parse_criteria("LO2")

        crit = _MPResterLegacy.parse_criteria("POPO2")
        assert "P2O3" in crit["pretty_formula"]["$in"]

    def test_include_user_agent(self):
        pytest.skip(
            "this test started failing with 'pymatgen.ext.matproj.MPRestError: REST query "
            "returned with error status code 403. Content: b'error code: 1020'"
        )
        headers = self.rester.session.headers
        assert "user-agent" in headers, "Include user-agent header by default"
        match = re.match(
            r"pymatgen/(\d+)\.(\d+)\.(\d+)\.?(\d+)? \(Python/(\d+)\.(\d)+\.(\d+) ([^\/]*)/([^\)]*)\)",
            headers["user-agent"],
        )
        assert match is not None, f"Unexpected user-agent value {headers['user-agent']}"
        self.rester = _MPResterLegacy(include_user_agent=False)
        assert "user-agent" not in self.rester.session.headers, "user-agent header unwanted"

    def test_database_version(self):
        with _MPResterLegacy(notify_db_version=True) as mpr:
            db_version = mpr.get_database_version()

        assert isinstance(db_version, str)
        yaml = YAML()
        with open(MP_LOG_FILE) as file:
            dct = yaml.load(file)

        assert dct["MAPI_DB_VERSION"]["LAST_ACCESSED"] == db_version
        assert isinstance(dct["MAPI_DB_VERSION"]["LOG"][db_version], int)

    def test_pourbaix_heavy(self):
        entries = self.rester.get_pourbaix_entries(["Na", "Ca", "Nd", "Y", "Ho", "F"])
        _ = PourbaixDiagram(entries, nproc=4, filter_solids=False)

    def test_pourbaix_mpr_pipeline(self):
        data = self.rester.get_pourbaix_entries(["Zn"])
        pbx = PourbaixDiagram(data, filter_solids=True, conc_dict={"Zn": 1e-8})
        pbx.find_stable_entry(10, 0)

        data = self.rester.get_pourbaix_entries(["Ag", "Te"])
        pbx = PourbaixDiagram(data, filter_solids=True, conc_dict={"Ag": 1e-8, "Te": 1e-8})
        assert len(pbx.stable_entries) == 30
        test_entry = pbx.find_stable_entry(8, 2)
        assert sorted(test_entry.entry_id) == ["ion-10", "mp-996958"]

        # Test against ion sets with multiple equivalent ions (Bi-V regression)
        entries = self.rester.get_pourbaix_entries(["Bi", "V"])
        pbx = PourbaixDiagram(entries, filter_solids=True, conc_dict={"Bi": 1e-8, "V": 1e-8})
        assert all("Bi" in entry.composition and "V" in entry.composition for entry in pbx.all_entries)

    @patch.dict(SETTINGS, {"PMG_MAPI_KEY": "foobar"})
    def test_api_key_is_none(self):
        # https://github.com/materialsproject/pymatgen/pull/3004
        with _MPResterLegacy(None) as mpr:
            assert mpr.api_key == "foobar"

        with _MPResterLegacy(api_key=None) as mpr:
            assert mpr.api_key == "foobar"


@pytest.mark.skipif(
    skip_mprester_tests or (not len(PMG_MAPI_KEY) > 20),
    reason="PMG_MAPI_KEY environment variable not set or MP API is down.",
)
class TestMPResterNewBasic(PymatgenTest):
    def setUp(self):
        self.rester = _MPResterBasic()

    def test_get_summary(self):
        docs = self.rester.get_summary({"formula": "Fe2O3"})
        assert len(docs) > 3

        mid = "mp-19770"
        doc = self.rester.get_summary_by_material_id(mid)
        assert doc["formula_pretty"] == "Fe2O3"

        doc = self.rester.summary.search(material_ids="mp-19770,mp-19017", _fields="formula_pretty,energy_above_hull")
        assert len(doc) == 2
        assert len(doc[0]) == 2
        assert doc[0]["energy_above_hull"] >= 0
        assert doc[1]["energy_above_hull"] >= 0

        # dos = self.rester.get_dos_by_material_id(mid)
        #
        # assert isinstance(dos, CompleteDos)

    def test_get_all_materials_ids_doc(self):
        mids = self.rester.get_material_ids("Al2O3")
        random.shuffle(mids)
        doc = self.rester.get_doc(mids.pop(0))
        assert doc["formula_pretty"] == "Al2O3"

    #
    # def test_get_xas_data(self):
    #     # Test getting XAS data
    #     data = self.rester.get_xas_data("mp-19017", "Li")
    #     assert data["mid_and_el"] == "mp-19017,Li"
    #     assert data["spectrum"]["x"][0] == approx(55.178)
    #     assert data["spectrum"]["y"][0] == approx(0.0164634)
    #
    # def test_get_data(self):
    #     props = {
    #         "energy",
    #         "energy_per_atom",
    #         "formation_energy_per_atom",
    #         "nsites",
    #         "unit_cell_formula",
    #         "pretty_formula",
    #         "is_hubbard",
    #         "elements",
    #         "nelements",
    #         "e_above_hull",
    #         "hubbards",
    #         "is_compatible",
    #         "task_ids",
    #         "density",
    #         "icsd_ids",
    #         "total_magnetization",
    #     }
    #     mp_id = "mp-1143"
    #     vals = requests.get(f"http://legacy.materialsproject.org/materials/{mp_id}/json/", timeout=600)
    #     expected_vals = vals.json()
    #
    #     for prop in props:
    #         if prop not in [
    #             "hubbards",
    #             "unit_cell_formula",
    #             "elements",
    #             "icsd_ids",
    #             "task_ids",
    #         ]:
    #             val = self.rester.get_data(mp_id, prop=prop)[0][prop]
    #             if prop in ["energy", "energy_per_atom"]:
    #                 prop = "final_" + prop
    #             assert expected_vals[prop] == approx(val), f"Failed with property {prop}"
    #         elif prop in ["elements", "icsd_ids", "task_ids"]:
    #             upstream_vals = set(self.rester.get_data(mp_id, prop=prop)[0][prop])
    #             assert set(expected_vals[prop]) <= upstream_vals
    #         else:
    #             assert expected_vals[prop] == self.rester.get_data(mp_id, prop=prop)[0][prop]
    #
    #     props = ["structure", "initial_structure", "final_structure", "entry"]
    #     for prop in props:
    #         obj = self.rester.get_data(mp_id, prop=prop)[0][prop]
    #         if prop.endswith("structure"):
    #             assert isinstance(obj, Structure)
    #         elif prop == "entry":
    #             obj = self.rester.get_data(mp_id, prop=prop)[0][prop]
    #             assert isinstance(obj, ComputedEntry)
    #
    #     # Test chemsys search
    #     data = self.rester.get_data("Fe-Li-O", prop="unit_cell_formula")
    #     assert len(data) > 1
    #     elements = {Element("Li"), Element("Fe"), Element("O")}
    #     for d in data:
    #         assert set(Composition(d["unit_cell_formula"]).elements).issubset(elements)
    #
    #     with pytest.raises(MPRestError, match="REST query returned with error status code 404"):
    #         self.rester.get_data("Fe2O3", "badmethod")
    #
    # def test_get_materials_id_from_task_id(self):
    #     assert self.rester.get_materials_id_from_task_id("mp-540081") == "mp-19017"
    #
    # def test_get_materials_id_references(self):
    #     mpr = _MPResterLegacy()
    #     data = mpr.get_materials_id_references("mp-123")
    #     assert len(data) > 1000

    def test_get_entries_and_in_chemsys(self):
        # One large system test.
        syms = ["Li", "Fe", "O", "P", "Mn"]

        # Small test.
        syms2 = "Fe-Li-O"
        entries = self.rester.get_entries_in_chemsys(syms)
        entries2 = self.rester.get_entries(syms2)
        elements = {Element(sym) for sym in syms}
        for entry in entries:
            assert isinstance(entry, ComputedEntry)
            assert set(entry.elements).issubset(elements)

        assert len(entries) > 1000

        for entry in entries2:
            assert isinstance(entry, ComputedEntry)
            assert set(entry.elements).issubset(elements)
        assert len(entries2) < 1000

        e1 = {i.entry_id for i in entries}
        e2 = {i.entry_id for i in entries2}
        assert e1.issuperset(e2)

    def test_get_structure_by_material_id(self):
        s1 = self.rester.get_structure_by_material_id("mp-1")
        assert s1.formula == "Cs1"

    def test_get_entry_by_material_id(self):
        entry = self.rester.get_entry_by_material_id("mp-19017")
        assert isinstance(entry, ComputedEntry)
        assert entry.reduced_formula == "LiFePO4"

        with pytest.raises(IndexError, match="list index out of range"):
            self.rester.get_entry_by_material_id("mp-2022")  # "mp-2022" does not exist

    # def test_query(self):
    #     criteria = {"elements": {"$in": ["Li", "Na", "K"], "$all": ["O"]}}
    #     props = ["pretty_formula", "energy"]
    #     data = self.rester.query(criteria=criteria, properties=props, chunk_size=0)
    #     assert len(data) > 6
    #     data = self.rester.query(criteria="*2O", properties=props, chunk_size=0)
    #     assert len(data) >= 52
    #     assert "Li2O" in (d["pretty_formula"] for d in data)

    # def test_get_exp_thermo_data(self):
    #     data = self.rester.get_exp_thermo_data("Fe2O3")
    #     assert len(data) > 0
    #     for d in data:
    #         assert d.formula == "Fe2O3"
    #
    # def test_get_dos_by_id(self):
    #     dos = self.rester.get_dos_by_material_id("mp-2254")
    #     assert isinstance(dos, CompleteDos)

    # def test_get_bandstructure_by_material_id(self):
    #     bs = self.rester.get_bandstructure_by_material_id("mp-2254")
    #     assert isinstance(bs, BandStructureSymmLine)
    #     bs_unif = self.rester.get_bandstructure_by_material_id("mp-2254", line_mode=False)
    #     assert isinstance(bs_unif, BandStructure)
    #     assert not isinstance(bs_unif, BandStructureSymmLine)
    #
    def test_get_phonon_data_by_material_id(self):
        bs = self.rester.get_phonon_bandstructure_by_material_id("mp-661")
        assert isinstance(bs, PhononBandStructureSymmLine)
        dos = self.rester.get_phonon_dos_by_material_id("mp-661")
        assert isinstance(dos, CompletePhononDos)

    #     ddb_str = self.rester.get_phonon_ddb_by_material_id("mp-661")
    #     assert isinstance(ddb_str, str)

    def test_get_structures(self):
        structs = self.rester.get_structures("Mn3O4")
        assert len(structs) > 0

    # def test_get_entries(self):
    #     entries = self.rester.get_entries("TiO2")
    #     assert len(entries) > 1
    #     for entry in entries:
    #         assert entry.reduced_formula == "TiO2"
    #
    #     entries = self.rester.get_entries("TiO2", inc_structure=True)
    #     assert len(entries) > 1
    #     for entry in entries:
    #         assert entry.structure.reduced_formula == "TiO2"

    # # all_entries = self.rester.get_entries("Fe", compatible_only=False)
    # # entries = self.rester.get_entries("Fe", compatible_only=True)
    # # assert len(entries) < len(all_entries)
    # entries = self.rester.get_entries("Fe", compatible_only=True, property_data=["cif"])
    # assert "cif" in entries[0].data
    #
    # for entry in self.rester.get_entries("CdO2", inc_structure=False):
    #     assert entry.data["oxide_type"] is not None
    #
    # # test if it will retrieve the conventional unit cell of Ni
    # entry = self.rester.get_entry_by_material_id("mp-23", inc_structure=True, conventional_unit_cell=True)
    # Ni = entry.structure
    # assert Ni.lattice.a == Ni.lattice.b
    # assert Ni.lattice.a == Ni.lattice.c
    # assert Ni.lattice.alpha == 90
    # assert Ni.lattice.beta == 90
    # assert Ni.lattice.gamma == 90
    #
    # # Ensure energy per atom is same
    # primNi = self.rester.get_entry_by_material_id("mp-23", inc_structure=True, conventional_unit_cell=False)
    # assert primNi.energy_per_atom == entry.energy_per_atom
    #
    # Ni = self.rester.get_structure_by_material_id("mp-23", conventional_unit_cell=True)
    # assert Ni.lattice.a == Ni.lattice.b
    # assert Ni.lattice.a == Ni.lattice.c
    # assert Ni.lattice.alpha == 90
    # assert Ni.lattice.beta == 90
    # assert Ni.lattice.gamma == 90
    #
    # # Test case where convs are different from initial and final
    # # th = self.rester.get_structure_by_material_id("mp-37", conventional_unit_cell=True)
    # # th_entry = self.rester.get_entry_by_material_id("mp-37", inc_structure=True, conventional_unit_cell=True)
    # # th_entry_initial = self.rester.get_entry_by_material_id(
    # #     "mp-37", inc_structure="initial", conventional_unit_cell=True
    # # )
    # # assert th == th_entry.structure
    # # assert len(th_entry.structure) == 4
    # # assert len(th_entry_initial.structure) == 2
    #
    # # Test if the polymorphs of Fe are properly sorted
    # # by e_above_hull when sort_by_e_above_hull=True
    # Fe_entries = self.rester.get_entries("Fe", sort_by_e_above_hull=True)
    # assert Fe_entries[0].data["e_above_hull"] == 0

    #
    # def test_get_exp_entry(self):
    #     entry = self.rester.get_exp_entry("Fe2O3")
    #     assert entry.energy == -825.5
    #
    # def test_get_stability(self):
    #     entries = self.rester.get_entries_in_chemsys(["Fe", "O"])
    #     modified_entries = [
    #         ComputedEntry(
    #             entry.composition,
    #             entry.uncorrected_energy + 0.01,
    #             parameters=entry.parameters,
    #             entry_id=f"mod_{entry.entry_id}",
    #         )
    #         for entry in entries
    #         if entry.reduced_formula == "Fe2O3"
    #     ]
    #     rester_ehulls = self.rester.get_stability(modified_entries)
    #     all_entries = entries + modified_entries
    #     compat = MaterialsProject2020Compatibility()
    #     all_entries = compat.process_entries(all_entries)
    #     pd = PhaseDiagram(all_entries)
    #     for entry in all_entries:
    #         if str(entry.entry_id).startswith("mod"):
    #             for dct in rester_ehulls:
    #                 if dct["entry_id"] == entry.entry_id:
    #                     data = dct
    #                     break
    #             assert pd.get_e_above_hull(entry) == approx(data["e_above_hull"])

    # def test_get_reaction(self):
    #     rxn = self.rester.get_reaction(["Li", "O"], ["Li2O"])
    #     assert "Li2O" in rxn["Experimental_references"]

    # def test_get_substrates(self):
    #     substrate_data = self.rester.get_substrates("mp-123", 5, [1, 0, 0])
    #     substrates = [sub_dict["sub_id"] for sub_dict in substrate_data]
    #     assert "mp-2534" in substrates

    # def test_get_surface_data(self):
    #     data = self.rester.get_surface_data("mp-126")  # Pt
    #     one_surf = self.rester.get_surface_data("mp-129", miller_index=[-2, -3, 1])
    #     assert one_surf["surface_energy"] == approx(2.99156963)
    #     assert_allclose(one_surf["miller_index"], [3, 2, 1])
    #     assert "surfaces" in data
    #     surfaces = data["surfaces"]
    #     assert len(surfaces) > 0
    #     surface = surfaces.pop()
    #     assert "miller_index" in surface
    #     assert "surface_energy" in surface
    #     assert "is_reconstructed" in surface
    #     data_inc = self.rester.get_surface_data("mp-126", inc_structures=True)
    #     assert "structure" in data_inc["surfaces"][0]
    #
    # def test_get_wulff_shape(self):
    #     ws = self.rester.get_wulff_shape("mp-126")
    #     assert isinstance(ws, WulffShape)
    #
    # def test_get_cohesive_energy(self):
    #     ecoh = self.rester.get_cohesive_energy("mp-13")
    #     assert ecoh, 5.04543279

    # def test_get_gb_data(self):
    #     mo_gbs = self.rester.get_gb_data(chemsys="Mo")
    #     assert len(mo_gbs) == 10
    #     mo_gbs_s5 = self.rester.get_gb_data(pretty_formula="Mo", sigma=5)
    #     assert len(mo_gbs_s5) == 3
    #     mo_s3_112 = self.rester.get_gb_data(
    #         material_id="mp-129",
    #         sigma=3,
    #         gb_plane=[1, -1, -2],
    #         include_work_of_separation=True,
    #     )
    #     assert len(mo_s3_112) == 1
    #     gb_f = mo_s3_112[0]["final_structure"]
    #     assert_allclose(gb_f.rotation_axis, [1, 1, 0])
    #     assert gb_f.rotation_angle == approx(109.47122)
    #     assert mo_s3_112[0]["gb_energy"] == approx(0.47965, rel=1e-4)
    #     assert mo_s3_112[0]["work_of_separation"] == approx(6.318144)
    #     assert "Mo24" in gb_f.formula
    #     hcp_s7 = self.rester.get_gb_data(material_id="mp-87", gb_plane=[0, 0, 0, 1], include_work_of_separation=True)
    #     assert hcp_s7[0]["gb_energy"] == approx(1.1206, rel=1e-4)
    #     assert hcp_s7[0]["work_of_separation"] == approx(2.4706, rel=1e-4)

    # def test_get_interface_reactions(self):
    #     kinks = self.rester.get_interface_reactions("LiCoO2", "Li3PS4")
    #     assert len(kinks) > 0
    #     kink = kinks[0]
    #     assert "energy" in kink
    #     assert "ratio_atomic" in kink
    #     assert "rxn" in kink
    #     assert isinstance(kink["rxn"], Reaction)
    #     kinks_open_O = self.rester.get_interface_reactions("LiCoO2", "Li3PS4", open_el="O", relative_mu=-1)
    #     assert len(kinks_open_O) > 0
    #     with pytest.warns(
    #         UserWarning,
    #         match="The reactant MnO9 has no matching entry with negative formation energy, "
    #         "instead convex hull energy for this composition will be used for reaction energy calculation.",
    #     ):
    #         self.rester.get_interface_reactions("LiCoO2", "MnO9")

    #
    # def test_pourbaix_heavy(self):
    #     entries = self.rester.get_pourbaix_entries(["Na", "Ca", "Nd", "Y", "Ho", "F"])
    #     _ = PourbaixDiagram(entries, nproc=4, filter_solids=False)
    #
    # def test_pourbaix_mpr_pipeline(self):
    #     data = self.rester.get_pourbaix_entries(["Zn"])
    #     pbx = PourbaixDiagram(data, filter_solids=True, conc_dict={"Zn": 1e-8})
    #     pbx.find_stable_entry(10, 0)
    #
    #     data = self.rester.get_pourbaix_entries(["Ag", "Te"])
    #     pbx = PourbaixDiagram(data, filter_solids=True, conc_dict={"Ag": 1e-8, "Te": 1e-8})
    #     assert len(pbx.stable_entries) == 30
    #     test_entry = pbx.find_stable_entry(8, 2)
    #     assert sorted(test_entry.entry_id) == ["ion-10", "mp-996958"]
    #
    #     # Test against ion sets with multiple equivalent ions (Bi-V regression)
    #     entries = self.rester.get_pourbaix_entries(["Bi", "V"])
    #     pbx = PourbaixDiagram(entries, filter_solids=True, conc_dict={"Bi": 1e-8, "V": 1e-8})
    #     assert all("Bi" in entry.composition and "V" in entry.composition for entry in pbx.all_entries)

    def test_parity_with_mp_api(self):
        try:
            from mp_api.client import MPRester as MpApi
        except Exception:
            pytest.skip("mp_api.client.MPRester cannot be imported for this test.")
        mpr_mp_api = MpApi(PMG_MAPI_KEY)
        # Test summary
        mp_data = mpr_mp_api.summary.search(formula="Al2O3")
        pmg_data = self.rester.get_summary({"formula": "Al2O3"})
        assert len(mp_data) == len(pmg_data)

        # Test get_entries
        chemsys = ["Li", "Fe", "O"]
        mp_entries = mpr_mp_api.get_entries_in_chemsys(chemsys)
        pmg_entries = self.rester.get_entries_in_chemsys(chemsys)

        assert len(mp_entries) == len(pmg_entries)
