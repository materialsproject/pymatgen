from __future__ import annotations

import numpy as np
import pytest
import requests
from pytest import approx

from pymatgen.core import SETTINGS, Element
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.ext.matproj import MPRester
from pymatgen.phonon.bandstructure import PhononBandStructureSymmLine
from pymatgen.phonon.dos import CompletePhononDos
from pymatgen.util.testing import MatSciTest

PMG_MAPI_KEY: str = SETTINGS.get("PMG_MAPI_KEY", "")

MP_URL = "https://api.materialsproject.org"

# Skip all MPRester tests if some downstream problem on the website, mp-api or whatever.
try:
    skip_mprester_tests = requests.get(MP_URL, timeout=60).status_code != 200
except (ModuleNotFoundError, ImportError, requests.exceptions.ConnectionError):
    skip_mprester_tests = True

if skip_mprester_tests:
    pytest.skip("MP API is down", allow_module_level=True)


@pytest.mark.skipif(
    not PMG_MAPI_KEY,
    reason="PMG_MAPI_KEY environment variable not set.",
)
class TestMPRester(MatSciTest):
    def setup_method(self):
        self.rester = MPRester()

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

    def test_get_all_materials_ids_doc(self):
        mids = self.rester.get_material_ids("Al2O3")
        np.random.default_rng().shuffle(mids)
        doc = self.rester.get_doc(mids.pop(0))
        assert doc["formula_pretty"] == "Al2O3"

    def test_get_entries_and_in_chemsys(self):
        # One large system test.
        syms = ["Li", "Fe", "O", "P", "Mn"]

        # Small test.
        syms2 = "Fe-Li-O"
        entries = self.rester.get_entries_in_chemsys(syms)

        entries2 = self.rester.get_entries(syms2, property_data=["band_gap"])
        elements = {Element(sym) for sym in syms}
        for entry in entries:
            assert isinstance(entry, ComputedEntry)
            assert set(entry.elements).issubset(elements)

        assert len(entries) > 1000

        # This gets everything in Li-Fe-O and within this subsystem. It should have more entries than
        # get_entries("Li-Fe-O"), which just gets only the ternary compounds.
        entries3 = self.rester.get_entries_in_chemsys(["Fe", "Li", "O"])
        assert len(entries3) > len(entries2)
        for entry in entries2:
            assert isinstance(entry, ComputedEntry)
            assert set(entry.elements).issubset(elements)
            assert "band_gap" in entry.data
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

    def test_get_phonon_data_by_material_id(self):
        bs = self.rester.get_phonon_bandstructure_by_material_id("mp-661")
        assert isinstance(bs, PhononBandStructureSymmLine)
        dos = self.rester.get_phonon_dos_by_material_id("mp-661")
        assert isinstance(dos, CompletePhononDos)

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

        all_entries = self.rester.get_entries("Fe", compatible_only=False)
        entries = self.rester.get_entries("Fe", compatible_only=True)
        assert len(entries) < len(all_entries)

        for entry in self.rester.get_entries("CdO2", inc_structure=False):
            assert entry.data["oxide_type"] is not None

        # test if it will retrieve the conventional unit cell of Ni
        entry = self.rester.get_entry_by_material_id("mp-23", inc_structure=True, conventional_unit_cell=True)
        Ni = entry.structure
        assert Ni.lattice.a == Ni.lattice.b
        assert Ni.lattice.a == Ni.lattice.c
        assert Ni.lattice.alpha == approx(60)
        assert Ni.lattice.beta == approx(60)
        assert Ni.lattice.gamma == approx(60)

        # Ensure energy per atom is same
        primNi = self.rester.get_entry_by_material_id("mp-23", inc_structure=True, conventional_unit_cell=False)
        assert primNi.energy_per_atom == entry.energy_per_atom

        Ni = self.rester.get_structure_by_material_id("mp-23", conventional_unit_cell=True)
        assert Ni.lattice.a == Ni.lattice.b
        assert Ni.lattice.a == Ni.lattice.c
        assert Ni.lattice.alpha == approx(90)
        assert Ni.lattice.beta == approx(90)
        assert Ni.lattice.gamma == approx(90)

    def test_api_parity(self):
        docs = [
            "summary",
            "core",
            "elasticity",
            "phonon",
            "eos",
            "similarity",
            "xas",
            "electronic_structure",
            "robocrys",
            "synthesis",
            "magnetism",
            "oxidation_states",
            "provenance",
            "alloys",
            "chemenv",
            "bonds",
            "dielectric",
        ]

        for doc in docs:
            # We should have Al2O3 data for these properties.
            data = self.rester.materials.__getattribute__(doc).search(material_ids="mp-1143")
            assert len(data) > 0, f"No Al2O3 data returned for {doc}"
            data = self.rester.__getattribute__(doc).search(material_ids="mp-1143")
            assert len(data) > 0, f"No Al2O3 data returned for {doc}"

        data = self.rester.materials.substrates.search(sub_id="mp-1143")
        assert len(data) > 0, "No substrate data returned."

        data = self.rester.materials.tasks.search(task_ids="mp-1143")
        assert len(data) > 0, "No tasks data returned."

        docs = ["surface_properties", "alloys"]

        for doc in docs:
            data = self.rester.materials.__getattribute__(doc).search(material_ids="mp-135")
            assert len(data) > 0, f"No Li data returned for {doc}"

        data = self.rester.materials.insertion_electrodes.search(formula="LiFePO4")
        assert len(data) > 0, "No insertion electrode data returned."

        # TODO: Test these docs: "grain_boundaries", "piezoelectric", "conversion_electrodes"
