from __future__ import annotations

import copy
import json
import unittest
from collections import defaultdict

import pytest
from monty.json import MontyDecoder
from pytest import approx

from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
from pymatgen.entries.computed_entries import (
    CompositionEnergyAdjustment,
    ComputedEntry,
    ComputedStructureEntry,
    ConstantEnergyAdjustment,
    EnergyAdjustment,
    GibbsComputedStructureEntry,
    ManualEnergyAdjustment,
    TemperatureEnergyAdjustment,
)
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.util.testing import TEST_FILES_DIR

filepath = f"{TEST_FILES_DIR}/vasprun.xml"
vasp_run = Vasprun(filepath)


def test_energy_adjustment():
    ea = EnergyAdjustment(10)
    assert ea.name == "Manual adjustment"
    assert not ea.cls
    ea_dct = ea.as_dict()
    ea2 = EnergyAdjustment.from_dict(ea_dct)
    assert str(ea_dct) == str(ea2.as_dict())


def test_energy_adjustment_repr():
    comp_cls = MaterialsProject2020Compatibility()
    cls_name = type(comp_cls).__name__
    for cls, label in ((None, "unknown"), (comp_cls, cls_name), ({"@class": cls_name}, cls_name)):
        ea = EnergyAdjustment(10, cls=cls)
        assert (
            repr(ea) == "EnergyAdjustment(name='Manual adjustment', value=10.0, uncertainty=nan, description=, "
            f"generated_by='{label}')"
        )


def test_manual_energy_adjustment():
    ea = ManualEnergyAdjustment(10)
    assert ea.name == "Manual energy adjustment"
    assert ea.value == 10
    assert ea.explain == "Manual energy adjustment (10.000 eV)"
    ea_dct = ea.as_dict()
    ea2 = ManualEnergyAdjustment.from_dict(ea_dct)
    assert str(ea_dct) == str(ea2.as_dict())


def test_constant_energy_adjustment():
    ea = ConstantEnergyAdjustment(8)
    assert ea.name == "Constant energy adjustment"
    assert ea.value == 8
    assert ea.explain == "Constant energy adjustment (8.000 eV)"
    ea_dct = ea.as_dict()
    ea2 = ConstantEnergyAdjustment.from_dict(ea_dct)
    assert str(ea_dct) == str(ea2.as_dict())


def test_composition_energy_adjustment():
    ea = CompositionEnergyAdjustment(2, 2, uncertainty_per_atom=0, name="H")
    assert ea.name == "H"
    assert ea.value == 4
    assert ea.explain == "Composition-based energy adjustment (2.000 eV/atom x 2 atoms)"
    ea_dct = ea.as_dict()
    ea2 = CompositionEnergyAdjustment.from_dict(ea_dct)
    assert str(ea_dct) == str(ea2.as_dict())


def test_temp_energy_adjustment():
    ea = TemperatureEnergyAdjustment(-0.1, 298, 5, uncertainty_per_deg=0, name="entropy")
    assert ea.name == "entropy"
    assert ea.value == -0.1 * 298 * 5
    assert ea.n_atoms == 5
    assert ea.temp == 298
    assert ea.explain == "Temperature-based energy adjustment (-0.1000 eV/K/atom x 298 K x 5 atoms)"
    ea_dct = ea.as_dict()
    ea2 = TemperatureEnergyAdjustment.from_dict(ea_dct)
    assert str(ea_dct) == str(ea2.as_dict())


class TestComputedEntry(unittest.TestCase):
    def setUp(self):
        self.entry = ComputedEntry(
            vasp_run.final_structure.composition,
            vasp_run.final_energy,
            parameters=vasp_run.incar,
        )
        self.entry2 = ComputedEntry({"Fe": 2, "O": 3}, 2.3)
        self.entry3 = ComputedEntry("Fe2O3", 2.3)
        self.entry4 = ComputedEntry("Fe2O3", 2.3, entry_id=1)
        self.entry5 = ComputedEntry("Fe6O9", 6.9)
        ea = ConstantEnergyAdjustment(-5, name="Dummy adjustment")
        self.entry6 = ComputedEntry("Fe6O9", 6.9, correction=-10)
        self.entry7 = ComputedEntry("Fe6O9", 6.9, energy_adjustments=[ea])

    def test_energy(self):
        assert self.entry.energy == approx(-269.38319884)
        self.entry.correction = 1.0
        assert self.entry.energy == approx(-268.38319884)
        assert self.entry3.energy_per_atom == approx(2.3 / 5)

    def test_composition(self):
        assert self.entry.composition.reduced_formula == "LiFe4(PO4)4"
        assert self.entry2.composition.reduced_formula == "Fe2O3"
        assert self.entry5.composition.reduced_formula == "Fe2O3"
        assert self.entry5.composition.get_reduced_formula_and_factor()[1] == 3

    def test_per_atom_props(self):
        entry = ComputedEntry("Fe6O9", 6.9)
        entry.energy_adjustments.append(CompositionEnergyAdjustment(-0.5, 9, uncertainty_per_atom=0.1, name="O"))
        assert entry.energy == approx(2.4)
        assert entry.energy_per_atom == approx(2.4 / 15)
        assert entry.uncorrected_energy == approx(6.9)
        assert entry.uncorrected_energy_per_atom == approx(6.9 / 15)
        assert entry.correction == approx(-4.5)
        assert entry.correction_per_atom == approx(-4.5 / 15)
        assert entry.correction_uncertainty == approx(0.9)
        assert entry.correction_uncertainty_per_atom == approx(0.9 / 15)

    def test_normalize(self):
        entry = ComputedEntry("Fe6O9", 6.9, correction=1)
        entry_formula = entry.normalize()
        assert entry_formula.composition.formula == "Fe2 O3"
        assert entry_formula.uncorrected_energy == approx(6.9 / 3)
        assert entry_formula.correction == approx(1 / 3)
        assert entry_formula.energy * 3 == approx(6.9 + 1)
        assert entry_formula.energy_adjustments[0].value == approx(1 / 3)
        entry_atom = entry.normalize("atom")
        assert entry_atom.composition.formula == "Fe0.4 O0.6"
        assert entry_atom.uncorrected_energy == approx(6.9 / 15)
        assert entry_atom.correction == approx(1 / 15)
        assert entry_atom.energy * 15 == approx(6.9 + 1)
        assert entry_atom.energy_adjustments[0].value == approx(1 / 15)

    def test_normalize_energy_adjustments(self):
        ene_adj_list = [
            ManualEnergyAdjustment(5),
            ConstantEnergyAdjustment(5),
            CompositionEnergyAdjustment(1, 5, uncertainty_per_atom=0, name="Na"),
            TemperatureEnergyAdjustment(0.005, 100, 10, uncertainty_per_deg=0),
        ]
        entry = ComputedEntry("Na5Cl5", 6.9, energy_adjustments=ene_adj_list)
        assert entry.correction == 20
        normed_entry = entry.normalize()
        assert normed_entry.correction == 4
        for ea in normed_entry.energy_adjustments:
            assert ea.value == 1

    def test_to_from_dict(self):
        dct = self.entry.as_dict()
        entry = ComputedEntry.from_dict(dct)
        assert self.entry == entry
        assert entry.energy == approx(-269.38319884)

    def test_to_from_dict_with_adjustment(self):
        """Legacy case where adjustment was provided manually."""
        dct = self.entry6.as_dict()
        entry = ComputedEntry.from_dict(dct)
        assert entry.uncorrected_energy == approx(6.9)
        assert entry.energy_adjustments[0].value == self.entry6.energy_adjustments[0].value

    def test_to_from_dict_with_adjustment_2(self):
        """Modern case where correction was provided manually."""
        dct = self.entry7.as_dict()
        entry = ComputedEntry.from_dict(dct)
        assert entry.uncorrected_energy == approx(6.9)
        assert entry.energy_adjustments[0].value == self.entry7.energy_adjustments[0].value

    def test_to_from_dict_with_adjustment_3(self):
        """
        Legacy case where the entry was serialized before the energy_adjustment
        attribute was part of ComputedEntry.
        """
        # same as entry6
        dct = {
            "@module": "pymatgen.entries.computed_entries",
            "@class": "ComputedEntry",
            "energy": 6.9,
            "composition": defaultdict(float, {"Fe": 6.0, "O": 9.0}),
            "parameters": {},
            "data": {},
            "entry_id": None,
            "correction": -10,
        }
        entry = ComputedEntry.from_dict(dct)
        assert entry.uncorrected_energy == approx(6.9)
        assert entry.correction == approx(-10)
        assert len(entry.energy_adjustments) == 1

    def test_conflicting_correction_adjustment(self):
        """
        Should raise a ValueError if a user tries to manually set both the correction
        and energy_adjustment, even if the values match.
        """
        ea = ConstantEnergyAdjustment(-10, name="Dummy adjustment")
        with pytest.raises(ValueError, match="Argument conflict!"):
            ComputedEntry("Fe6O9", 6.9, correction=-10, energy_adjustments=[ea])

    def test_entry_id(self):
        assert self.entry4.entry_id == 1
        assert self.entry2.entry_id is None

    def test_str(self):
        assert str(self.entry) is not None

    def test_sulfide_energy(self):
        self.entry = ComputedEntry("BaS", -10.21249155)
        assert self.entry.energy == approx(-10.21249155)
        assert self.entry.energy_per_atom == approx(-10.21249155 / 2)
        self.entry.correction = 1.0
        assert self.entry.energy == approx(-9.21249155)

    def test_is_element(self):
        entry = ComputedEntry("Fe3", 2.3)
        assert entry.is_element

    def test_copy(self):
        for entry in (self.entry, self.entry2, self.entry3, self.entry4):
            copy = entry.copy()
            assert entry == copy
            assert str(entry) == str(copy)


class TestComputedStructureEntry(unittest.TestCase):
    def setUp(self):
        self.entry = ComputedStructureEntry(vasp_run.final_structure, vasp_run.final_energy, parameters=vasp_run.incar)

    def test_energy(self):
        assert self.entry.energy == approx(-269.38319884)
        self.entry.correction = 1.0
        assert self.entry.energy == approx(-268.38319884)

    def test_composition(self):
        assert self.entry.composition.reduced_formula == "LiFe4(PO4)4"

    def test_to_from_dict(self):
        dct = self.entry.as_dict()
        entry = ComputedStructureEntry.from_dict(dct)
        assert self.entry == entry
        assert entry.energy == approx(-269.38319884)

    def test_str(self):
        assert str(self.entry) is not None

    def test_to_from_dict_structure_with_adjustment_3(self):
        """
        Legacy case where the structure entry was serialized before the energy_adjustment
        attribute was part of ComputedEntry.
        """
        # ComputedStructureEntry for Oxygen, mp-12957, as of April 2020
        # with an arbitrary 1 eV correction added
        dct = {
            "@module": "pymatgen.entries.computed_entries",
            "@class": "ComputedStructureEntry",
            "energy": -39.42116819,
            "composition": defaultdict(float, {"O": 8.0}),
            "parameters": {
                "run_type": "GGA",
                "is_hubbard": False,
                "pseudo_potential": {
                    "functional": "PBE",
                    "labels": ["O"],
                    "pot_type": "paw",
                },
                "hubbards": {},
                "potcar_symbols": ["PBE O"],
                "oxide_type": "None",
            },
            "data": {"oxide_type": "None"},
            "entry_id": "mp-12957",
            "correction": 1,
            "structure": {
                "@module": "pymatgen.core.structure",
                "@class": "Structure",
                "charge": None,
                "lattice": {
                    "matrix": [
                        [-1.7795583, 0.0, 3.86158265],
                        [4.17564656, -3.03266995, -0.01184798],
                        [4.17564656, 3.03266995, -0.01184798],
                    ],
                    "a": 4.251899376264673,
                    "b": 5.160741380296335,
                    "c": 5.160741380296335,
                    "alpha": 71.97975354157973,
                    "beta": 109.9211782454931,
                    "gamma": 109.9211782454931,
                    "volume": 97.67332322031668,
                },
                "sites": [
                    {
                        "species": [{"element": "O", "occu": 1}],
                        "abc": [0.8531272, 0.15466029, 0.15466029],
                        "xyz": [
                            -0.22657617390155504,
                            -1.750215367360042e-17,
                            3.2907563697176516,
                        ],
                        "label": "O",
                        "properties": {"magmom": 0.002},
                    },
                    {
                        "species": [{"element": "O", "occu": 1}],
                        "abc": [0.84038763, 0.71790132, 0.21754949],
                        "xyz": [
                            2.410593174641884,
                            -1.5174019592685084,
                            3.234143088794756,
                        ],
                        "label": "O",
                        "properties": {"magmom": -0.002},
                    },
                    {
                        "species": [{"element": "O", "occu": 1}],
                        "abc": [0.17255465, 0.21942628, 0.21942628],
                        "xyz": [
                            1.5254221229000986,
                            -2.121360826524921e-18,
                            0.6611345262629937,
                        ],
                        "label": "O",
                        "properties": {"magmom": 0.002},
                    },
                    {
                        "species": [{"element": "O", "occu": 1}],
                        "abc": [0.15961237, 0.78245051, 0.28209968],
                        "xyz": [
                            4.161145821004675,
                            -1.5173989265985586,
                            0.6037435893572642,
                        ],
                        "label": "O",
                        "properties": {"magmom": -0.002},
                    },
                    {
                        "species": [{"element": "O", "occu": 1}],
                        "abc": [0.84038763, 0.21754949, 0.71790132],
                        "xyz": [
                            2.410593174641884,
                            1.5174019592685082,
                            3.234143088794756,
                        ],
                        "label": "O",
                        "properties": {"magmom": -0.002},
                    },
                    {
                        "species": [{"element": "O", "occu": 1}],
                        "abc": [0.82744535, 0.78057372, 0.78057372],
                        "xyz": [
                            5.046312697099901,
                            -1.3574974398403584e-16,
                            3.176752163737006,
                        ],
                        "label": "O",
                        "properties": {"magmom": 0.002},
                    },
                    {
                        "species": [{"element": "O", "occu": 1}],
                        "abc": [0.15961237, 0.28209968, 0.78245051],
                        "xyz": [
                            4.161145821004675,
                            1.5173989265985584,
                            0.6037435893572642,
                        ],
                        "label": "O",
                        "properties": {"magmom": -0.002},
                    },
                    {
                        "species": [{"element": "O", "occu": 1}],
                        "abc": [0.1468728, 0.84533971, 0.84533971],
                        "xyz": [
                            6.798310993901555,
                            -1.7769364890338579e-16,
                            0.5471303202823484,
                        ],
                        "label": "O",
                        "properties": {"magmom": 0.002},
                    },
                ],
            },
        }
        entry = ComputedEntry.from_dict(dct)
        assert entry.uncorrected_energy == approx(-39.42116819)
        assert entry.energy == approx(-38.42116819)
        assert entry.correction == approx(1)
        assert len(entry.energy_adjustments) == 1

    def test_copy(self):
        copy = self.entry.copy()
        assert copy == self.entry
        assert copy is not self.entry
        assert copy.structure is not self.entry.structure
        assert str(copy) == str(self.entry)

    def test_eq(self):
        copy1 = self.entry.copy()
        copy2 = self.entry.copy()
        copy3 = self.entry.copy()

        # if no entry ids but same energy/comps --> equal
        assert copy1 == copy2

        # check different entry_ids --> not equal
        copy1.entry_id = "mp_1"
        copy2.entry_id = "mp_2"
        assert copy1 != copy2

        # same entry_id and energy --> equal
        copy3.entry_id = "mp_1"
        assert copy3 == copy1

        # check different energy adjustments (but same id) --> not equal
        copy3.energy_adjustments.append(ConstantEnergyAdjustment(0.1))
        assert copy3 != copy1


class TestGibbsComputedStructureEntry(unittest.TestCase):
    def setUp(self):
        self.temps = [300, 600, 900, 1200, 1500, 1800]
        self.struct = vasp_run.final_structure
        self.num_atoms = self.struct.composition.num_atoms
        self.entries_with_temps = {
            temp: GibbsComputedStructureEntry(
                self.struct,
                -2.436,
                temp=temp,
                gibbs_model="SISSO",
                parameters=vasp_run.incar,
                entry_id="test",
            )
            for temp in self.temps
        }

        with open(f"{TEST_FILES_DIR}/Mn-O_entries.json") as f:
            data = json.load(f)
        with open(f"{TEST_FILES_DIR}/structure_CO2.json") as f:
            self.co2_struct = MontyDecoder().process_decoded(json.load(f))

        self.mp_entries = [MontyDecoder().process_decoded(d) for d in data]

    def test_gf_sisso(self):
        energies = {
            300: -56.21273010866969,
            600: -51.52997063074788,
            900: -47.29888391585979,
            1200: -42.942338738866304,
            1500: -37.793417248809774,
            1800: -32.32513382051749,
        }
        for t in self.temps:
            assert self.entries_with_temps[t].energy == approx(energies[t])

    def test_interpolation(self):
        temp = 450
        entry = GibbsComputedStructureEntry(self.struct, -2.436, temp=temp)
        assert entry.energy == approx(-53.7243542548528)

    def test_expt_gas_entry(self):
        co2_entry = GibbsComputedStructureEntry(self.co2_struct, 0, temp=900)
        assert co2_entry.energy == approx(-16.406560223724014)
        assert co2_entry.energy_per_atom == approx(-1.3672133519770011)

    def test_from_entries(self):
        gibbs_entries = GibbsComputedStructureEntry.from_entries(self.mp_entries)
        assert gibbs_entries is not None

    def test_from_pd(self):
        pd = PhaseDiagram(self.mp_entries)
        gibbs_entries = GibbsComputedStructureEntry.from_pd(pd)
        assert gibbs_entries is not None

    def test_to_from_dict(self):
        test_entry = self.entries_with_temps[300]
        dct = test_entry.as_dict()
        entry = GibbsComputedStructureEntry.from_dict(dct)
        assert test_entry == entry
        assert entry.energy == approx(test_entry.energy)

    def test_str(self):
        assert str(self.entries_with_temps[300]) is not None

    def test_normalize(self):
        for entry in self.entries_with_temps.values():
            entry = copy.deepcopy(entry)
            normed_entry = entry.normalize(mode="atom")
            assert entry.uncorrected_energy == approx(normed_entry.uncorrected_energy * self.num_atoms)
