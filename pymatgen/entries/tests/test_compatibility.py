# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import annotations

import copy
import json
import os
import unittest
import warnings
from collections import defaultdict
from math import sqrt
from pathlib import Path

import pytest
from monty.json import MontyDecoder

from pymatgen.core.composition import Composition
from pymatgen.core.lattice import Lattice
from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Structure
from pymatgen.entries.compatibility import (
    MU_H2O,
    AqueousCorrection,
    Compatibility,
    CompatibilityError,
    MaterialsProject2020Compatibility,
    MaterialsProjectAqueousCompatibility,
    MaterialsProjectCompatibility,
    MITAqueousCompatibility,
    MITCompatibility,
)
from pymatgen.entries.computed_entries import (
    ComputedEntry,
    ComputedStructureEntry,
    ConstantEnergyAdjustment,
)
from pymatgen.util.testing import PymatgenTest


class CorrectionSpecificityTest(unittest.TestCase):
    """
    Make sure corrections are only applied to GGA or GGA+U entries
    """

    def setUp(self):
        warnings.simplefilter("ignore")
        self.entry1 = ComputedEntry(
            "Fe2O3",
            -1,
            0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 5.3, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        self.entry2 = ComputedEntry(
            "FeS",
            -1,
            0.0,
            parameters={
                "is_hubbard": False,
                "run_type": "GGA",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE S 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )

        self.entry3 = ComputedEntry(
            "Fe2O3",
            -1,
            0.0,
            parameters={
                "is_hubbard": False,
                "run_type": "R2SCAN",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        self.compat = MaterialsProjectCompatibility(check_potcar_hash=False)

    def test_correction_specificity(self):
        processed = self.compat.process_entries([self.entry1, self.entry2, self.entry3])

        assert len(processed) == 2

        assert self.entry1.correction != 0
        assert self.entry2.correction != 0
        assert self.entry3.correction == 0.0


# abstract Compatibility tests
class DummyCompatibility(Compatibility):
    """
    Dummy class to test abstract Compatibility interface
    """

    def get_adjustments(self, entry):
        return [ConstantEnergyAdjustment(-10, name="Dummy adjustment")]


def test_process_entries_return_type():
    """
    process_entries should accept single entries or a list, and always return a list
    """
    entry = ComputedEntry("Fe2O3", -2)
    compat = DummyCompatibility()

    assert isinstance(compat.process_entries(entry), list)
    assert isinstance(compat.process_entries([entry]), list)


def test_no_duplicate_corrections():
    """
    Compatibility should never apply the same correction twice
    """
    entry = ComputedEntry("Fe2O3", -2)
    compat = DummyCompatibility()

    assert entry.correction == 0
    compat.process_entries(entry)
    assert entry.correction == -10
    compat.process_entries(entry)
    assert entry.correction == -10
    compat.process_entries(entry, clean=True)
    assert entry.correction == -10


def test_clean_arg():
    """
    clean=False should preserve existing corrections, clean=True should delete
    them before processing
    """
    entry = ComputedEntry("Fe2O3", -2, correction=-4)
    compat = DummyCompatibility()

    assert entry.correction == -4
    compat.process_entries(entry, clean=False)
    assert entry.correction == -14
    compat.process_entries(entry)
    assert entry.correction == -10


def test_energy_adjustment_normalize():
    """
    Both manual and automatically generated energy adjustments should be scaled
    by the normalize method
    """
    entry = ComputedEntry("Fe4O6", -2, correction=-4)
    entry = entry.normalize()
    for ea in entry.energy_adjustments:
        if "Manual" in ea.name:
            assert ea.value == -2

    compat = DummyCompatibility()
    entry = ComputedEntry("Fe4O6", -2, correction=-4)
    entry = compat.process_entries(entry)[0]
    entry = entry.normalize()
    for ea in entry.energy_adjustments:
        if "Dummy" in ea.name:
            assert ea.value == -5


def test_overlapping_adjustments():
    """
    Compatibility should raise a CompatibilityError if there is already a
    correction with the same name, but a different value, and process_entries
    should skip that entry.
    """
    ea = ConstantEnergyAdjustment(-5, name="Dummy adjustment")
    entry = ComputedEntry("Fe2O3", -2, energy_adjustments=[ea])
    compat = DummyCompatibility()

    assert entry.correction == -5

    # in case of a collision between EnergyAdjustment, check for a UserWarning
    with pytest.warns(UserWarning, match="already has an energy adjustment called Dummy"):
        processed = compat.process_entries(entry, clean=False)

    assert len(processed) == 0


class MaterialsProjectCompatibilityTest(unittest.TestCase):
    def setUp(self):
        warnings.simplefilter("ignore")
        self.entry1 = ComputedEntry(
            "Fe2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 5.3, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )

        self.entry_sulfide = ComputedEntry(
            "FeS",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": False,
                "run_type": "GGA",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE S 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )

        self.entry4 = ComputedEntry(
            "H8",
            -27.1,
            correction=0.0,
            parameters={
                "run_type": "LDA",
                "is_hubbard": False,
                "pseudo_potential": {
                    "functional": "PBE",
                    "labels": ["H"],
                    "pot_type": "paw",
                },
                "hubbards": {},
                "potcar_symbols": ["PBE H"],
                "oxide_type": "None",
            },
        )

        self.entry2 = ComputedEntry(
            "Fe3O4",
            -2,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 5.3, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        self.entry3 = ComputedEntry(
            "FeO",
            -2,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 4.3, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )

        self.compat = MaterialsProjectCompatibility(check_potcar_hash=False)
        self.ggacompat = MaterialsProjectCompatibility("GGA", check_potcar_hash=False)

    def tearDown(self):
        warnings.simplefilter("default")

    def test_process_entry(self):
        # Correct parameters
        assert self.compat.process_entry(self.entry1) is not None
        assert self.ggacompat.process_entry(self.entry1) is None

        # Correct parameters
        entry = ComputedEntry(
            "Fe2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": False,
                "hubbards": {},
                "run_type": "GGA",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        assert self.compat.process_entry(entry) is None
        assert self.ggacompat.process_entry(entry) is not None

        entry = ComputedEntry(
            "Fe2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 5.3, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        assert self.compat.process_entry(entry) is not None

    def test_correction_values(self):
        # test_corrections
        assert self.compat.process_entry(self.entry1).correction == pytest.approx(-2.733 * 2 - 0.70229 * 3)

        entry = ComputedEntry(
            "FeF3",
            -2,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 5.3, "F": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE F 08Apr2002",
                        "hash": "180141c33d032bfbfff30b3bea9d23dd",
                    },
                ],
            },
        )
        assert self.compat.process_entry(entry) is not None

        # Check actual correction
        assert self.compat.process_entry(entry).correction == pytest.approx(-2.733)

        assert self.compat.process_entry(self.entry_sulfide).correction == pytest.approx(-0.66346)

    def test_U_values(self):
        # Wrong U value
        entry = ComputedEntry(
            "Fe2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 5.2, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        assert self.compat.process_entry(entry) is None

        # GGA run of U
        entry = ComputedEntry(
            "Fe2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": False,
                "hubbards": None,
                "run_type": "GGA",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        assert self.compat.process_entry(entry) is None

        # GGA+U run of non-U
        entry = ComputedEntry(
            "Al2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Al": 5.3, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Al 06Sep2000",
                        "hash": "805c888bbd2793e462311f6a20d873d9",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        assert self.compat.process_entry(entry) is None

        # Materials project should not have a U for sulfides
        entry = ComputedEntry(
            "FeS2",
            -2,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 5.3, "S": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE S 08Apr2002",
                        "hash": "f7f8e4a74a6cbb8d63e41f4373b54df2",
                    },
                ],
            },
        )
        assert self.compat.process_entry(entry) is None

    def test_wrong_psp(self):
        # Wrong psp
        entry = ComputedEntry(
            "Fe2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 5.3, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe 06Sep2000",
                        "hash": "9530da8244e4dac17580869b4adab115",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        assert self.compat.process_entry(entry) is None

    def test_element_processing(self):
        entry = ComputedEntry(
            "O",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": False,
                "hubbards": {},
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    }
                ],
                "run_type": "GGA",
            },
        )
        entry = self.compat.process_entry(entry)
        #        self.assertEqual(entry.entry_id, -8)
        assert entry.energy == pytest.approx(-1)
        assert self.ggacompat.process_entry(entry).energy == pytest.approx(-1)

    def test_get_explanation_dict(self):
        compat = MaterialsProjectCompatibility(check_potcar_hash=False)
        entry = ComputedEntry(
            "Fe2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 5.3, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        d = compat.get_explanation_dict(entry)
        assert d["corrections"][0]["name"] == "MPRelaxSet Potcar Correction"

    def test_get_corrections_dict(self):
        compat = MaterialsProjectCompatibility(check_potcar_hash=False)
        ggacompat = MaterialsProjectCompatibility("GGA", check_potcar_hash=False)

        # Correct parameters
        entry = ComputedEntry(
            "Fe2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 5.3, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        c = compat.get_corrections_dict(entry)[0]
        assert c["MP Anion Correction"] == pytest.approx(-2.10687)
        assert c["MP Advanced Correction"] == pytest.approx(-5.466)

        entry.parameters["is_hubbard"] = False
        del entry.parameters["hubbards"]
        c = ggacompat.get_corrections_dict(entry)[0]
        assert "MP Advanced Correction" not in c

    def test_process_entries(self):
        entries = self.compat.process_entries([self.entry1, self.entry2, self.entry3, self.entry4])
        assert len(entries) == 2

    def test_msonable(self):
        compat_dict = self.compat.as_dict()
        decoder = MontyDecoder()
        temp_compat = decoder.process_decoded(compat_dict)
        assert isinstance(temp_compat, MaterialsProjectCompatibility)


class MaterialsProjectCompatibility2020Test(unittest.TestCase):
    def setUp(self):
        warnings.simplefilter("ignore")
        self.entry1 = ComputedEntry(
            "Fe2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 5.3, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )

        self.entry_sulfide = ComputedEntry(
            "FeS",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": False,
                "run_type": "GGA",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE S 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )

        self.entry2 = ComputedEntry(
            "Fe3O4",
            -2,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 5.3, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        self.entry3 = ComputedEntry(
            "FeO",
            -2,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 4.3, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )

        self.compat = MaterialsProject2020Compatibility(check_potcar_hash=False)
        self.ggacompat = MaterialsProject2020Compatibility("GGA", check_potcar_hash=False)

    def tearDown(self):
        warnings.simplefilter("default")

    def test_process_entry(self):
        # Correct parameters
        assert self.compat.process_entry(self.entry1) is not None
        assert self.ggacompat.process_entry(self.entry1) is None

        # Correct parameters
        entry = ComputedEntry(
            "Fe2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": False,
                "hubbards": {},
                "run_type": "GGA",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        assert self.compat.process_entry(entry) is None
        assert self.ggacompat.process_entry(entry) is not None

        entry = ComputedEntry(
            "Fe2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 5.3, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        assert self.compat.process_entry(entry) is not None

    def test_oxi_state_guess(self):
        # An entry where Composition.oxi_state_guesses will return an empty list
        entry_blank = ComputedEntry(
            "Ga3Te",
            -12.1900,
            correction=0.0,
            parameters={
                "run_type": "GGA",
                "is_hubbard": False,
                "pseudo_potential": {"functional": "PBE", "labels": ["Ga_d", "Te"], "pot_type": "paw"},
                "hubbards": {},
                "potcar_symbols": ["PBE Ga_d", "PBE Te"],
                "oxide_type": "None",
            },
        )

        # An entry where one anion will only be corrected if oxidation_states is populated
        entry_oxi = ComputedEntry(
            "Mo2Cl8O",
            -173.0655,
            correction=0.0,
            parameters={
                "run_type": "GGA+U",
                "is_hubbard": True,
                "pseudo_potential": {"functional": "PBE", "labels": ["Mo_pv", "Cl", "O"], "pot_type": "paw"},
                "hubbards": {"Mo": 4.38, "Cl": 0.0, "O": 0.0},
                "potcar_symbols": ["PBE Mo_pv", "PBE Cl", "PBE O"],
                "oxide_type": "oxide",
            },
        )

        # An entry that should receive multiple anion corrections if oxidation
        # states are populated
        entry_multi_anion = ComputedEntry(
            "C8N4Cl4",
            -87.69656726,
            correction=0.0,
            parameters={
                "run_type": "GGA",
                "is_hubbard": False,
                "pseudo_potential": {"functional": "PBE", "labels": ["C", "N", "Cl"], "pot_type": "paw"},
                "hubbards": {},
                "potcar_symbols": ["PBE C", "PBE N", "PBE Cl"],
                "oxide_type": "None",
            },
        )

        with pytest.warns(UserWarning, match="Failed to guess oxidation state"):
            e1 = self.compat.process_entry(entry_blank)
            assert e1.correction == pytest.approx(-0.422)

        e2 = self.compat.process_entry(entry_oxi)
        assert e2.correction == pytest.approx(-0.687 + -3.202 * 2 + -0.614 * 8)

        e3 = self.compat.process_entry(entry_multi_anion)
        assert e3.correction == pytest.approx(-0.361 * 4 + -0.614 * 4)

    def test_correction_values(self):
        # test_corrections
        assert self.compat.process_entry(self.entry1).correction == pytest.approx(-2.256 * 2 - 0.687 * 3)

        entry = ComputedEntry(
            "FeF3",
            -2,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 5.3, "F": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE F 08Apr2002",
                        "hash": "180141c33d032bfbfff30b3bea9d23dd",
                    },
                ],
            },
        )
        assert self.compat.process_entry(entry) is not None

        # Check actual correction
        assert self.compat.process_entry(entry).correction == pytest.approx(-0.462 * 3 + -2.256)

        assert self.compat.process_entry(self.entry_sulfide).correction == pytest.approx(-0.503)

    def test_oxdiation_by_electronegativity(self):
        # make sure anion corrections are only applied when the element has
        # a negative oxidation state (e.g., correct CaSi but not SiO2 for Si)
        # as determined by electronegativity (i.e., the data.oxidation_states key is absent)

        entry1 = ComputedEntry.from_dict(
            {
                "@module": "pymatgen.entries.computed_entries",
                "@class": "ComputedEntry",
                "energy": -17.01015622,
                "composition": defaultdict(float, {"Si": 2.0, "Ca": 2.0}),
                "energy_adjustments": [],
                "parameters": {
                    "run_type": "GGA",
                    "is_hubbard": False,
                    "pseudo_potential": {
                        "functional": "PBE",
                        "labels": ["Ca_sv", "Si"],
                        "pot_type": "paw",
                    },
                    "hubbards": {},
                    "potcar_symbols": ["PBE Ca_sv", "PBE Si"],
                    "oxide_type": "None",
                },
                "data": {"oxide_type": "None"},
                "entry_id": "mp-1563",
                "correction": 0.0,
            }
        )

        entry2 = ComputedEntry.from_dict(
            {
                "@module": "pymatgen.entries.computed_entries",
                "@class": "ComputedEntry",
                "energy": -47.49120119,
                "composition": defaultdict(float, {"Si": 2.0, "O": 4.0}),
                "energy_adjustments": [],
                "parameters": {
                    "run_type": "GGA",
                    "is_hubbard": False,
                    "pseudo_potential": {
                        "functional": "PBE",
                        "labels": ["Si", "O"],
                        "pot_type": "paw",
                    },
                    "hubbards": {},
                    "potcar_symbols": ["PBE Si", "PBE O"],
                    "oxide_type": "oxide",
                },
                "data": {"oxide_type": "oxide"},
                "entry_id": "mp-546794",
                "correction": 0.0,
            }
        )

        # CaSi; only correction should be Si
        assert self.compat.process_entry(entry1).correction == pytest.approx(0.071 * 2)

        # SiO2; only corrections should be oxide
        assert self.compat.process_entry(entry2).correction == pytest.approx(-0.687 * 4)

    def test_oxdiation(self):
        # make sure anion corrections are only applied when the element has
        # a negative oxidation state (e.g., correct CaSi but not SiO2 for Si)
        # as determined by the data.oxidation_states key

        entry1 = ComputedEntry.from_dict(
            {
                "@module": "pymatgen.entries.computed_entries",
                "@class": "ComputedEntry",
                "energy": -17.01015622,
                "composition": defaultdict(float, {"Si": 2.0, "Ca": 2.0}),
                "energy_adjustments": [],
                "parameters": {
                    "run_type": "GGA",
                    "is_hubbard": False,
                    "pseudo_potential": {
                        "functional": "PBE",
                        "labels": ["Ca_sv", "Si"],
                        "pot_type": "paw",
                    },
                    "hubbards": {},
                    "potcar_symbols": ["PBE Ca_sv", "PBE Si"],
                    "oxide_type": "None",
                },
                "data": {
                    "oxide_type": "None",
                    "oxidation_states": {"Ca": 2.0, "Si": -2.0},
                },
                "entry_id": "mp-1563",
                "correction": 0.0,
            }
        )

        entry2 = ComputedEntry.from_dict(
            {
                "@module": "pymatgen.entries.computed_entries",
                "@class": "ComputedEntry",
                "energy": -47.49120119,
                "composition": defaultdict(float, {"Si": 2.0, "O": 4.0}),
                "energy_adjustments": [],
                "parameters": {
                    "run_type": "GGA",
                    "is_hubbard": False,
                    "pseudo_potential": {
                        "functional": "PBE",
                        "labels": ["Si", "O"],
                        "pot_type": "paw",
                    },
                    "hubbards": {},
                    "potcar_symbols": ["PBE Si", "PBE O"],
                    "oxide_type": "oxide",
                },
                "data": {
                    "oxide_type": "oxide",
                    "oxidation_states": {"Si": 4.0, "O": -2.0},
                },
                "entry_id": "mp-546794",
                "correction": 0.0,
            }
        )

        # CaSi; only correction should be Si
        assert self.compat.process_entry(entry1).correction == pytest.approx(0.071 * 2)

        # SiO2; only corrections should be oxide
        assert self.compat.process_entry(entry2).correction == pytest.approx(-0.687 * 4)

    def test_U_values(self):
        # Wrong U value
        entry = ComputedEntry(
            "Fe2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 5.2, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        assert self.compat.process_entry(entry) is None

        # GGA run of U
        entry = ComputedEntry(
            "Fe2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": False,
                "hubbards": None,
                "run_type": "GGA",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        assert self.compat.process_entry(entry) is None

        # GGA+U run of non-U
        entry = ComputedEntry(
            "Al2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Al": 5.3, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Al 06Sep2000",
                        "hash": "805c888bbd2793e462311f6a20d873d9",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        assert self.compat.process_entry(entry) is None

        # Materials project should not have a U for sulfides
        entry = ComputedEntry(
            "FeS2",
            -2,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 5.3, "S": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE S 08Apr2002",
                        "hash": "f7f8e4a74a6cbb8d63e41f4373b54df2",
                    },
                ],
            },
        )
        assert self.compat.process_entry(entry) is None

    def test_wrong_psp(self):
        # Wrong psp
        entry = ComputedEntry(
            "Fe2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 5.3, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe 06Sep2000",
                        "hash": "9530da8244e4dac17580869b4adab115",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        assert self.compat.process_entry(entry) is None

    def test_element_processing(self):
        entry = ComputedEntry(
            "O",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": False,
                "hubbards": {},
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    }
                ],
                "run_type": "GGA",
            },
        )
        entry = self.compat.process_entry(entry)
        assert entry.energy == pytest.approx(-1)
        assert self.ggacompat.process_entry(entry).energy == pytest.approx(-1)

    def test_get_explanation_dict(self):
        compat = MaterialsProjectCompatibility(check_potcar_hash=False)
        entry = ComputedEntry(
            "Fe2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 5.3, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        d = compat.get_explanation_dict(entry)
        assert d["corrections"][0]["name"] == "MPRelaxSet Potcar Correction"

    def test_energy_adjustments(self):
        compat = MaterialsProject2020Compatibility(check_potcar_hash=False)
        ggacompat = MaterialsProject2020Compatibility("GGA", check_potcar_hash=False)

        # Fe 4 Co 2 O 8 (Fe2CoO4)
        entry = {
            "@module": "pymatgen.entries.computed_entries",
            "@class": "ComputedEntry",
            "energy": -91.94962744,
            "composition": defaultdict(float, {"Fe": 4.0, "Co": 2.0, "O": 8.0}),
            "energy_adjustments": [],
            "parameters": {
                "run_type": "GGA+U",
                "is_hubbard": True,
                "pseudo_potential": {
                    "functional": "PBE",
                    "labels": ["Fe_pv", "Co", "O"],
                    "pot_type": "paw",
                },
                "hubbards": {"Fe": 5.3, "Co": 3.32, "O": 0.0},
                "potcar_symbols": ["PBE Fe_pv", "PBE Co", "PBE O"],
                "oxide_type": "oxide",
            },
            "data": {"oxide_type": "oxide"},
            "entry_id": "mp-753222",
            "correction": 0,
        }
        entry = ComputedEntry.from_dict(entry)

        c = compat.process_entry(entry)
        assert "MP2020 anion correction (oxide)" in [ea.name for ea in c.energy_adjustments]
        assert "MP2020 GGA/GGA+U mixing correction (Fe)" in [ea.name for ea in c.energy_adjustments]
        assert "MP2020 GGA/GGA+U mixing correction (Co)" in [ea.name for ea in c.energy_adjustments]

        for ea in c.energy_adjustments:
            if ea.name == "MP2020 GGA/GGA+U mixing correction (Fe)":
                assert ea.value == pytest.approx(-2.256 * 4)
                assert ea.uncertainty == pytest.approx(0.0101 * 4)
            elif ea.name == "MP2020 GGA/GGA+U mixing correction (Co)":
                assert ea.value == pytest.approx(-1.638 * 2)
                assert ea.uncertainty == pytest.approx(0.006 * 2)
            elif ea.name == "MP2020 anion correction (oxide)":
                assert ea.value == pytest.approx(-0.687 * 8)
                assert ea.uncertainty == pytest.approx(0.002 * 8)

        entry.parameters["is_hubbard"] = False
        del entry.parameters["hubbards"]
        c = ggacompat.process_entry(entry)
        assert "MP2020 GGA/GGA+U mixing correction" not in [ea.name for ea in c.energy_adjustments]

    def test_process_entries(self):
        entries = self.compat.process_entries([self.entry1, self.entry2, self.entry3])
        assert len(entries) == 2

    def test_config_file(self):
        config_file = Path(PymatgenTest.TEST_FILES_DIR / "MP2020Compatibility_alternate.yaml")
        compat = MaterialsProject2020Compatibility(config_file=config_file)
        entry = compat.process_entry(self.entry1)
        for ea in entry.energy_adjustments:
            if ea.name == "MP2020 GGA/GGA+U mixing correction (Fe)":
                assert ea.value == pytest.approx(-0.224 * 2)

    def test_msonable(self):
        compat_dict = self.compat.as_dict()
        decoder = MontyDecoder()
        temp_compat = decoder.process_decoded(compat_dict)
        assert isinstance(temp_compat, MaterialsProject2020Compatibility)

    def test_processing_entries_inplace(self):
        # load two entries in GGA_GGA_U_R2SCAN thermo type
        entriesJson = Path(PymatgenTest.TEST_FILES_DIR / "entries_thermo_type_GGA_GGA_U_R2SCAN.json")
        with open(entriesJson) as file:
            entries = json.load(file, cls=MontyDecoder)
        # check whether the compatibility scheme can keep input entries unchanged
        entries_copy = copy.deepcopy(entries)
        self.compat.process_entries(entries, inplace=False)
        assert all(e.correction == e_copy.correction for e, e_copy in zip(entries, entries_copy))


class MITCompatibilityTest(unittest.TestCase):
    def tearDown(self):
        warnings.simplefilter("default")

    def setUp(self):
        warnings.simplefilter("ignore")
        self.compat = MITCompatibility(check_potcar_hash=True)
        self.ggacompat = MITCompatibility("GGA", check_potcar_hash=True)
        self.entry_O = ComputedEntry(
            "Fe2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 4.0, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe 06Sep2000",
                        "hash": "9530da8244e4dac17580869b4adab115",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )

        self.entry_F = ComputedEntry(
            "FeF3",
            -2,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 4.0, "F": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe 06Sep2000",
                        "hash": "9530da8244e4dac17580869b4adab115",
                    },
                    {
                        "titel": "PAW_PBE F 08Apr2002",
                        "hash": "180141c33d032bfbfff30b3bea9d23dd",
                    },
                ],
            },
        )
        self.entry_S = ComputedEntry(
            "FeS2",
            -2,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 1.9, "S": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe 06Sep2000",
                        "hash": "9530da8244e4dac17580869b4adab115",
                    },
                    {
                        "titel": "PAW_PBE S 08Apr2002",
                        "hash": "d368db6899d8839859bbee4811a42a88",
                    },
                ],
            },
        )

    def test_process_entry(self):
        # Correct parameters
        assert self.compat.process_entry(self.entry_O) is not None
        assert self.compat.process_entry(self.entry_F) is not None

    def test_correction_value(self):
        # Check actual correction
        assert self.compat.process_entry(self.entry_O).correction == pytest.approx(-1.723 * 2 - 0.66975 * 3)
        assert self.compat.process_entry(self.entry_F).correction == pytest.approx(-1.723)
        assert self.compat.process_entry(self.entry_S).correction == pytest.approx(-1.113)

    def test_U_value(self):
        # MIT should have a U value for Fe containing sulfides
        assert self.compat.process_entry(self.entry_S) is not None

        # MIT should not have a U value for Ni containing sulfides
        entry = ComputedEntry(
            "NiS2",
            -2,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Ni": 1.9, "S": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Ni 06Sep2000",
                        "hash": "653f5772e68b2c7fd87ffd1086c0d710",
                    },
                    {
                        "titel": "PAW_PBE S 08Apr2002",
                        "hash": "d368db6899d8839859bbee4811a42a88",
                    },
                ],
            },
        )

        assert self.compat.process_entry(entry) is None

        entry = ComputedEntry(
            "NiS2",
            -2,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": None,
                "run_type": "GGA",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Ni 06Sep2000",
                        "hash": "653f5772e68b2c7fd87ffd1086c0d710",
                    },
                    {
                        "titel": "PAW_PBE S 08Apr2002",
                        "hash": "d368db6899d8839859bbee4811a42a88",
                    },
                ],
            },
        )

        assert self.ggacompat.process_entry(entry) is not None

    def test_wrong_U_value(self):
        # Wrong U value
        entry = ComputedEntry(
            "Fe2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 5.2, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe 06Sep2000",
                        "hash": "9530da8244e4dac17580869b4adab115",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )

        assert self.compat.process_entry(entry) is None

        # GGA run
        entry = ComputedEntry(
            "Fe2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": False,
                "hubbards": None,
                "run_type": "GGA",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe 06Sep2000",
                        "hash": "9530da8244e4dac17580869b4adab115",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        assert self.compat.process_entry(entry) is None
        assert self.ggacompat.process_entry(entry) is not None

    def test_wrong_psp(self):
        # Wrong psp
        entry = ComputedEntry(
            "Fe2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 4.0, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        assert self.compat.process_entry(entry) is None

    def test_element_processing(self):
        # Testing processing of elements.
        entry = ComputedEntry(
            "O",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": False,
                "hubbards": {},
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    }
                ],
                "run_type": "GGA",
            },
        )
        entry = self.compat.process_entry(entry)
        assert entry.energy == pytest.approx(-1)

    def test_same_potcar_symbol(self):
        # Same symbol different hash thus a different potcar
        # Correct Hash Correct Symbol
        entry = ComputedEntry(
            "Fe2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 4.0, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe 06Sep2000",
                        "hash": "9530da8244e4dac17580869b4adab115",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        # Incorrect Hash Correct Symbol
        entry2 = ComputedEntry(
            "Fe2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 4.0, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {"titel": "PAW_PBE Fe 06Sep2000", "hash": "DifferentHash"},
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )

        compat = MITCompatibility()
        assert len(compat.process_entries([entry, entry2])) == 2
        assert len(self.compat.process_entries([entry, entry2])) == 1

    def test_revert_to_symbols(self):
        # Test that you can revert to potcar_symbols if potcar_spec is not present
        compat = MITCompatibility()
        entry = ComputedEntry(
            "Fe2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 4.0, "O": 0},
                "run_type": "GGA+U",
                "potcar_symbols": ["PAW_PBE Fe 06Sep2000", "PAW_PBE O 08Apr2002"],
            },
        )

        assert compat.process_entry(entry) is not None
        # raise if check_potcar_hash is set
        with pytest.raises(ValueError):
            self.compat.process_entry(entry)

    def test_potcar_doenst_match_structure(self):
        compat = MITCompatibility()
        entry = ComputedEntry(
            "Li2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 4.0, "O": 0},
                "run_type": "GGA+U",
                "potcar_symbols": ["PAW_PBE Fe_pv 06Sep2000", "PAW_PBE O 08Apr2002"],
            },
        )

        assert compat.process_entry(entry) is None

    def test_potcar_spec_is_none(self):
        compat = MITCompatibility(check_potcar_hash=True)
        entry = ComputedEntry(
            "Li2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 4.0, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [None, None],
            },
        )

        assert compat.process_entry(entry) is None

    def test_get_explanation_dict(self):
        compat = MITCompatibility(check_potcar_hash=False)
        entry = ComputedEntry(
            "Fe2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 4.0, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe 06Sep2000",
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )
        d = compat.get_explanation_dict(entry)
        assert d["corrections"][0]["name"] == "MITRelaxSet Potcar Correction"

    def test_msonable(self):
        compat_dict = self.compat.as_dict()
        decoder = MontyDecoder()
        temp_compat = decoder.process_decoded(compat_dict)
        assert isinstance(temp_compat, MITCompatibility)


class OxideTypeCorrectionTest(unittest.TestCase):
    def setUp(self):
        self.compat = MITCompatibility(check_potcar_hash=True)

    def test_no_struct_compat(self):
        lio2_entry_nostruct = ComputedEntry(
            Composition("Li2O4"),
            -3,
            data={"oxide_type": "superoxide"},
            parameters={
                "is_hubbard": False,
                "hubbards": None,
                "run_type": "GGA",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Li 17Jan2003",
                        "hash": "65e83282d1707ec078c1012afbd05be8",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )

        lio2_entry_corrected = self.compat.process_entry(lio2_entry_nostruct)
        assert lio2_entry_corrected.energy == pytest.approx(-3 - 0.13893 * 4)

    def test_process_entry_superoxide(self):
        el_li = Element("Li")
        el_o = Element("O")
        latt = Lattice([[3.985034, 0.0, 0.0], [0.0, 4.881506, 0.0], [0.0, 0.0, 2.959824]])
        elts = [el_li, el_li, el_o, el_o, el_o, el_o]
        coords = []
        coords.append([0.500000, 0.500000, 0.500000])
        coords.append([0.0, 0.0, 0.0])
        coords.append([0.632568, 0.085090, 0.500000])
        coords.append([0.367432, 0.914910, 0.500000])
        coords.append([0.132568, 0.414910, 0.000000])
        coords.append([0.867432, 0.585090, 0.000000])
        struct = Structure(latt, elts, coords)
        lio2_entry = ComputedStructureEntry(
            struct,
            -3,
            parameters={
                "is_hubbard": False,
                "hubbards": None,
                "run_type": "GGA",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Li 17Jan2003",
                        "hash": "65e83282d1707ec078c1012afbd05be8",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )

        lio2_entry_corrected = self.compat.process_entry(lio2_entry)
        assert lio2_entry_corrected.energy == pytest.approx(-3 - 0.13893 * 4)

    def test_process_entry_peroxide(self):
        latt = Lattice.from_parameters(3.159597, 3.159572, 7.685205, 89.999884, 89.999674, 60.000510)
        el_li = Element("Li")
        el_o = Element("O")
        elts = [el_li, el_li, el_li, el_li, el_o, el_o, el_o, el_o]
        coords = [
            [0.666656, 0.666705, 0.750001],
            [0.333342, 0.333378, 0.250001],
            [0.000001, 0.000041, 0.500001],
            [0.000001, 0.000021, 0.000001],
            [0.333347, 0.333332, 0.649191],
            [0.333322, 0.333353, 0.850803],
            [0.666666, 0.666686, 0.350813],
            [0.666665, 0.666684, 0.149189],
        ]
        struct = Structure(latt, elts, coords)
        li2o2_entry = ComputedStructureEntry(
            struct,
            -3,
            parameters={
                "is_hubbard": False,
                "hubbards": None,
                "run_type": "GGA",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Li 17Jan2003",
                        "hash": "65e83282d1707ec078c1012afbd05be8",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )

        li2o2_entry_corrected = self.compat.process_entry(li2o2_entry)
        assert li2o2_entry_corrected.energy == pytest.approx(-3 - 0.44317 * 4)

    def test_process_entry_ozonide(self):
        el_li = Element("Li")
        el_o = Element("O")
        elts = [el_li, el_o, el_o, el_o]
        latt = Lattice.from_parameters(3.999911, 3.999911, 3.999911, 133.847504, 102.228244, 95.477342)
        coords = [
            [0.513004, 0.513004, 1.000000],
            [0.017616, 0.017616, 0.000000],
            [0.649993, 0.874790, 0.775203],
            [0.099587, 0.874790, 0.224797],
        ]
        struct = Structure(latt, elts, coords)
        lio3_entry = ComputedStructureEntry(
            struct,
            -3,
            parameters={
                "is_hubbard": False,
                "hubbards": None,
                "run_type": "GGA",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Li 17Jan2003",
                        "hash": "65e83282d1707ec078c1012afbd05be8",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )

        lio3_entry_corrected = self.compat.process_entry(lio3_entry)
        assert lio3_entry_corrected.energy == pytest.approx(-3.0)

    def test_process_entry_oxide(self):
        el_li = Element("Li")
        el_o = Element("O")
        elts = [el_li, el_li, el_o]
        latt = Lattice.from_parameters(3.278, 3.278, 3.278, 60, 60, 60)
        coords = [[0.25, 0.25, 0.25], [0.75, 0.75, 0.75], [0.0, 0.0, 0.0]]
        struct = Structure(latt, elts, coords)
        li2o_entry = ComputedStructureEntry(
            struct,
            -3,
            parameters={
                "is_hubbard": False,
                "hubbards": None,
                "run_type": "GGA",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Li 17Jan2003",
                        "hash": "65e83282d1707ec078c1012afbd05be8",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )

        li2o_entry_corrected = self.compat.process_entry(li2o_entry)
        assert li2o_entry_corrected.energy == pytest.approx(-3.0 - 0.66975)


class SulfideTypeCorrection2020Test(unittest.TestCase):
    def setUp(self):
        self.compat = MaterialsProject2020Compatibility(check_potcar_hash=False)

    def test_struct_no_struct(self):
        # Processing an Entry should produce the same correction whether or not
        # that entry has a Structure attached to it.

        # Na2S2, entry mp-2400, with and without structure
        from collections import defaultdict

        entry_struct_as_dict = {
            "@module": "pymatgen.entries.computed_entries",
            "@class": "ComputedStructureEntry",
            "energy": -28.42580746,
            "composition": defaultdict(float, {"Na": 4.0, "S": 4.0}),
            "correction": 0,
            "parameters": {
                "run_type": "GGA",
                "is_hubbard": False,
                "pseudo_potential": {
                    "functional": "PBE",
                    "labels": ["Na_pv", "S"],
                    "pot_type": "paw",
                },
                "hubbards": {},
                "potcar_symbols": ["PBE Na_pv", "PBE S"],
                "oxide_type": "None",
            },
            "data": {"oxide_type": "None"},
            "entry_id": "mp-2400",
            "structure": {
                "@module": "pymatgen.core.structure",
                "@class": "Structure",
                "charge": None,
                "lattice": {
                    "matrix": [
                        [4.5143094, 0.0, 0.0],
                        [-2.2571547, 3.90950662, 0.0],
                        [0.0, 0.0, 10.28414905],
                    ],
                    "a": 4.5143094,
                    "b": 4.514309399183436,
                    "c": 10.28414905,
                    "alpha": 90.0,
                    "beta": 90.0,
                    "gamma": 120.00000000598358,
                    "volume": 181.50209256783256,
                },
                "sites": [
                    {
                        "species": [{"element": "Na", "occu": 1}],
                        "abc": [0.0, 0.0, 0.0],
                        "xyz": [0.0, 0.0, 0.0],
                        "label": "Na",
                        "properties": {"magmom": 0.0},
                    },
                    {
                        "species": [{"element": "Na", "occu": 1}],
                        "abc": [0.0, 0.0, 0.5],
                        "xyz": [0.0, 0.0, 5.142074525],
                        "label": "Na",
                        "properties": {"magmom": 0.0},
                    },
                    {
                        "species": [{"element": "Na", "occu": 1}],
                        "abc": [0.33333333, 0.66666667, 0.25],
                        "xyz": [
                            -2.2571547075855847e-08,
                            2.6063377596983557,
                            2.5710372625,
                        ],
                        "label": "Na",
                        "properties": {"magmom": 0.0},
                    },
                    {
                        "species": [{"element": "Na", "occu": 1}],
                        "abc": [0.66666667, 0.33333333, 0.75],
                        "xyz": [2.2571547225715474, 1.3031688603016447, 7.7131117875],
                        "label": "Na",
                        "properties": {"magmom": 0.0},
                    },
                    {
                        "species": [{"element": "S", "occu": 1}],
                        "abc": [0.33333333, 0.66666667, 0.644551],
                        "xyz": [
                            -2.2571547075855847e-08,
                            2.6063377596983557,
                            6.62865855432655,
                        ],
                        "label": "S",
                        "properties": {"magmom": 0.0},
                    },
                    {
                        "species": [{"element": "S", "occu": 1}],
                        "abc": [0.66666667, 0.33333333, 0.144551],
                        "xyz": [
                            2.2571547225715474,
                            1.3031688603016447,
                            1.4865840293265502,
                        ],
                        "label": "S",
                        "properties": {"magmom": 0.0},
                    },
                    {
                        "species": [{"element": "S", "occu": 1}],
                        "abc": [0.66666667, 0.33333333, 0.355449],
                        "xyz": [
                            2.2571547225715474,
                            1.3031688603016447,
                            3.65549049567345,
                        ],
                        "label": "S",
                        "properties": {"magmom": 0.0},
                    },
                    {
                        "species": [{"element": "S", "occu": 1}],
                        "abc": [0.33333333, 0.66666667, 0.855449],
                        "xyz": [
                            -2.2571547075855847e-08,
                            2.6063377596983557,
                            8.79756502067345,
                        ],
                        "label": "S",
                        "properties": {"magmom": 0.0},
                    },
                ],
            },
        }

        entry_no_struct_as_dict = {
            "@module": "pymatgen.entries.computed_entries",
            "@class": "ComputedEntry",
            "energy": -28.42580746,
            "composition": defaultdict(float, {"Na": 4.0, "S": 4.0}),
            "correction": 0,
            "parameters": {
                "run_type": "GGA",
                "is_hubbard": False,
                "pseudo_potential": {
                    "functional": "PBE",
                    "labels": ["Na_pv", "S"],
                    "pot_type": "paw",
                },
                "hubbards": {},
                "potcar_symbols": ["PBE Na_pv", "PBE S"],
                "oxide_type": "None",
            },
            "data": {"oxide_type": "None"},
            "entry_id": "mp-2400",
        }

        na2s2_entry_struct = ComputedStructureEntry.from_dict(entry_struct_as_dict)
        na2s2_entry_nostruct = ComputedEntry.from_dict(entry_no_struct_as_dict)

        struct_corrected = self.compat.process_entry(na2s2_entry_struct)
        nostruct_corrected = self.compat.process_entry(na2s2_entry_nostruct)

        assert struct_corrected.correction == pytest.approx(nostruct_corrected.correction)


class OxideTypeCorrectionNoPeroxideCorrTest(unittest.TestCase):
    def setUp(self):
        self.compat = MITCompatibility(correct_peroxide=False)

    def test_oxide_energy_corr(self):
        el_li = Element("Li")
        el_o = Element("O")
        elts = [el_li, el_li, el_o]
        latt = Lattice.from_parameters(3.278, 3.278, 3.278, 60, 60, 60)
        coords = [[0.25, 0.25, 0.25], [0.75, 0.75, 0.75], [0.0, 0.0, 0.0]]
        struct = Structure(latt, elts, coords)
        li2o_entry = ComputedStructureEntry(
            struct,
            -3,
            parameters={
                "is_hubbard": False,
                "hubbards": None,
                "run_type": "GGA",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Li 17Jan2003",
                        "hash": "65e83282d1707ec078c1012afbd05be8",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )

        li2o_entry_corrected = self.compat.process_entry(li2o_entry)
        assert li2o_entry_corrected.energy == pytest.approx(-3.0 - 0.66975)

    def test_peroxide_energy_corr(self):
        latt = Lattice.from_parameters(3.159597, 3.159572, 7.685205, 89.999884, 89.999674, 60.000510)
        el_li = Element("Li")
        el_o = Element("O")
        elts = [el_li, el_li, el_li, el_li, el_o, el_o, el_o, el_o]
        coords = [
            [0.666656, 0.666705, 0.750001],
            [0.333342, 0.333378, 0.250001],
            [0.000001, 0.000041, 0.500001],
            [0.000001, 0.000021, 0.000001],
            [0.333347, 0.333332, 0.649191],
            [0.333322, 0.333353, 0.850803],
            [0.666666, 0.666686, 0.350813],
            [0.666665, 0.666684, 0.149189],
        ]
        struct = Structure(latt, elts, coords)
        li2o2_entry = ComputedStructureEntry(
            struct,
            -3,
            parameters={
                "is_hubbard": False,
                "hubbards": None,
                "run_type": "GGA",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Li 17Jan2003",
                        "hash": "65e83282d1707ec078c1012afbd05be8",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )

        li2o2_entry_corrected = self.compat.process_entry(li2o2_entry)
        with pytest.raises(AssertionError):
            self.assertAlmostEqual(*(li2o2_entry_corrected.energy, -3 - 0.44317 * 4, 4))
        assert li2o2_entry_corrected.energy == pytest.approx(-3 - 0.66975 * 4)

    def test_ozonide(self):
        el_li = Element("Li")
        el_o = Element("O")
        elts = [el_li, el_o, el_o, el_o]
        latt = Lattice.from_parameters(3.999911, 3.999911, 3.999911, 133.847504, 102.228244, 95.477342)
        coords = [
            [0.513004, 0.513004, 1.000000],
            [0.017616, 0.017616, 0.000000],
            [0.649993, 0.874790, 0.775203],
            [0.099587, 0.874790, 0.224797],
        ]
        struct = Structure(latt, elts, coords)
        lio3_entry = ComputedStructureEntry(
            struct,
            -3,
            parameters={
                "is_hubbard": False,
                "hubbards": None,
                "run_type": "GGA",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Li 17Jan2003",
                        "hash": "65e83282d1707ec078c1012afbd05be8",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )

        lio3_entry_corrected = self.compat.process_entry(lio3_entry)
        assert lio3_entry_corrected.energy == pytest.approx(-3.0 - 3 * 0.66975)


class TestMaterialsProjectAqueousCompatibility:
    """
    Test MaterialsProjectAqueousCompatibility

    -x- formation energy of H2O should always be -2.458 eV/H2O
    -x- H2 energy should always be the same value
    -x- H2O energy should always be the same value
    -x- Should get warnings if you init without all energy args
    -x- Should get CompatibilityError if you get_entry without all energy args
    -x- energy args should auto-populate from entries passed to process_entries
    -x- check compound entropies appropriately added
    -x- check hydrate adjustment appropriately applied

    Notes:
        Argument values from MaterialsProjectCompatibility as of April 2020:
            corrected DFT energy of H2O = -15.5875 eV/H2O (mp-697111) or -5.195 eV/atom
            corrected DFT energy of O2 = -4.9276 eV/atom (mp-12957)
            total energy corrections applied to H2O (eV/H2O) -0.70229 eV/H2O or -0.234 eV/atom
    """

    def test_h_h2o_energy_with_args_single(self):
        compat = MaterialsProjectAqueousCompatibility(
            o2_energy=-4.9276,
            h2o_energy=-5,
            h2o_adjustments=-0.234,
            solid_compat=None,
        )

        h2o_entry_1 = ComputedEntry(Composition("H2O"), -15)  # -5 eV/atom
        h2o_entry_2 = ComputedEntry(Composition("H4O2"), -6)  # -1 eV/atom
        h2_entry_1 = ComputedEntry(Composition("H8"), -100)  # -12.5 eV/atom
        h2_entry_2 = ComputedEntry(Composition("H2"), -16)  # -8 eV/atom

        for entry in [h2o_entry_1, h2o_entry_2]:
            compat.process_entries(entry)

        for entry in [h2_entry_1, h2_entry_2]:
            with pytest.warns(UserWarning, match="Processing single H2 entries"):
                compat.process_entries(entry)

        # the corrections should set the energy of any H2 polymorph the same, because
        # we have only processed one entry at time. Energy differences of H2O
        # polymorphs should be preserved.
        assert h2o_entry_2.energy_per_atom == pytest.approx(h2o_entry_1.energy_per_atom + 4)
        assert h2_entry_2.energy_per_atom == pytest.approx(h2_entry_1.energy_per_atom)

        o2_entry_1 = ComputedEntry(Composition("O2"), -4.9276 * 2)
        o2_entry_1 = compat.process_entries(o2_entry_1)[0]

        h2o_form_e = 3 * h2o_entry_1.energy_per_atom - 2 * h2_entry_2.energy_per_atom - o2_entry_1.energy_per_atom
        assert h2o_form_e == pytest.approx(MU_H2O)

    def test_h_h2o_energy_with_args_multi(self):
        compat = MaterialsProjectAqueousCompatibility(
            o2_energy=-4.9276,
            h2o_energy=-5,
            h2o_adjustments=-0.234,
            solid_compat=None,
        )

        h2o_entry_1 = ComputedEntry(Composition("H2O"), -15)  # -5 eV/atom
        h2o_entry_2 = ComputedEntry(Composition("H4O2"), -6)  # -1 eV/atom
        h2_entry_1 = ComputedEntry(Composition("H8"), -100)  # -12.5 eV/atom
        h2_entry_2 = ComputedEntry(Composition("H2"), -16)  # -8 eV/atom

        compat.process_entries([h2o_entry_1, h2o_entry_2, h2_entry_1, h2_entry_2])

        # Energy differences of H2O and H2 polymorphs should be preserved.
        assert h2o_entry_2.energy_per_atom == pytest.approx(h2o_entry_1.energy_per_atom + 4)
        assert h2_entry_2.energy_per_atom == pytest.approx(h2_entry_1.energy_per_atom + 4.5)

        o2_entry_1 = ComputedEntry(Composition("O2"), -4.9276 * 2)
        o2_entry_1 = compat.process_entries(o2_entry_1)[0]

        h2o_form_e = 3 * h2o_entry_1.energy_per_atom - 2 * h2_entry_1.energy_per_atom - o2_entry_1.energy_per_atom
        assert h2o_form_e == pytest.approx(MU_H2O)

    def test_h_h2o_energy_no_args(self):
        with pytest.warns(UserWarning, match="You did not provide the required O2 and H2O energies."):
            compat = MaterialsProjectAqueousCompatibility(solid_compat=None)

        h2o_entry_1 = ComputedEntry(Composition("H2O"), (-5.195 + 0.234) * 3, correction=-0.234 * 3)  # -5.195 eV/atom
        h2o_entry_2 = ComputedEntry(Composition("H4O2"), -6)  # -1 eV/atom
        h2_entry_1 = ComputedEntry(Composition("H8"), -100)  # -12.5 eV/atom``
        h2_entry_2 = ComputedEntry(Composition("H2"), -16)  # -8 eV/atom
        o2_entry_1 = ComputedEntry(Composition("O2"), -4.9276 * 2)

        with pytest.raises(CompatibilityError, match="Either specify the energies as arguments to "):
            compat.get_adjustments(h2_entry_1)

        compat.process_entries([h2o_entry_1, h2o_entry_2, h2_entry_1, h2_entry_2, o2_entry_1])

        assert compat.o2_energy == -4.9276
        assert compat.h2o_energy == -5.195
        assert compat.h2o_adjustments == -0.234

        # the corrections should preserve the difference in energy among H2O and H2 polymorphs
        assert h2o_entry_2.energy_per_atom == pytest.approx(h2o_entry_1.energy_per_atom + 4.195)
        assert h2_entry_2.energy_per_atom == pytest.approx(h2_entry_1.energy_per_atom + 4.5)

        # the water formation energy, calculated from the lowest energy polymorphs,
        # should equal the experimental value
        h2o_form_e = 3 * h2o_entry_1.energy_per_atom - 2 * h2_entry_1.energy_per_atom - o2_entry_1.energy_per_atom
        assert h2o_form_e == pytest.approx(MU_H2O)

    def test_compound_entropy(self):
        compat = MaterialsProjectAqueousCompatibility(
            o2_energy=-10, h2o_energy=-20, h2o_adjustments=-0.5, solid_compat=None
        )

        o2_entry_1 = ComputedEntry(Composition("O2"), -4.9276 * 2)

        initial_energy = o2_entry_1.energy_per_atom
        o2_entry_1 = compat.process_entries(o2_entry_1)[0]
        processed_energy = o2_entry_1.energy_per_atom

        assert initial_energy - processed_energy == pytest.approx(compat.cpd_entropies["O2"])

    def test_hydrate_adjustment(self):
        compat = MaterialsProjectAqueousCompatibility(
            o2_energy=-10, h2o_energy=-20, h2o_adjustments=-0.5, solid_compat=None
        )

        hydrate_entry = ComputedEntry(Composition("FeH4O2"), -10)

        initial_energy = hydrate_entry.energy
        hydrate_entry = compat.process_entries(hydrate_entry)[0]
        processed_energy = hydrate_entry.energy

        assert initial_energy - processed_energy == pytest.approx(2 * (compat.h2o_adjustments * 3 + MU_H2O))

    def test_processing_entries_inplace(self):
        h2o_entry = ComputedEntry(Composition("H2O"), (-5.195 + 0.234) * 3, correction=-0.234 * 3)  # -5.195 eV/atom
        o2_entry = ComputedEntry(Composition("O2"), -4.9276 * 2)
        # check that compatibility scheme does not change input entries themselves
        entries = [h2o_entry, o2_entry]
        entries_copy = copy.deepcopy(entries)
        MaterialsProjectAqueousCompatibility().process_entries(entries, inplace=False)
        assert all(e.correction == e_copy.correction for e, e_copy in zip(entries, entries_copy))


class AqueousCorrectionTest(unittest.TestCase):
    def setUp(self):
        module_dir = os.path.dirname(os.path.abspath(__file__))
        fp = os.path.join(module_dir, os.path.pardir, "MITCompatibility.yaml")
        self.corr = AqueousCorrection(fp)

    def test_compound_energy(self):
        O2_entry = self.corr.correct_entry(
            ComputedEntry(Composition("O2"), -4.9355 * 2, parameters={"run_type": "GGA"})
        )
        H2_entry = self.corr.correct_entry(ComputedEntry(Composition("H2"), 3, parameters={"run_type": "GGA"}))
        H2O_entry = self.corr.correct_entry(ComputedEntry(Composition("H2O"), 3, parameters={"run_type": "GGA"}))
        H2O_formation_energy = H2O_entry.energy - (H2_entry.energy + O2_entry.energy / 2.0)
        assert H2O_formation_energy == pytest.approx(-2.46)

        entry = ComputedEntry(Composition("H2O"), -16, parameters={"run_type": "GGA"})
        entry = self.corr.correct_entry(entry)
        assert entry.energy == pytest.approx(-14.916)

        entry = ComputedEntry(Composition("H2O"), -24, parameters={"run_type": "GGA"})
        entry = self.corr.correct_entry(entry)
        assert entry.energy == pytest.approx(-14.916)

        entry = ComputedEntry(Composition("Cl"), -24, parameters={"run_type": "GGA"})
        entry = self.corr.correct_entry(entry)
        assert entry.energy == pytest.approx(-24.344373)


class MITAqueousCompatibilityTest(unittest.TestCase):
    def setUp(self):
        self.compat = MITCompatibility(check_potcar_hash=True)
        self.aqcompat = MITAqueousCompatibility(check_potcar_hash=True)
        module_dir = os.path.dirname(os.path.abspath(__file__))
        fp = os.path.join(module_dir, os.path.pardir, "MITCompatibility.yaml")
        self.aqcorr = AqueousCorrection(fp)

    def test_aqueous_compat(self):
        el_li = Element("Li")
        el_o = Element("O")
        el_h = Element("H")
        latt = Lattice.from_parameters(3.565276, 3.565276, 4.384277, 90.000000, 90.000000, 90.000000)
        elts = [el_h, el_h, el_li, el_li, el_o, el_o]
        coords = [
            [0.000000, 0.500000, 0.413969],
            [0.500000, 0.000000, 0.586031],
            [0.000000, 0.000000, 0.000000],
            [0.500000, 0.500000, 0.000000],
            [0.000000, 0.500000, 0.192672],
            [0.500000, 0.000000, 0.807328],
        ]
        struct = Structure(latt, elts, coords)
        lioh_entry = ComputedStructureEntry(
            struct,
            -3,
            parameters={
                "is_hubbard": False,
                "hubbards": None,
                "run_type": "GGA",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Li 17Jan2003",
                        "hash": "65e83282d1707ec078c1012afbd05be8",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                    {
                        "titel": "PAW_PBE H 15Jun2001",
                        "hash": "bb43c666e3d36577264afe07669e9582",
                    },
                ],
            },
        )
        lioh_entry_compat = self.compat.process_entry(lioh_entry)
        lioh_entry_compat_aqcorr = self.aqcorr.correct_entry(lioh_entry_compat)
        lioh_entry_aqcompat = self.aqcompat.process_entry(lioh_entry)
        assert lioh_entry_compat_aqcorr.energy == pytest.approx(lioh_entry_aqcompat.energy)

    def test_potcar_doenst_match_structure(self):
        compat = MITCompatibility()
        el_li = Element("Li")
        el_o = Element("O")
        el_h = Element("H")
        latt = Lattice.from_parameters(3.565276, 3.565276, 4.384277, 90.000000, 90.000000, 90.000000)
        elts = [el_h, el_h, el_li, el_li, el_o, el_o]
        coords = [
            [0.000000, 0.500000, 0.413969],
            [0.500000, 0.000000, 0.586031],
            [0.000000, 0.000000, 0.000000],
            [0.500000, 0.500000, 0.000000],
            [0.000000, 0.500000, 0.192672],
            [0.500000, 0.000000, 0.807328],
        ]
        struct = Structure(latt, elts, coords)

        lioh_entry = ComputedStructureEntry(
            struct,
            -3,
            parameters={
                "is_hubbard": False,
                "hubbards": None,
                "run_type": "GGA",
                "potcar_symbols": ["PAW_PBE Fe 17Jan2003", "PAW_PBE O 08Apr2002", "PAW_PBE H 15Jun2001"],
            },
        )

        assert compat.process_entry(lioh_entry) is None

    def test_msonable(self):
        compat_dict = self.aqcompat.as_dict()
        decoder = MontyDecoder()
        temp_compat = decoder.process_decoded(compat_dict)
        assert isinstance(temp_compat, MITAqueousCompatibility)

    def test_dont_error_on_weird_elements(self):
        entry = ComputedEntry(
            "AmSi",
            -1,
            correction=0.0,
            parameters={
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Am 08May2007",
                        "hash": "ed5eebd8a143e35a0c19e9f8a2c42a93",
                    },
                    {
                        "titel": "PAW_PBE Si 05Jan2001",
                        "hash": "b2b0ea6feb62e7cde209616683b8f7f5",
                    },
                ]
            },
        )
        assert self.compat.process_entry(entry) is None


class CorrectionErrors2020CompatibilityTest(unittest.TestCase):
    def setUp(self):
        warnings.simplefilter("ignore")
        self.compat = MaterialsProject2020Compatibility()

        self.entry1 = ComputedEntry(
            "Fe2O3",
            -1,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 5.3, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )

        self.entry_sulfide = ComputedEntry(
            "FeS",
            -1,
            0.0,
            parameters={
                "is_hubbard": False,
                "run_type": "GGA",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE S 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )

        self.entry2 = ComputedEntry(
            "Fe3O4",
            -2,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 5.3, "O": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE O 08Apr2002",
                        "hash": "7a25bc5b9a5393f46600a4939d357982",
                    },
                ],
            },
        )

        self.entry_fluoride = ComputedEntry(
            "FeF3",
            -2,
            correction=0.0,
            parameters={
                "is_hubbard": True,
                "hubbards": {"Fe": 5.3, "F": 0},
                "run_type": "GGA+U",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Fe_pv 06Sep2000",
                        "hash": "994537de5c4122b7f1b77fb604476db4",
                    },
                    {
                        "titel": "PAW_PBE F 08Apr2002",
                        "hash": "180141c33d032bfbfff30b3bea9d23dd",
                    },
                ],
            },
        )

        self.entry_hydride = ComputedEntry(
            "LiH",
            -2,
            correction=0.0,
            parameters={
                "is_hubbard": False,
                "run_type": "GGA",
                "potcar_spec": [
                    {
                        "titel": "PAW_PBE Li_sv 10Sep2004",
                        "hash": "8245d7383d7556214082aa40a887cd96",
                    },
                    {
                        "titel": "PAW_PBE H 15Jun2001",
                        "hash": "bb43c666e3d36577264afe07669e9582",
                    },
                ],
            },
        )

    def tearDown(self):
        warnings.simplefilter("default")

    def test_errors(self):
        entry1_corrected = self.compat.process_entry(self.entry1)
        assert entry1_corrected.correction_uncertainty == pytest.approx(sqrt((2 * 0.0101) ** 2 + (3 * 0.002) ** 2))

        entry2_corrected = self.compat.process_entry(self.entry2)
        assert entry2_corrected.correction_uncertainty == pytest.approx(sqrt((3 * 0.0101) ** 2 + (4 * 0.002) ** 2))

        entry_sulfide_corrected = self.compat.process_entry(self.entry_sulfide)
        assert entry_sulfide_corrected.correction_uncertainty == pytest.approx(0.0093)

        entry_fluoride_corrected = self.compat.process_entry(self.entry_fluoride)
        assert entry_fluoride_corrected.correction_uncertainty == pytest.approx(sqrt((3 * 0.0026) ** 2 + 0.0101**2))

        entry_hydride_corrected = self.compat.process_entry(self.entry_hydride)
        assert entry_hydride_corrected.correction_uncertainty == pytest.approx(0.0013)


if __name__ == "__main__":
    unittest.main()
