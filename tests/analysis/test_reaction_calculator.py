from __future__ import annotations

import unittest
import warnings
from collections import defaultdict
from math import isnan

import numpy as np
import pytest
from pytest import approx

from pymatgen.analysis.reaction_calculator import BalancedReaction, ComputedReaction, Reaction, ReactionError
from pymatgen.core.composition import Composition
from pymatgen.entries.computed_entries import ComputedEntry


class TestReaction(unittest.TestCase):
    def setUp(self):
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_init(self):
        reactants = [Composition("Fe"), Composition("O2")]
        products = [Composition("Fe2O3")]
        rxn = Reaction(reactants, products)
        assert str(rxn) == "2 Fe + 1.5 O2 -> Fe2O3"
        assert rxn.normalized_repr == "4 Fe + 3 O2 -> 2 Fe2O3"

        d = rxn.as_dict()
        rxn = Reaction.from_dict(d)
        repr, factor = rxn.normalized_repr_and_factor()
        assert repr == "4 Fe + 3 O2 -> 2 Fe2O3"
        assert factor == approx(2)

        reactants = [Composition("FePO4"), Composition("Mn")]
        products = [Composition("FePO4"), Composition("Xe")]
        rxn = Reaction(reactants, products)
        assert str(rxn) == "FePO4 -> FePO4"

        products = [Composition("Ti2 O4"), Composition("O1")]
        reactants = [Composition("Ti1 O2")]
        rxn = Reaction(reactants, products)
        assert str(rxn) == "2 TiO2 -> 2 TiO2"

        reactants = [Composition("FePO4"), Composition("Li")]
        products = [Composition("LiFePO4")]
        rxn = Reaction(reactants, products)
        assert str(rxn) == "FePO4 + Li -> LiFePO4"

        reactants = [Composition("MgO")]
        products = [Composition("MgO")]

        rxn = Reaction(reactants, products)
        assert str(rxn) == "MgO -> MgO"

        reactants = [Composition("Mg")]
        products = [Composition("Mg")]

        rxn = Reaction(reactants, products)
        assert str(rxn) == "Mg -> Mg"

        reactants = [Composition("FePO4"), Composition("LiPO3")]
        products = [Composition("LiFeP2O7")]

        rxn = Reaction(reactants, products)
        assert str(rxn) == "FePO4 + LiPO3 -> LiFeP2O7"

        reactants = [Composition("Na"), Composition("K2O")]
        products = [Composition("Na2O"), Composition("K")]
        rxn = Reaction(reactants, products)
        assert str(rxn) == "2 Na + K2O -> Na2O + 2 K"

        # Test for an old bug which has a problem when excess product is
        # defined.
        products = [Composition("FePO4"), Composition("O")]
        reactants = [Composition("FePO4")]
        rxn = Reaction(reactants, products)

        assert str(rxn) == "FePO4 -> FePO4"

        products = list(map(Composition, ["LiCrO2", "La8Ti8O12", "O2"]))
        reactants = [Composition("LiLa3Ti3CrO12")]
        rxn = Reaction(reactants, products)
        assert str(rxn) == "LiLa3Ti3CrO12 -> LiCrO2 + 1.5 La2Ti2O3 + 2.75 O2"

    def test_rank(self):
        reactants = [Composition("La2Zr2O7"), Composition("LiCoO2")]
        products = [
            Composition("La2O3"),
            Composition("Co2O3"),
            Composition("Li2ZrO3"),
            Composition("Li2O"),
        ]

        assert str(Reaction(reactants, products)) == "La2Zr2O7 + 6 LiCoO2 -> La2O3 + 3 Co2O3 + 2 Li2ZrO3 + Li2O"

        reactants = [Composition("La2O3"), Composition("Co2O3"), Composition("Li2ZrO3")]
        products = [
            Composition("Li2O"),
            Composition("La2Zr2O7"),
            Composition("Li3CoO3"),
        ]
        assert (
            str(Reaction(reactants, products)) == "La2O3 + 0.3333 Co2O3 + 2 Li2ZrO3 -> Li2O + La2Zr2O7 + 0.6667 Li3CoO3"
        )

        reactants = [Composition("La2O3"), Composition("Co2O3"), Composition("Li2ZrO3")]
        products = [
            Composition("Xe"),
            Composition("Li2O"),
            Composition("La2Zr2O7"),
            Composition("Li3CoO3"),
        ]
        assert (
            str(Reaction(reactants, products)) == "La2O3 + 0.3333 Co2O3 + 2 Li2ZrO3 -> Li2O + La2Zr2O7 + 0.6667 Li3CoO3"
        )

        reactants = [Composition("La2O3"), Composition("Co2O3"), Composition("Li2ZrO3")]
        products = [
            Composition("Xe"),
            Composition("Li2O"),
            Composition("La2Zr2O7"),
            Composition("Li3CoO3"),
            Composition("XeNe"),
        ]
        assert (
            str(Reaction(reactants, products)) == "La2O3 + 0.3333 Co2O3 + 2 Li2ZrO3 -> Li2O + La2Zr2O7 + 0.6667 Li3CoO3"
        )

        reactants = [Composition("LiCoO2")]
        products = [
            Composition("La2O3"),
            Composition("Co2O3"),
            Composition("Li2O1"),
            Composition("Li1F1"),
            Composition("Co1F3"),
        ]
        assert str(Reaction(reactants, products)) == "1.667 LiCoO2 + 0.3333 CoF3 -> Co2O3 + 0.3333 Li2O + LiF"

        # this test can fail because of numerical rank calculation issues
        reactants = [Composition("LiCoO2"), Composition("Li2O1")]
        products = [Composition("ZrF4"), Composition("Co2O3")]
        assert str(Reaction(reactants, products)) == "2 LiCoO2 -> Li2O + Co2O3"

    def test_singular_case(self):
        rxn = Reaction(
            [Composition("XeMn"), Composition("Li")],
            [Composition("S"), Composition("LiS2"), Composition("FeCl")],
        )
        assert str(rxn) == "Li + 2 S -> LiS2"

    def test_overdetermined(self):
        with pytest.raises(ReactionError, match="Reaction cannot be balanced"):
            Reaction([Composition("Li")], [Composition("LiO2")])

    def test_scientific_notation(self):
        products = [Composition("FePO3.9999"), Composition("O2")]
        reactants = [Composition("FePO4")]
        rxn = Reaction(reactants, products)
        assert str(rxn) == "FePO4 -> Fe1P1O3.9999 + 5e-05 O2"
        assert rxn == Reaction.from_str(str(rxn))

        rxn2 = Reaction.from_str("FePO4 + 20 CO -> 1e1 O2 + Fe1P1O4 + 20 C")
        assert str(rxn2) == "20 CO -> 20 C + 10 O2"

    def test_equals(self):
        reactants = [Composition("Fe"), Composition("O2")]
        products = [Composition("Fe2O3")]
        rxn = Reaction(reactants, products)
        reactants = [Composition("O2"), Composition("Fe")]
        products = [Composition("Fe2O3")]
        rxn2 = Reaction(reactants, products)
        assert rxn == rxn2

    def test_normalize_to(self):
        products = [Composition("Fe"), Composition("O2")]
        reactants = [Composition("Fe2O3")]
        rxn = Reaction(reactants, products)
        rxn.normalize_to(Composition("Fe"), 3)
        assert str(rxn) == "1.5 Fe2O3 -> 3 Fe + 2.25 O2"

    def test_calculate_energy(self):
        reactants = [Composition("MgO"), Composition("Al2O3")]
        products = [Composition("MgAl2O4")]
        energies = {
            Composition("MgO"): -0.1,
            Composition("Al2O3"): -0.2,
            Composition("MgAl2O4"): -0.5,
        }
        rxn = Reaction(reactants, products)
        assert str(rxn) == "MgO + Al2O3 -> MgAl2O4"
        assert rxn.normalized_repr == "MgO + Al2O3 -> MgAl2O4"
        assert rxn.calculate_energy(energies) == approx(-0.2, abs=1e-5)

    def test_as_entry(self):
        reactants = [Composition("MgO"), Composition("Al2O3")]
        products = [Composition("MgAl2O4")]
        energies = {
            Composition("MgO"): -0.1,
            Composition("Al2O3"): -0.2,
            Composition("MgAl2O4"): -0.5,
        }
        rxn = Reaction(reactants, products)
        entry = rxn.as_entry(energies)
        assert entry.name == "MgO + Al2O3 -> MgAl2O4"
        assert entry.energy == approx(-0.2, abs=1e-5)

        products = [Composition("Fe"), Composition("O2")]
        reactants = [Composition("Fe2O3")]
        rxn = Reaction(reactants, products)
        energies = {
            Composition("Fe"): 0,
            Composition("O2"): 0,
            Composition("Fe2O3"): 0.5,
        }
        entry = rxn.as_entry(energies)
        assert entry.composition == Composition("Fe1 O1.5")
        assert entry.energy == approx(-0.25, abs=1e-5)

    def test_products_reactants(self):
        reactants = [
            Composition("Li3Fe2(PO4)3"),
            Composition("Fe2O3"),
            Composition("O2"),
        ]
        products = [Composition("LiFePO4")]
        energies = {
            Composition("Li3Fe2(PO4)3"): -0.1,
            Composition("Fe2O3"): -0.2,
            Composition("O2"): -0.2,
            Composition("LiFePO4"): -0.5,
        }
        rxn = Reaction(reactants, products)

        assert Composition("O2") in rxn.products, "O not in products!"
        assert Composition("Li3Fe2(PO4)3") in rxn.reactants, "Li3Fe2(PO4)4 not in reactants!"
        assert str(rxn) == "0.3333 Li3Fe2(PO4)3 + 0.1667 Fe2O3 -> 0.25 O2 + LiFePO4"
        assert rxn.normalized_repr == "4 Li3Fe2(PO4)3 + 2 Fe2O3 -> 3 O2 + 12 LiFePO4"
        assert rxn.calculate_energy(energies) == approx(-0.48333333, abs=1e-5)

    def test_to_from_dict(self):
        reactants = [Composition("Fe"), Composition("O2")]
        products = [Composition("Fe2O3")]
        rxn = Reaction(reactants, products)
        d = rxn.as_dict()
        rxn = Reaction.from_dict(d)
        assert rxn.normalized_repr == "4 Fe + 3 O2 -> 2 Fe2O3"

    def test_underdetermined(self):
        reactants = [Composition("Fe"), Composition("O2")]
        products = [Composition("Fe"), Composition("O2")]
        rxn = Reaction(reactants, products)
        assert str(rxn) == "Fe + O2 -> Fe + O2"

        reactants = [
            Composition("Fe"),
            Composition("O2"),
            Composition("Na"),
            Composition("Li"),
            Composition("Cl"),
        ]
        products = [Composition("FeO2"), Composition("NaCl"), Composition("Li2Cl2")]
        rxn = Reaction(reactants, products)
        assert str(rxn) == "Fe + O2 + Na + 2 Li + 1.5 Cl2 -> FeO2 + NaCl + 2 LiCl"

        reactants = [
            Composition("Fe"),
            Composition("Na"),
            Composition("Li2O"),
            Composition("Cl"),
        ]
        products = [
            Composition("LiCl"),
            Composition("Na2O"),
            Composition("Xe"),
            Composition("FeCl"),
            Composition("Mn"),
        ]
        rxn = Reaction(reactants, products)
        # this can't normalize to 1 LiCl + 1 Na2O (not enough O), so chooses LiCl and FeCl
        assert str(rxn) == "Fe + Na + 0.5 Li2O + Cl2 -> LiCl + 0.5 Na2O + FeCl"

    def test_underdetermined_reactants(self):
        reactants = [Composition("Li"), Composition("Cl"), Composition("Cl")]
        products = [Composition("LiCl")]
        rxn = Reaction(reactants, products)
        assert str(rxn) == "Li + 0.25 Cl2 + 0.25 Cl2 -> LiCl"

        reactants = [Composition("LiMnCl3"), Composition("LiCl"), Composition("MnCl2")]
        products = [Composition("Li2MnCl4")]
        rxn = Reaction(reactants, products)
        assert str(rxn) == "LiMnCl3 + 3 LiCl + MnCl2 -> 2 Li2MnCl4"


class TestBalancedReaction(unittest.TestCase):
    def setUp(self):
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_init(self):
        rct = {Composition("K2SO4"): 3, Composition("Na2S"): 1, Composition("Li"): 24}
        prod = {Composition("KNaS"): 2, Composition("K2S"): 2, Composition("Li2O"): 12}
        rxn = BalancedReaction(rct, prod)
        assert str(rxn) is not None

        # Test unbalanced exception
        rct = {Composition("K2SO4"): 1, Composition("Na2S"): 1, Composition("Li"): 24}
        prod = {Composition("KNaS"): 2, Composition("K2S"): 2, Composition("Li2O"): 12}
        with pytest.raises(ReactionError, match="Reaction is unbalanced"):
            BalancedReaction(rct, prod)

    def test_to_from_dict(self):
        rct = {Composition("K2SO4"): 3, Composition("Na2S"): 1, Composition("Li"): 24}
        prod = {Composition("KNaS"): 2, Composition("K2S"): 2, Composition("Li2O"): 12}
        rxn = BalancedReaction(rct, prod)
        d = rxn.as_dict()
        new_rxn = BalancedReaction.from_dict(d)
        for comp in new_rxn.all_comp:
            assert new_rxn.get_coeff(comp) == rxn.get_coeff(comp)

    def test_from_str(self):
        rxn = BalancedReaction({Composition("Li"): 4, Composition("O2"): 1}, {Composition("Li2O"): 2})
        assert rxn == BalancedReaction.from_str("4 Li + O2 -> 2Li2O")

        rxn = BalancedReaction(
            {Composition("Li(NiO2)3"): 1},
            {
                Composition("O2"): 0.5,
                Composition("Li(NiO2)2"): 1,
                Composition("NiO"): 1,
            },
        )

        assert rxn == BalancedReaction.from_str("1.000 Li(NiO2)3 -> 0.500 O2 + 1.000 Li(NiO2)2 + 1.000 NiO")

    def test_remove_spectator_species(self):
        rxn = BalancedReaction(
            {Composition("Li"): 4, Composition("O2"): 1, Composition("Na"): 1},
            {Composition("Li2O"): 2, Composition("Na"): 1},
        )

        assert Composition("Na") not in rxn.all_comp


class TestComputedReaction(unittest.TestCase):
    def setUp(self):
        d = [
            {
                "correction": 0,
                "data": {},
                "energy": -108.56492362,
                "parameters": {},
                "composition": {"Li": 54},
            },
            {
                "correction": 0,
                "data": {},
                "energy": -577.94689128,
                "parameters": {},
                "composition": {"O": 32, "Li": 64},
            },
            {
                "correction": 0,
                "data": {},
                "energy": -17.02844794,
                "parameters": {},
                "composition": {"O": 2},
            },
            {
                "correction": 0,
                "data": {},
                "energy": -959.64693323,
                "parameters": {},
                "composition": {"O": 72, "Li": 72},
            },
        ]
        entries = []
        for e in d:
            entries.append(ComputedEntry.from_dict(e))
        rcts = list(filter(lambda e: e.composition.reduced_formula in ["Li", "O2"], entries))
        prods = list(filter(lambda e: e.composition.reduced_formula == "Li2O2", entries))

        self.rxn = ComputedReaction(rcts, prods)

    def test_calculated_reaction_energy(self):
        assert self.rxn.calculated_reaction_energy == approx(-5.60748821935)

    def test_calculated_reaction_energy_uncertainty(self):
        d = [
            {
                "correction": 0,
                "data": {},
                "energy": -108.56492362,
                "parameters": {},
                "composition": {"Li": 54},
            },
            {
                "correction": 0,
                "data": {},
                "energy": -17.02844794,
                "parameters": {},
                "composition": {"O": 2},
            },
            {
                "@module": "pymatgen.entries.computed_entries",
                "@class": "ComputedEntry",
                "energy": -38.76889738,
                "composition": defaultdict(float, {"Li": 4, "O": 4}),
                "energy_adjustments": [
                    {
                        "@module": "pymatgen.entries.computed_entries",
                        "@class": "ConstantEnergyAdjustment",
                        "@version": "2020.6.8",
                        "value": -1.864,
                        "uncertainty": 0.0744,
                        "name": "MP2020 Composition Correction",
                        "cls": {
                            "@module": "pymatgen.entries.compatibility",
                            "@class": "MaterialsProject2020Compatibility",
                            "@version": "2020.6.8",
                            "compat_type": "Advanced",
                            "correct_peroxide": True,
                            "check_potcar_hash": False,
                        },
                        "description": "Constant energy adjustment (-1.864 eV)",
                    }
                ],
                "parameters": {
                    "run_type": "GGA",
                    "is_hubbard": False,
                    "pseudo_potential": {
                        "functional": "PBE",
                        "labels": ["Li_sv", "O"],
                        "pot_type": "paw",
                    },
                    "hubbards": {},
                    "potcar_symbols": ["PBE Li_sv", "PBE O"],
                    "oxide_type": "peroxide",
                },
                "data": {"oxide_type": "peroxide"},
                "entry_id": "mp-841",
                "correction": -1.864,
            },
        ]
        entries = []
        for e in d:
            entries.append(ComputedEntry.from_dict(e))
        rcts = list(filter(lambda e: e.composition.reduced_formula in ["Li", "O2"], entries))
        prods = list(filter(lambda e: e.composition.reduced_formula == "Li2O2", entries))

        rxn_with_uncertainty = ComputedReaction(rcts, prods)
        assert rxn_with_uncertainty.calculated_reaction_energy_uncertainty == approx(0.5 * 0.0744)

    def test_calculated_reaction_energy_uncertainty_for_no_uncertainty(self):
        # test that reaction_energy_uncertainty property doesn't cause errors
        # when products/reactants have no uncertainties
        assert self.rxn.calculated_reaction_energy_uncertainty == 0

    def test_calculated_reaction_energy_uncertainty_for_nan(self):
        # test that reaction_energy_uncertainty property is nan when the uncertainty
        # for any product/reactant is nan
        d = [
            {
                "correction": 0,
                "data": {},
                "energy": -108.56492362,
                "parameters": {},
                "composition": {"Li": 54},
            },
            {
                "correction": 0,
                "data": {},
                "energy": -17.02844794,
                "parameters": {},
                "composition": {"O": 2},
            },
            {
                "@module": "pymatgen.entries.computed_entries",
                "@class": "ComputedEntry",
                "energy": -38.76889738,
                "composition": defaultdict(float, {"Li": 4, "O": 4}),
                "energy_adjustments": [
                    {
                        "@module": "pymatgen.entries.computed_entries",
                        "@class": "ConstantEnergyAdjustment",
                        "@version": "2020.6.8",
                        "value": -1.864,
                        "uncertainty": np.nan,
                        "name": "MP2020 Composition Correction",
                        "cls": {
                            "@module": "pymatgen.entries.compatibility",
                            "@class": "MaterialsProject2020Compatibility",
                            "@version": "2020.6.8",
                            "compat_type": "Advanced",
                            "correct_peroxide": True,
                            "check_potcar_hash": False,
                        },
                        "description": "Constant energy adjustment (-1.864 eV)",
                    }
                ],
                "parameters": {
                    "run_type": "GGA",
                    "is_hubbard": False,
                    "pseudo_potential": {
                        "functional": "PBE",
                        "labels": ["Li_sv", "O"],
                        "pot_type": "paw",
                    },
                    "hubbards": {},
                    "potcar_symbols": ["PBE Li_sv", "PBE O"],
                    "oxide_type": "peroxide",
                },
                "data": {"oxide_type": "peroxide"},
                "entry_id": "mp-841",
                "correction": -1.864,
            },
        ]
        entries = []
        for e in d:
            entries.append(ComputedEntry.from_dict(e))
        rcts = list(filter(lambda e: e.composition.reduced_formula in ["Li", "O2"], entries))
        prods = list(filter(lambda e: e.composition.reduced_formula == "Li2O2", entries))

        rxn_with_uncertainty = ComputedReaction(rcts, prods)
        assert isnan(rxn_with_uncertainty.calculated_reaction_energy_uncertainty)

    def test_init(self):
        assert str(self.rxn) == "2 Li + O2 -> Li2O2"

    def test_to_from_dict(self):
        d = self.rxn.as_dict()
        new_rxn = ComputedReaction.from_dict(d)
        assert str(new_rxn) == "2 Li + O2 -> Li2O2"

    def test_all_entries(self):
        for c, e in zip(self.rxn.coeffs, self.rxn.all_entries):
            if c > 0:
                assert e.composition.reduced_formula == "Li2O2"
                assert e.energy == approx(-959.64693323)
