from __future__ import annotations

import unittest
import warnings

import numpy as np
import pytest
from matplotlib.figure import Figure as mpl_figure
from pandas import DataFrame
from plotly.graph_objects import Figure as plotly_figure
from scipy.spatial import ConvexHull

from pymatgen.analysis.interface_reactions import (
    GrandPotentialInterfacialReactivity,
    InterfacialReactivity,
)
from pymatgen.analysis.phase_diagram import GrandPotentialPhaseDiagram, PhaseDiagram
from pymatgen.analysis.reaction_calculator import Reaction
from pymatgen.core.composition import Composition, Element
from pymatgen.entries.computed_entries import ComputedEntry


class InterfaceReactionTest(unittest.TestCase):
    def setUp(self):
        self.entries = [
            ComputedEntry(Composition("Li"), 0),
            ComputedEntry(Composition("Mn"), 0),
            ComputedEntry(Composition("O2"), 0),
            ComputedEntry(Composition("MnO2"), -10),
            ComputedEntry(Composition("Mn2O4"), -60),
            ComputedEntry(Composition("MnO3"), 20),
            ComputedEntry(Composition("Li2O"), -10),
            ComputedEntry(Composition("Li2O2"), -8),
            ComputedEntry(Composition("LiMnO2"), -30),
        ]
        self.pd = PhaseDiagram(self.entries)

        chempots = {Element("Li"): -3}
        self.gpd = GrandPotentialPhaseDiagram(self.entries, chempots)

        ir_0 = InterfacialReactivity(
            c1=Composition("O2"),
            c2=Composition("Mn"),
            pd=self.pd,
            norm=False,
            use_hull_energy=False,
        )
        ir_1 = GrandPotentialInterfacialReactivity(
            c1=Composition("MnO2"),
            c2=Composition("Mn"),
            grand_pd=self.gpd,
            pd_non_grand=self.pd,
            norm=False,
            include_no_mixing_energy=True,
            use_hull_energy=False,
        )
        ir_2 = GrandPotentialInterfacialReactivity(
            c1=Composition("Mn"),
            c2=Composition("O2"),
            grand_pd=self.gpd,
            pd_non_grand=self.pd,
            norm=True,
            include_no_mixing_energy=True,
            use_hull_energy=False,
        )
        ir_3 = GrandPotentialInterfacialReactivity(
            c1=Composition("Li2O"),
            c2=Composition("Mn"),
            grand_pd=self.gpd,
            norm=False,
            include_no_mixing_energy=True,
            pd_non_grand=self.pd,
            use_hull_energy=False,
        )
        ir_4 = GrandPotentialInterfacialReactivity(
            c1=Composition("Mn"),
            c2=Composition("O2"),
            grand_pd=self.gpd,
            norm=True,
            include_no_mixing_energy=False,
            pd_non_grand=self.pd,
            use_hull_energy=False,
        )
        ir_5 = GrandPotentialInterfacialReactivity(
            c1=Composition("Mn"),
            c2=Composition("Li2O"),
            grand_pd=self.gpd,
            pd_non_grand=self.pd,
            norm=True,
            include_no_mixing_energy=True,
            use_hull_energy=False,
        )
        ir_6 = InterfacialReactivity(
            c1=Composition("Li2O2"),
            c2=Composition("Li"),
            pd=self.pd,
            norm=False,
            use_hull_energy=True,
        )
        ir_7 = InterfacialReactivity(
            c1=Composition("Li2O2"),
            c2=Composition("Li"),
            pd=self.pd,
            norm=False,
            use_hull_energy=False,
        )
        ir_8 = GrandPotentialInterfacialReactivity(
            c1=Composition("Li2O2"),
            c2=Composition("MnO2"),
            grand_pd=self.gpd,
            pd_non_grand=self.pd,
            norm=False,
            include_no_mixing_energy=False,
            use_hull_energy=True,
        )
        ir_9 = GrandPotentialInterfacialReactivity(
            c1=Composition("Li2O2"),
            c2=Composition("MnO2"),
            grand_pd=self.gpd,
            pd_non_grand=self.pd,
            norm=False,
            include_no_mixing_energy=False,
            use_hull_energy=False,
        )
        ir_10 = InterfacialReactivity(
            Composition("O2"),
            Composition("Mn"),
            pd=self.pd,
            norm=True,
            use_hull_energy=False,
        )
        ir_11 = GrandPotentialInterfacialReactivity(
            Composition("Li2O2"),
            Composition("Li2O2"),
            grand_pd=self.gpd,
            norm=True,
            include_no_mixing_energy=True,
            pd_non_grand=self.pd,
            use_hull_energy=False,
        )
        ir_12 = InterfacialReactivity(
            Composition("Li2O2"),
            Composition("Li2O2"),
            pd=self.pd,
            norm=True,
            use_hull_energy=False,
        )
        with pytest.raises(Exception) as context1:
            _ = InterfacialReactivity(Composition("Li2O2"), Composition("Li"), pd=self.gpd, norm=True)
            assert (
                str(context1.exception) == "Please use the GrandPotentialInterfacialReactivity "
                "class for interfacial reactions with open elements!"
            )
        with pytest.raises(Exception) as context2:
            _ = GrandPotentialInterfacialReactivity(
                Composition("O2"),
                Composition("Mn"),
                grand_pd=self.gpd,
                pd_non_grand=None,
                norm=False,
                include_no_mixing_energy=True,
            )
            assert str(context2.exception) == "Please provide non-grand phase diagram to compute no_mixing_energy!"

        self.ir = [ir_0, ir_1, ir_2, ir_3, ir_4, ir_5, ir_6, ir_7, ir_8, ir_9, ir_10, ir_11, ir_12]

    def test_get_entry_energy(self):
        comp = Composition("MnO3")
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            energy = InterfacialReactivity._get_entry_energy(self.pd, comp)
            assert len(w) == 1
            assert (
                "The reactant MnO3 has no matching entry with"
                " negative formation energy, instead convex "
                "hull energy for this composition will be used"
                " for reaction energy calculation." in str(w[-1].message)
            )
        test1 = np.isclose(energy, -30, atol=1e-03)
        assert test1, f"_get_entry_energy: energy for {comp.reduced_formula} is wrong!"
        # Test normal functionality
        comp = Composition("MnO2")
        test2 = np.isclose(InterfacialReactivity._get_entry_energy(self.pd, comp), -30, atol=1e-03)
        assert test2, f"_get_entry_energy: energy for {comp.reduced_formula} is wrong!"

    def test_get_grand_potential(self):
        comp = Composition("LiMnO2")
        # Test non-normalized case
        test1 = np.isclose(self.ir[1]._get_grand_potential(comp), -27, atol=1e-03)
        assert test1, "_get_grand_potential: Non-normalized case gets error!"

        # Test normalized case
        test2 = np.isclose(self.ir[2]._get_grand_potential(comp), -9, atol=1e-03)
        assert test2, "_get_grand_potential: Normalized case gets error!"

        comp2 = Composition("Li2O2")
        # Test use_hull_energy option.
        test3 = np.isclose(self.ir[8]._get_grand_potential(comp2), -4, atol=1e-03)
        assert test3, "_get_grand_potential: get hull energy gets error!"

        test4 = np.isclose(self.ir[9]._get_grand_potential(comp2), -2, atol=1e-03)
        assert test4, f"_get_grand_potential: gets error for {comp2.reduced_formula}!"

    def test_get_energy(self):
        test1 = np.isclose(self.ir[0]._get_energy(0.5), -15, atol=1e-03)
        assert test1, "_get_energy: phase diagram gets error!"

        test2 = np.isclose(self.ir[3]._get_energy(0.6666666), -7.333333, atol=1e-03)
        assert test2, "_get_energy: grand canonical phase diagram gets error!"

        test3 = np.isclose(self.ir[6]._get_energy(0.3333333), -3.333333, atol=1e-03)
        assert test3, "_get_energy: convex hull energy gets error. "

        test4 = np.isclose(self.ir[7]._get_energy(0.3333333), -4, atol=1e-03)
        assert test4, "_get_energy: gets error. "

    def test_get_reaction(self):
        assert (
            str(self.ir[0]._get_reaction(0.5)) == "0.5 Mn + 0.5 O2 -> 0.5 MnO2"
        ), "_get_reaction: reaction not involving chempots species gets error!"
        assert (
            str(self.ir[3]._get_reaction(0.666666)) == "0.5 Li2O + 0.5 Mn -> Li + 0.25 MnO2 + 0.25 Mn"
        ), "_get_reaction: reaction involving chempots species gets error!"

    def test_get_get_elmt_amt_in_rxt(self):
        rxt1 = Reaction(
            [Composition("Mn"), Composition("O2"), Composition("Li")],
            [Composition("LiMnO2")],
        )
        test1 = np.isclose(self.ir[2]._get_elmt_amt_in_rxn(rxt1), 3)
        assert test1, "_get_get_elmt_amt_in_rxt: gpd elements amounts gets error!"

        rxt2 = rxt1
        rxt2.normalize_to(Composition("Li"), 0.5)
        test2 = np.isclose(self.ir[2]._get_elmt_amt_in_rxn(rxt2), 1.5)
        assert test2, "_get_get_elmt_amt_in_rxt: gpd elements amounts gets error!"

        rxt3 = Reaction([Composition("O2"), Composition("Li")], [Composition("Li2O")])
        # Li is not counted
        test3 = np.isclose(self.ir[2]._get_elmt_amt_in_rxn(rxt3), 1)
        assert test3, "_get_get_elmt_amt_in_rxt: gpd elements amounts gets error!"

        # Li is counted
        test4 = np.isclose(self.ir[6]._get_elmt_amt_in_rxn(rxt3), 3)
        assert test4, "_get_get_elmt_amt_in_rxt: pd elements amounts gets error!"

    def test_convert(self):
        test_array = [(0.5, 1, 3), (0.4, 2, 3), (0, 1, 9), (1, 2, 7)]
        result = [InterfacialReactivity._convert(x, f1, f2) for x, f1, f2 in test_array]
        answer = [0.75, 0.5, 0, 1]
        assert np.allclose(result, answer), f"_convert: conversion gets error! {answer} expected, but gets {result}"

    def test_reverse_convert(self):
        test_array = [(0.5, 1, 3), (0.4, 2, 3), (0, 1, 9), (1, 2, 7)]
        result = [InterfacialReactivity._reverse_convert(x, f1, f2) for x, f1, f2 in test_array]
        answer = [0.25, 0.3076923, 0, 1]
        assert np.allclose(result, answer), f"_convert: conversion gets error! {answer} expected, but gets {result}"

    def test_products_property(self):
        test1 = sorted(self.ir[0].products) == sorted(["MnO2", "O2", "Mn"])
        assert test1, "decomposition products gets error for reaction not involving chempots species!"

        test2 = sorted(self.ir[3].products) == sorted(["Li", "MnO2", "Mn", "Li2O"])
        assert test2, "decomposition products gets error for reaction involving chempots species!"

    def test_get_kinks(self):
        def test_get_kinks_helper(
            ir,
            index_expect,
            x_kink_expect,
            energy_kink_expect,
            react_kink_expect,
            energy_per_rxt_kink_expect,
        ):
            kinks = ir.get_kinks()
            index = [i[0] for i in kinks]
            x_kink = [i[1] for i in kinks]
            energy_kink = [i[2] for i in kinks]
            react_kink = [str(i[3]) for i in kinks]
            energy_per_rxt_kink = [i[4] for i in kinks]

            test1 = index == index_expect
            assert test1, "get_kinks:index gets error!"

            test2 = np.allclose(x_kink, x_kink_expect)
            assert test2, "get_kinks:x kinks gets error!"

            test3 = np.allclose(energy_kink, energy_kink_expect)
            assert test3, "get_kinks:energy kinks gets error!"

            # Testing reaction strings are hard,
            # as species could be arranged in random order.
            test4 = len(react_kink) == len(react_kink_expect)
            assert test4, (
                f"get_kinks: reaction kinks gets error for {ir.c1_original.reduced_formula} and "
                f"{ir.c2_original.reduced_formula} reaction!"
            )

            test5 = np.allclose(energy_per_rxt_kink, energy_per_rxt_kink_expect)
            assert test5, "get_kinks: energy_per_rxt_kinks gets error!"

        test_get_kinks_helper(
            self.ir[0],
            [1, 2, 3],
            [0, 0.5, 1],
            [0, -15, 0],
            ["Mn -> Mn", "0.5 O2 + 0.5 Mn -> 0.5 MnO2", "O2 -> O2"],
            [0, -15 * InterfacialReactivity.EV_TO_KJ_PER_MOL, 0],
        )
        test_get_kinks_helper(
            self.ir[10],
            [1, 2, 3],
            [0, 0.66667, 1],
            [0, -10, 0],
            ["Mn -> Mn", "0.5 O2 + 0.5 Mn -> 0.5 MnO2", "O2 -> O2"],
            [0, -15 * InterfacialReactivity.EV_TO_KJ_PER_MOL, 0],
        )
        test_get_kinks_helper(
            self.ir[11],
            [1, 2],
            [0, 1],
            [-3, -3],
            ["Li2O2 + 2 Li -> 2 Li2O", "Li2O2 + 2 Li -> 2 Li2O"],
            [-6 * InterfacialReactivity.EV_TO_KJ_PER_MOL] * 2,
        )
        test_get_kinks_helper(
            self.ir[12],
            [1, 2],
            [0, 1],
            [-0.5, -0.5],
            ["Li2O2 -> Li2O + 0.5 O2", "Li2O2 -> Li2O + 0.5 O2"],
            [-2 * InterfacialReactivity.EV_TO_KJ_PER_MOL] * 2,
        )

    def test_convexity(self):
        def test_convexity_helper(ir):
            lst = list(ir.get_kinks())
            x_kink = [i[1] for i in lst]
            energy_kink = [i[2] for i in lst]
            points = list(zip(x_kink, energy_kink))
            if len(points) >= 3:
                # To test convexity of the plot, construct convex hull from
                # the kinks and make sure
                # 1. all points are below the end points
                # 2. all points are on the convex hull.
                relative_vectors_1 = [(x - x_kink[0], e - energy_kink[0]) for x, e in points]
                relative_vectors_2 = [(x - x_kink[-1], e - energy_kink[-1]) for x, e in points]
                relative_vectors = zip(relative_vectors_1, relative_vectors_2)
                positions = [np.cross(v1, v2) for v1, v2 in relative_vectors]
                test1 = np.all(np.array(positions) <= 0)

                hull = ConvexHull(points)
                test2 = len(hull.vertices) == len(points)
                assert test1 and test2, "Error: Generating non-convex plot!"

        for ir in self.ir:
            test_convexity_helper(ir)

    def test_get_original_composition_ratio(self):
        # expected reaction1: 0.5 O2 + 0.5 Mn -> 0.5 MnO2
        reaction1 = self.ir[0]._get_reaction(0.5)
        test1 = np.isclose(self.ir[0]._get_original_composition_ratio(reaction1), 0.5)
        assert test1, "_get_original_composition_ratio: reaction not involving chempots species gets error!"

        #  expected reaction2: 0.5 Mn + 0.5 Li2O -> Li + 0.25 MnO2 + 0.25 Mn
        reaction2 = self.ir[3]._get_reaction(0.666666)
        test2 = np.isclose(self.ir[3]._get_original_composition_ratio(reaction2), 0.5)
        assert test2, "_get_original_composition_ratio: reaction involving chempots species gets error!"

    def test_get_critical_original_kink_ratio(self):
        test1 = np.allclose(self.ir[0].get_critical_original_kink_ratio(), [0, 0.5, 1])
        assert test1, "get_critical_original_kink_ratio: gets error!"
        test2 = np.allclose(self.ir[10].get_critical_original_kink_ratio(), [0, 0.5, 1])
        assert test2, "get_critical_original_kink_ratio: gets error!"
        test3 = np.allclose(self.ir[11].get_critical_original_kink_ratio(), [0, 1])
        assert test3, "get_critical_original_kink_ratio: gets error!"
        test4 = np.allclose(self.ir[2].get_critical_original_kink_ratio(), [0, 0.5, 1])
        assert test4, "get_critical_original_kink_ratio: gets error!"
        test5 = np.allclose(self.ir[3].get_critical_original_kink_ratio(), [0, 0.66666, 1])
        assert test5, "get_critical_original_kink_ratio: gets error!"

    def test_labels(self):
        d_pymg = self.ir[0].labels
        d_test = {
            1: "x= 0.0 energy in eV/atom = 0.0 Mn -> Mn",
            2: "x= 0.5 energy in eV/atom = -15.0 0.5 Mn + 0.5 O2 -> 0.5 MnO2",
            3: "x= 1.0 energy in eV/atom = 0.0 O2 -> O2",
        }

        assert d_pymg == d_test, (
            "labels:label does not match for interfacial system "
            f"with {self.ir[0].c1_original.reduced_formula} and {self.ir[0].c2_original.reduced_formula}."
        )

    def test_plot(self):
        for i in self.ir:
            fig = i.plot(backend="matplotlib")
            assert fig, isinstance(fig, mpl_figure)

            fig = i.plot(backend="plotly")
            assert isinstance(fig, plotly_figure)

    def test_get_dataframe(self):
        for i in self.ir:
            df = i.get_dataframe()
            assert isinstance(df, DataFrame)

    def test_minimum(self):
        answer = [
            (0.5, -15),
            (0, 0),
            (0.3333333, -10),
            (0.6666666, -7.333333),
            (0.3333333, -7.333333),
            (0.1428571, -7.333333),
            (0.3333333, -3.333333),
            (0.3333333, -4.0),
        ]
        for i, j in zip(self.ir, answer):
            assert np.allclose(i.minimum, j), (
                f"minimum: the system with {i.c1_original.reduced_formula} and {i.c2_original.reduced_formula} "
                f"gets error!{j} expected, but gets {i.minimum}"
            )

    def test_get_no_mixing_energy(self):
        answer = [
            [("MnO2 (eV/f.u.)", 0.0), ("Mn (eV/f.u.)", 0.0)],
            [("Mn (eV/atom)", 0.0), ("O2 (eV/atom)", -4.0)],
            [("Li2O (eV/f.u.)", 0.0), ("Mn (eV/f.u.)", 0.0)],
            [("Mn (eV/atom)", 0.0), ("O2 (eV/atom)", -4.0)],
            [("Mn (eV/atom)", 0.0), ("Li2O (eV/atom)", 0.0)],
        ]

        def name_lst(lst):
            return lst[0][0], lst[1][0]

        def energy_lst(lst):
            return lst[0][1], lst[1][1]

        result_info = [i.get_no_mixing_energy() for i in self.ir if i.grand]
        print(result_info)
        for i, j in zip(result_info, answer):
            assert name_lst(i) == name_lst(
                j
            ), f"get_no_mixing_energy: names get error, {name_lst(j)} expected but gets {name_lst(i)}"
            assert np.allclose(energy_lst(i), energy_lst(j)), (
                "get_no_mixing_energy: "
                "no_mixing energies get error, "
                f"{energy_lst(j)} expected but gets {energy_lst(i)}"
            )

    def test_get_chempot_correction(self):
        # test data from fig. 6 in ref:
        # Prediction of A2BX4 metal-chalcogenide compounds via
        # first-principles thermodynamics, PHYSICAL REVIEW B 86, 014109 (2012)

        # test pressure effect.
        actual = InterfacialReactivity.get_chempot_correction("O", 298.15, 100e5)
        expect = 0.05916
        assert np.isclose(
            actual, expect, atol=1e-2
        ), f"get_chempot_correction gets error, {expect} expected but gets {actual}"
        # test temperature effect.
        actual_2 = InterfacialReactivity.get_chempot_correction("O", 1000, 1e5)
        expect_2 = -0.82352
        assert np.isclose(
            actual_2, expect_2, atol=1e-2
        ), f"get_chempot_correction gets error, {expect_2} expected but gets {actual_2}"

        actual_3 = InterfacialReactivity.get_chempot_correction("O", 500, 1e5)
        expect_3 = -0.223
        assert np.isclose(
            actual_3, expect_3, atol=1e-2
        ), f"get_chempot_correction gets error, {expect_3} expected but gets {actual_3}"
        # test mixed effect.
        actual_4 = InterfacialReactivity.get_chempot_correction("O", 1000, 1e-25)
        expect_4 = -3.800
        assert np.isclose(
            actual_4, expect_4, atol=1e-2
        ), f"get_chempot_correction gets error, {expect_4} expected but gets {actual_4}"
        actual_5 = InterfacialReactivity.get_chempot_correction("O", 1250, 1e-25)
        expect_5 = -4.86
        assert np.isclose(
            actual_5, expect_5, atol=1e-2
        ), f"get_chempot_correction gets error, {expect_5} expected but gets {actual_5}"
        actual_6 = InterfacialReactivity.get_chempot_correction("O", 1500, 1e-25)
        expect_6 = -5.928
        assert np.isclose(
            actual_6, expect_6, atol=1e-2
        ), f"get_chempot_correction gets error, {expect_6} expected but gets {actual_6}"
        actual_7 = InterfacialReactivity.get_chempot_correction("O", 1000, 1e-15)
        expect_7 = -2.808
        assert np.isclose(
            actual_7, expect_7, atol=1e-2
        ), f"get_chempot_correction gets error, {expect_7} expected but gets {actual_7}"
        # test non-gas phase.
        actual_8 = InterfacialReactivity.get_chempot_correction("Li", 1000, 1e15)
        expect_8 = 0
        assert np.isclose(
            actual_8, expect_8, atol=1e-5
        ), f"get_chempot_correction gets error, {expect_8} expected but gets {actual_8}"


if __name__ == "__main__":
    unittest.main()
