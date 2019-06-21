from __future__ import division, unicode_literals

import warnings
import unittest as unittest
import numpy as np

from scipy.spatial import ConvexHull
from pymatgen import Composition
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.analysis.phase_diagram import PhaseDiagram, \
    GrandPotentialPhaseDiagram
from pymatgen.analysis.reaction_calculator import Reaction
from pymatgen.analysis.interface_reactions import InterfacialReactivity


class InterfaceReactionTest(unittest.TestCase):
    def setUp(self):
        self.entries = [ComputedEntry(Composition('Li'), 0),
                        ComputedEntry(Composition('Mn'), 0),
                        ComputedEntry(Composition('O2'), 0),
                        ComputedEntry(Composition('MnO2'), -10),
                        ComputedEntry(Composition('Mn2O4'), -60),
                        ComputedEntry(Composition('MnO3'), 20),
                        ComputedEntry(Composition('Li2O'), -10),
                        ComputedEntry(Composition('Li2O2'), -8),
                        ComputedEntry(Composition('LiMnO2'), -30)
                        ]
        self.pd = PhaseDiagram(self.entries)
        chempots = {'Li': -3}
        self.gpd = GrandPotentialPhaseDiagram(self.entries, chempots)
        self.ir = []
        # ir[0]
        self.ir.append(
            InterfacialReactivity(Composition('O2'), Composition('Mn'),
                                  self.pd, norm=0, include_no_mixing_energy=0,
                                  pd_non_grand=None, use_hull_energy=False))
        # ir[1]
        self.ir.append(
            InterfacialReactivity(Composition('MnO2'), Composition('Mn'),
                                  self.gpd, norm=0, include_no_mixing_energy=1,
                                  pd_non_grand=self.pd, use_hull_energy=False))
        # ir[2]
        self.ir.append(
            InterfacialReactivity(Composition('Mn'), Composition('O2'),
                                  self.gpd, norm=1, include_no_mixing_energy=1,
                                  pd_non_grand=self.pd, use_hull_energy=False))
        # ir[3]
        self.ir.append(
            InterfacialReactivity(Composition('Li2O'), Composition('Mn'),
                                  self.gpd, norm=0, include_no_mixing_energy=1,
                                  pd_non_grand=self.pd, use_hull_energy=False))
        # ir[4]
        self.ir.append(
            InterfacialReactivity(Composition('Mn'), Composition('O2'),
                                  self.gpd, norm=1, include_no_mixing_energy=0,
                                  pd_non_grand=self.pd, use_hull_energy=False))
        # ir[5]
        self.ir.append(
            InterfacialReactivity(Composition('Mn'), Composition('Li2O'),
                                  self.gpd, norm=1, include_no_mixing_energy=1,
                                  pd_non_grand=self.pd, use_hull_energy=False))
        # ir[6]
        self.ir.append(
            InterfacialReactivity(Composition('Li2O2'), Composition('Li'),
                                  self.pd, norm=0, include_no_mixing_energy=0,
                                  pd_non_grand=None, use_hull_energy=True))
        # ir[7]
        self.ir.append(
            InterfacialReactivity(Composition('Li2O2'), Composition('Li'),
                                  self.pd, norm=0, include_no_mixing_energy=0,
                                  pd_non_grand=None, use_hull_energy=False))
        # ir[8]
        self.ir.append(
            InterfacialReactivity(Composition('Li2O2'), Composition('MnO2'),
                                  self.gpd, norm=0, include_no_mixing_energy=0,
                                  pd_non_grand=self.pd, use_hull_energy=True))
        # ir[9]
        self.ir.append(
            InterfacialReactivity(Composition('Li2O2'), Composition('MnO2'),
                                  self.gpd, norm=0, include_no_mixing_energy=0,
                                  pd_non_grand=self.pd, use_hull_energy=False))
        # ir[10]
        self.ir.append(
            InterfacialReactivity(Composition('O2'), Composition('Mn'),
                                  self.pd, norm=1, include_no_mixing_energy=0,
                                  pd_non_grand=None, use_hull_energy=False))
        # ir[11]
        self.ir.append(
            InterfacialReactivity(Composition('Li2O2'), Composition('Li2O2'),
                                  self.gpd, norm=1, include_no_mixing_energy=1,
                                  pd_non_grand=self.pd, use_hull_energy=False))
        # ir[12]
        self.ir.append(
            InterfacialReactivity(Composition('Li2O2'), Composition('Li2O2'),
                                  self.pd, norm=1, include_no_mixing_energy=0,
                                  pd_non_grand=None, use_hull_energy=False))

        with self.assertRaises(Exception) as context1:
            self.ir.append(
                InterfacialReactivity(Composition('Li2O2'), Composition('Li'),
                                      self.pd, norm=1,
                                      include_no_mixing_energy=1,
                                      pd_non_grand=None))
        self.assertTrue(
            'Please provide grand phase diagram '
            'to compute no_mixing_energy!' == str(context1.exception))

        with self.assertRaises(Exception) as context2:
            self.ir.append(
                InterfacialReactivity(Composition('O2'), Composition('Mn'),
                                      self.gpd, norm=0,
                                      include_no_mixing_energy=1,
                                      pd_non_grand=None))
        self.assertTrue(
            'Please provide non-grand phase diagram '
            'to compute no_mixing_energy!' == str(context2.exception))

    def test_get_entry_energy(self):
        # Test warning message.
        comp = Composition('MnO3')
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            energy = InterfacialReactivity._get_entry_energy(self.pd, comp)
            self.assertTrue(len(w) == 1)
            self.assertTrue("The reactant MnO3 has no matching entry with"
                            " negative formation energy, instead convex "
                            "hull energy for this composition will be used"
                            " for reaction energy calculation."
                            in str(w[-1].message))
        test1 = np.isclose(energy, -30, atol=1e-03)
        self.assertTrue(test1,
                        '_get_entry_energy: energy for {} is wrong!'.format(
                            comp.reduced_formula))
        # Test normal functionality
        comp = Composition('MnO2')
        test2 = np.isclose(
            InterfacialReactivity._get_entry_energy(self.pd, comp),
            -30, atol=1e-03)
        self.assertTrue(test2,
                        '_get_entry_energy: energy for {} is wrong!'.format(
                            comp.reduced_formula))

    def test_get_grand_potential(self):
        comp = Composition('LiMnO2')
        # Test non-normalized case
        test1 = np.isclose(self.ir[1]._get_grand_potential(comp), -27,
                           atol=1e-03)
        self.assertTrue(test1,
                        '_get_grand_potential: '
                        'Non-normalized case gets error!')

        # Test normalized case
        test2 = np.isclose(self.ir[2]._get_grand_potential(comp), -9,
                           atol=1e-03)
        self.assertTrue(test2,
                        '_get_grand_potential: '
                        'Normalized case gets error!')

        comp2 = Composition('Li2O2')
        # Test use_hull_energy option.
        test3 = np.isclose(self.ir[8]._get_grand_potential(comp2), -4,
                           atol=1e-03)
        self.assertTrue(test3,
                        '_get_grand_potential: '
                        'get hull energy gets error!')

        test4 = np.isclose(self.ir[9]._get_grand_potential(comp2), -2,
                           atol=1e-03)
        self.assertTrue(test4,
                        '_get_grand_potential: '
                        'gets error for {}!'.format(comp2.reduced_formula))

    def test_get_energy(self):
        test1 = (np.isclose(self.ir[0]._get_energy(0.5), -15, atol=1e-03))
        self.assertTrue(test1, '_get_energy: phase diagram gets error!')

        test2 = (
            np.isclose(self.ir[3]._get_energy(0.6666666),
                       -7.333333, atol=1e-03))
        self.assertTrue(test2,
                        '_get_energy: '
                        'grand canonical phase diagram gets error!')

        test3 = (
            np.isclose(self.ir[6]._get_energy(0.3333333),
                       -3.333333, atol=1e-03))
        self.assertTrue(test3,
                        '_get_energy: convex hull energy gets error. ')

        test4 = (
            np.isclose(self.ir[7]._get_energy(0.3333333),
                       -4, atol=1e-03))
        self.assertTrue(test4,
                        '_get_energy: gets error. ')

    def test_get_reaction(self):
        test1 = str(self.ir[0]._get_reaction(0.5)) == '0.5 O2 + 0.5 Mn -> ' \
                                                      '0.5 MnO2'
        self.assertTrue(test1,
                        '_get_reaction: '
                        'reaction not involving chempots species gets error!')

        test2 = str(self.ir[3]._get_reaction(0.666666)) \
            == '0.5 Mn + 0.5 Li2O -> Li + 0.25 MnO2 + 0.25 Mn' \
            or str(self.ir[3]._get_reaction(0.666666)) \
            == '0.5 Mn + 0.5  Li2O -> Li + 0.25 Mn + 0.25 MnO2'
        self.assertTrue(test2,
                        '_get_reaction: '
                        'reaction involving chempots species gets error!')

    def test_get_get_elmt_amt_in_rxt(self):
        rxt1 = Reaction(
            [Composition('Mn'), Composition('O2'), Composition('Li')],
            [Composition('LiMnO2')])
        test1 = np.isclose(self.ir[2]._get_elmt_amt_in_rxt(rxt1), 3)
        self.assertTrue(test1,
                        '_get_get_elmt_amt_in_rxt: '
                        'gpd elements amounts gets error!')

        rxt2 = rxt1
        rxt2.normalize_to(Composition('Li'), 0.5)
        test2 = np.isclose(self.ir[2]._get_elmt_amt_in_rxt(rxt2), 1.5)
        self.assertTrue(test2,
                        '_get_get_elmt_amt_in_rxt: '
                        'gpd elements amounts gets error!')

        rxt3 = Reaction([Composition('O2'), Composition('Li')],
                        [Composition('Li2O')])
        # Li is not counted
        test3 = np.isclose(self.ir[2]._get_elmt_amt_in_rxt(rxt3), 1)
        self.assertTrue(test3,
                        '_get_get_elmt_amt_in_rxt: '
                        'gpd elements amounts gets error!')

        # Li is counted
        test4 = np.isclose(self.ir[6]._get_elmt_amt_in_rxt(rxt3), 3)
        self.assertTrue(test4,
                        '_get_get_elmt_amt_in_rxt: '
                        'pd elements amounts gets error!')

    def test_convert(self):
        test_array = [(0.5, 1, 3), (0.4, 2, 3), (0, 1, 9), (1, 2, 7)]
        result = [InterfacialReactivity._convert(x, f1, f2)
                  for x, f1, f2 in test_array]
        answer = [0.75, 0.5, 0, 1]
        self.assertTrue(np.allclose(result, answer),
                        '_convert: conversion gets error! {0} expected,'
                        ' but gets {1}'.format(answer, result))

    def test_reverse_convert(self):
        test_array = [(0.5, 1, 3), (0.4, 2, 3), (0, 1, 9), (1, 2, 7)]
        result = [InterfacialReactivity._reverse_convert(x, f1, f2)
                  for x, f1, f2 in test_array]
        answer = [0.25, 0.3076923, 0, 1]
        self.assertTrue(np.allclose(result, answer),
                        '_convert: conversion gets error! {0} expected,'
                        ' but gets {1}'.format(answer, result))

    def test_get_products(self):
        test1 = sorted(self.ir[0].get_products()) == sorted(
            ['MnO2', 'O2', 'Mn'])
        self.assertTrue(test1,
                        'get_products: decomposition products gets error '
                        'for reaction not involving chempots species!')

        test2 = sorted(self.ir[3].get_products()) == sorted(
            ['Li', 'MnO2', 'Mn', 'Li2O'])
        self.assertTrue(test2,
                        'get_decomp: decomposition products gets error '
                        'for reaction involving chempots species!')

    def test_get_kinks(self):
        def test_get_kinks_helper(ir, index_expect,
                                  x_kink_expect, energy_kink_expect,
                                  react_kink_expect,
                                  energy_per_rxt_kink_expect):
            lst = list(ir.get_kinks())
            index = [i[0] for i in lst]
            x_kink = [i[1] for i in lst]
            energy_kink = [i[2] for i in lst]
            react_kink = [str(i[3]) for i in lst]
            energy_per_rxt_kink = [i[4] for i in lst]
            test1 = index == index_expect
            self.assertTrue(test1, 'get_kinks:index gets error!')

            test2 = np.allclose(x_kink, x_kink_expect)
            self.assertTrue(test2, 'get_kinks:x kinks gets error!')

            test3 = np.allclose(energy_kink, energy_kink_expect)
            self.assertTrue(test3, 'get_kinks:energy kinks gets error!')

            # Testing reaction strings are hard,
            # as species could be arranged in random order.
            test4 = len(react_kink) == len(react_kink_expect)
            self.assertTrue(test4,
                            'get_kinks: reaction kinks '
                            'gets error for {0} and {1} reaction!'.format(
                                ir.c1_original.reduced_formula,
                                ir.c2_original.reduced_formula))

            test5 = np.allclose(energy_per_rxt_kink,
                                energy_per_rxt_kink_expect)
            self.assertTrue(test5,
                            'get_kinks: energy_per_rxt_kinks gets error!')

        test_get_kinks_helper(self.ir[0], [1, 2, 3], [0, 0.5, 1],
                              [0, -15, 0],
                              ['Mn -> Mn', '0.5 O2 + 0.5 Mn -> 0.5 MnO2',
                               'O2 -> O2'],
                              [0,
                               -15 * InterfacialReactivity.EV_TO_KJ_PER_MOL,
                               0])
        test_get_kinks_helper(self.ir[10], [1, 2, 3], [0, 0.66667, 1],
                              [0, -10, 0],
                              ['Mn -> Mn', '0.5 O2 + 0.5 Mn -> 0.5 MnO2',
                               'O2 -> O2'],
                              [0,
                               -15 * InterfacialReactivity.EV_TO_KJ_PER_MOL,
                               0])
        test_get_kinks_helper(self.ir[11], [1, 2], [0, 1], [-3, -3],
                              ['Li2O2 + 2 Li -> 2 Li2O',
                               'Li2O2 + 2 Li -> 2 Li2O'],
                              [-6 * InterfacialReactivity.EV_TO_KJ_PER_MOL] *
                              2)
        test_get_kinks_helper(self.ir[12], [1, 2], [0, 1], [-0.5, -0.5],
                              ['Li2O2 -> Li2O + 0.5 O2',
                               'Li2O2 -> Li2O + 0.5 O2'],
                              [-2 * InterfacialReactivity.EV_TO_KJ_PER_MOL] *
                              2)

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
                relative_vectors_1 = [(x - x_kink[0], e - energy_kink[0])
                                      for x, e in points]
                relative_vectors_2 = [(x - x_kink[-1], e - energy_kink[-1])
                                      for x, e in points]
                relative_vectors = zip(relative_vectors_1, relative_vectors_2)
                positions = [np.cross(v1, v2) for v1, v2 in relative_vectors]
                test1 = np.all(np.array(positions) <= 0)

                hull = ConvexHull(points)
                test2 = len(hull.vertices) == len(points)
                self.assertTrue(test1 and test2,
                                'Error: Generating non-convex plot!')

        test_convexity_helper(self.ir[0])
        test_convexity_helper(self.ir[1])
        test_convexity_helper(self.ir[2])
        test_convexity_helper(self.ir[3])
        test_convexity_helper(self.ir[4])
        test_convexity_helper(self.ir[5])
        test_convexity_helper(self.ir[6])
        test_convexity_helper(self.ir[7])
        test_convexity_helper(self.ir[8])
        test_convexity_helper(self.ir[9])
        test_convexity_helper(self.ir[10])
        test_convexity_helper(self.ir[11])
        test_convexity_helper(self.ir[12])

    def test_get_original_composition_ratio(self):
        # expected reaction1: 0.5 O2 + 0.5 Mn -> 0.5 MnO2
        reaction1 = self.ir[0]._get_reaction(0.5)
        test1 = np.isclose(self.ir[0]._get_original_composition_ratio(
            reaction1), 0.5)
        self.assertTrue(test1,
                        '_get_original_composition_ratio: '
                        'reaction not involving chempots species gets error!')

        #  expected reaction2: 0.5 Mn + 0.5 Li2O -> Li + 0.25 MnO2 + 0.25 Mn
        reaction2 = self.ir[3]._get_reaction(0.666666)
        test2 = np.isclose(self.ir[3]._get_original_composition_ratio(
                reaction2), 0.5)
        self.assertTrue(test2,
                        '_get_original_composition_ratio: '
                        'reaction involving chempots species gets error!')

    def test_get_critical_original_kink_ratio(self):
        test1 = np.allclose(self.ir[0].get_critical_original_kink_ratio(),
                            [0, 0.5, 1])
        self.assertTrue(test1, 'get_critical_original_kink_ratio:'
                               ' gets error!')
        test2 = np.allclose(self.ir[10].get_critical_original_kink_ratio(),
                            [0, 0.5, 1])
        self.assertTrue(test2, 'get_critical_original_kink_ratio:'
                               ' gets error!')
        test3 = np.allclose(self.ir[11].get_critical_original_kink_ratio(),
                            [0, 1])
        self.assertTrue(test3, 'get_critical_original_kink_ratio:'
                               ' gets error!')
        test4 = np.allclose(self.ir[2].get_critical_original_kink_ratio(),
                            [0, 0.5, 1])
        self.assertTrue(test4, 'get_critical_original_kink_ratio:'
                               ' gets error!')
        test5 = np.allclose(self.ir[3].get_critical_original_kink_ratio(),
                            [0, 0.66666, 1])
        self.assertTrue(test5, 'get_critical_original_kink_ratio:'
                               ' gets error!')

    def test_labels(self):
        ir = self.ir[0]
        dict = ir.labels()
        test1 = dict == {1: 'x= 0.0 energy in eV/atom = 0.0 Mn -> Mn',
                         2: 'x= 0.5 energy in eV/atom = -15.0 0.5 O2 + 0.5 '
                            'Mn -> 0.5 MnO2',
                         3: 'x= 1.0 energy in eV/atom = 0.0 O2 -> O2'}
        self.assertTrue(test1,
                        'labels:label does not match for interfacial system '
                        'with {0} and {1}.'.format(
                            ir.c1_original.reduced_formula,
                            ir.c2_original.reduced_formula))

    def test_plot(self):
        # Test plot is hard. Here just to call the plot function to see if any
        #  error occurs.
        for i in self.ir:
            i.plot()

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
            self.assertTrue(np.allclose(i.minimum(), j),
                            'minimum: the system with {0} and {1} '
                            'gets error!{2} expected, but gets {3}'.format(
                                i.c1_original.reduced_formula,
                                i.c2_original.reduced_formula, str(j),
                                str(i.minimum())))

    def test_get_no_mixing_energy(self):
        with self.assertRaises(Exception) as context1:
            self.ir[0].get_no_mixing_energy()
        self.assertTrue(
            'Please provide grand potential phase diagram'
            ' for computing no_mixing_energy!' == str(context1.exception))
        answer = [
            [(u'MnO2 (eV/f.u.)', 0.0), (u'Mn (eV/f.u.)', 0.0)],
            [(u'Mn (eV/atom)', 0.0), (u'O2 (eV/atom)', -4.0)],
            [(u'Li2O (eV/f.u.)', 0.0), (u'Mn (eV/f.u.)', 0.0)],
            [(u'Mn (eV/atom)', 0.0), (u'O2 (eV/atom)', -4.0)],
            [(u'Mn (eV/atom)', 0.0), (u'Li2O (eV/atom)', 0.0)]
        ]

        def name_lst(lst):
            return (lst[0][0], lst[1][0])

        def energy_lst(lst):
            return (lst[0][1], lst[1][1])

        result_info = [i.get_no_mixing_energy() for i in self.ir if i.grand]
        for i, j in zip(result_info, answer):
            self.assertTrue(name_lst(i) == name_lst(j),
                            'get_no_mixing_energy: names get error,'
                            ' {0} expected but gets {1}'.format(
                                name_lst(j), name_lst(i)))
            self.assertTrue(np.allclose(energy_lst(i), energy_lst(j)),
                            'get_no_mixing_energy: '
                            'no_mixing energies get error, '
                            '{0} expected but gets {1}'.format(
                                energy_lst(j), energy_lst(i)))

    def test_get_chempot_correction(self):
        # test data from fig. 6 in ref:
        # Prediction of A2BX4 metal-chalcogenide compounds via
        # first-principles thermodynamics, PHYSICAL REVIEW B 86, 014109 (2012)

        # test pressure effect.
        actual = InterfacialReactivity.get_chempot_correction("O", 298.15,
                                                              100E5)
        expect = 0.05916
        self.assertTrue(np.isclose(actual, expect, atol=1E-2),
                        "get_chempot_correction gets "
                        "error, {0} expected but gets {1}".format(expect,
                                                                  actual))
        # test temperature effect.
        actual_2 = InterfacialReactivity.get_chempot_correction("O", 1000,
                                                                1E5)
        expect_2 = -0.82352
        self.assertTrue(np.isclose(actual_2, expect_2, atol=1E-2),
                        "get_chempot_correction gets "
                        "error, {0} expected but gets {1}".format(expect_2,
                                                                  actual_2))

        actual_3 = InterfacialReactivity.get_chempot_correction("O", 500,
                                                                1E5)
        expect_3 = -0.223
        self.assertTrue(np.isclose(actual_3, expect_3, atol=1E-2),
                        "get_chempot_correction gets "
                        "error, {0} expected but gets {1}".format(expect_3,
                                                                  actual_3))
        # test mixed effect.
        actual_4 = InterfacialReactivity.get_chempot_correction("O", 1000,
                                                                1E-25)
        expect_4 = -3.800
        self.assertTrue(np.isclose(actual_4, expect_4, atol=1E-2),
                        "get_chempot_correction gets "
                        "error, {0} expected but gets {1}".format(expect_4,
                                                                  actual_4))
        actual_5 = InterfacialReactivity.get_chempot_correction("O", 1250,
                                                                1E-25)
        expect_5 = -4.86
        self.assertTrue(np.isclose(actual_5, expect_5, atol=1E-2),
                        "get_chempot_correction gets "
                        "error, {0} expected but gets {1}".format(expect_5,
                                                                  actual_5))
        actual_6 = InterfacialReactivity.get_chempot_correction("O", 1500,
                                                                1E-25)
        expect_6 = -5.928
        self.assertTrue(np.isclose(actual_6, expect_6, atol=1E-2),
                        "get_chempot_correction gets "
                        "error, {0} expected but gets {1}".format(expect_6,
                                                                  actual_6))
        actual_7 = InterfacialReactivity.get_chempot_correction("O", 1000,
                                                                1E-15)
        expect_7 = -2.808
        self.assertTrue(np.isclose(actual_7, expect_7, atol=1E-2),
                        "get_chempot_correction gets "
                        "error, {0} expected but gets {1}".format(expect_7,
                                                                  actual_7))
        # test non-gas phase.
        actual_8 = InterfacialReactivity.get_chempot_correction("Li", 1000,
                                                                1E15)
        expect_8 = 0
        self.assertTrue(np.isclose(actual_8, expect_8, atol=1E-5),
                        "get_chempot_correction gets "
                        "error, {0} expected but gets {1}".format(expect_8,
                                                                  actual_8))


if __name__ == '__main__':
    unittest.main()
