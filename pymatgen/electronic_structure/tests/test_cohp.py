# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals
import unittest
import json
import os

from pymatgen.electronic_structure.cohp import CompleteCohp, Cohp, IcohpValue, IcohpCollection
from pymatgen.electronic_structure.core import Spin, Orbital
from pymatgen.util.testing import PymatgenTest

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "cohp")


class CohpTest(unittest.TestCase):
    def setUp(self):
        with open(os.path.join(test_dir, "cohp.json"), "r") as f:
            self.cohp = Cohp.from_dict(json.load(f))
        self.cohp_only = Cohp(self.cohp.efermi,
                              self.cohp.energies,
                              self.cohp.cohp)
        with open(os.path.join(test_dir, "coop.json"), "r") as f:
            self.coop = Cohp.from_dict(json.load(f))

    def test_as_from_dict(self):
        with open(os.path.join(test_dir, "cohp.json"), "r") as f:
            cohp_dict = json.load(f)
        self.assertEqual(self.cohp.as_dict(), cohp_dict)

    def test_attributes(self):
        self.assertEqual(len(self.cohp.energies), 301)
        self.assertEqual(self.cohp.efermi, 9.75576)
        self.assertEqual(self.coop.efermi, 5.90043)
        self.assertFalse(self.cohp.are_coops)
        self.assertTrue(self.coop.are_coops)

    def test_get_icohp(self):
        self.assertEqual(self.cohp.get_icohp(),
                         self.cohp.get_cohp(integrated=True))
        self.assertEqual(None, self.cohp_only.get_icohp())

    def test_get_interpolated_value(self):
        # icohp_ef are the ICHOP(Ef) values taken from
        # the ICOHPLIST.lobster file.
        icohp_ef_dict = {Spin.up: -0.10218, Spin.down: -0.19701}
        icoop_ef_dict = {Spin.up: 0.24714}
        icohp_ef = self.cohp.get_interpolated_value(self.cohp.efermi,
                                                    integrated=True)
        icoop_ef = self.coop.get_interpolated_value(self.coop.efermi,
                                                    integrated=True)
        self.assertAlmostEqual(icohp_ef_dict, icohp_ef)
        self.assertAlmostEqual(icoop_ef_dict, icoop_ef)
        with self.assertRaises(ValueError):
            self.cohp_only.get_interpolated_value(5.0, integrated=True)

    def test_str(self):
        with open(os.path.join(test_dir, "cohp.str"), "rt") as f:
            str_cohp = f.read()
        with open(os.path.join(test_dir, "coop.str"), "rt") as f:
            str_coop = f.read()
        self.assertEqual(self.cohp.__str__(), str_cohp)
        self.assertEqual(self.coop.__str__(), str_coop)


class IcohpValueTest(unittest.TestCase):
    def setUp(self):
        # without spin polarization
        label = "1"
        atom1 = "K1"
        atom2 = "F2"
        length = "2.3"
        translation = [-1, 0, 0]
        num = 1
        icohp = {Spin.up: -2.0}
        are_coops = False
        self.icohpvalue = IcohpValue(label=label, atom1=atom1, atom2=atom2, length=length, translation=translation,
                                     num=num, icohp=icohp, are_coops=are_coops)

        label_sp = "1"
        atom1_sp = "K1"
        atom2_sp = "F2"
        length_sp = "2.3"
        translation_sp = [-1, 0, 0]
        num_sp = 1
        icohp_sp = {Spin.up: -1.1, Spin.down: -1.0}
        are_coops_sp = False
        self.icohpvalue_sp = IcohpValue(label=label_sp, atom1=atom1_sp, atom2=atom2_sp, length=length_sp,
                                        translation=translation_sp, num=num_sp, icohp=icohp_sp, are_coops=are_coops_sp)

    def test_attributes(self):
        # without spin polarization
        self.assertEqual(self.icohpvalue_sp.num_bonds, 1)
        self.assertEqual(self.icohpvalue_sp.are_coops, False)
        self.assertEqual(self.icohpvalue_sp.is_spin_polarized, True)
        self.assertDictEqual(self.icohpvalue.icohp, {Spin.up: -2.0})

        # with spin polarization
        self.assertEqual(self.icohpvalue_sp.num_bonds, 1)
        self.assertEqual(self.icohpvalue_sp.are_coops, False)
        self.assertEqual(self.icohpvalue_sp.is_spin_polarized, True)
        self.assertDictEqual(self.icohpvalue_sp.icohp, {Spin.up: -1.1, Spin.down: -1.0})

    def test_icohpvalue(self):
        # without spin polarization
        self.assertEqual(self.icohpvalue.icohpvalue(spin=Spin.up), -2.0)

        # with spin polarization
        self.assertEqual(self.icohpvalue_sp.icohpvalue(spin=Spin.up), -1.1)
        self.assertEqual(self.icohpvalue_sp.icohpvalue(spin=Spin.down), -1.0)

    def test_summed_icohp(self):
        # without spin polarization
        self.assertEqual(self.icohpvalue.summed_icohp, -2.0)

        # with spin polarization
        self.assertEqual(self.icohpvalue_sp.summed_icohp, -2.1)


class CombinedIcohpTest(unittest.TestCase):
    def setUp(self):
        # without spin polarization:
        are_coops = False
        is_spin_polarized = False
        list_atom2 = ['K2', 'K2', 'K2', 'K2', 'K2', 'K2']
        list_icohp = [{Spin.up: -0.40075}, {Spin.up: -0.40074}, {Spin.up: -0.40079}, {Spin.up: -0.40079},
                      {Spin.up: -0.40074}, {Spin.up: -0.40075}]
        list_icoop = [{Spin.up: 0.02342}, {Spin.up: 0.02342}, {Spin.up: 0.02343}, {Spin.up: 0.02343},
                      {Spin.up: 0.02342}, {Spin.up: 0.02342}]
        list_labels = ['1', '2', '3', '4', '5', '6']
        list_length = [2.71199, 2.71199, 2.71199, 2.71199, 2.71199, 2.71199]
        list_num = [1, 1, 1, 1, 1, 1]
        list_atom1 = ['F1', 'F1', 'F1', 'F1', 'F1', 'F1']
        list_translation = [[0, -1, -1], [-1, 0, -1], [0, 0, -1], [-1, -1, 0], [0, -1, 0], [-1, 0, 0]]
        self.icohpcollection_KF = IcohpCollection(is_spin_polarized=is_spin_polarized, are_coops=are_coops,
                                                  list_labels=list_labels, list_atom1=list_atom1, list_atom2=list_atom2,
                                                  list_length=list_length, list_translation=list_translation,
                                                  list_num=list_num, list_icohp=list_icohp)

        self.icoopcollection_KF = IcohpCollection(is_spin_polarized=is_spin_polarized, are_coops=True,
                                                  list_labels=list_labels, list_atom1=list_atom1, list_atom2=list_atom2,
                                                  list_length=list_length, list_translation=list_translation,
                                                  list_num=list_num, list_icohp=list_icoop)

        # with spin polarization:
        list_atom2_sp = ['Fe7', 'Fe9']
        list_labels_sp = ['1', '2']
        list_translation_sp = [[0, 0, 0], [0, 0, 0]]
        list_length_sp = [2.83189, 2.45249]
        list_atom1_sp = ['Fe8', 'Fe8']
        is_spin_polarized_sp = True
        are_coops_sp = False
        list_num_sp = [2, 1]
        list_icohp_sp = [{Spin.up: -0.10218, Spin.down: -0.19701}, {Spin.up: -0.28485, Spin.down: -0.58279}]
        list_icoop_sp = [{Spin.up: -0.11389, Spin.down: -0.20828}, {Spin.up: -0.04087, Spin.down: -0.05756}]

        self.icohpcollection_Fe = IcohpCollection(is_spin_polarized=is_spin_polarized_sp, are_coops=are_coops_sp,
                                                  list_labels=list_labels_sp, list_atom1=list_atom1_sp,
                                                  list_atom2=list_atom2_sp, list_length=list_length_sp,
                                                  list_translation=list_translation_sp, list_num=list_num_sp,
                                                  list_icohp=list_icohp_sp)
        self.icoopcollection_Fe = IcohpCollection(is_spin_polarized=is_spin_polarized_sp, are_coops=True,
                                                  list_labels=list_labels_sp, list_atom1=list_atom1_sp,
                                                  list_atom2=list_atom2_sp, list_length=list_length_sp,
                                                  list_translation=list_translation_sp, list_num=list_num_sp,
                                                  list_icohp=list_icoop_sp)

    def test_get_icohp_by_label(self):
        # without spin polarization

        # ICOHPs
        self.assertEqual(self.icohpcollection_KF.get_icohp_by_label("1"), -0.40075)
        self.assertEqual(self.icohpcollection_KF.get_icohp_by_label("2"), -0.40074)
        self.assertEqual(self.icohpcollection_KF.get_icohp_by_label("3"), -0.40079)
        self.assertEqual(self.icohpcollection_KF.get_icohp_by_label("4"), -0.40079)
        self.assertEqual(self.icohpcollection_KF.get_icohp_by_label("5"), -0.40074)
        self.assertEqual(self.icohpcollection_KF.get_icohp_by_label("6"), -0.40075)

        # with spin polarization
        # summed spin
        # ICOHPs
        self.assertEqual(self.icohpcollection_Fe.get_icohp_by_label("1"), -0.10218 - 0.19701)
        self.assertEqual(self.icohpcollection_Fe.get_icohp_by_label("2"), -0.28485 - 0.58279)

        # Spin up
        # ICOHPs
        self.assertEqual(self.icohpcollection_Fe.get_icohp_by_label("1", summed_spin_channels=False), -0.10218)
        self.assertEqual(self.icohpcollection_Fe.get_icohp_by_label("2", summed_spin_channels=False), -0.28485)

        # Spin down
        # ICOHPs
        self.assertEqual(self.icohpcollection_Fe.get_icohp_by_label("1", summed_spin_channels=False, spin=Spin.down),
                         -0.19701)
        self.assertEqual(self.icohpcollection_Fe.get_icohp_by_label("2", summed_spin_channels=False, spin=Spin.down),
                         -0.58279)

    def test_get_summed_icohp_by_labellist(self):
        # without spin polarization
        self.assertAlmostEqual(
            self.icohpcollection_KF.get_summed_icohp_by_labellist(["1", "2", "3", "4", "5", "6"], divisor=6.0),
            -0.40076)

        # with spin polarization
        sum1 = (-0.10218 - 0.19701 - 0.28485 - 0.58279) / 2.0
        sum2 = (-0.10218 - 0.28485) / 2.0
        sum3 = (-0.19701 - 0.58279) / 2.0
        self.assertAlmostEqual(self.icohpcollection_Fe.get_summed_icohp_by_labellist(["1", "2"], divisor=2.0), sum1)
        self.assertAlmostEqual(
            self.icohpcollection_Fe.get_summed_icohp_by_labellist(["1", "2"], summed_spin_channels=False, divisor=2.0),
            sum2)
        self.assertAlmostEqual(
            self.icohpcollection_Fe.get_summed_icohp_by_labellist(["1", "2"], summed_spin_channels=False,
                                                                  spin=Spin.down, divisor=2.0), sum3)

    def test_get_icohp_dict_by_bondlengths(self):
        # without spin polarization
        icohpvalue = {}
        icohpvalue["1"] = {'@module': 'pymatgen.electronic_structure.cohp', 'num': 1, 'length': 2.71199,
                           'icohp': {Spin.up: -0.40075},
                           'are_coops': False, 'label': '1', 'atom2': 'K2', '@class': 'IcohpValue', 'atom1': 'F1',
                           'translation': [0, -1, -1]}
        icohpvalue["2"] = {'@module': 'pymatgen.electronic_structure.cohp', 'num': 1, 'length': 2.71199,
                           'icohp': {Spin.up: -0.40074},
                           'are_coops': False, 'label': '2', 'atom2': 'K2', '@class': 'IcohpValue', 'atom1': 'F1',
                           'translation': [-1, 0, -1]}
        icohpvalue["3"] = {'@module': 'pymatgen.electronic_structure.cohp', 'num': 1, 'length': 2.71199,
                           'icohp': {Spin.up: -0.40079},
                           'are_coops': False, 'label': '3', 'atom2': 'K2', '@class': 'IcohpValue', 'atom1': 'F1',
                           'translation': [0, 0, -1]}
        icohpvalue["4"] = {'@module': 'pymatgen.electronic_structure.cohp', 'num': 1, 'length': 2.71199,
                           'icohp': {Spin.up: -0.40079},
                           'are_coops': False, 'label': '4', 'atom2': 'K2', '@class': 'IcohpValue', 'atom1': 'F1',
                           'translation': [-1, -1, 0]}
        icohpvalue["5"] = {'@module': 'pymatgen.electronic_structure.cohp', 'num': 1, 'length': 2.71199,
                           'icohp': {Spin.up: -0.40074},
                           'are_coops': False, 'label': '5', 'atom2': 'K2', '@class': 'IcohpValue', 'atom1': 'F1',
                           'translation': [0, -1, 0]}
        icohpvalue["6"] = {'@module': 'pymatgen.electronic_structure.cohp', 'num': 1, 'length': 2.71199,
                           'icohp': {Spin.up: -0.40075},
                           'are_coops': False, 'label': '6', 'atom2': 'K2', '@class': 'IcohpValue', 'atom1': 'F1',
                           'translation': [-1, 0, 0]}

        dict_KF = self.icohpcollection_KF.get_icohp_dict_by_bondlengths(minbondlength=0.0, maxbondlength=8.0)
        for key, value in sorted(dict_KF.items()):
            self.assertDictEqual(value.as_dict(), icohpvalue[key])

        self.assertDictEqual({}, self.icohpcollection_KF.get_icohp_dict_by_bondlengths(minbondlength=0.0,
                                                                                       maxbondlength=1.0))

        # with spin polarization
        icohpvalue_spin = {}
        icohpvalue_spin["1"] = {'num': 2, 'atom2': 'Fe7', 'translation': [0, 0, 0],
                                '@module': 'pymatgen.electronic_structure.cohp',
                                'are_coops': False, 'atom1': 'Fe8',
                                'label': '1', 'length': 2.83189, '@class': 'IcohpValue',
                                'icohp': {Spin.up: -0.10218, Spin.down: -0.19701}}
        icohpvalue_spin["2"] = {'num': 1, 'atom2': 'Fe9', 'translation': [0, 0, 0],
                                '@module': 'pymatgen.electronic_structure.cohp',
                                'are_coops': False, 'atom1': 'Fe8',
                                'label': '2', 'length': 2.45249, '@class': 'IcohpValue',
                                'icohp': {Spin.up: -0.28485, Spin.down: -0.58279}}

        dict_Fe = self.icohpcollection_Fe.get_icohp_dict_by_bondlengths(minbondlength=0.0, maxbondlength=8.0)
        for key, value in sorted(dict_Fe.items()):
            self.assertDictEqual(value.as_dict(), icohpvalue_spin[key])

        dict_Fe2 = self.icohpcollection_Fe.get_icohp_dict_by_bondlengths(minbondlength=2.5, maxbondlength=2.9)
        self.assertEqual(len(dict_Fe2), 1)
        for key, value in sorted(dict_Fe2.items()):
            self.assertDictEqual(value.as_dict(), icohpvalue_spin[key])

    def test_get_icohp_dict_of_certain_site(self):
        # without spin polarization
        icohpvalue = {}
        icohpvalue["1"] = {'translation': [0, -1, -1], 'are_coops': False,
                           '@module': 'pymatgen.electronic_structure.cohp', 'length': 2.71199,
                           'atom2': 'K2',
                           '@class': 'IcohpValue', 'atom1': 'F1', 'num': 1, 'label': '1', 'icohp': {Spin.up: -0.40075}}
        icohpvalue["2"] = {'translation': [-1, 0, -1], 'are_coops': False,
                           '@module': 'pymatgen.electronic_structure.cohp', 'length': 2.71199,
                           'atom2': 'K2',
                           '@class': 'IcohpValue', 'atom1': 'F1', 'num': 1, 'label': '2', 'icohp': {Spin.up: -0.40074}}
        icohpvalue["3"] = {'translation': [0, 0, -1], 'are_coops': False,
                           '@module': 'pymatgen.electronic_structure.cohp', 'length': 2.71199,
                           'atom2': 'K2',
                           '@class': 'IcohpValue', 'atom1': 'F1', 'num': 1, 'label': '3', 'icohp': {Spin.up: -0.40079}}
        icohpvalue["4"] = {'translation': [-1, -1, 0], 'are_coops': False,
                           '@module': 'pymatgen.electronic_structure.cohp', 'length': 2.71199,
                           'atom2': 'K2',
                           '@class': 'IcohpValue', 'atom1': 'F1', 'num': 1, 'label': '4', 'icohp': {Spin.up: -0.40079}}
        icohpvalue["5"] = {'translation': [0, -1, 0], 'are_coops': False,
                           '@module': 'pymatgen.electronic_structure.cohp', 'length': 2.71199,
                           'atom2': 'K2',
                           '@class': 'IcohpValue', 'atom1': 'F1', 'num': 1, 'label': '5', 'icohp': {Spin.up: -0.40074}}
        icohpvalue["6"] = {'translation': [-1, 0, 0], 'are_coops': False,
                           '@module': 'pymatgen.electronic_structure.cohp', 'length': 2.71199,
                           'atom2': 'K2',
                           '@class': 'IcohpValue', 'atom1': 'F1', 'num': 1, 'label': '6', 'icohp': {Spin.up: -0.40075}}

        dict_KF = self.icohpcollection_KF.get_icohp_dict_of_certain_site(site=0)

        for key, value in sorted(dict_KF.items()):
            self.assertDictEqual(value.as_dict(), icohpvalue[key])

        # compare number of results dependent on minsummedicohp, maxsummedicohp,minbondlength, maxbondlength
        dict_KF_2 = self.icohpcollection_KF.get_icohp_dict_of_certain_site(site=0, minsummedicohp=None,
                                                                           maxsummedicohp=-0.0, minbondlength=0.0,
                                                                           maxbondlength=8.0)
        dict_KF_3 = self.icohpcollection_KF.get_icohp_dict_of_certain_site(site=0, minsummedicohp=None,
                                                                           maxsummedicohp=-0.5, minbondlength=0.0,
                                                                           maxbondlength=8.0)
        dict_KF_4 = self.icohpcollection_KF.get_icohp_dict_of_certain_site(site=0, minsummedicohp=0.0,
                                                                           maxsummedicohp=None, minbondlength=0.0,
                                                                           maxbondlength=8.0)
        dict_KF_5 = self.icohpcollection_KF.get_icohp_dict_of_certain_site(site=0, minsummedicohp=None,
                                                                           maxsummedicohp=None, minbondlength=0.0,
                                                                           maxbondlength=2.0)
        dict_KF_6 = self.icohpcollection_KF.get_icohp_dict_of_certain_site(site=0, minsummedicohp=None,
                                                                           maxsummedicohp=None, minbondlength=3.0,
                                                                           maxbondlength=8.0)

        self.assertEqual(len(dict_KF_2), 6)
        self.assertEqual(len(dict_KF_3), 0)
        self.assertEqual(len(dict_KF_4), 0)
        self.assertEqual(len(dict_KF_5), 0)
        self.assertEqual(len(dict_KF_6), 0)

        # spin polarization

        dict_Fe = self.icohpcollection_Fe.get_icohp_dict_of_certain_site(site=0)
        self.assertEqual(len(dict_Fe), 0)

        # Fe8
        dict_Fe2 = self.icohpcollection_Fe.get_icohp_dict_of_certain_site(site=7)
        self.assertEqual(len(dict_Fe2), 2)
        # Test the values

        icohplist_Fe = {}
        icohplist_Fe["1"] = {'are_coops': False, 'translation': [0, 0, 0],
                             'icohp': {Spin.down: -0.19701, Spin.up: -0.10218}, 'length': 2.83189,
                             '@module': 'pymatgen.electronic_structure.cohp', 'atom1': 'Fe8', 'atom2': 'Fe7',
                             'label': '1',
                             '@class': 'IcohpValue', 'num': 2}
        icohplist_Fe["2"] = {'are_coops': False, 'translation': [0, 0, 0],
                             'icohp': {Spin.down: -0.58279, Spin.up: -0.28485}, 'length': 2.45249,
                             '@module': 'pymatgen.electronic_structure.cohp', 'atom1': 'Fe8', 'atom2': 'Fe9',
                             'label': '2',
                             '@class': 'IcohpValue', 'num': 1}

        for key, value in sorted(dict_Fe2.items()):
            self.assertEqual(value.as_dict(), icohplist_Fe[key])

        # Fe9
        dict_Fe3 = self.icohpcollection_Fe.get_icohp_dict_of_certain_site(site=8)
        self.assertEqual(len(dict_Fe3), 1)

        # compare number of results dependent on minsummedicohp, maxsummedicohp,minbondlength, maxbondlength
        # Fe8
        dict_Fe4 = self.icohpcollection_Fe.get_icohp_dict_of_certain_site(site=7, minsummedicohp=-0.3,
                                                                          maxsummedicohp=None, minbondlength=0.0,
                                                                          maxbondlength=8.0)
        self.assertEqual(len(dict_Fe4), 1)
        values = []
        for key, value in dict_Fe4.items():
            values.append(value)
        self.assertDictEqual(values[0].as_dict(), icohplist_Fe["1"])

        dict_Fe5 = self.icohpcollection_Fe.get_icohp_dict_of_certain_site(site=7, minsummedicohp=None,
                                                                          maxsummedicohp=-0.3, minbondlength=0.0,
                                                                          maxbondlength=8.0)
        self.assertEqual(len(dict_Fe5), 1)
        values = []
        for key, value in dict_Fe5.items():
            values.append(value)
        self.assertDictEqual(values[0].as_dict(), icohplist_Fe["2"])

        dict_Fe6 = self.icohpcollection_Fe.get_icohp_dict_of_certain_site(site=7, minsummedicohp=None,
                                                                          maxsummedicohp=None, minbondlength=0.0,
                                                                          maxbondlength=2.5)

        self.assertEqual(len(dict_Fe6), 1)
        values = []
        for key, value in dict_Fe6.items():
            values.append(value)
        self.assertDictEqual(values[0].as_dict(), icohplist_Fe["2"])

        dict_Fe7 = self.icohpcollection_Fe.get_icohp_dict_of_certain_site(site=7, minsummedicohp=None,
                                                                          maxsummedicohp=None, minbondlength=2.5,
                                                                          maxbondlength=8.0)
        self.assertEqual(len(dict_Fe7), 1)
        values = []
        for key, value in dict_Fe7.items():
            values.append(value)
        self.assertDictEqual(values[0].as_dict(), icohplist_Fe["1"])

    def test_extremum_icohpvalue(self):
        # without spin polarization
        # ICOHPs
        self.assertEqual(self.icohpcollection_KF.extremum_icohpvalue(), -0.40079)
        # ICOOPs
        self.assertEqual(self.icoopcollection_KF.extremum_icohpvalue(), 0.02343)
        # with spin polarization
        # summed spin
        # ICOHPs
        self.assertEqual(self.icohpcollection_Fe.extremum_icohpvalue(), -0.86764)
        self.assertAlmostEqual(self.icoopcollection_Fe.extremum_icohpvalue(), -0.09842999999999999)
        # ICOOPs
        # spin up
        # ICOHPs
        self.assertEqual(self.icohpcollection_Fe.extremum_icohpvalue(summed_spin_channels=False), -0.28485)
        # ICOOPs
        self.assertEqual(self.icoopcollection_Fe.extremum_icohpvalue(summed_spin_channels=False), -0.04087)
        # spin down
        # ICOHPs
        self.assertEqual(self.icohpcollection_Fe.extremum_icohpvalue(summed_spin_channels=False, spin=Spin.down),
                         -0.58279)
        # ICOOPs
        self.assertEqual(self.icoopcollection_Fe.extremum_icohpvalue(summed_spin_channels=False, spin=Spin.down),
                         -0.05756)


class CompleteCohpTest(PymatgenTest):
    def setUp(self):
        filepath = os.path.join(test_dir, "complete_cohp_lobster.json")
        with open(filepath, "r") as f:
            self.cohp_lobster_dict = CompleteCohp.from_dict(json.load(f))
        filepath = os.path.join(test_dir, "complete_coop_lobster.json")
        with open(filepath, "r") as f:
            self.coop_lobster_dict = CompleteCohp.from_dict(json.load(f))
        filepath = os.path.join(test_dir, "complete_cohp_lmto.json")
        with open(filepath, "r") as f:
            self.cohp_lmto_dict = CompleteCohp.from_dict(json.load(f))
        filepath = os.path.join(test_dir, "complete_cohp_orbitalwise.json")
        with open(filepath, "r") as f:
            self.cohp_orb_dict = CompleteCohp.from_dict(json.load(f))
        # Lobster 3.0
        filepath = os.path.join(test_dir, "complete_cohp_forb.json")
        with open(filepath, "r") as f:
            self.cohp_lobster_forb_dict = CompleteCohp.from_dict(json.load(f))

            # Lobster 2.0
        filepath = os.path.join(test_dir, "COPL.BiSe")
        structure = os.path.join(test_dir, "CTRL.BiSe")
        self.cohp_lmto = CompleteCohp.from_file("lmto", filename=filepath,
                                                structure_file=structure)
        filepath = os.path.join(test_dir, "COHPCAR.lobster")
        structure = os.path.join(test_dir, "POSCAR")
        self.cohp_lobster = CompleteCohp.from_file("lobster",
                                                   filename=filepath,
                                                   structure_file=structure)
        filepath = os.path.join(test_dir, "COOPCAR.lobster.BiSe")
        structure = os.path.join(test_dir, "POSCAR.BiSe")
        self.coop_lobster = CompleteCohp.from_file("lobster",
                                                   filename=filepath,
                                                   structure_file=structure,
                                                   are_coops=True)
        filepath = os.path.join(test_dir, "COHPCAR.lobster.orbitalwise")
        structure = os.path.join(test_dir, "POSCAR.orbitalwise")
        self.cohp_orb = CompleteCohp.from_file("lobster",
                                               filename=filepath,
                                               structure_file=structure)
        filepath = os.path.join(test_dir, "COHPCAR.lobster.notot.orbitalwise")
        self.cohp_notot = CompleteCohp.from_file("lobster",
                                                 filename=filepath,
                                                 structure_file=structure)
        # Lobster 3.0
        filepath = os.path.join(test_dir, "COHPCAR.lobster.Na2UO4")
        structure = os.path.join(test_dir, "POSCAR.Na2UO4")
        self.cohp_lobster_forb = CompleteCohp.from_file("lobster", filename=filepath, structure_file=structure)

    def test_attiributes(self):
        self.assertFalse(self.cohp_lobster.are_coops)
        self.assertFalse(self.cohp_lobster_dict.are_coops)
        self.assertFalse(self.cohp_lmto.are_coops)
        self.assertFalse(self.cohp_lmto_dict.are_coops)
        self.assertTrue(self.coop_lobster.are_coops)
        self.assertTrue(self.coop_lobster_dict.are_coops)
        self.assertFalse(self.cohp_lobster_forb.are_coops)
        self.assertFalse(self.cohp_lobster_forb_dict.are_coops)

        self.assertEqual(len(self.cohp_lobster.energies), 301)
        self.assertEqual(len(self.cohp_lmto.energies), 801)
        self.assertEqual(len(self.coop_lobster.energies), 241)
        self.assertEqual(len(self.cohp_lobster_forb.energies), 7)

        self.assertEqual(self.cohp_lobster.efermi, 9.75576)
        self.assertEqual(self.cohp_lmto.efermi, -2.3433)
        self.assertEqual(self.coop_lobster.efermi, 5.90043)
        self.assertEqual(self.cohp_lobster_forb.efermi, 4.12875)

    def test_dict(self):
        # The json files are dict representations of the COHPs from the LMTO
        # and LOBSTER calculations and should thus be the same.

        self.assertEqual(self.cohp_lobster.as_dict(),
                         self.cohp_lobster_dict.as_dict())
        self.assertEqual(self.cohp_orb.as_dict(),
                         self.cohp_orb_dict.as_dict())
        # Lobster 3.0, including f orbitals
        self.assertEqual(self.cohp_lobster_forb.as_dict(),
                         self.cohp_lobster_forb_dict.as_dict())

        # Testing the LMTO dicts will be more involved. Since the average
        # is calculated and not read, there may be differences in rounding
        # with a very small number of matrix elements, which would cause the
        # test to fail
        for key in ["COHP", "ICOHP"]:
            self.assertArrayAlmostEqual(
                self.cohp_lmto.as_dict()[key]["average"]["1"],
                self.cohp_lmto_dict.as_dict()[key]["average"]["1"], 5)
        for key in self.cohp_lmto.as_dict():
            if key not in ["COHP", "ICOHP"]:
                self.assertEqual(self.cohp_lmto.as_dict()[key],
                                 self.cohp_lmto_dict.as_dict()[key])
            else:
                for bond in self.cohp_lmto.as_dict()[key]:
                    if bond != "average":
                        self.assertEqual(self.cohp_lmto.as_dict()[key][bond],
                                         self.cohp_lmto_dict.as_dict()[key][bond])

    def test_icohp_values(self):
        # icohp_ef are the ICHOP(Ef) values taken from
        # the ICOHPLIST.lobster file.
        icohp_ef_dict = {"1": {Spin.up: -0.10218, Spin.down: -0.19701},
                         "2": {Spin.up: -0.28485, Spin.down: -0.58279}}
        all_cohps_lobster = self.cohp_lobster.all_cohps
        for bond in icohp_ef_dict:
            icohp_ef = all_cohps_lobster[bond].get_interpolated_value(
                self.cohp_lobster.efermi, integrated=True)
            self.assertEqual(icohp_ef_dict[bond], icohp_ef)

        icoop_ef_dict = {"1": {Spin.up: 0.14245},
                         "2": {Spin.up: -0.04118},
                         "3": {Spin.up: 0.14245},
                         "4": {Spin.up: -0.04118},
                         "5": {Spin.up: -0.03516},
                         "6": {Spin.up: 0.10745},
                         "7": {Spin.up: -0.03516},
                         "8": {Spin.up: 0.10745},
                         "9": {Spin.up: -0.12395},
                         "10": {Spin.up: 0.24714},
                         "11": {Spin.up: -0.12395}}
        all_coops_lobster = self.coop_lobster.all_cohps
        for bond in icoop_ef_dict:
            icoop_ef = all_coops_lobster[bond].get_interpolated_value(
                self.coop_lobster.efermi, integrated=True)
            self.assertEqual(icoop_ef_dict[bond], icoop_ef)

    def test_orbital_resolved_cohp(self):
        # When read from a COHPCAR file, total COHPs are calculated from
        # the orbital-resolved COHPs if the total is missing. This may be
        # case for LOBSTER version 2.2.0 and earlier due to a bug with the
        # cohpgenerator keyword. The calculated total should be approximately
        # the total COHP calculated by LOBSTER. Due to numerical errors in
        # the LOBSTER calculation, the precision is not very high though.

        self.assertArrayAlmostEqual(
            self.cohp_orb.all_cohps["1"].cohp[Spin.up],
            self.cohp_notot.all_cohps["1"].cohp[Spin.up], decimal=3)
        self.assertArrayAlmostEqual(
            self.cohp_orb.all_cohps["1"].icohp[Spin.up],
            self.cohp_notot.all_cohps["1"].icohp[Spin.up], decimal=3)

        # Tests different methods for getting orbital-resolved COHPs
        ref = self.cohp_orb.orb_res_cohp["1"]["4s-4px"]
        cohp_label = self.cohp_orb.get_orbital_resolved_cohp("1",
                                                             "4s-4px")
        self.assertEqual(cohp_label.cohp, ref["COHP"])
        self.assertEqual(cohp_label.icohp, ref["ICOHP"])
        orbitals = [[Orbital.s, Orbital.px], ["s", "px"], [0, 3]]
        cohps = [self.cohp_orb.get_orbital_resolved_cohp("1",
                                                         [[4, orb[0]], [4, orb[1]]]) for orb in orbitals]
        # print(cohps)
        for cohp in cohps:
            self.assertEqual(cohp.as_dict(), cohp_label.as_dict())


if __name__ == "__main__":
    unittest.main()
