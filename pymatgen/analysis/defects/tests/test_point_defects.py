# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals
import os
import unittest
import random
import matplotlib
matplotlib.use("pdf")
from pymatgen.util.testing import PymatgenTest
from pymatgen.analysis.local_env import ValenceIonicRadiusEvaluator
from pymatgen.analysis.defects.point_defects import *
from pymatgen.core import PeriodicSite
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.core.periodic_table import Element
from pymatgen.io.vasp.outputs import Chgcar
from pymatgen.analysis.bond_valence import BVAnalyzer
from monty.os.path import which
from pymatgen.io.cif import CifParser

try:
    import zeo
except ImportError:
    zeo = None

try:
    from skimage.feature import peak_local_max

    peak_local_max_found = True
except ImportError:
    peak_local_max_found = False

gulp_present = which('gulp')
test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        'test_files')


@unittest.skipIf(not zeo, "zeo not present.")
class VacancyTest(PymatgenTest):
    def setUp(self):
        """
        Setup MgO rocksalt structure for testing Vacancy
        """
        mgo_latt = [[4.212, 0, 0], [0, 4.212, 0], [0, 0, 4.212]]
        mgo_specie = ["Mg"] * 4 + ["O"] * 4
        mgo_frac_cord = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
                         [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5], [0.5, 0.5, 0.5]]
        self._mgo_uc = Structure(mgo_latt, mgo_specie, mgo_frac_cord, True,
                                 True)

        bv = BVAnalyzer()
        self._mgo_uc = bv.get_oxi_state_decorated_structure(self._mgo_uc)
        self._mgo_val_rad_eval = ValenceIonicRadiusEvaluator(self._mgo_uc)
        self._mgo_val = self._mgo_val_rad_eval.valences
        self._mgo_rad = self._mgo_val_rad_eval.radii
        self._mgo_vac = Vacancy(self._mgo_uc, self._mgo_val, self._mgo_rad)

    def test_defectsite_count(self):
        self.assertTrue(self._mgo_vac.defectsite_count() == 2,
                        "Vacancy count wrong")

    def test_enumerate_defectsites(self):
        """
        The vacancy sites should be the lattice sites.
        And there should be only two unique vacancy sites for MgO.
        """
        uniq_sites = []
        uniq_sites.append(self._mgo_uc.sites[3])
        uniq_sites.append(self._mgo_uc.sites[7])
        uniq_def_sites = self._mgo_vac.enumerate_defectsites()
        # Match uniq_sites iwth uniq_def_sites
        # self.assertTrue(len(uniq_def_sites) == 2, "Vacancy init failed")
        # mgo_spg = Spacegroup(int_number=225)
        # self.assertTrue(mgo_spg.are_symmetrically_equivalent(uniq_sites,
        #                uniq_def_sites),  "Vacancy init failed")

    def test_get_defectsite_index(self):
        for i in range(self._mgo_vac.defectsite_count()):
            self.assertTrue(self._mgo_vac.get_defectsite_structure_index(i) <
                            len(self._mgo_uc.sites),
                            "Defect site index beyond range")

    def test_gt_defectsite_coordination_number(self):
        for i in range(self._mgo_vac.defectsite_count()):
            self.assertTrue(
                round(self._mgo_vac.get_defectsite_coordination_number(
                    i)) == 6.0, "Wrong coordination number")

    def test_get_defectsite_coordinated_elements(self):
        for i in range(self._mgo_vac.defectsite_count()):
            site_index = self._mgo_vac.get_defectsite_structure_index(i)
            site_el = self._mgo_uc[site_index].species_string
            self.assertTrue(
                site_el not in self._mgo_vac.get_coordinated_elements(
                    i), "Coordinated elements are wrong")

    def test_get_defectsite_effective_charge(self):
        for i in range(self._mgo_vac.defectsite_count()):
            site_index = self._mgo_vac.get_defectsite_structure_index(i)
            site_el = self._mgo_uc[site_index].species_and_occu
            eff_charge = self._mgo_vac.get_defectsite_effective_charge(i)
            if site_el["Mg2+"] == 1:
                self.assertEqual(eff_charge, -2)
            if site_el["O2-"] == 1:
                self.assertEqual(eff_charge, 2)

    def test_get_coordinatedsites_min_max_charge(self):
        for i in range(self._mgo_vac.defectsite_count()):
            min_chrg, max_chrg = self._mgo_vac.get_coordsites_min_max_charge(i)
            self.assertEqual(min_chrg, max_chrg)

    def test_make_supercells_with_defects(self):
        scaling_matrix = [2, 2, 2]
        vac_specie = ['Mg']
        vac_scs = self._mgo_vac.make_supercells_with_defects(
            scaling_matrix, vac_specie)
        expected_structure_formulae = ["Mg32 O32", "Mg32 O31", "Mg31 O32"]
        self.assertEqual(len(vac_scs), 2)
        for sc in vac_scs:
            self.assertIn(sc.formula, expected_structure_formulae)

        vac_scs = self._mgo_vac.make_supercells_with_defects(scaling_matrix)
        expected_structure_formulae = ["Mg32 O32", "Mg32 O31", "Mg31 O32"]
        self.assertEqual(len(vac_scs), 3)
        for sc in vac_scs:
            self.assertIn(sc.formula, expected_structure_formulae)

        vac_scs = self._mgo_vac.make_supercells_with_defects(
            scaling_matrix, limit_return_structures=1)
        expected_structure_formulae = ["Mg32 O32", "Mg32 O31", "Mg31 O32"]
        self.assertEqual(len(vac_scs), 2)
        for sc in vac_scs:
            self.assertIn(sc.formula, expected_structure_formulae)


@unittest.skipIf(not gulp_present, "gulp not present.")
class VacancyFormationEnergyTest(PymatgenTest):
    def setUp(self):
        mgo_latt = [[4.212, 0, 0], [0, 4.212, 0], [0, 0, 4.212]]
        mgo_specie = ["Mg"] * 4 + ["O"] * 4
        mgo_frac_cord = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
                         [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5], [0.5, 0.5, 0.5]]
        self.mgo_uc = Structure(mgo_latt, mgo_specie, mgo_frac_cord, True, True)
        mgo_valrad_eval = ValenceIonicRadiusEvaluator(self.mgo_uc)
        val = mgo_valrad_eval.valences
        rad = mgo_valrad_eval.radii
        self.mgo_vac = Vacancy(self.mgo_uc, val, rad)
        self.mgo_vfe = VacancyFormationEnergy(self.mgo_vac)

        # This test doesn't pass!
        # def test_get_energy(self):
        #    for i in range(len(self.mgo_vac.enumerate_defectsites())):
        #        vfe = self.mgo_vfe.get_energy(i)
        #        print(vfe)
        #        self.assertIsInstance(vfe, float)


@unittest.skipIf(not zeo, "zeo not present.")
class InterstitialTest(PymatgenTest):
    def setUp(self):
        """
        Setup MgO rocksalt structure for testing Interstitial
        """
        mgo_latt = [[4.212, 0, 0], [0, 4.212, 0], [0, 0, 4.212]]
        mgo_specie = ["Mg"] * 4 + ["O"] * 4
        mgo_frac_cord = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
                         [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5], [0.5, 0.5, 0.5]]
        self._mgo_uc = Structure(mgo_latt, mgo_specie, mgo_frac_cord, True,
                                 True)
        mgo_val_rad_eval = ValenceIonicRadiusEvaluator(self._mgo_uc)
        self._mgo_val = mgo_val_rad_eval.valences
        self._mgo_rad = mgo_val_rad_eval.radii
        self._mgo_interstitial = Interstitial(
            self._mgo_uc, self._mgo_val, self._mgo_rad
        )

    def test_enumerate_defectsites(self):
        """
        The interstitial sites should be within the lattice
        """
        uniq_def_sites = self._mgo_interstitial.enumerate_defectsites()
        self.assertTrue(len(uniq_def_sites) == 2, "Interstitial init failed")

    def test_defectsite_count(self):
        print(self._mgo_interstitial.defectsite_count())
        self.assertTrue(self._mgo_interstitial.defectsite_count() == 2,
                        "Vacancy count wrong")

    def test_get_defectsite_coordination_number(self):
        for i in range(self._mgo_interstitial.defectsite_count()):
            print(self._mgo_interstitial.get_defectsite_coordination_number(
                i))

    def test_get_coordinated_sites(self):
        for i in range(self._mgo_interstitial.defectsite_count()):
            print(self._mgo_interstitial.get_coordinated_sites(
                i))

    def test_get_coordsites_charge_sum(self):
        for i in range(self._mgo_interstitial.defectsite_count()):
            print(self._mgo_interstitial.get_coordsites_charge_sum(
                i))

    def test_get_defectsite_coordinated_elements(self):
        struct_el = self._mgo_uc.composition.elements
        for i in range(self._mgo_interstitial.defectsite_count()):
            for el in self._mgo_interstitial.get_coordinated_elements(i):
                self.assertTrue(
                    Element(el) in struct_el, "Coordinated elements are wrong"
                )

    def test_get_radius(self):
        for i in range(self._mgo_interstitial.defectsite_count()):
            rad = self._mgo_interstitial.get_radius(i)
            print(rad)
            self.assertTrue(rad, float)


@unittest.skipIf(not zeo, "zeo not present.")
class InterstitialVoronoiFaceCenterTest(unittest.TestCase):
    def setUp(self):
        """
        Setup MgO rocksalt structure for testing Interstitial
        """
        mgo_latt = [[4.212, 0, 0], [0, 4.212, 0], [0, 0, 4.212]]
        mgo_specie = ["Mg"] * 4 + ["O"] * 4
        mgo_frac_cord = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
                         [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5], [0.5, 0.5, 0.5]]
        self._mgo_uc = Structure(mgo_latt, mgo_specie, mgo_frac_cord, True,
                                 True)
        mgo_val_rad_eval = ValenceIonicRadiusEvaluator(self._mgo_uc)
        self._mgo_val = mgo_val_rad_eval.valences
        self._mgo_rad = mgo_val_rad_eval.radii
        self._mgo_interstitial = Interstitial(
            self._mgo_uc, self._mgo_val, self._mgo_rad,
            site_type='voronoi_facecenter'
        )

    def test_enumerate_defectsites(self):
        """
        The interstitial sites should be within the lattice
        """
        uniq_def_sites = self._mgo_interstitial.enumerate_defectsites()
        print("Length of voronoi face centers", len(uniq_def_sites))
        self.assertTrue(len(uniq_def_sites) == 2, "Defect site count wrong")

    def test_defectsite_count(self):
        print(self._mgo_interstitial.defectsite_count())
        self.assertTrue(self._mgo_interstitial.defectsite_count() == 2,
                        "Vacancy count wrong")

    def test_get_defectsite_coordination_number(self):
        for i in range(self._mgo_interstitial.defectsite_count()):
            coord_no = \
                self._mgo_interstitial.get_defectsite_coordination_number(
                    i)
            self.assertTrue(isinstance(coord_no, float))

    def test_get_coordinated_sites(self):
        for i in range(self._mgo_interstitial.defectsite_count()):
            print(self._mgo_interstitial.get_coordinated_sites(
                i))

    def test_get_coordsites_charge_sum(self):
        for i in range(self._mgo_interstitial.defectsite_count()):
            print(self._mgo_interstitial.get_coordsites_charge_sum(
                i))

    def test_get_defectsite_coordinated_elements(self):
        struct_el = self._mgo_uc.composition.elements
        for i in range(self._mgo_interstitial.defectsite_count()):
            for el in self._mgo_interstitial.get_coordinated_elements(i):
                self.assertTrue(
                    Element(el) in struct_el, "Coordinated elements are wrong"
                )

    def test_get_radius(self):
        for i in range(self._mgo_interstitial.defectsite_count()):
            rad = self._mgo_interstitial.get_radius(i)
            self.assertAlmostEqual(rad, 0.0)


@unittest.skipIf(not zeo, "zeo not present.")
class InterstitialHighAccuracyTest(unittest.TestCase):
    def setUp(self):
        """
        Setup MgO rocksalt structure for testing Interstitial
        """
        mgo_latt = [[4.212, 0, 0], [0, 4.212, 0], [0, 0, 4.212]]
        mgo_specie = ["Mg"] * 4 + ["O"] * 4
        mgo_frac_cord = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
                         [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5], [0.5, 0.5, 0.5]]
        self._mgo_uc = Structure(mgo_latt, mgo_specie, mgo_frac_cord, True,
                                 True)
        mgo_val_rad_eval = ValenceIonicRadiusEvaluator(self._mgo_uc)
        self._mgo_val = mgo_val_rad_eval.valences
        self._mgo_rad = mgo_val_rad_eval.radii
        self._mgo_interstitial = Interstitial(
            self._mgo_uc, self._mgo_val, self._mgo_rad, accuracy='High'
        )

    def test_enumerate_defectsites(self):
        """
        The interstitial sites should be within the lattice
        """
        uniq_def_sites = self._mgo_interstitial.enumerate_defectsites()
        for site in uniq_def_sites:
            self.assertIsInstance(site, PeriodicSite,
                                  "Returned objects are not sites")

            # print len(uniq_def_sites)
            # self.assertTrue(len(uniq_def_sites) == 2, "Interstitial init
            # failed")

    def test_defectsite_count(self):
        self.assertIsNotNone(self._mgo_interstitial.defectsite_count(),
                             "Interstitial count wrong")

    def test_get_defectsite_coordination_number(self):
        for i in range(self._mgo_interstitial.defectsite_count()):
            print(self._mgo_interstitial.get_defectsite_coordination_number(
                i))

    def test_get_coordinated_sites(self):
        for i in range(self._mgo_interstitial.defectsite_count()):
            print(self._mgo_interstitial.get_coordinated_sites(
                i))

    def test_get_coordsites_charge_sum(self):
        for i in range(self._mgo_interstitial.defectsite_count()):
            print(self._mgo_interstitial.get_coordsites_charge_sum(
                i))

    def test_get_defectsite_coordinated_elements(self):
        struct_el = self._mgo_uc.composition.elements
        for i in range(self._mgo_interstitial.defectsite_count()):
            for el in self._mgo_interstitial.get_coordinated_elements(i):
                self.assertTrue(
                    Element(el) in struct_el, "Coordinated elements are wrong"
                )

    def test_get_radius(self):
        for i in range(self._mgo_interstitial.defectsite_count()):
            rad = self._mgo_interstitial.get_radius(i)
            print(rad)
            self.assertTrue(rad, float)


"""
Some of the tests are nearly useless. Better tests are needed
"""


@unittest.skipIf(not (gulp_present and zeo), "gulp or zeo not present.")
class InterstitialAnalyzerTest(unittest.TestCase):
    def setUp(self):
        mgo_latt = [[4.212, 0, 0], [0, 4.212, 0], [0, 0, 4.212]]
        mgo_specie = ["Mg"] * 4 + ["O"] * 4
        mgo_frac_cord = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
                         [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5], [0.5, 0.5, 0.5]]
        self.mgo_uc = Structure(mgo_latt, mgo_specie, mgo_frac_cord, True, True)
        mgo_valrad_eval = ValenceIonicRadiusEvaluator(self.mgo_uc)
        val = mgo_valrad_eval.valences
        rad = mgo_valrad_eval.radii
        self.mgo_val = val
        self.mgo_rad = rad
        self.mgo_inter = Interstitial(self.mgo_uc, val, rad)
        self.mgo_ia = InterstitialAnalyzer(self.mgo_inter, 'Mg', 2)

    def test_get_relaxedenergy(self):
        for i in range(len(self.mgo_inter.enumerate_defectsites())):
            ife = self.mgo_ia.get_energy(i, True)
            site_coords = self.mgo_inter.get_defectsite(i).coords
            site_radius = self.mgo_inter.get_radius(i)
            # print(i, site_coords, site_radius, ife)
            self.assertIsInstance(ife, float)

    def test_get_norelaxedenergy(self):
        for i in range(self.mgo_inter.defectsite_count()):
            ife = self.mgo_ia.get_energy(i, False)
            site_coords = self.mgo_inter.get_defectsite(i).coords
            site_radius = self.mgo_inter.get_radius(i)
            print(i, site_coords, site_radius, ife)
            self.assertIsInstance(ife, float)

    def test_get_percentage_volume_change(self):
        for i in range(self.mgo_inter.defectsite_count()):
            del_vol = self.mgo_ia.get_percentage_volume_change(i)
            # print(i, del_vol)

    def test_get_percentage_lattice_parameter_change(self):
        for i in range(self.mgo_inter.defectsite_count()):
            del_lat = self.mgo_ia.get_percentage_lattice_parameter_change(i)
            # print(i, del_lat)

    def test_get_percentage_bond_distance_change(self):
        for i in range(self.mgo_inter.defectsite_count()):
            del_bd = self.mgo_ia.get_percentage_bond_distance_change(i)
            # print(i, del_bd)

    def test_relaxed_structure_match(self):
        for i in range(self.mgo_inter.defectsite_count()):
            for j in range(self.mgo_inter.defectsite_count()):
                match = self.mgo_ia.relaxed_structure_match(i, j)
                # print(i, j, match)
                if i == j:
                    self.assertTrue(match)


@unittest.skipIf(not (gulp_present and zeo), "gulp or zeo not present.")
class InterstitialStructureRelaxerTest(unittest.TestCase):
    def setUp(self):
        mgo_latt = [[4.212, 0, 0], [0, 4.212, 0], [0, 0, 4.212]]
        mgo_specie = ["Mg"] * 4 + ["O"] * 4
        mgo_frac_cord = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
                         [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5], [0.5, 0.5, 0.5]]
        self.mgo_uc = Structure(mgo_latt, mgo_specie, mgo_frac_cord, True, True)
        mgo_valrad_eval = ValenceIonicRadiusEvaluator(self.mgo_uc)
        val = mgo_valrad_eval.valences
        rad = mgo_valrad_eval.radii
        self.mgo_val = val
        self.mgo_rad = rad
        self.mgo_inter = Interstitial(self.mgo_uc, val, rad)
        self.isr = InterstitialStructureRelaxer(self.mgo_inter, 'Mg', 2)

    def test_relaxed_structure_match(self):
        for i in range(self.mgo_inter.defectsite_count()):
            for j in range(self.mgo_inter.defectsite_count()):
                match = self.isr.relaxed_structure_match(i, j)
                # print i, j, match
                if i == j:
                    self.assertTrue(match)

    def test_relaxed_energy_match(self):
        for i in range(self.mgo_inter.defectsite_count()):
            for j in range(self.mgo_inter.defectsite_count()):
                match = self.isr.relaxed_energy_match(i, j)
                # print i, j, match
                if i == j:
                    self.assertTrue(match)

    def test_get_relaxed_structure(self):
        for i in range(self.mgo_inter.defectsite_count()):
            relax_struct = self.isr.get_relaxed_structure(i)
            self.assertIsInstance(relax_struct, Structure)

    def test_get_relaxed_energy(self):
        for i in range(self.mgo_inter.defectsite_count()):
            energy = self.isr.get_relaxed_energy(i)
            self.assertIsInstance(energy, float)

    def test_get_relaxed_interstitial(self):
        ri = self.isr.get_relaxed_interstitial()
        self.assertIsInstance(ri, RelaxedInterstitial)


@unittest.skipIf(not (gulp_present and zeo), "gulp or zeo not present.")
class RelaxedInsterstitialTest(unittest.TestCase):
    def setUp(self):
        mgo_latt = [[4.212, 0, 0], [0, 4.212, 0], [0, 0, 4.212]]
        mgo_specie = ["Mg"] * 4 + ["O"] * 4
        mgo_frac_cord = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
                         [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5], [0.5, 0.5, 0.5]]
        self.mgo_uc = Structure(mgo_latt, mgo_specie, mgo_frac_cord, True, True)
        mgo_valrad_eval = ValenceIonicRadiusEvaluator(self.mgo_uc)
        val = mgo_valrad_eval.valences
        rad = mgo_valrad_eval.radii
        self.mgo_val = val
        self.mgo_rad = rad
        self.mgo_inter = Interstitial(self.mgo_uc, val, rad)
        isr = InterstitialStructureRelaxer(self.mgo_inter, 'Mg', 2)
        self.ri = isr.get_relaxed_interstitial()

    def test_formation_energy(self):
        for i in range(self.mgo_inter.defectsite_count()):
            ife = self.ri.formation_energy(i)
            self.assertIsInstance(ife, float)
            print("ife", ife)

    def test_get_percentage_volume_change(self):
        for i in range(self.mgo_inter.defectsite_count()):
            del_vol = self.ri.get_percentage_volume_change(i)
            self.assertIsInstance(del_vol, float)
            print("del_vol", del_vol)

    def test_get_percentage_lattice_parameter_change(self):
        for i in range(self.mgo_inter.defectsite_count()):
            del_lat = self.ri.get_percentage_lattice_parameter_change(i)
            self.assertNotEqual(del_lat['a'], 0)
            self.assertNotEqual(del_lat['b'], 0)
            self.assertNotEqual(del_lat['c'], 0)

    def test_get_percentage_bond_distance_change(self):
        for i in range(self.mgo_inter.defectsite_count()):
            del_bd = self.ri.get_percentage_bond_distance_change(i)
            # self.assertIsInstance(del_bd, float)
            # print del_bd


class StructureMotifInterstitialTest(PymatgenTest):
    def setUp(self):
        self.silicon = Structure(
            Lattice.from_lengths_and_angles(
                [5.47, 5.47, 5.47],
                [90.0, 90.0, 90.0]),
            ["Si", "Si", "Si", "Si", "Si", "Si", "Si", "Si"],
            [[0.000000, 0.000000, 0.500000],
             [0.750000, 0.750000, 0.750000],
             [0.000000, 0.500000, 1.000000],
             [0.750000, 0.250000, 0.250000],
             [0.500000, 0.000000, 1.000000],
             [0.250000, 0.750000, 0.250000],
             [0.500000, 0.500000, 0.500000],
             [0.250000, 0.250000, 0.750000]],
            validate_proximity=False, to_unit_cell=False,
            coords_are_cartesian=False, site_properties=None)
        self.smi = StructureMotifInterstitial(self.silicon, "Si",
                                              motif_types=["tet", "oct"],
                                              op_threshs=[0.3, 0.5],
                                              dl=0.4, doverlap=1.0,
                                              facmaxdl=1.01)
        self.diamond = Structure(
            Lattice([[2.189, 0, 1.264], [0.73, 2.064, 1.264],
                     [0, 0, 2.528]]), ["C0+", "C0+"], [[2.554, 1.806, 4.423],
                                                       [0.365, 0.258, 0.632]],
            validate_proximity=False,
            to_unit_cell=False, coords_are_cartesian=True,
            site_properties=None)
        self.nacl = Structure(
            Lattice([[3.485, 0, 2.012], [1.162, 3.286, 2.012],
                     [0, 0, 4.025]]), ["Na1+", "Cl1-"], [[0, 0, 0],
                                                         [2.324, 1.643, 4.025]],
            validate_proximity=False,
            to_unit_cell=False, coords_are_cartesian=True,
            site_properties=None)
        self.cscl = Structure(
            Lattice([[4.209, 0, 0], [0, 4.209, 0], [0, 0, 4.209]]),
            ["Cl1-", "Cs1+"], [[2.105, 2.105, 2.105], [0, 0, 0]],
            validate_proximity=False, to_unit_cell=False,
            coords_are_cartesian=True, site_properties=None)
        self.square_pyramid = Structure(
            Lattice([[100, 0, 0], [0, 100, 0], [0, 0, 100]]),
            ["C", "C", "C", "C", "C", "C"], [
                [0, 0, 0], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], \
                [0, 0, 1]], validate_proximity=False, to_unit_cell=False,
            coords_are_cartesian=True, site_properties=None)
        self.trigonal_bipyramid = Structure(
            Lattice([[100, 0, 0], [0, 100, 0], [0, 0, 100]]),
            ["P", "Cl", "Cl", "Cl", "Cl", "Cl"], [
                [0, 0, 0], [0, 0, 2.14], [0, 2.02, 0], [1.74937, -1.01, 0], \
                [-1.74937, -1.01, 0], [0, 0, -2.14]], validate_proximity=False,
            to_unit_cell=False, coords_are_cartesian=True,
            site_properties=None)

    def test_all(self):
        self.assertIsInstance(self.smi, StructureMotifInterstitial)

        self.assertEqual(len(self.smi.enumerate_defectsites()), 1)
        self.assertIsInstance(self.smi.enumerate_defectsites()[0], PeriodicSite)
        self.assertEqual("Si",
                         self.smi.enumerate_defectsites()[0].species_string)

        self.assertEqual("tet", self.smi.get_motif_type(0))

        elem_cn_dict = self.smi.get_coordinating_elements_cns(0)
        self.assertEqual(len(list(elem_cn_dict.keys())), 1)
        self.assertEqual(list(elem_cn_dict.keys())[0], "Si")
        self.assertEqual(elem_cn_dict["Si"], 4)

        structs = self.smi.make_supercells_with_defects(np.array([1, 1, 1]))
        self.assertEqual(len(structs), 2)
        self.assertIsInstance(structs[0], Structure)

    def tearDown(self):
        del self.smi
        del self.silicon
        del self.diamond
        del self.nacl
        del self.cscl


class TopographyAnalyzerTest(unittest.TestCase):
    def setUp(self):
        feo4 = Structure.from_file(os.path.join(test_dir, "LiFePO4.cif"))
        feo4.remove_species(["Li"])
        feo4.remove_oxidation_states()
        self.feo4 = feo4

    def test_topography_analyzer(self):
        # check interstitial sites for FePO4 using Voronoi Tessellation
        vor_feo4 = TopographyAnalyzer(self.feo4,
                                      framework_ions=["O"],
                                      cations=["P", "Fe"], check_volume=False)
        vor_feo4.cluster_nodes(tol=1.2)
        vor_feo4.remove_collisions(1.2)
        s_feo4 = vor_feo4.get_structure_with_nodes()
        sites_feo4 = np.array(
            [s_feo4[i].frac_coords for i in range(len(s_feo4))
             if s_feo4[i].species_string == "X0+"])

        # check total number of vnodes
        self.assertAlmostEqual(len(vor_feo4.vnodes), 24)

        # check four sites that match Li sites in LiFePO4(mp-19017)
        site_predicted = [[0, 0, 0], [0.5, 0.5, 0.5], [0.5, 0, 0.5],
                          [0, 0.5, 0]]
        for i in range(0, 4):
            is_site_matched = False
            for site in sites_feo4:
                distance = s_feo4.lattice. \
                    get_distance_and_image(site, site_predicted[i])
                if distance[0] < 0.01:
                    is_site_matched = True
                else:
                    continue
            self.assertTrue(is_site_matched)


@unittest.skipIf(not peak_local_max_found,
                 "skimage.feature.peak_local_max module not present.")
class ChgDenAnalyzerTest(unittest.TestCase):
    def setUp(self):
        # This is a CHGCAR_sum file with reduced grid size
        chgcar_path = os.path.join(test_dir, "CHGCAR.FePO4")
        chg_FePO4 = Chgcar.from_file(chgcar_path)
        self.chgcar_path = chgcar_path
        self.chg_FePO4 = chg_FePO4
        self.ca_FePO4 = ChargeDensityAnalyzer(chg_FePO4)
        self.s_LiFePO4 = Structure.from_file(
            os.path.join(test_dir, "LiFePO4.cif"))

    def test_get_local_extrema(self):
        ca = ChargeDensityAnalyzer.from_file(self.chgcar_path)
        threshold_frac = random.random()
        threshold_abs_min = random.randrange(2, 14)
        threshold_abs_max = random.randrange(27e2, 28e4)

        # Minima test
        full_list_min = self.ca_FePO4.get_local_extrema(
            find_min=True, threshold_frac=1.0)
        frac_list_min_frac = self.ca_FePO4.get_local_extrema(
            find_min=True, threshold_frac=threshold_frac)
        frac_list_min_abs = self.ca_FePO4.get_local_extrema(
            find_min=True, threshold_abs=threshold_abs_min)

        self.assertAlmostEqual(len(full_list_min) * threshold_frac,
                               len(frac_list_min_frac), delta=1)

        ca.get_local_extrema(find_min=True)
        df_expected = ca.extrema_df[
            ca.extrema_df["Charge Density"] <= threshold_abs_min]
        self.assertEqual(len(frac_list_min_abs), len(df_expected))

        # Maxima test
        full_list_max = self.ca_FePO4.get_local_extrema(
            find_min=False, threshold_frac=1.0)
        frac_list_max = self.ca_FePO4.get_local_extrema(
            find_min=False, threshold_frac=threshold_frac)
        frac_list_max_abs = self.ca_FePO4.get_local_extrema(
            find_min=False, threshold_abs=threshold_abs_max)

        self.assertAlmostEqual(len(full_list_max) * threshold_frac,
                               len(frac_list_max), delta=1)

        # Local maxima should finds all center of atoms
        self.assertEqual(len(self.ca_FePO4.structure), len(full_list_max))

        ca.get_local_extrema(find_min=False)
        df_expected = ca.extrema_df[
            ca.extrema_df["Charge Density"] >= threshold_abs_max]
        self.assertEqual(len(frac_list_max_abs), len(df_expected))

    def test_remove_collisions(self):
        ca = ChargeDensityAnalyzer(self.chg_FePO4)
        ca.get_local_extrema(threshold_frac=0)
        ca.remove_collisions()  # should not trigger error
        self.assertEqual(len(ca.extrema_df), 0)

        self.ca_FePO4.get_local_extrema(
            find_min=False, threshold_frac=1.0)
        self.ca_FePO4.remove_collisions(min_dist=0.5)
        self.assertEqual(len(self.ca_FePO4.extrema_df), 0)

    def test_cluster_nodes(self):
        ca = ChargeDensityAnalyzer(self.chg_FePO4)
        ca.get_local_extrema()
        ca.cluster_nodes(tol=20)
        self.assertEqual(len(ca.extrema_df), 1)

    def test_get_structure_with_nodes(self):
        s_FePO4 = self.ca_FePO4.get_structure_with_nodes(find_min=True)

        sites_predicted = np.array(
            [self.s_LiFePO4[i].frac_coords for i in range(len(self.s_LiFePO4))
             if self.s_LiFePO4[i].species_string == "Li"])
        sites_guess = np.array(
            [s_FePO4[i].frac_coords for i in range(len(s_FePO4))
             if s_FePO4[i].species_string == "X0+"])
        distances = s_FePO4.lattice.get_all_distances(sites_predicted,
                                                      sites_guess).flatten()
        distances = [d for d in distances if d < 0.1]
        self.assertEqual(len(distances), len(sites_predicted))

    def test_from_file(self):
        ca = ChargeDensityAnalyzer.from_file(self.chgcar_path)
        assert isinstance(ca, ChargeDensityAnalyzer)


if __name__ == "__main__":
    unittest.main()
    # suite = unittest.TestLoader().loadTestsFromTestCase(
    # ValenceIonicRadiusEvaluatorTest)
    # suite = unittest.TestLoader().loadTestsFromTestCase(InterstitialTest)
    # suite = unittest.TestLoader().loadTestsFromTestCase(VacancyTest)
    # suite = unittest.TestLoader().loadTestsFromTestCase(
    # VacancyFormationEnergyTest)
    # suite = unittest.TestLoader().loadTestsFromTestCase(
    # InterstitialAnalyzerTest)
    # unittest.TextTestRunner(verbosity=3).run(suite)