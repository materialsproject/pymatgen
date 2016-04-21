# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest2 as unittest
import sys

from pymatgen.analysis.defects.point_defects import *
from pymatgen.core.structure import Structure
from pymatgen.core.periodic_table import Element
from pymatgen.analysis.bond_valence import BVAnalyzer
from monty.os.path import which
from pymatgen.io.cif import CifParser

try:
    import zeo
except ImportError:
    zeo = None

gulp_present = which('gulp')
test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        'test_files')


class ValenceIonicRadiusEvaluatorTest(unittest.TestCase):
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
        self._mgo_valrad_evaluator = ValenceIonicRadiusEvaluator(self._mgo_uc)
        # self._si = Cssr.from_file("../../../../test_files/Si.cssr").structure
        # self._ci_valrad_evaluator = ValenceIonicRadiusEvaluator(self._si)

    def test_valences_ionic_structure(self):
        valence_dict = self._mgo_valrad_evaluator.valences
        for val in list(valence_dict.values()):
            self.assertTrue(val in {2, -2})

    def test_radii_ionic_structure(self):
        radii_dict = self._mgo_valrad_evaluator.radii
        for rad in list(radii_dict.values()):
            self.assertTrue(rad in {0.86, 1.26})


class ValenceIonicRadiusEvaluatorMultiOxiTest(unittest.TestCase):
    def setUp(self):
        """
        Setup Fe3O4  structure for testing multiple oxidation states
        """
        cif_ob = CifParser(os.path.join(test_dir, "Fe3O4.cif"))
        self._struct = cif_ob.get_structures()[0]
        self._valrad_evaluator = ValenceIonicRadiusEvaluator(self._struct)
        self._length = len(self._struct.sites)

    def test_valences_ionic_structure(self):
        valence_set = set(self._valrad_evaluator.valences.values())
        self.assertEqual(valence_set, {2, 3, -2})

    def test_radii_ionic_structure(self):
        radii_set = set(self._valrad_evaluator.radii.values())
        self.assertEqual(len(radii_set), 3)
        self.assertEqual(radii_set, {0.72, 0.75, 1.26})


@unittest.skipIf(not zeo, "zeo not present.")
class VacancyTest(unittest.TestCase):
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

    @unittest.skip("deprecated")
    def test_get_volume(self):
        for i in range(self._mgo_vac.defectsite_count()):
            vol = self._mgo_vac.get_volume(i)
            # Once the zeo++ is properly working, make sure vol is +ve
            self.assertIsInstance(vol, float)

    @unittest.skip("deprecated")
    def test_get_surface_area(self):
        for i in range(self._mgo_vac.defectsite_count()):
            sa = self._mgo_vac.get_surface_area(i)
            # Once the zeo++ is properly working, make sure vol is +ve
            self.assertIsInstance(sa, float)


@unittest.skipIf(not gulp_present, "gulp not present.")
class VacancyFormationEnergyTest(unittest.TestCase):
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

    def test_get_energy(self):
        for i in range(len(self.mgo_vac.enumerate_defectsites())):
            vfe = self.mgo_vfe.get_energy(i)
            print(vfe)
            self.assertIsInstance(vfe, float)


@unittest.skipIf(not zeo, "zeo not present.")
class InterstitialTest(unittest.TestCase):
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
