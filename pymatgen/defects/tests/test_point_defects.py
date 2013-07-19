#!/usr/bin/python

import unittest
import sys

import numpy as np

from pymatgen.defects.point_defects import *
from pymatgen.core.structure import  Structure
from pymatgen.core.lattice import Lattice
from pymatgen.core.periodic_table import  Element
from pymatgen.symmetry.finder import SymmetryFinder 
from pymatgen.symmetry.spacegroup import Spacegroup


class StructureWithValenceIonicRadiusTest(unittest.TestCase):
    def setUp(self):
        """
        Setup MgO rocksalt structure for testing Vacancy
        """
        mgo_latt = [[4.212, 0, 0], [0, 4.212, 0], [0, 0, 4.212]]
        mgo_specie = ["Mg"]*4 + ["O"]*4
        mgo_frac_cord = [[0,0,0], [0.5,0.5,0], [0.5,0,0.5], [0, 0.5,0.5],
                         [0.5,0,0], [0,0.5,0], [0,0,0.5], [0.5,0.5,0.5]]
        self._mgo_uc = Structure(mgo_latt, mgo_specie, mgo_frac_cord, True, True)
        self._mgo_withvalrad = StructWithValenceIonicRadius(self._mgo_uc)
    
    def test_valences(self):
        valence_dict = self._mgo_withvalrad.valences
        for val in valence_dict.values():
            self.assertTrue(val in  {2,-2})

    def test_radii(self):
        rad_dict = self._mgo_withvalrad.radii
        for rad in rad_dict.values():
            self.assertTrue(rad in {0.86,1.26})


class VacancyTest(unittest.TestCase):
    def setUp(self):
        """
        Setup MgO rocksalt structure for testing Vacancy
        """
        mgo_latt = [[4.212, 0, 0], [0, 4.212, 0], [0, 0, 4.212]]
        mgo_specie = ["Mg"]*4 +  ["O"]*4
        mgo_frac_cord = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
                         [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5], [0.5, 0.5, 0.5]]
        self._mgo_uc1 = Structure(mgo_latt, mgo_specie, mgo_frac_cord, True, True)
        self._mgo_uc = StructWithValenceIonicRadius(self._mgo_uc1)
        self._mgo_vac = Vacancy(self._mgo_uc)
    
    def test_defectsite_count(self):
        self.assertTrue(self._mgo_vac.defectsite_count() == 2, 
                "Vacancy count wrong")

    def test_enumerate_defectsites(self):
        """
        The vacancy sites should be the lattice sites. 
        And there should be only two unique vacancy sites for MgO.
        """
        uniq_sites = []
        uniq_sites.append(self._mgo_uc1.sites[3])
        uniq_sites.append(self._mgo_uc1.sites[7])
        uniq_def_sites = self._mgo_vac.enumerate_defectsites()
        #Match uniq_sites iwth uniq_def_sites
        #self.assertTrue(len(uniq_def_sites) == 2, "Vacancy init failed")
        #mgo_spg = Spacegroup(int_number=225)
        #self.assertTrue(mgo_spg.are_symmetrically_equivalent(uniq_sites,
        #                uniq_def_sites),  "Vacancy init failed")

    def test_get_defectsite_index(self):
        for i in range(self._mgo_vac.defectsite_count()):
            self.assertTrue(self._mgo_vac.get_defectsite_structure_index(i) < 
                    len(self._mgo_uc1.sites), "Defect site index beyond range")

    def test_gt_defectsite_coordination_number(self):
        for i in range(self._mgo_vac.defectsite_count()):
            self.assertTrue(round(self._mgo_vac.get_defectsite_coordination_number(
                i))==6.0, "Wrong coordination number")

    def test_get_defectsite_coordinated_elements(self):
        for i in range(self._mgo_vac.defectsite_count()):
            site_index = self._mgo_vac.get_defectsite_structure_index(i)
            site_el = self._mgo_uc1[site_index].species_string
            self.assertTrue(
                    site_el not in self._mgo_vac.get_coordinated_elements(
                        i), "Coordinated elements are wrong")

    def test_get_defectsite_effective_charge(self):
        for i in range(self._mgo_vac.defectsite_count()):
            site_index = self._mgo_vac.get_defectsite_structure_index(i)
            site_el = self._mgo_uc1[site_index].species_and_occu
            eff_charge = self._mgo_vac.get_defectsite_effective_charge(i)
            if site_el["Mg"] == 1:
                self.assertEqual(eff_charge, -2)
            if site_el["O"] == 1:
                self.assertEqual(eff_charge, 2)

    def test_get_coordinatedsites_min_max_charge(self):
        for i in range(self._mgo_vac.defectsite_count()):
            min_chrg, max_chrg = self._mgo_vac.get_coordsites_min_max_charge(i)
            self.assertEqual(min_chrg, max_chrg)

    def test_get_volume(self):
        for i in range(self._mgo_vac.defectsite_count()):
            vol = self._mgo_vac.get_volume(i)
            #Once the zeo++ is properly working, make sure vol is +ve
            self.assertIsInstance(vol, float)

    def test_get_surface_area(self):
        for i in range(self._mgo_vac.defectsite_count()):
            sa = self._mgo_vac.get_surface_area(i)
            #Once the zeo++ is properly working, make sure vol is +ve
            self.assertIsInstance(sa, float)
        
class VacancyFormationEnergyTest(unittest.TestCase):
    def setUp(self):
        mgo_latt = [[4.212, 0, 0], [0, 4.212, 0], [0, 0, 4.212]]
        mgo_specie = ["Mg"]*4 +  ["O"]*4
        mgo_frac_cord = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
                         [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5], [0.5, 0.5, 0.5]]
        self.mgo_uc = Structure(mgo_latt, mgo_specie, mgo_frac_cord, True, True)
        mgo_uc_valrad = StructWithValenceIonicRadius(self.mgo_uc)
        self.mgo_vac = Vacancy(mgo_uc_valrad)
        self.mgo_vfe = VacancyFormationEnergy(self.mgo_vac)
        
    def test_get_energy(self):
        for i in range(len(self.mgo_vac.enumerate_defectsites())):
            vfe = self.mgo_vfe.get_energy(i)
            print vfe
            self.assertIsInstance(vfe, float)

class InterstitialTest(unittest.TestCase):
    def setUp(self):
        """
        Setup MgO rocksalt structure for testing Interstitial
        """
        mgo_latt = [[4.212, 0, 0], [0, 4.212, 0], [0, 0, 4.212]]
        mgo_specie = ["Mg"]*4 +  ["O"]*4
        mgo_frac_cord = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
                         [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5], [0.5, 0.5, 0.5]]
        self._mgo_uc1 = Structure(mgo_latt, mgo_specie, mgo_frac_cord, True, True)
        self._mgo_uc = StructWithValenceIonicRadius(self._mgo_uc1)
        self._mgo_interstitial = Interstitial(self._mgo_uc)
    
    def test_enumerate_defectsites(self):
        """
        The interstitial sites should be within the lattice
        """
        uniq_def_sites = self._mgo_interstitial.enumerate_defectsites()
        self.assertTrue(len(uniq_def_sites) == 2, "Interstitial init failed")
        #mgo_spg = Spacegroup(int_number=225)
        #self.assertTrue(mgo_spg.are_symmetrically_equivalent(uniq_sites,
        #                uniq_def_sites),  "Vacancy init failed")

    def test_defectsite_count(self):
        self.assertTrue(self._mgo_interstitial.defectsite_count() == 2, 
                "Vacancy count wrong")
        
    def test_get_defectsite_coordination_number(self):
        for i in range(self._mgo_interstitial.defectsite_count()):
            print >>sys.stderr, self._mgo_interstitial.get_defectsite_coordination_number(i)

    def test_get_coordsites_charge_sum(self):
        for i in range(self._mgo_interstitial.defectsite_count()):
            print >>sys.stderr, self._mgo_interstitial.get_coordsites_charge_sum(i)

    def test_get_defectsite_coordinated_elements(self):
        struct_el = self._mgo_uc1.composition.elements
        for i in range(self._mgo_interstitial.defectsite_count()):
            for el in self._mgo_interstitial.get_coordinated_elements(i):
                self.assertTrue(
                        Element(el) in struct_el, "Coordinated elements are wrong"
                        )

    def test_get_radius(self):
        for i in range(self._mgo_interstitial.defectsite_count()):
            rad = self._mgo_interstitial.get_radius(i)
            print >> sys.stderr, rad
            self.assertTrue(rad, float)

class InterstitialFormationEnergyTest(unittest.TestCase):
    def setUp(self):
        mgo_latt = [[4.212, 0, 0], [0, 4.212, 0], [0, 0, 4.212]]
        mgo_specie = ["Mg"]*4 +  ["O"]*4
        mgo_frac_cord = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
                         [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5], [0.5, 0.5, 0.5]]
        self.mgo_uc = Structure(mgo_latt, mgo_specie, mgo_frac_cord, True, True)
        mgo_uc_valrad = StructWithValenceIonicRadius(self.mgo_uc)
        self.mgo_inter = Interstitial(mgo_uc_valrad)
        self.mgo_ife = InterstitialFormationEnergy(self.mgo_inter)
        
    def test_get_relaxedenergy(self):
        for i in range(len(self.mgo_inter.enumerate_defectsites())):
            ife = self.mgo_ife.get_energy(i, "Mg", 2, 0.4, True)
            site_coords = self.mgo_inter.get_defectsite(i).coords
            site_radius = self.mgo_inter.get_radius(i)
            print i, site_coords, site_radius, ife
            self.assertIsInstance(ife, float)
    def test_get_norelaxedenergy(self):
        for i in range(len(self.mgo_inter.enumerate_defectsites())):
            ife = self.mgo_ife.get_energy(i, "Mg", 2, 1, False)
            site_coords = self.mgo_inter.get_defectsite(i).coords
            site_radius = self.mgo_inter.get_radius(i)
            print i, site_coords, site_radius, ife
            self.assertIsInstance(ife, float)


if __name__ == "__main__":
    #unittest.main()
    #suite = unittest.TestLoader().loadTestsFromTestCase(StructureWithValenceIonicRadiusTest)
    #suite = unittest.TestLoader().loadTestsFromTestCase(InterstitialTest)
    #suite = unittest.TestLoader().loadTestsFromTestCase(VacancyTest)
    suite = unittest.TestLoader().loadTestsFromTestCase(VacancyFormationEnergyTest)
    suite = unittest.TestLoader().loadTestsFromTestCase(InterstitialFormationEnergyTest)
    unittest.TextTestRunner(verbosity=3).run(suite)
