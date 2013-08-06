#!/usr/bin/python

import unittest
import sys

import numpy as np

from pymatgen.analysis.defects.point_defects import *
from pymatgen.core.structure import  Structure
from pymatgen.core.lattice import Lattice
from pymatgen.core.periodic_table import  Element
from pymatgen.symmetry.finder import SymmetryFinder 
from pymatgen.symmetry.spacegroup import Spacegroup


class ValenceIonicRadiusEvaluatorTest(unittest.TestCase):
    def setUp(self):
        """
        Setup MgO rocksalt structure for testing Vacancy
        """
        mgo_latt = [[4.212, 0, 0], [0, 4.212, 0], [0, 0, 4.212]]
        mgo_specie = ["Mg"]*4 + ["O"]*4
        mgo_frac_cord = [[0,0,0], [0.5,0.5,0], [0.5,0,0.5], [0, 0.5,0.5],
                         [0.5,0,0], [0,0.5,0], [0,0,0.5], [0.5,0.5,0.5]]
        self._mgo_uc = Structure(mgo_latt, mgo_specie, mgo_frac_cord, True, True)
        self._mgo_valrad_evaluator = ValenceIonicRadiusEvaluator(self._mgo_uc)
    
    def test_valences(self):
        valence_dict = self._mgo_valrad_evaluator.valences
        for val in valence_dict.values():
            self.assertTrue(val in  {2,-2})

    def test_radii(self):
        rad_dict = self._mgo_valrad_evaluator.radii
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
        self._mgo_uc = Structure(mgo_latt, mgo_specie, mgo_frac_cord, True, True)
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
        #Match uniq_sites iwth uniq_def_sites
        #self.assertTrue(len(uniq_def_sites) == 2, "Vacancy init failed")
        #mgo_spg = Spacegroup(int_number=225)
        #self.assertTrue(mgo_spg.are_symmetrically_equivalent(uniq_sites,
        #                uniq_def_sites),  "Vacancy init failed")

    def test_get_defectsite_index(self):
        for i in range(self._mgo_vac.defectsite_count()):
            self.assertTrue(self._mgo_vac.get_defectsite_structure_index(i) < 
                    len(self._mgo_uc.sites), "Defect site index beyond range")

    def test_gt_defectsite_coordination_number(self):
        for i in range(self._mgo_vac.defectsite_count()):
            self.assertTrue(round(self._mgo_vac.get_defectsite_coordination_number(
                i))==6.0, "Wrong coordination number")

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
        mgo_valrad_eval = ValenceIonicRadiusEvaluator(self.mgo_uc)
        val = mgo_valrad_eval.valences
        rad = mgo_valrad_eval.radii
        self.mgo_vac = Vacancy(self.mgo_uc, val, rad)
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
        self._mgo_uc = Structure(mgo_latt, mgo_specie, mgo_frac_cord, True, True)
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
        struct_el = self._mgo_uc.composition.elements
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


class InterstitialAnalyzerTest(unittest.TestCase):
    def setUp(self):
        mgo_latt = [[4.212, 0, 0], [0, 4.212, 0], [0, 0, 4.212]]
        mgo_specie = ["Mg"]*4 +  ["O"]*4
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
            print i, site_coords, site_radius, ife
            self.assertIsInstance(ife, float)

    def test_get_norelaxedenergy(self):
        for i in range(self.mgo_inter.defectsite_count()):
            ife = self.mgo_ia.get_energy(i, False)
            site_coords = self.mgo_inter.get_defectsite(i).coords
            site_radius = self.mgo_inter.get_radius(i)
            print i, site_coords, site_radius, ife
            self.assertIsInstance(ife, float)

    def test_get_percentage_volume_change(self):
        for i in range(self.mgo_inter.defectsite_count()):
            del_vol = self.mgo_ia.get_percentage_volume_change(i)
            print i, del_vol

    def test_get_percentage_lattice_parameter_change(self):
        for i in range(self.mgo_inter.defectsite_count()):
            del_lat = self.mgo_ia.get_percentage_lattice_parameter_change(i)
            print i, del_lat

    def test_get_percentage_bond_distance_change(self):
        for i in range(self.mgo_inter.defectsite_count()):
            del_bd = self.mgo_ia.get_percentage_bond_distance_change(i)
            print i, del_bd

    def test_relaxed_structure_match(self):
        for i in range(self.mgo_inter.defectsite_count()):
            for j in range(self.mgo_inter.defectsite_count()):
                match = self.mgo_ia.relaxed_structure_match(i,j)
                print i, j, match
                if i == j:
                    self.assertTrue(match)


class InterstitialStructureRelaxerTest(unittest.TestCase):
    def setUp(self):
        mgo_latt = [[4.212, 0, 0], [0, 4.212, 0], [0, 0, 4.212]]
        mgo_specie = ["Mg"]*4 +  ["O"]*4
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
                match = self.isr.relaxed_structure_match(i,j)
                #print i, j, match
                if i == j:
                    self.assertTrue(match)

    def test_relaxed_energy_match(self):
        for i in range(self.mgo_inter.defectsite_count()):
            for j in range(self.mgo_inter.defectsite_count()):
                match = self.isr.relaxed_energy_match(i,j)
                #print i, j, match
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


class RelaxedInsterstitialTest(unittest.TestCase):
    def setUp(self):
        mgo_latt = [[4.212, 0, 0], [0, 4.212, 0], [0, 0, 4.212]]
        mgo_specie = ["Mg"]*4 +  ["O"]*4
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
            print "ife", ife

    def test_get_percentage_volume_change(self):
        for i in range(self.mgo_inter.defectsite_count()):
            del_vol = self.ri.get_percentage_volume_change(i)
            self.assertIsInstance(del_vol, float)
            print "del_vol", del_vol

    def test_get_percentage_lattice_parameter_change(self):
        for i in range(self.mgo_inter.defectsite_count()):
            del_lat = self.ri.get_percentage_lattice_parameter_change(i)
            self.assertNotEqual(del_lat['a'], 0)
            self.assertNotEqual(del_lat['b'], 0)
            self.assertNotEqual(del_lat['c'], 0)
            print "del_lat", del_lat

    def test_get_percentage_bond_distance_change(self):
        for i in range(self.mgo_inter.defectsite_count()):
            del_bd = self.ri.get_percentage_bond_distance_change(i)
            #self.assertIsInstance(del_bd, float)
            #print del_bd


if __name__ == "__main__":
    unittest.main()
    #suite = unittest.TestLoader().loadTestsFromTestCase(ValenceIonicRadiusEvaluatorTest)
    #suite = unittest.TestLoader().loadTestsFromTestCase(InterstitialTest)
    #suite = unittest.TestLoader().loadTestsFromTestCase(VacancyTest)
    #suite = unittest.TestLoader().loadTestsFromTestCase(VacancyFormationEnergyTest)
    #suite = unittest.TestLoader().loadTestsFromTestCase(InterstitialAnalyzerTest)
    #unittest.TextTestRunner(verbosity=3).run(suite)
