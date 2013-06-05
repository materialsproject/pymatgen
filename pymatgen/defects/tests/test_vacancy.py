#!/usr/bin/python

import unittest
from pymatgen.defects.point_defects import Vacancy
from pymatgen.core.structure import  Structure
from pymatgen.core.lattice import Lattice
from pymatgen.symmetry.finder import SymmetryFinder 
from pymatgen.symmetry.spacegroup import Spacegroup
import numpy as np

class VacancyTest(unittest.TestCase):
    def setUp(self):
        """
        Setup MgO rocksalt structure for testing Vacancy
        """
        mgo_latt = [[4.212, 0, 0], [0, 4.212, 0], [0, 0, 4.212]]
        mgo_specie = ["Mg2+"]*4 +  ["O2-"]*4
        mgo_frac_cord = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
                         [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5], [0.5, 0.5, 0.5]]
        self.mgo_uc = Structure(mgo_latt, mgo_specie, mgo_frac_cord, True, True)
    
    def test_enumerate_uniq_defectsites(self):
        """
        The vacancy sites should be the lattice sites. 
        And there should be only two unique vacancy sites for MgO.
        """
        mgo_uc1 = self.mgo_uc.copy()
        uniq_sites = []
        uniq_sites.append(mgo_uc1.sites[3])
        uniq_sites.append(mgo_uc1.sites[7])
        mgo_vac = Vacancy(mgo_uc1)
        uniq_def_sites = mgo_vac.enumerate_uniq_defectsites()
        self.assertTrue(len(uniq_def_sites) == 2, "Vacancy init failed")
        #mgo_spg = Spacegroup(int_number=225)
        #self.assertTrue(mgo_spg.are_symmetrically_equivalent(uniq_sites,
        #                uniq_def_sites),  "Vacancy init failed")
        
class InterstitialTest(unittest.TestCase):
    def setUp(self):
        """
        Setup MgO rocksalt structure for testing Interstitial
        """
        mgo_latt = [[4.212, 0, 0], [0, 4.212, 0], [0, 0, 4.212]]
        mgo_specie = ["Mg2+"]*4 +  ["O2-"]*4
        mgo_frac_cord = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
                         [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5], [0.5, 0.5, 0.5]]
        self.mgo_uc = Structure(mgo_latt, mgo_specie, mgo_frac_cord, True, True)
    
    def test_enumerate_uniq_defectsites(self):
        """
        The interstitial sites should be within the lattice
        """
        mgo_uc1 = self.mgo_uc.copy()
        uniq_sites = []
        uniq_sites.append(mgo_uc1.sites[3])
        uniq_sites.append(mgo_uc1.sites[7])
        mgo_vac = Vacancy(mgo_uc1)
        uniq_def_sites = mgo_vac.enumerate_uniq_defectsites()
        self.assertTrue(len(uniq_def_sites) == 2, "Vacancy init failed")
        #mgo_spg = Spacegroup(int_number=225)
        #self.assertTrue(mgo_spg.are_symmetrically_equivalent(uniq_sites,
        #                uniq_def_sites),  "Vacancy init failed")
        


if __name__ == "__main__":
    unittest.main()
