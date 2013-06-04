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
        Generate MgO rocksalt structure and compare the vacancies
        """
        mgo_latt = [[4.212, 0, 0], [0, 4.212, 0], [0, 0, 4.212]]
        mgo_specie = ["Mg2+"]*4 +  ["O2-"]*4
        mgo_frac_cord = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
                         [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5], [0.5, 0.5, 0.5]]
        self.mgo_uc = Structure(mgo_latt, mgo_specie, mgo_frac_cord, True, True)
        #vac1 = mgo_uc.remove(0)
        #vac2 = mgo_uc.remove(4)
        #mgo_symmfinder = SymmetryFinder(mgo_uc)
    
    def test_init(self):
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
