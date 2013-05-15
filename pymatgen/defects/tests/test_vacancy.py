#!/usr/bin/python

import unittest
from pymatgen.defects.point_defects import Vacancy
from pymatgen.core.structure import  Structure
from pymatgen.core.lattice import Lattice
from pymatgen.symmetry.finder import SymmetryFinder 
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
        mgo_uc = Structure(mgo_latt, mgo_specie, mgo_frac_cord, True, True)
        vac1 = mgo_uc.remove(0)
        vac2 = mgo_uc.remove(4)
        mgo_symmfinder = SymmetryFinder(mgo_uc)
        
        


