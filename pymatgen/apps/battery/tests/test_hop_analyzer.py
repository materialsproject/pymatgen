# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
Created on Jan 25, 2012
"""

__author__ = "Jimmy Shen"
__copyright__ = "Copyright 2019, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Jimmy Shen"
__email__ = "jmmshn@lbl.gov"
__date__ = "April 1, 2019"

import unittest
import os
import json

from monty.serialization import loadfn
from pymatgen.io.vasp import Chgcar
from pymatgen.apps.battery.hop_analyzer import MigrationPathAnalyzer
from pymatgen.analysis.structure_matcher import StructureMatcher, ElementComparator


test_dir = os.path.join(
    os.path.dirname(__file__), "..", "..", "..", "..", 'test_files')


class HopTest(unittest.TestCase):
    def setUp(self):
        self.test_ents_MOF = loadfn(test_dir + '/Mn6O5F7_cat_migration.json')
        self.aeccar_MOF = Chgcar.from_file(test_dir + '/AECCAR_Mn6O5F7.vasp')
        self.mpa_MOF = MigrationPathAnalyzer(
            base_entry=self.test_ents_MOF['ent_base'],
            single_cat_entries=self.test_ents_MOF['one_cation'],
            base_aeccar=self.aeccar_MOF)
        self.rough_sm = StructureMatcher(
            comparator=ElementComparator(),
            primitive_cell=False,
            ltol=0.5,
            stol=0.5,
            angle_tol=7)

    def test_get_all_sym_sites(self):
        """
        Once we apply all the symmetry operations to get the new positions the structures
        with the Li in different positions should still match
        (using sloppier tolerences)
        """
        for ent in self.mpa_MOF.translated_single_cat_entries:
            li_sites = self.mpa_MOF.get_all_sym_sites(ent).sites
            s0 = self.mpa_MOF.base_entry.structure.copy()
            s0.insert(0, 'Li', li_sites[0].frac_coords)
            for isite in li_sites[1:]:
                s1 = self.mpa_MOF.base_entry.structure.copy()
                s1.insert(0, 'Li', isite.frac_coords)
                self.assertTrue(self.rough_sm.fit(s0, s1))

    def test_get_all_sym_sites(self):
        """
        Once we apply all the symmetry operations to get the new positions the structures
        with the Li in different positions should still match
        (using sloppier tolerences)
        """
        for ent in self.mpa_MOF.translated_single_cat_entries:
            li_sites = self.mpa_MOF.get_all_sym_sites(ent).sites
            s0 = self.mpa_MOF.base_entry.structure.copy()
            s0.insert(0, 'Li', li_sites[0].frac_coords)
            for isite in li_sites[1:]:
                s1 = self.mpa_MOF.base_entry.structure.copy()
                s1.insert(0, 'Li', isite.frac_coords)
                self.assertTrue(self.rough_sm.fit(s0, s1))

if __name__ == '__main__':
    unittest.main()
