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
            base_struct_entry=self.test_ents_MOF['ent_base'],
            single_cat_entries=self.test_ents_MOF['one_cation'],
            base_aeccar=self.aeccar_MOF)
        self.rough_sm = StructureMatcher(
            comparator=ElementComparator(),
            primitive_cell=False,
            ltol=0.5,
            stol=0.5,
            angle_tol=7)
        self.mpa_MOF.get_full_sites()

    def test_get_all_sym_sites(self):
        """
        Once we apply all the symmetry operations to get the new positions the structures
        with the Li in different positions should still match
        (using sloppier tolerences)
        """
        for ent in self.mpa_MOF.translated_single_cat_entries:
            li_sites = self.mpa_MOF.get_all_sym_sites(ent).sites
            s0 = self.mpa_MOF.base_struct_entry.structure.copy()
            s0.insert(0, 'Li', li_sites[0].frac_coords)
            for isite in li_sites[1:]:
                s1 = self.mpa_MOF.base_struct_entry.structure.copy()
                s1.insert(0, 'Li', isite.frac_coords)
                self.assertTrue(self.rough_sm.fit(s0, s1))

    def test_get_full_sites(self):
        """
        For Mn6O5F7 there should be 8 symmetry inequivalent final relaxed positions
        """
        self.assertEqual(len(self.mpa_MOF.full_sites), 8)

    def test_integration(self):
        """
        Sanity check: for a long enough diagonaly hop, if we turn the radius of the tube way up, it should cover the entire unit cell
        """
        self.mpa_MOF.populate_unique_edges_df()
        self.mpa_MOF._setup_grids()
        self.mpa_MOF._tube_radius = 1000
        total_chg_per_vol = self.mpa_MOF.base_aeccar.data['total'].sum(
        ) / self.mpa_MOF.base_aeccar.ngridpts / self.mpa_MOF.base_aeccar.structure.volume
        self.assertAlmostEqual(
            self.mpa_MOF._get_chg_between_sites_tube(0), total_chg_per_vol)

    def test_unique_edges(self):
        """
        check that the edgelist and unique_edges has the right number of entries
        """
        self.mpa_MOF.populate_unique_edges_df()

        self.assertEqual(len(self.mpa_MOF.edgelist), 87)
        self.assertEqual(len(self.mpa_MOF.unique_edges), 16)

        self.mpa_MOF.get_chg_values_for_unique_hops()
        self.assertAlmostEqual(
            self.mpa_MOF.unique_edges['chg_total_tube'].values.min(),
            0.00033953800319920594)
        self.assertAlmostEqual(
            self.mpa_MOF.unique_edges['chg_total_tube'].values.max(),
            0.08161937859517454)

    def test_from_aeccar(self):
        self.mpa_MOF_chgonly = MigrationPathAnalyzer.from_aeccar(base_aeccar=self.aeccar_MOF)
        self.assertEqual(len(self.mpa_MOF_chgonly.full_sites), 16)

if __name__ == '__main__':
    unittest.main()
