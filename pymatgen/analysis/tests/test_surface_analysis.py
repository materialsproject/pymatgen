import unittest
import os
import random
import glob
import json

import numpy as np

from pymatgen.analysis.surface_analysis import SurfaceEnergyCalculator, \
    SurfaceEnergyPlotter, Composition
from pymatgen.util.testing import PymatgenTest
from pymatgen.entries.computed_entries import ComputedStructureEntry

__author__ = "Richard Tran"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Richard Tran"
__email__ = "rit001@eng.ucsd.edu"
__date__ = "Aug 24, 2017"


def get_path(path_str):
    cwd = os.path.abspath(os.path.dirname(__file__))
    path = os.path.join(cwd, "..", "..", "..", "test_files",
                        "surface_tests", path_str)
    return path

class SurfaceEnergyCalculatorTest(PymatgenTest):

    def setUp(self):

        self.entry_dict = get_entry_dict(os.path.join(get_path(""),
                                                      "Cu_entries.txt"))
        with open(os.path.join(get_path(""), 'Cu_ucell_entry.txt')) as ucell_entry:
            ucell_entry = json.loads(ucell_entry.read())
        ucell_entry = ComputedStructureEntry.from_dict(ucell_entry)
        self.Cu_analyzer = SurfaceEnergyCalculator(ucell_entry)

    def test_gamma_calculator(self):

        # Test case for clean Cu surface

        # make sure we've loaded all our files correctly
        self.assertEqual(len(self.entry_dict.keys()), 13)
        all_se = []
        for hkl, entries in self.entry_dict.items():
            u1, u2 = -1, 0
            se1 = self.Cu_analyzer.calculate_gamma_at_u(list(entries.keys())[0], u_ref=u1)
            se2 = self.Cu_analyzer.calculate_gamma_at_u(list(entries.keys())[0], u_ref=u2)
            # For a stoichiometric system, we expect surface
            # energy to be independent of chemical potential
            self.assertEqual(se1, se2)
            all_se.append(se1)

        clean111_entry = list(self.entry_dict[(1,1,1)].keys())[0]
        # The (111) facet should be the most stable
        self.assertEqual(min(all_se),
                         self.Cu_analyzer.calculate_gamma_at_u( \
                             clean111_entry, 0))

        # Get the coefficients for surfacce energy. For a clean
        # stoichiometric system, the third term should be the surface energy
        # and the adsorption and nonstoichiometric terms should be 0
        b1, b2, b3 = self.Cu_analyzer.surface_energy_coefficients(clean111_entry)
        self.assertEqual(b1, b2)
        self.assertEqual(b1, 0)
        self.assertEqual(b3, min(all_se))

    def test_adsorption_quantities(self):

        # monolayer
        # gibbs_binding_energy
        # Nads_in_slab
        # Nsurfs_ads_in_slab

        print()

    def test_solve_2_linear_eqns(self):

        # For clean nonstoichiometric system, the two equations should
        # be parallel because the surface energy is a constant
        clean111_entry = list(self.entry_dict[(1, 1, 1)].keys())[0]
        clean100_entry = list(self.entry_dict[(1, 0, 0)].keys())[0]
        c111 = self.Cu_analyzer.surface_energy_coefficients(clean111_entry)
        c100 = self.Cu_analyzer.surface_energy_coefficients(clean100_entry)
        soln = self.Cu_analyzer.solve_2_linear_eqns(c111, c100)
        self.assertFalse(soln[0])
        self.assertEqual(soln[1], soln[0])

class SurfaceEnergyPlotterTest(PymatgenTest):

    def setUp(self):

        entry_dict = get_entry_dict(os.path.join(get_path(""),
                                                 "Cu_entries.txt"))
        with open(os.path.join(get_path(""), 'Cu_ucell_entry.txt')) as ucell_entry:
            ucell_entry = json.loads(ucell_entry.read())
        ucell_entry = ComputedStructureEntry.from_dict(ucell_entry)
        calculator = SurfaceEnergyCalculator(ucell_entry)
        self.Cu_analyzer = SurfaceEnergyPlotter(entry_dict, calculator, [-1, 0])

    def test_max_adsorption_chempot_range(self):

        print()

    def test_chempot_range_adsorption(self):

        print()

    def test_wulff_shape_from_chempot(self):

        # Test if it generates a Wulff shape, test if
        # all the facets for Cu wulff shape are inside.
        Cu_wulff = self.Cu_analyzer.wulff_shape_from_chempot()
        area_frac_dict = Cu_wulff.area_fraction_dict
        facets_hkl = [(1,1,1), (3,3,1), (3,1,0), (1,0,0),
                      (3,1,1), (2,1,0), (2,2,1)]
        for hkl in area_frac_dict.keys():
            if hkl in facets_hkl:
                self.assertNotEqual(area_frac_dict[hkl], 0)
            else:
                self.assertEqual(area_frac_dict[hkl], 0)

    def test_return_stable_slab_entry_at_u(self):

        print()

    def test_area_frac_vs_chempot_plot(self):

        print()

    def test_chempot_vs_gamma_plot_one(self):

        print()

    def test_chempot_vs_gamma_clean(self):

        print()

    def test_chempot_vs_gamma_facet(self):

        print()

    def test_chempot_plot_addons(self):

        print()

    def test_create_slab_label(self):

        print()

    def test_color_palette_dict(self):

        print()

    def test_stable_u_range_dict(self):

        print()

    def test_get_clean_ads_entry_pair(self):

        print()



if __name__ == "__main__":
    unittest.main()

def get_entry_dict(filename):
    # helper to generate an entry_dict

    entry_dict = {}
    with open(filename) as entries:
        entries = json.loads(entries.read())
    for k in entries.keys():
        n = k[25:]
        miller_index = []
        for i, s in enumerate(n):
            if i > 2:
                break
            miller_index.append(int(s))
        hkl = tuple(miller_index)
        if hkl not in entry_dict.keys():
            entry_dict[hkl] = {}
        entry_dict[hkl][ComputedStructureEntry.from_dict(entries[k])] = []

    return entry_dict