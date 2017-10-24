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
        with open(os.path.join(get_path(""), 'ucell_entries.txt')) as ucell_entries:
            ucell_entries = json.loads(ucell_entries.read())
        self.ucell_entries = ucell_entries

        Cu_ucell_entry = ComputedStructureEntry.from_dict(self.ucell_entries["Cu"])
        self.Cu_analyzer = SurfaceEnergyCalculator(Cu_ucell_entry)

        # entry_dict for the adsorption case, O adsorption on Ni, Rh and Pt
        metals_O_entry_dict = {"Ni": {(1, 1, 1): {}, (1, 0, 0): {}},
                               "Pt": {(1, 1, 1): {}},
                               "Rh": {(1, 0, 0): {}}
                               }

        with open(os.path.join(get_path(""), "csentries_slabs.json")) as entries:
            entries = json.loads(entries.read())
        for k in entries.keys():
            entry = ComputedStructureEntry.from_dict(entries[k])
            for el in metals_O_entry_dict.keys():
                if el in k:
                    if "111" in k:
                        metals_O_entry_dict[el][(1, 1, 1)][entry] = []
                    if "110" in k:
                        metals_O_entry_dict[el][(1, 1, 0)][entry] = []
                    if "100" in k:
                        metals_O_entry_dict[el][(1, 0, 0)][entry] = []

        with open(os.path.join(get_path(""), "csentries_o_ads.json")) as entries:
            entries = json.loads(entries.read())
        for k in entries.keys():
            entry = ComputedStructureEntry.from_dict(entries[k])
            for el in metals_O_entry_dict.keys():
                if el in k:
                    if "111" in k:
                        clean = list(metals_O_entry_dict[el][(1, 1, 1)].keys())[0]
                        metals_O_entry_dict[el][(1, 1, 1)][clean] = [entry]
                    if "110" in k:
                        clean = list(metals_O_entry_dict[el][(1, 1, 0)].keys())[0]
                        metals_O_entry_dict[el][(1, 1, 0)][clean] = [entry]
                    if "100" in k:
                        clean = list(metals_O_entry_dict[el][(1, 0, 0)].keys())[0]
                        metals_O_entry_dict[el][(1, 0, 0)][clean] = [entry]

        self.metals_O_entry_dict = metals_O_entry_dict

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

        # Test cases for getting adsorption related quantities for a 1/4
        # monolalyer adsorption of O on the low MMI surfaces of Pt, Ni and Rh

        with open(os.path.join(get_path(""), 'isolated_O_entry.txt')) as isolated_O_entry:
            isolated_O_entry = json.loads(isolated_O_entry.read())
        O = ComputedStructureEntry.from_dict(isolated_O_entry)

        for el in self.metals_O_entry_dict.keys():
            el_ucell = ComputedStructureEntry.from_dict(self.ucell_entries[el])
            se_calc = SurfaceEnergyCalculator(el_ucell, adsorbate_entry=O)
            for hkl in self.metals_O_entry_dict[el].keys():
                for clean in self.metals_O_entry_dict[el][hkl]:
                    for ads in self.metals_O_entry_dict[el][hkl][clean]:
                        # Determine the correct number of monolayers.
                        ml = se_calc.get_monolayer(ads, clean)
                        self.assertEqual(int(round(ml)), 4)
                        # Determine the correct number of adsorbates
                        Nads = se_calc.Nads_in_slab(ads)
                        self.assertEqual(Nads, 1)
                        # Determine the correct number of surfaces with an adsorbate
                        Nsurf = se_calc.Nsurfs_ads_in_slab(ads)
                        self.assertEqual(Nsurf, 1)
                        # Determine the correct binding energy
                        gbind = (ads.energy - ml*clean.energy)/Nads - O.energy_per_atom
                        self.assertEqual(gbind, se_calc.gibbs_binding_energy(ads, clean))
                        # Determine the correction Gibbs adsorption energy
                        eads = Nads * gbind
                        self.assertEqual(eads, se_calc.gibbs_binding_energy(ads, clean, eads=True))

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
        with open(os.path.join(get_path(""), 'ucell_entries.txt')) as ucell_entries:
            ucell_entries = json.loads(ucell_entries.read())

        Cu_ucell_entry = ComputedStructureEntry.from_dict(ucell_entries["Cu"])
        calculator = SurfaceEnergyCalculator(Cu_ucell_entry)
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