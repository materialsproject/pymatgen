import unittest
import os
import random
import glob
import json

import numpy as np
from sympy import Number

from pymatgen.analysis.surface_analysis import SurfaceEnergyCalculator, \
    Composition, SlabEntry#, SurfaceEnergyPlotter
from pymatgen.util.testing import PymatgenTest
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatgen import Structure, Lattice

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

        self.Cu_entry_dict = get_entry_dict(os.path.join(get_path(""),
                                                         "Cu_entries.txt"))
        self.assertEqual(len(self.Cu_entry_dict.keys()), 13)

        with open(os.path.join(get_path(""), 'ucell_entries.txt')) as ucell_entries:
            ucell_entries = json.loads(ucell_entries.read())
        self.ucell_entries = ucell_entries
        Cu_ucell_entry = ComputedStructureEntry.from_dict(self.ucell_entries["Cu"])
        self.Cu_calc = SurfaceEnergyCalculator(Cu_ucell_entry)
        self.metals_O_entry_dict = load_O_adsorption()
        with open(os.path.join(get_path(""),
                               'isolated_O_entry.txt')) as isolated_O_entry:
            isolated_O_entry = json.loads(isolated_O_entry.read())
        self.O = ComputedStructureEntry.from_dict(isolated_O_entry)

        # Load dummy MgO slab entries
        MgO_ucell_entry = ComputedStructureEntry.from_dict(self.ucell_entries["MgO"])
        Mg_ucell_entry = ComputedStructureEntry.from_dict(self.ucell_entries["Mg"])
        self.MgO_calc = SurfaceEnergyCalculator(MgO_ucell_entry,
                                                ref_entries=[Mg_ucell_entry])
        self.MgO_slab_entry_dict = get_entry_dict(os.path.join(get_path(""),
                                                               "MgO_slab_entries.txt"))

    def test_surface_energy_coefficients(self):
        # For a nonstoichiometric case, the cheimcal potentials do not
        # cancel out, they serve as a reservoir for any missing element
        for slab_entry in self.MgO_slab_entry_dict[(1,1,1)].keys():
            coeffs = self.MgO_calc.surface_energy_coefficients(slab_entry)
            coeff_keys = [str(k) for k in coeffs.keys()]
            self.assertEqual(tuple(coeff_keys), ("1", "Mg"))

        # For the case of a clean, stoichiometric slab, the the key to
        # the coefficient should be 1 (i.e. surface energy is a constant).
        for hkl in self.Cu_entry_dict.keys():
            slab_entry = list(self.Cu_entry_dict[hkl].keys())[0]
            coeffs = self.Cu_calc.surface_energy_coefficients(slab_entry)
            self.assertEqual(len(coeffs.keys()), 1)
            self.assertEqual(list(coeffs.keys())[0], Number(1))

    def test_gamma_calculator(self):

        # Test case for clean Cu surface

        all_se = []
        u1, u2 = -1, 0
        for hkl in self.Cu_entry_dict.keys():
            entry = list(self.Cu_entry_dict[hkl].keys())[0]
            se1 = self.Cu_calc.calculate_gamma(entry, u_default=u1)
            se2 = self.Cu_calc.calculate_gamma(entry, u_default=u2)
            # For a stoichiometric system, we expect surface
            # energy to be independent of chemical potential
            self.assertEqual(se1, se2)
            all_se.append(se1)

        # The (111) facet should be the most stable
        clean111_entry = list(self.Cu_entry_dict[(1,1,1)].keys())[0]
        se_Cu111 = self.Cu_calc.calculate_gamma(clean111_entry)
        self.assertEqual(min(all_se), se_Cu111)

        # Check surface eenrgy properly calculated
        Cu_entry = ComputedStructureEntry.from_dict(self.ucell_entries["Cu"])
        slab = clean111_entry.structure
        sa = clean111_entry.surface_area
        se111 = (clean111_entry.energy - len(slab)*Cu_entry.energy_per_atom)/(2*sa)
        self.assertEqual(se111, se_Cu111)

    def test_adsorption_quantities(self):

        # Test cases for getting adsorption related quantities for a 1/4
        # monolalyer adsorption of O on the low MMI surfaces of Pt, Ni and Rh

        for el in self.metals_O_entry_dict.keys():
            el_ucell = ComputedStructureEntry.from_dict(self.ucell_entries[el])
            se_calc = SurfaceEnergyCalculator(el_ucell, adsorbate_entries=[self.O])
            for hkl in self.metals_O_entry_dict[el].keys():
                for clean in self.metals_O_entry_dict[el][hkl]:
                    for ads in self.metals_O_entry_dict[el][hkl][clean]:
                        # Determine the correct number of monolayers.
                        ml = se_calc.get_unit_primitive_area(ads, clean)
                        self.assertEqual(int(round(ml)), 4)
                        # Determine the correct number of adsorbates
                        Nads = se_calc.Nads_in_slab(ads)
                        self.assertEqual(Nads, 1)
                        # Determine the correct number of surfaces with an adsorbate
                        Nsurf = se_calc.Nsurfs_ads_in_slab(ads)
                        self.assertEqual(Nsurf, 1)
                        # Determine the correct binding energy
                        gbind = (ads.energy - ml*clean.energy)/Nads - self.O.energy_per_atom
                        self.assertEqual(gbind, se_calc.gibbs_binding_energy(ads, clean))
                        # Determine the correction Gibbs adsorption energy
                        eads = Nads * gbind
                        self.assertEqual(eads, se_calc.gibbs_binding_energy(ads, clean, eads=True))

    def test_get_surface_equilibrium(self):

        # For clean stoichiometric system, the two equations should
        # be parallel because the surface energy is a constant. Then
        # get_surface_equilibrium should return None
        clean111_entry = list(self.Cu_entry_dict[(1, 1, 1)].keys())[0]
        clean100_entry = list(self.Cu_entry_dict[(1, 0, 0)].keys())[0]
        soln = self.Cu_calc.get_surface_equilibrium([clean111_entry, clean100_entry])
        self.assertFalse(soln)

        # For adsorbed system, we should find one intercept
        Pt_entries = self.metals_O_entry_dict["Pt"]
        clean = list(Pt_entries[(1,1,1)].keys())[0]
        ads = Pt_entries[(1,1,1)][clean][0]
        Pt_ucell_entry = ComputedStructureEntry.from_dict(self.ucell_entries["Pt"])
        Pt_analyzer = SurfaceEnergyCalculator(Pt_ucell_entry, adsorbate_entries=[self.O])
        soln = Pt_analyzer.get_surface_equilibrium([clean, ads])

        self.assertNotEqual(list(soln.values())[0], list(soln.values())[1])

        # Check if the number of parameters for adsorption are correct
        self.assertEqual(("O", "gamma"), tuple(soln.keys()))
        # Adsorbed systems have a b2=(-1*Nads) / (Nsurfs * Aads)
        self.assertEqual(soln["O"], -1/ads.surface_area)
        # Check the intercept
        gbind = Pt_analyzer.gibbs_binding_energy(ads, clean)
        b3 = cclean[2] + gbind*(1/(sa))
        self.assertEqual(cads[2], b3)


# class SurfaceEnergyPlotterTest(PymatgenTest):
#
#     def setUp(self):
#
#         entry_dict = get_entry_dict(os.path.join(get_path(""),
#                                                  "Cu_entries.txt"))
#         with open(os.path.join(get_path(""), 'ucell_entries.txt')) as ucell_entries:
#             ucell_entries = json.loads(ucell_entries.read())
#
#         Cu_ucell_entry = ComputedStructureEntry.from_dict(ucell_entries["Cu"])
#         calculator = SurfaceEnergyCalculator(Cu_ucell_entry)
#         self.Cu_analyzer = SurfaceEnergyPlotter(entry_dict, calculator, [-1, 0])
#
#         with open(os.path.join(get_path(""), 'isolated_O_entry.txt')) as isolated_O_entry:
#             isolated_O_entry = json.loads(isolated_O_entry.read())
#         self.O = ComputedStructureEntry.from_dict(isolated_O_entry)
#         self.metals_O_entry_dict = load_O_adsorption()
#         ucell_entry = ComputedStructureEntry.from_dict(ucell_entries["Pt"])
#         self.Pt_calculator = SurfaceEnergyCalculator(ucell_entry, adsorbate_entry=self.O)
#         self.Pt_analyzer = SurfaceEnergyPlotter(self.metals_O_entry_dict["Pt"],
#                                                 self.Pt_calculator, [-1, 0])
#         ucell_entry = ComputedStructureEntry.from_dict(ucell_entries["Ni"])
#         self.Ni_calculator = SurfaceEnergyCalculator(ucell_entry, adsorbate_entry=self.O)
#         self.Ni_analyzer = SurfaceEnergyPlotter(self.metals_O_entry_dict["Ni"],
#                                                 self.Ni_calculator, [-1, 0])
#         ucell_entry = ComputedStructureEntry.from_dict(ucell_entries["Rh"])
#         self.Rh_calculator = SurfaceEnergyCalculator(ucell_entry, adsorbate_entry=self.O)
#         self.Rh_analyzer = SurfaceEnergyPlotter(self.metals_O_entry_dict["Rh"],
#                                                 self.Rh_calculator, [-1, 0])
#         self.Oads_analyzer_dict = {"Pt": self.Pt_analyzer,
#                                    "Ni": self.Ni_analyzer,
#                                    "Rh": self.Rh_analyzer}
#
#     def test_max_adsorption_chempot_range(self):
#
#         # Tests if the range of chemical
#         # potential has been reasonably chosen
#
#         for el in self.Oads_analyzer_dict.keys():
#             analyzer = self.Oads_analyzer_dict[el]
#             # Is max chempot 0 or greater
#             chempot = analyzer.max_adsorption_chempot_range(0)
#             self.assertLessEqual(chempot[0], 0)
#
#             # Is the min chempot going to give us the clean gamma
#             if (1,1,1) in self.metals_O_entry_dict[el].keys():
#                 se = analyzer.return_stable_slab_entry_at_u((1,1,1,),
#                                                             u_ads=chempot[0])[1]
#                 clean = list(self.metals_O_entry_dict[el][(1,1,1)].keys())[0]
#                 se_clean = analyzer.se_calculator.surface_energy_coefficients(clean)[2]
#                 self.assertEqual(se, se_clean)
#                 se = analyzer.return_stable_slab_entry_at_u((1,1,1,),
#                                                             u_ads=chempot[1])[1]
#                 self.assertGreater(se, 0)
#
#     def test_wulff_shape_from_chempot(self):
#
#         # Test if it generates a Wulff shape, test if
#         # all the facets for Cu wulff shape are inside.
#         Cu_wulff = self.Cu_analyzer.wulff_shape_from_chempot()
#         area_frac_dict = Cu_wulff.area_fraction_dict
#         facets_hkl = [(1,1,1), (3,3,1), (3,1,0), (1,0,0),
#                       (3,1,1), (2,1,0), (2,2,1)]
#         for hkl in area_frac_dict.keys():
#             if hkl in facets_hkl:
#                 self.assertNotEqual(area_frac_dict[hkl], 0)
#             else:
#                 self.assertEqual(area_frac_dict[hkl], 0)
#
#         for el in self.Oads_analyzer_dict.keys():
#             # Test WulffShape for adsorbed surfaces
#             analyzer = self.Oads_analyzer_dict[el]
#             chempot = analyzer.max_adsorption_chempot_range(0)
#             wulff = analyzer.wulff_shape_from_chempot(u_ads=np.mean(chempot))
#             se = wulff.weighted_surface_energy
#
#     def test_create_slab_label(self):
#
#         for el in self.metals_O_entry_dict.keys():
#             analyzer = self.Oads_analyzer_dict[el]
#             for hkl in self.metals_O_entry_dict[el].keys():
#                 # Test WulffShape for adsorbed surfaces
#                 for clean in self.metals_O_entry_dict[el][hkl]:
#                     label = analyzer.create_slab_label(clean, miller_index=hkl)
#
#                     self.assertEqual(str(hkl), label)
#
#                     for ads in self.metals_O_entry_dict[el][hkl][clean]:
#                         label = analyzer.create_slab_label(clean, miller_index=hkl,
#                                                            ads_entry=ads)
#                         self.assertEqual(label, str(hkl)+"+O")
#
#     def test_color_palette_dict(self):
#
#         for el in self.metals_O_entry_dict.keys():
#             analyzer = self.Oads_analyzer_dict[el]
#             color_dict = analyzer.color_palette_dict()
#
#     def test_stable_u_range_dict(self):
#
#         # Test if it generates a Wulff shape, test if
#         # all the facets for Cu wulff shape are inside.
#         stable_u_range_dict = self.Cu_analyzer.stable_u_range_dict()
#         for entry in stable_u_range_dict.keys():
#             urange = stable_u_range_dict[entry]
#             self.assertEqual(-1, urange[0])
#             self.assertEqual(0, urange[1])
#
#         for el in self.Oads_analyzer_dict.keys():
#             # Test WulffShape for adsorbed surfaces
#             analyzer = self.Oads_analyzer_dict[el]
#             urange = analyzer.stable_u_range_dict(clean_only=False)
#             for entry in urange.keys():
#                 u = urange[entry]
#                 self.assertLess(u[0], u[-1])
#
#     def test_get_clean_ads_entry_pair(self):
#
#         for el in self.metals_O_entry_dict.keys():
#             # Test WulffShape for adsorbed surfaces
#             analyzer = self.Oads_analyzer_dict[el]
#             for hkl in self.metals_O_entry_dict[el].keys():
#                 for clean in self.metals_O_entry_dict[el][hkl].keys():
#                     for ads in self.metals_O_entry_dict[el][hkl][clean]:
#                         pair = analyzer.get_clean_ads_entry_pair(ads)
#                         c = analyzer.se_calculator.surface_energy_coefficients(
#                             pair[0], ads_slab_entry=pair[1])
#                         c1 = analyzer.se_calculator.surface_energy_coefficients(
#                             clean, ads_slab_entry=ads)
#                         self.assertEqual(tuple(c), tuple(c1))
#
#     # def test_monolayer_vs_BE(self):
#     #     for el in self.Oads_analyzer_dict.keys():
#     #         # Test WulffShape for adsorbed surfaces
#     #         analyzer = self.Oads_analyzer_dict[el]
#     #         plt = analyzer.monolayer_vs_BE()
#     #
#     # def test_area_frac_vs_chempot_plot(self):
#     #
#     #     for el in self.Oads_analyzer_dict.keys():
#     #         # Test WulffShape for adsorbed surfaces
#     #         analyzer = self.Oads_analyzer_dict[el]
#     #         plt = analyzer.area_frac_vs_chempot_plot(x_is_u_ads=True)
#     #
#     # def test_chempot_vs_gamma_clean(self):
#     #
#     #     plt = self.Cu_analyzer.chempot_vs_gamma_clean()
#     #     for el in self.Oads_analyzer_dict.keys():
#     #         # Test WulffShape for adsorbed surfaces
#     #         analyzer = self.Oads_analyzer_dict[el]
#     #         plt = analyzer.chempot_vs_gamma_clean(x_is_u_ads=True)
#     #
#     # def test_chempot_vs_gamma_facet(self):
#     #
#     #     for el in self.metals_O_entry_dict.keys():
#     #         for hkl in self.metals_O_entry_dict[el].keys():
#     #             # Test WulffShape for adsorbed surfaces
#     #             analyzer = self.Oads_analyzer_dict[el]
#     #             plt = analyzer.chempot_vs_gamma_facet(hkl)
#
#
# if __name__ == "__main__":
#     unittest.main()
#
#
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
        entry = ComputedStructureEntry.from_dict(entries[k])
        entry_dict[hkl][SlabEntry(entry, hkl, name=k)] = []

    return entry_dict


def load_O_adsorption():

    # Loads the dictionary for clean and O adsorbed Rh, Pt, and Ni entries

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
                    metals_O_entry_dict[el][(1, 1, 1)][SlabEntry(entry, (1,1,1), name=k)] = []
                if "110" in k:
                    metals_O_entry_dict[el][(1, 1, 0)][SlabEntry(entry, (1,1,0), name=k)] = []
                if "100" in k:
                    metals_O_entry_dict[el][(1, 0, 0)][SlabEntry(entry, (1,0,0), name=k)] = []

    with open(os.path.join(get_path(""), "csentries_o_ads.json")) as entries:
        entries = json.loads(entries.read())
    for k in entries.keys():
        entry = ComputedStructureEntry.from_dict(entries[k])
        for el in metals_O_entry_dict.keys():
            if el in k:
                if "111" in k:
                    clean = list(metals_O_entry_dict[el][(1, 1, 1)].keys())[0]
                    metals_O_entry_dict[el][(1, 1, 1)][clean] = [SlabEntry(entry, (1,1,1), name=k, adsorbate="O")]
                if "110" in k:
                    clean = list(metals_O_entry_dict[el][(1, 1, 0)].keys())[0]
                    metals_O_entry_dict[el][(1, 1, 0)][clean] = [SlabEntry(entry, (1,1,0), name=k, adsorbate="O")]
                if "100" in k:
                    clean = list(metals_O_entry_dict[el][(1, 0, 0)].keys())[0]
                    metals_O_entry_dict[el][(1, 0, 0)][clean] = [SlabEntry(entry, (1,0,0), name=k, adsorbate="O")]

    return metals_O_entry_dict