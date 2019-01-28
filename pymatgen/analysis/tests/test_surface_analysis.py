
import unittest
import os
import warnings
import json
from sympy import Number, Symbol

from pymatgen.analysis.surface_analysis import SlabEntry, SurfaceEnergyPlotter, \
    NanoscaleStability, WorkFunctionAnalyzer
from pymatgen.util.testing import PymatgenTest
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatgen.analysis.wulff import WulffShape

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


class SlabEntryTest(PymatgenTest):

    def setUp(self):

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

        with open(os.path.join(get_path(""), 'ucell_entries.txt')) as ucell_entries:
            ucell_entries = json.loads(ucell_entries.read())
        self.ucell_entries = ucell_entries

        # Load objects for O adsorption tests
        self.metals_O_entry_dict = load_O_adsorption()

        # Load objects for Cu test
        self.Cu_entry_dict = get_entry_dict(os.path.join(get_path(""),
                                                         "Cu_entries.txt"))
        self.assertEqual(len(self.Cu_entry_dict.keys()), 13)
        self.Cu_ucell_entry = ComputedStructureEntry.from_dict(self.ucell_entries["Cu"])

        # Load dummy MgO slab entries
        self.MgO_ucell_entry = ComputedStructureEntry.from_dict(self.ucell_entries["MgO"])
        self.Mg_ucell_entry = ComputedStructureEntry.from_dict(self.ucell_entries["Mg"])
        self.MgO_slab_entry_dict = get_entry_dict(os.path.join(get_path(""),
                                                               "MgO_slab_entries.txt"))

    def test_properties(self):
        # Test cases for getting adsorption related quantities for a 1/4
        # monolalyer adsorption of O on the low MMI surfaces of Pt, Ni and Rh

        for el in self.metals_O_entry_dict.keys():
            el_ucell = ComputedStructureEntry.from_dict(self.ucell_entries[el])
            for hkl in self.metals_O_entry_dict[el].keys():
                for clean in self.metals_O_entry_dict[el][hkl].keys():
                    for ads in self.metals_O_entry_dict[el][hkl][clean]:
                        ml = ads.get_unit_primitive_area
                        self.assertAlmostEqual(ml, 4, 2)
                        self.assertAlmostEqual(ads.get_monolayer, 1/4, 2)
                        Nads = ads.Nads_in_slab
                        self.assertEqual(Nads, 1)
                        self.assertEqual(ads.Nsurfs_ads_in_slab, 1)

                        # Determine the correct binding energy
                        with open(os.path.join(get_path(""),
                                               'isolated_O_entry.txt')) as isolated_O_entry:
                            isolated_O_entry = json.loads(isolated_O_entry.read())
                        O = ComputedStructureEntry.from_dict(isolated_O_entry)
                        gbind = (ads.energy - ml*clean.energy)/Nads - O.energy_per_atom
                        self.assertEqual(gbind, ads.gibbs_binding_energy())
                        # Determine the correction Gibbs adsorption energy
                        eads = Nads * gbind
                        self.assertEqual(eads, ads.gibbs_binding_energy(eads=True))
                        se = ads.surface_energy(el_ucell)
                        self.assertAlmostEqual(se.as_coefficients_dict()[Symbol("delu_O")],
                                               (-1/2)*ads.surface_area**(-1))

    def test_create_slab_label(self):

        for el in self.metals_O_entry_dict.keys():
            for hkl in self.metals_O_entry_dict[el].keys():
                # Test WulffShape for adsorbed surfaces
                for clean in self.metals_O_entry_dict[el][hkl]:
                    label = clean.create_slab_label
                    comp = str(clean.composition.reduced_composition)
                    self.assertEqual(str(hkl)+" %s" %(comp), label)

                    for ads in self.metals_O_entry_dict[el][hkl][clean]:
                        label = ads.create_slab_label
                        self.assertEqual(label, str(hkl)+" %s+O, 0.250 ML" %(comp))

    def test_surface_energy(self):
        # For a nonstoichiometric case, the cheimcal potentials do not
        # cancel out, they serve as a reservoir for any missing atoms
        for slab_entry in self.MgO_slab_entry_dict[(1,1,1)].keys():
            se = slab_entry.surface_energy(self.MgO_ucell_entry,
                                           ref_entries=[self.Mg_ucell_entry])
            self.assertEqual(tuple(se.as_coefficients_dict().keys()),
                             (Number(1), Symbol("delu_Mg")))

        # For the case of a clean, stoichiometric slab, the surface energy
        # should be constant (i.e. surface energy is a constant).
        all_se = []
        ECu = self.Cu_ucell_entry.energy_per_atom
        for hkl in self.Cu_entry_dict.keys():
            slab_entry = list(self.Cu_entry_dict[hkl].keys())[0]
            se = slab_entry.surface_energy(self.Cu_ucell_entry)
            all_se.append(se)
            # Manually calculate surface energy
            manual_se = (slab_entry.energy - \
                         ECu *len(slab_entry.structure))/(2*slab_entry.surface_area)
            self.assertArrayAlmostEqual(float(se), manual_se, 10)

        # The (111) facet should be the most stable
        clean111_entry = list(self.Cu_entry_dict[(1,1,1)].keys())[0]
        se_Cu111 = clean111_entry.surface_energy(self.Cu_ucell_entry)
        self.assertEqual(min(all_se), se_Cu111)

    def test_cleaned_up_slab(self):
        # The cleaned up slab should have the same reduced formula as a clean slab
        for el in self.metals_O_entry_dict.keys():
            for hkl in self.metals_O_entry_dict[el].keys():
                for clean in self.metals_O_entry_dict[el][hkl].keys():
                    for ads in self.metals_O_entry_dict[el][hkl][clean]:
                        s = ads.cleaned_up_slab
                        self.assertEqual(s.composition.reduced_composition,
                                         clean.composition.reduced_composition)


class SurfaceEnergyPlotterTest(PymatgenTest):

    def setUp(self):

        entry_dict = get_entry_dict(os.path.join(get_path(""),
                                                 "Cu_entries.txt"))
        self.Cu_entry_dict = entry_dict
        with open(os.path.join(get_path(""), 'ucell_entries.txt')) as ucell_entries:
            ucell_entries = json.loads(ucell_entries.read())

        self.Cu_ucell_entry = ComputedStructureEntry.from_dict(ucell_entries["Cu"])
        self.Cu_analyzer = SurfaceEnergyPlotter(entry_dict, self.Cu_ucell_entry)

        self.metals_O_entry_dict = load_O_adsorption()
        ucell_entry = ComputedStructureEntry.from_dict(ucell_entries["Pt"])
        self.Pt_analyzer = SurfaceEnergyPlotter(self.metals_O_entry_dict["Pt"],
                                                ucell_entry)
        ucell_entry = ComputedStructureEntry.from_dict(ucell_entries["Ni"])
        self.Ni_analyzer = SurfaceEnergyPlotter(self.metals_O_entry_dict["Ni"],
                                                ucell_entry)
        ucell_entry = ComputedStructureEntry.from_dict(ucell_entries["Rh"])
        self.Rh_analyzer = SurfaceEnergyPlotter(self.metals_O_entry_dict["Rh"],
                                                ucell_entry)
        self.Oads_analyzer_dict = {"Pt": self.Pt_analyzer,
                                   "Ni": self.Ni_analyzer,
                                   "Rh": self.Rh_analyzer}

    def test_get_stable_entry_at_u(self):

        for el in self.Oads_analyzer_dict.keys():
            plotter = self.Oads_analyzer_dict[el]
            for hkl in plotter.all_slab_entries.keys():
                # Test that the surface energy is clean for specific range of chempot
                entry1, gamma1 = \
                    plotter.get_stable_entry_at_u(hkl, delu_dict={Symbol("delu_O"): -7})
                entry2, gamma2 = \
                    plotter.get_stable_entry_at_u(hkl, delu_dict={Symbol("delu_O"): -6})
                self.assertEqual(gamma1, gamma2)
                self.assertEqual(entry1.label, entry2.label)

                # Now test that for a high chempot, adsorption
                # occurs and gamma is not equal to clean gamma
                entry3, gamma3 = \
                    plotter.get_stable_entry_at_u(hkl, delu_dict={Symbol("delu_O"): -1})
                self.assertNotEqual(entry3.label, entry2.label)
                self.assertNotEqual(gamma3, gamma2)

                # For any chempot greater than -6, surface energy should vary
                # but the configuration should remain the same
                entry4, gamma4 = \
                    plotter.get_stable_entry_at_u(hkl, delu_dict={Symbol("delu_O"): 0})
                self.assertEqual(entry3.label, entry4.label)
                self.assertNotEqual(gamma3, gamma4)

    def test_wulff_from_chempot(self):

        # Test if it generates a Wulff shape, test if
        # all the facets for Cu wulff shape are inside.
        Cu_wulff = self.Cu_analyzer.wulff_from_chempot()
        area_frac_dict = Cu_wulff.area_fraction_dict
        facets_hkl = [(1,1,1), (3,3,1), (3,1,0), (1,0,0),
                      (3,1,1), (2,1,0), (2,2,1)]
        for hkl in area_frac_dict.keys():
            if hkl in facets_hkl:
                self.assertNotEqual(area_frac_dict[hkl], 0)
            else:
                self.assertEqual(area_frac_dict[hkl], 0)

        for el in self.Oads_analyzer_dict.keys():
            # Test WulffShape for adsorbed surfaces
            analyzer = self.Oads_analyzer_dict[el]
            # chempot = analyzer.max_adsorption_chempot_range(0)
            wulff = analyzer.wulff_from_chempot(delu_default=-6)
            se = wulff.weighted_surface_energy

        # Test if a different Wulff shape is generated
        # for Ni when adsorption comes into play
        wulff_neg7 = self.Oads_analyzer_dict["Ni"].wulff_from_chempot(delu_default=-7)
        wulff_neg6 = self.Oads_analyzer_dict["Ni"].wulff_from_chempot(delu_default=-6)
        self.assertEqual(wulff_neg7.weighted_surface_energy,
                         wulff_neg6.weighted_surface_energy)
        wulff_neg55 = self.Oads_analyzer_dict["Ni"].wulff_from_chempot(delu_default=-5.5)
        self.assertNotEqual(wulff_neg55.weighted_surface_energy,
                            wulff_neg6.weighted_surface_energy)
        wulff_neg525 = self.Oads_analyzer_dict["Ni"].wulff_from_chempot(delu_default=-5.25)
        self.assertNotEqual(wulff_neg55.weighted_surface_energy,
                            wulff_neg525.weighted_surface_energy)

    def test_color_palette_dict(self):

        for el in self.metals_O_entry_dict.keys():
            analyzer = self.Oads_analyzer_dict[el]
            color_dict = analyzer.color_palette_dict()
            for hkl in self.metals_O_entry_dict[el].keys():
                for clean in self.metals_O_entry_dict[el][hkl].keys():
                    color = color_dict[clean]
                    for ads in self.metals_O_entry_dict[el][hkl][clean]:
                        color = color_dict[ads]

    def test_get_surface_equilibrium(self):
        # For clean stoichiometric system, the two equations should
        # be parallel because the surface energy is a constant. Then
        # get_surface_equilibrium should return None
        clean111_entry = list(self.Cu_entry_dict[(1, 1, 1)].keys())[0]
        clean100_entry = list(self.Cu_entry_dict[(1, 0, 0)].keys())[0]
        soln = self.Cu_analyzer.get_surface_equilibrium([clean111_entry,
                                                         clean100_entry])
        self.assertFalse(soln)

        # For adsorbed system, we should find one intercept
        Pt_entries = self.metals_O_entry_dict["Pt"]
        clean = list(Pt_entries[(1, 1, 1)].keys())[0]
        ads = Pt_entries[(1, 1, 1)][clean][0]
        Pt_analyzer = self.Oads_analyzer_dict["Pt"]
        soln = Pt_analyzer.get_surface_equilibrium([clean, ads])

        self.assertNotEqual(list(soln.values())[0], list(soln.values())[1])

        # Check if the number of parameters for adsorption are correct
        self.assertEqual((Symbol("delu_O"), Symbol("gamma")), tuple(soln.keys()))
        # Adsorbed systems have a b2=(-1*Nads) / (Nsurfs * Aads)
        se = ads.surface_energy(Pt_analyzer.ucell_entry, Pt_analyzer.ref_entries)
        self.assertAlmostEqual(se.as_coefficients_dict()[Symbol("delu_O")],
                               -1 / (2 * ads.surface_area))

    def test_stable_u_range_dict(self):
        for el in self.Oads_analyzer_dict.keys():
            analyzer = self.Oads_analyzer_dict[el]

        stable_u_range = analyzer.stable_u_range_dict([-1,0], Symbol("delu_O"),
                                                      no_doped=False)
        all_u = []
        for entry in stable_u_range.keys():
            all_u.extend(stable_u_range[entry])
        self.assertGreater(len(all_u), 1)

    def test_entry_dict_from_list(self):

        # Plug in a list of entries to see if it works
        all_Pt_slab_entries = []
        Pt_entries = self.Pt_analyzer.all_slab_entries
        for hkl in Pt_entries.keys():
            for clean in Pt_entries[hkl].keys():
                all_Pt_slab_entries.append(clean)
                all_Pt_slab_entries.extend(Pt_entries[hkl][clean])

        a = SurfaceEnergyPlotter(all_Pt_slab_entries,
                                 self.Pt_analyzer.ucell_entry)
        self.assertEqual(type(a).__name__, "SurfaceEnergyPlotter")

    # def test_monolayer_vs_BE(self):
    #     for el in self.Oads_analyzer_dict.keys():
    #         # Test WulffShape for adsorbed surfaces
    #         analyzer = self.Oads_analyzer_dict[el]
    #         plt = analyzer.monolayer_vs_BE()
    #
    # def test_area_frac_vs_chempot_plot(self):
    #
    #     for el in self.Oads_analyzer_dict.keys():
    #         # Test WulffShape for adsorbed surfaces
    #         analyzer = self.Oads_analyzer_dict[el]
    #         plt = analyzer.area_frac_vs_chempot_plot(x_is_u_ads=True)
    #
    # def test_chempot_vs_gamma_clean(self):
    #
    #     plt = self.Cu_analyzer.chempot_vs_gamma_clean()
    #     for el in self.Oads_analyzer_dict.keys():
    #         # Test WulffShape for adsorbed surfaces
    #         analyzer = self.Oads_analyzer_dict[el]
    #         plt = analyzer.chempot_vs_gamma_clean(x_is_u_ads=True)
    #
    # def test_chempot_vs_gamma_facet(self):
    #
    #     for el in self.metals_O_entry_dict.keys():
    #         for hkl in self.metals_O_entry_dict[el].keys():
    #             # Test WulffShape for adsorbed surfaces
    #             analyzer = self.Oads_analyzer_dict[el]
    #             plt = analyzer.chempot_vs_gamma_facet(hkl)
    # def test_surface_chempot_range_map(self):
    #
    #     for el in self.metals_O_entry_dict.keys():
    #         for hkl in self.metals_O_entry_dict[el].keys():
    #             # Test WulffShape for adsorbed surfaces
    #             analyzer = self.Oads_analyzer_dict[el]
    #             plt = analyzer.chempot_vs_gamma_facet(hkl)


class WorkfunctionAnalyzerTest(PymatgenTest):

    def setUp(self):

        self.kwargs = {"poscar_filename": get_path("CONTCAR.relax1.gz"),
                       "locpot_filename": get_path("LOCPOT.gz"),
                       "outcar_filename": get_path("OUTCAR.relax1.gz")}
        self.wf_analyzer = WorkFunctionAnalyzer.from_files(**self.kwargs)

    def test_attributes(self):
        wf_analyzer_shift = WorkFunctionAnalyzer.from_files(shift=0.25, **self.kwargs)
        self.assertEqual("%.1f" %(self.wf_analyzer.ave_bulk_p),
                         "%.1f" %(wf_analyzer_shift.ave_bulk_p))

    # def test_plt(self):
    #     plt = self.wf_analyzer.get_locpot_along_slab_plot()
    #     self.assertEqual(type(plt).__name__, "module")

    def test_is_converged(self):
        self.assertTrue(self.wf_analyzer.is_converged())


class NanoscaleStabilityTest(PymatgenTest):

    def setUp(self):

        # Load all entries
        La_hcp_entry_dict = get_entry_dict(os.path.join(get_path(""),
                                                        "La_hcp_entries.txt"))
        La_fcc_entry_dict = get_entry_dict(os.path.join(get_path(""),
                                                        "La_fcc_entries.txt"))
        with open(os.path.join(get_path(""), 'ucell_entries.txt')) as ucell_entries:
            ucell_entries = json.loads(ucell_entries.read())
        La_hcp_ucell_entry = ComputedStructureEntry.from_dict(ucell_entries["La_hcp"])
        La_fcc_ucell_entry = ComputedStructureEntry.from_dict(ucell_entries["La_fcc"])

        # Set up the NanoscaleStabilityClass
        self.La_hcp_analyzer = SurfaceEnergyPlotter(La_hcp_entry_dict,
                                                    La_hcp_ucell_entry)
        self.La_fcc_analyzer = SurfaceEnergyPlotter(La_fcc_entry_dict,
                                                    La_fcc_ucell_entry)
        self.nanoscale_stability = NanoscaleStability([self.La_fcc_analyzer,
                                                       self.La_hcp_analyzer])

    def test_stability_at_r(self):

        # Check that we have a different polymorph that is
        # stable below or above the equilibrium particle size
        r = self.nanoscale_stability.solve_equilibrium_point(self.La_hcp_analyzer,
                                                             self.La_fcc_analyzer)*10

        # hcp phase of La particle should be the stable
        # polymorph above the equilibrium radius
        hcp_wulff = self.La_hcp_analyzer.wulff_from_chempot()
        bulk = self.La_hcp_analyzer.ucell_entry
        ghcp, rhcp = self.nanoscale_stability.wulff_gform_and_r(hcp_wulff, bulk, r+10,
                                                                from_sphere_area=True)
        fcc_wulff = self.La_fcc_analyzer.wulff_from_chempot()
        bulk = self.La_fcc_analyzer.ucell_entry
        gfcc, rfcc = self.nanoscale_stability.wulff_gform_and_r(fcc_wulff, bulk, r+10,
                                                                from_sphere_area=True)
        self.assertGreater(gfcc, ghcp)

        # fcc phase of La particle should be the stable
        # polymorph below the equilibrium radius
        hcp_wulff = self.La_hcp_analyzer.wulff_from_chempot()
        bulk = self.La_hcp_analyzer.ucell_entry
        ghcp, rhcp = self.nanoscale_stability.wulff_gform_and_r(hcp_wulff, bulk, r-10,
                                                                from_sphere_area=True)
        fcc_wulff = self.La_fcc_analyzer.wulff_from_chempot()
        bulk = self.La_fcc_analyzer.ucell_entry
        gfcc, rfcc = self.nanoscale_stability.wulff_gform_and_r(fcc_wulff, bulk, r-10,
                                                                from_sphere_area=True)
        self.assertLess(gfcc, ghcp)

    def test_scaled_wulff(self):
        # Ensure for a given radius, the effective radius
        # of the Wulff shape is the same (correctly scaled)
        hcp_wulff = self.La_hcp_analyzer.wulff_from_chempot()
        fcc_wulff = self.La_fcc_analyzer.wulff_from_chempot()
        w1 = self.nanoscale_stability.scaled_wulff(hcp_wulff, 10)
        w2 = self.nanoscale_stability.scaled_wulff(fcc_wulff, 10)
        self.assertAlmostEqual(w1.effective_radius, w2.effective_radius)
        self.assertAlmostEqual(w1.effective_radius, 10)
        self.assertAlmostEqual(10, w2.effective_radius)


def get_entry_dict(filename):
    # helper to generate an entry_dict

    entry_dict = {}
    with open(filename) as entries:
        entries = json.loads(entries.read())
    for k in entries.keys():
        n = k[25:]
        miller_index = []
        for i, s in enumerate(n):
            if s == "_":
                break
            if s == "-":
                continue
            t = int(s)
            if n[i - 1] == "-":
                t *= -1
            miller_index.append(t)
        hkl = tuple(miller_index)
        if hkl not in entry_dict.keys():
            entry_dict[hkl] = {}
        entry = ComputedStructureEntry.from_dict(entries[k])
        entry_dict[hkl][SlabEntry(entry.structure, entry.energy, hkl, label=k)] = []

    return entry_dict


def load_O_adsorption():

    # Loads the dictionary for clean and O adsorbed Rh, Pt, and Ni entries

    # Load the adsorbate as an entry
    with open(os.path.join(get_path(""),
                           'isolated_O_entry.txt')) as isolated_O_entry:
        isolated_O_entry = json.loads(isolated_O_entry.read())
    O = ComputedStructureEntry.from_dict(isolated_O_entry)

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
                    clean = SlabEntry(entry.structure, entry.energy,
                                      (1,1,1), label=k+"_clean")
                    metals_O_entry_dict[el][(1, 1, 1)][clean] = []
                if "110" in k:
                    clean = SlabEntry(entry.structure, entry.energy,
                                      (1, 1, 0), label=k + "_clean")
                    metals_O_entry_dict[el][(1, 1, 0)][clean] = []
                if "100" in k:
                    clean = SlabEntry(entry.structure, entry.energy,
                                      (1,0,0), label=k+"_clean")
                    metals_O_entry_dict[el][(1, 0, 0)][clean] = []

    with open(os.path.join(get_path(""), "csentries_o_ads.json")) as entries:
        entries = json.loads(entries.read())
    for k in entries.keys():
        entry = ComputedStructureEntry.from_dict(entries[k])
        for el in metals_O_entry_dict.keys():
            if el in k:
                if "111" in k:
                    clean = list(metals_O_entry_dict[el][(1, 1, 1)].keys())[0]
                    ads = SlabEntry(entry.structure, entry.energy, (1,1,1),
                                    label=k+"_O", adsorbates=[O], clean_entry=clean)
                    metals_O_entry_dict[el][(1, 1, 1)][clean] = [ads]
                if "110" in k:
                    clean = list(metals_O_entry_dict[el][(1, 1, 0)].keys())[0]
                    ads = SlabEntry(entry.structure, entry.energy, (1,1,0),
                                    label=k+"_O", adsorbates=[O], clean_entry=clean)
                    metals_O_entry_dict[el][(1, 1, 0)][clean] = [ads]
                if "100" in k:
                    clean = list(metals_O_entry_dict[el][(1, 0, 0)].keys())[0]
                    ads = SlabEntry(entry.structure, entry.energy, (1,0,0),
                                    label=k+"_O", adsorbates=[O], clean_entry=clean)
                    metals_O_entry_dict[el][(1, 0, 0)][clean] = [ads]

    return metals_O_entry_dict


if __name__ == "__main__":
    unittest.main()