from __future__ import annotations

import json

import matplotlib.pyplot as plt
import numpy as np
import orjson
import pytest
from monty.json import MontyEncoder
from numpy.testing import assert_allclose
from pytest import approx
from sympy import Number, Symbol

from pymatgen.analysis.surface_analysis import (
    NanoscaleStability,
    SlabEntry,
    SurfaceEnergyPlotter,
    WorkFunctionAnalyzer,
    sub_chempots,
)
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatgen.util.testing import MatSciTest
from tests.testing import TEST_FILES_DIR

__author__ = "Richard Tran"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Richard Tran"
__email__ = "rit001@eng.ucsd.edu"
__date__ = "Aug 24, 2017"


TEST_DIR = f"{TEST_FILES_DIR}/surfaces"


class TestSlabEntry(MatSciTest):
    def setup_method(self):
        with open(f"{TEST_DIR}/ucell_entries.txt", encoding="utf-8") as file:
            ucell_entries = orjson.loads(file.read())
        self.ucell_entries = ucell_entries

        # Load objects for O adsorption tests
        self.metals_O_entry_dict = load_O_adsorption()

        # Load objects for Cu test
        self.Cu_entry_dict = get_entry_dict(f"{TEST_DIR}/Cu_entries.txt")
        assert len(self.Cu_entry_dict) == 13
        self.Cu_ucell_entry = ComputedStructureEntry.from_dict(self.ucell_entries["Cu"])

        # Load dummy MgO slab entries
        self.MgO_ucell_entry = ComputedStructureEntry.from_dict(self.ucell_entries["MgO"])
        self.Mg_ucell_entry = ComputedStructureEntry.from_dict(self.ucell_entries["Mg"])
        self.MgO_slab_entry_dict = get_entry_dict(f"{TEST_DIR}/MgO_slab_entries.txt")

    def test_properties(self):
        # Test cases for getting adsorption related quantities for a 1/4
        # monolayer adsorption of O on the low MMI surfaces of Pt, Ni and Rh

        for el, val in self.metals_O_entry_dict.items():
            el_ucell = ComputedStructureEntry.from_dict(self.ucell_entries[el])
            for hkl in val:
                for clean in self.metals_O_entry_dict[el][hkl]:
                    for ads in self.metals_O_entry_dict[el][hkl][clean]:
                        ml = ads.get_unit_primitive_area
                        assert ml == approx(4, abs=1e-2)
                        assert ads.get_monolayer == approx(1 / 4, abs=1e-2)
                        Nads = ads.Nads_in_slab
                        assert Nads == 1
                        assert ads.Nsurfs_ads_in_slab == 1

                        # Determine the correct binding energy
                        with open(f"{TEST_DIR}/isolated_O_entry.txt", encoding="utf-8") as txt_file:
                            isolated_O_entry = orjson.loads(txt_file.read())
                        O_cse = ComputedStructureEntry.from_dict(isolated_O_entry)
                        g_bind = (ads.energy - ml * clean.energy) / Nads - O_cse.energy_per_atom
                        assert g_bind == ads.gibbs_binding_energy()
                        # Determine the correction Gibbs adsorption energy
                        eads = Nads * g_bind
                        assert eads == ads.gibbs_binding_energy(eads=True)
                        se = ads.surface_energy(el_ucell)
                        assert se.as_coefficients_dict()[Symbol("delu_O")] == approx(
                            (-1 / 2) * ads.surface_area ** (-1)
                        )

    def test_create_slab_label(self):
        for el, val in self.metals_O_entry_dict.items():
            for hkl in val:
                # Test WulffShape for adsorbed surfaces
                for clean in self.metals_O_entry_dict[el][hkl]:
                    label = clean.create_slab_label
                    comp = str(clean.composition.reduced_composition)
                    assert f"{hkl} {comp}" == label

                    for ads in self.metals_O_entry_dict[el][hkl][clean]:
                        label = ads.create_slab_label
                        assert label == f"{hkl} {comp}+O, 0.250 ML"

    def test_surface_energy(self):
        # For a non-stoichiometric case, the chemical potentials do not
        # cancel out, they serve as a reservoir for any missing atoms
        for slab_entry in self.MgO_slab_entry_dict[1, 1, 1]:
            se = slab_entry.surface_energy(self.MgO_ucell_entry, ref_entries=[self.Mg_ucell_entry])
            assert tuple(se.as_coefficients_dict()) == (Number(1), Symbol("delu_Mg"))

        # For the case of a clean, stoichiometric slab, the surface energy
        # should be constant (i.e. surface energy is a constant).
        all_se = []
        ECu = self.Cu_ucell_entry.energy_per_atom
        for val in self.Cu_entry_dict.values():
            slab_entry = next(iter(val))
            se = slab_entry.surface_energy(self.Cu_ucell_entry)
            all_se.append(se)
            # Manually calculate surface energy
            manual_se = (slab_entry.energy - ECu * len(slab_entry.structure)) / (2 * slab_entry.surface_area)
            assert_allclose(float(se), manual_se, 10)

        # The (111) facet should be the most stable
        clean111_entry = next(iter(self.Cu_entry_dict[1, 1, 1]))
        se_Cu111 = clean111_entry.surface_energy(self.Cu_ucell_entry)
        assert min(all_se) == se_Cu111

    def test_cleaned_up_slab(self):
        # The cleaned up slab should have the same reduced formula as a clean slab
        for el, val in self.metals_O_entry_dict.items():
            for hkl in val:
                for clean in self.metals_O_entry_dict[el][hkl]:
                    for ads in self.metals_O_entry_dict[el][hkl][clean]:
                        slab_clean = ads.cleaned_up_slab
                        assert slab_clean.composition.reduced_composition == clean.composition.reduced_composition

    def test_as_dict_from_dict(self):
        clean = next(iter(self.metals_O_entry_dict["Pt"][1, 1, 1]))
        ads = self.metals_O_entry_dict["Pt"][1, 1, 1][clean][0]
        dct = ads.as_dict()
        assert dct["@class"] == "SlabEntry"
        assert dct["miller_index"] == (1, 1, 1)
        assert dct["clean_entry"] is clean

        # In-memory round trip
        ads2 = SlabEntry.from_dict(dct)
        assert ads2.energy == approx(ads.energy)
        assert ads2.structure == ads.structure
        assert ads2.clean_entry is clean
        assert ads2.Nads_in_slab == ads.Nads_in_slab

        # Round trip through JSON
        dct = json.loads(json.dumps(ads.as_dict(), cls=MontyEncoder))
        ads3 = SlabEntry.from_dict(dct)
        assert isinstance(ads3.clean_entry, SlabEntry)
        assert ads3.clean_entry.energy == approx(clean.energy)
        assert tuple(ads3.miller_index) == (1, 1, 1)
        assert ads3.label == ads.label
        assert ads3.composition == ads.composition
        assert ads3.get_monolayer == approx(ads.get_monolayer)
        assert ads3.gibbs_binding_energy() == approx(ads.gibbs_binding_energy())

    def test_from_computed_structure_entry(self):
        cse = ComputedStructureEntry(self.Cu_ucell_entry.structure, -10.0)
        entry = SlabEntry.from_computed_structure_entry(cse, (1, 0, 0), label="Cu100", marker="--", color="r")
        assert entry.miller_index == (1, 0, 0)
        assert entry.label == "Cu100"
        assert entry.energy == approx(-10.0)
        assert entry.structure == cse.structure
        assert entry.adsorbates == []
        assert entry.clean_entry is None
        assert entry.mark == "--"
        assert entry.color == "r"

    def test_create_slab_label_from_data(self):
        clean = next(iter(self.Cu_entry_dict[1, 1, 1]))
        entry = SlabEntry(clean.structure, clean.energy, (1, 1, 1), data={"label": "custom"})
        assert entry.create_slab_label == "custom"

    def test_surface_energy_missing_ref_entries(self):
        slab_entry = next(iter(self.MgO_slab_entry_dict[1, 1, 1]))
        with pytest.warns(UserWarning, match="Elemental references missing"):
            se = slab_entry.surface_energy(self.MgO_ucell_entry)
        assert isinstance(se, float)


class TestSurfaceEnergyPlotter(MatSciTest):
    def setup_method(self):
        entry_dict = get_entry_dict(f"{TEST_DIR}/Cu_entries.txt")
        self.Cu_entry_dict = entry_dict
        with open(f"{TEST_DIR}/ucell_entries.txt", encoding="utf-8") as file:
            ucell_entries = orjson.loads(file.read())

        self.Cu_ucell_entry = ComputedStructureEntry.from_dict(ucell_entries["Cu"])
        self.Cu_analyzer = SurfaceEnergyPlotter(entry_dict, self.Cu_ucell_entry)

        self.metals_O_entry_dict = load_O_adsorption()
        ucell_entry = ComputedStructureEntry.from_dict(ucell_entries["Pt"])
        self.Pt_analyzer = SurfaceEnergyPlotter(self.metals_O_entry_dict["Pt"], ucell_entry)
        ucell_entry = ComputedStructureEntry.from_dict(ucell_entries["Ni"])
        self.Ni_analyzer = SurfaceEnergyPlotter(self.metals_O_entry_dict["Ni"], ucell_entry)
        ucell_entry = ComputedStructureEntry.from_dict(ucell_entries["Rh"])
        self.Rh_analyzer = SurfaceEnergyPlotter(self.metals_O_entry_dict["Rh"], ucell_entry)
        self.Oads_analyzer_dict = {
            "Pt": self.Pt_analyzer,
            "Ni": self.Ni_analyzer,
            "Rh": self.Rh_analyzer,
        }

    def test_get_stable_entry_at_u(self):
        for plotter in self.Oads_analyzer_dict.values():
            for hkl in plotter.all_slab_entries:
                # Test that the surface energy is clean for specific range of chempot
                entry1, gamma1 = plotter.get_stable_entry_at_u(hkl, delu_dict={Symbol("delu_O"): -7})
                entry2, gamma2 = plotter.get_stable_entry_at_u(hkl, delu_dict={Symbol("delu_O"): -6})
                assert gamma1 == gamma2
                assert entry1.label == entry2.label

                # Now test that for a high chempot, adsorption
                # occurs and gamma is not equal to clean gamma
                entry3, gamma3 = plotter.get_stable_entry_at_u(hkl, delu_dict={Symbol("delu_O"): -1})
                assert entry3.label != entry2.label
                assert gamma3 != gamma2

                # For any chempot greater than -6, surface energy should vary
                # but the configuration should remain the same
                entry4, gamma4 = plotter.get_stable_entry_at_u(hkl, delu_dict={Symbol("delu_O"): 0})
                assert entry3.label == entry4.label
                assert gamma3 != gamma4

    def test_wulff_from_chempot(self):
        # Test if it generates a Wulff shape, test if
        # all the facets for Cu wulff shape are inside.
        Cu_wulff = self.Cu_analyzer.wulff_from_chempot()
        area_frac_dict = Cu_wulff.area_fraction_dict
        facets_hkl = [
            (1, 1, 1),
            (3, 3, 1),
            (3, 1, 0),
            (1, 0, 0),
            (3, 1, 1),
            (2, 1, 0),
            (2, 2, 1),
        ]
        for hkl in area_frac_dict:
            if hkl in facets_hkl:
                assert area_frac_dict[hkl] != 0
            else:
                assert area_frac_dict[hkl] == 0

        for analyzer in self.Oads_analyzer_dict.values():
            # Test WulffShape for adsorbed surfaces
            # chempot = analyzer.max_adsorption_chempot_range(0)
            wulff = analyzer.wulff_from_chempot(delu_default=-6)
            _ = wulff.weighted_surface_energy

        # Test if a different Wulff shape is generated
        # for Ni when adsorption comes into play
        wulff_neg7 = self.Oads_analyzer_dict["Ni"].wulff_from_chempot(delu_default=-7)
        wulff_neg6 = self.Oads_analyzer_dict["Ni"].wulff_from_chempot(delu_default=-6)
        assert wulff_neg7.weighted_surface_energy == wulff_neg6.weighted_surface_energy
        wulff_neg55 = self.Oads_analyzer_dict["Ni"].wulff_from_chempot(delu_default=-5.5)
        assert wulff_neg55.weighted_surface_energy != wulff_neg6.weighted_surface_energy
        wulff_neg525 = self.Oads_analyzer_dict["Ni"].wulff_from_chempot(delu_default=-5.25)
        assert wulff_neg55.weighted_surface_energy != wulff_neg525.weighted_surface_energy

    def test_color_palette_dict(self):
        for el, val in self.metals_O_entry_dict.items():
            analyzer = self.Oads_analyzer_dict[el]
            color_dict = analyzer.color_palette_dict()
            for hkl in val:
                for clean in self.metals_O_entry_dict[el][hkl]:
                    _ = color_dict[clean]
                    for ads in self.metals_O_entry_dict[el][hkl][clean]:
                        _ = color_dict[ads]

    def test_get_surface_equilibrium(self):
        # For clean stoichiometric system, the two equations should
        # be parallel because the surface energy is a constant. Then
        # get_surface_equilibrium should return None
        clean111_entry = next(iter(self.Cu_entry_dict[1, 1, 1]))
        clean100_entry = next(iter(self.Cu_entry_dict[1, 0, 0]))
        soln = self.Cu_analyzer.get_surface_equilibrium([clean111_entry, clean100_entry])
        assert not soln

        # For adsorbed system, we should find one intercept
        Pt_entries = self.metals_O_entry_dict["Pt"]
        clean = next(iter(Pt_entries[1, 1, 1]))
        ads = Pt_entries[1, 1, 1][clean][0]
        Pt_analyzer = self.Oads_analyzer_dict["Pt"]
        soln = Pt_analyzer.get_surface_equilibrium([clean, ads])

        assert next(iter(soln.values())) != list(soln.values())[1]

        # Check if the number of parameters for adsorption are correct
        assert (Symbol("delu_O"), Symbol("gamma")) == tuple(soln)
        # Adsorbed systems have a b2=(-1*Nads) / (Nsurfs * Aads)
        se = ads.surface_energy(Pt_analyzer.ucell_entry, Pt_analyzer.ref_entries)
        assert se.as_coefficients_dict()[Symbol("delu_O")] == approx(-1 / (2 * ads.surface_area))

    def test_stable_u_range_dict(self):
        analyzer = list(self.Oads_analyzer_dict.values())[-1]

        stable_u_range = analyzer.stable_u_range_dict([-1, 0], Symbol("delu_O"), no_doped=False)
        all_u = []
        for entry in stable_u_range:
            all_u.extend(stable_u_range[entry])
        assert len(all_u) > 1

    def test_entry_dict_from_list(self):
        # Plug in a list of entries to see if it works
        all_Pt_slab_entries = []
        Pt_entries = self.Pt_analyzer.all_slab_entries
        for hkl in Pt_entries:
            for clean in Pt_entries[hkl]:
                all_Pt_slab_entries.append(clean)
                all_Pt_slab_entries.extend(Pt_entries[hkl][clean])

        surf_ene_plotter = SurfaceEnergyPlotter(all_Pt_slab_entries, self.Pt_analyzer.ucell_entry)
        assert surf_ene_plotter.list_of_chempots == self.Pt_analyzer.list_of_chempots

    def test_monolayer_vs_BE(self):
        for analyzer in self.Oads_analyzer_dict.values():
            ax = analyzer.monolayer_vs_BE()
            assert isinstance(ax, plt.Axes)
            assert ax.get_xlabel() == "O Coverage (ML)"
            assert ax.get_ylabel() == "Binding Energy (eV)"
            lines = {line.get_label(): line for line in ax.get_lines()}
            assert set(lines) == {str(hkl) for hkl in analyzer.all_slab_entries}
            for hkl, entries in analyzer.all_slab_entries.items():
                ads = next(iter(entries.values()))[0]
                assert_allclose(lines[str(hkl)].get_xdata(), [ads.get_monolayer])
                assert_allclose(lines[str(hkl)].get_ydata(), [ads.gibbs_binding_energy()])
            plt.close("all")

        ax = self.Ni_analyzer.monolayer_vs_BE(plot_eads=True)
        assert ax.get_ylabel() == "Adsorption Energy (eV)"
        ads = self.metals_O_entry_dict["Ni"][1, 1, 1][next(iter(self.metals_O_entry_dict["Ni"][1, 1, 1]))][0]
        line = next(line for line in ax.get_lines() if line.get_label() == "(1, 1, 1)")
        assert_allclose(line.get_ydata(), [ads.gibbs_binding_energy(eads=True)])

    def test_monolayer_vs_BE_skips_clean_only_facets(self):
        # A facet without adsorbed entries must be skipped rather than crashing
        Ni_entries = self.metals_O_entry_dict["Ni"]
        entries = {
            (1, 1, 1): Ni_entries[1, 1, 1],
            (1, 0, 0): {clean: [] for clean in Ni_entries[1, 0, 0]},
        }
        analyzer = SurfaceEnergyPlotter(entries, self.Ni_analyzer.ucell_entry)
        ax = analyzer.monolayer_vs_BE()
        assert [line.get_label() for line in ax.get_lines()] == ["(1, 1, 1)"]
        plt.close("all")

    def test_area_frac_vs_chempot_plot(self):
        delu_O = Symbol("delu_O")
        ax = self.Ni_analyzer.area_frac_vs_chempot_plot(delu_O, [0, -8], increments=5)
        assert isinstance(ax, plt.Axes)
        assert "Fractional area" in ax.get_ylabel()
        assert ax.get_xlabel() == r"Chemical potential $\Delta\mu_{O}$ (eV)"
        assert ax.get_ylim() == approx((0, 1))
        facet_lines = [line for line in ax.get_lines() if line.get_label() in ("(1, 1, 1)", "(1, 0, 0)")]
        assert len(facet_lines) == 2
        for line in facet_lines:
            assert_allclose(line.get_xdata(), [-8, -6, -4, -2, 0])
        # Area fractions of all facets sum to unity at each chempot
        assert_allclose(sum(line.get_ydata() for line in facet_lines), 1)

        # The x-axis endpoint should match the Wulff shape at that chempot
        wulff = self.Ni_analyzer.wulff_from_chempot(delu_dict={delu_O: -8})
        line111 = next(line for line in facet_lines if line.get_label() == "(1, 1, 1)")
        assert line111.get_ydata()[0] == approx(wulff.area_fraction_dict[1, 1, 1])

    @pytest.mark.parametrize("show_unstable", [True, False])
    @pytest.mark.parametrize("JPERM2", [True, False])
    def test_chempot_vs_gamma(self, show_unstable, JPERM2):
        delu_O = Symbol("delu_O")
        ax = self.Ni_analyzer.chempot_vs_gamma(delu_O, [0, -8], show_unstable=show_unstable, JPERM2=JPERM2)
        assert isinstance(ax, plt.Axes)
        unit = "J/$m^{2}$" if JPERM2 else r"eV/$\AA^{2}$"
        assert ax.get_ylabel() == f"Surface energy ({unit})"
        assert ax.get_xlabel() == r"Chemical potential $\Delta\mu_{O}$ (eV)"

        # 4 surface energy lines (clean and O-adsorbed for 2 facets),
        # followed by 2 vertical "--k" lines from chempot_plot_addons
        lines = ax.get_lines()
        assert len(lines) == 6
        scale = 16.0217656 if JPERM2 else 1
        clean111 = next(iter(self.metals_O_entry_dict["Ni"][1, 1, 1]))
        gamma111 = self.Ni_analyzer.surfe_dict[clean111] * scale
        clean_line = lines[0]
        assert_allclose(np.array(clean_line.get_ydata(), dtype=float), [gamma111, gamma111])
        if show_unstable:
            assert_allclose(np.array(clean_line.get_xdata(), dtype=float), [-8, 0])
        else:
            u_eq = self.Ni_analyzer.get_surface_equilibrium(
                [clean111, self.metals_O_entry_dict["Ni"][1, 1, 1][clean111][0]]
            )[delu_O]
            assert_allclose(np.array(clean_line.get_xdata(), dtype=float), [-8, float(u_eq)])
        plt.close("all")

    def test_chempot_vs_gamma_options(self):
        delu_O = Symbol("delu_O")
        ax = self.Ni_analyzer.chempot_vs_gamma(delu_O, [-8, 0], miller_index=(1, 1, 1), use_entry_labels=True)
        labels = [line.get_label() for line in ax.get_lines()[:2]]
        assert labels == ["Ni_111_clean", "Ni_111_O"]
        plt.close("all")

        _, ax = plt.subplots()
        ax_out = self.Ni_analyzer.chempot_vs_gamma(delu_O, [-8, 0], plt=ax, no_label=True, no_doped=True)
        assert ax_out is ax
        # Only the two clean slabs (no labels), plus the 2 vertical lines
        lines = ax.get_lines()
        assert len(lines) == 4
        assert all(line.get_label().startswith("_") for line in lines[:2])
        plt.close("all")

    def test_chempot_vs_gamma_plot_one(self):
        delu_O = Symbol("delu_O")
        clean = next(iter(self.metals_O_entry_dict["Ni"][1, 1, 1]))
        ads = self.metals_O_entry_dict["Ni"][1, 1, 1][clean][0]
        _, ax = plt.subplots()
        ax = self.Ni_analyzer.chempot_vs_gamma_plot_one(ax, ads, delu_O, [0, -2], label="ads", JPERM2=True)
        (line,) = ax.get_lines()
        assert line.get_label() == "ads"
        # Stoichiometric w.r.t. bulk once adsorbates are removed: solid line
        assert line.get_linestyle() == "-"
        coeffs = self.Ni_analyzer.as_coeffs_dict[ads]
        expected = [16.0217656 * (coeffs[1] + coeffs[delu_O] * u) for u in (-2, 0)]
        assert_allclose(line.get_xdata(), [-2, 0])
        assert_allclose(np.array(line.get_ydata(), dtype=float), np.array(expected, dtype=float))
        plt.close("all")

    @pytest.mark.parametrize("JPERM2", [True, False])
    def test_BE_vs_clean_SE(self, JPERM2):
        ax = self.Ni_analyzer.BE_vs_clean_SE({Symbol("delu_O"): -6}, JPERM2=JPERM2)
        assert isinstance(ax, plt.Axes)
        assert ax.get_ylabel() == "Binding Energy (eV)"
        assert ("J/m^2" in ax.get_xlabel()) is JPERM2
        scale = 16.0217656 if JPERM2 else 1

        points = [coll.get_offsets()[0] for coll in ax.collections]
        expected = [
            [self.Ni_analyzer.surfe_dict[clean] * scale, ads.gibbs_binding_energy()]
            for entries in self.metals_O_entry_dict["Ni"].values()
            for clean, ads_list in entries.items()
            for ads in ads_list
        ]
        assert_allclose(points, expected)
        assert [txt.get_text() for txt in ax.texts] == ["0.25", "0.25"]
        plt.close("all")

    def test_BE_vs_clean_SE_options(self):
        ax = self.Ni_analyzer.BE_vs_clean_SE({}, plot_eads=True, annotate_monolayer=False)
        assert ax.get_ylabel() == "Adsorption Energy (eV)"
        assert len(ax.texts) == 0
        plt.close("all")

    def test_surface_chempot_range_map(self):
        delu_O = Symbol("delu_O")
        ax = self.Ni_analyzer.surface_chempot_range_map(["Ni", "O"], (1, 1, 1), [[-1, 0], [-8, 0]], incr=3)
        assert isinstance(ax, plt.Axes)
        assert ax.get_xlim() == approx((-1, 0))
        assert ax.get_ylim() == approx((-8, 0))
        assert ax.get_xlabel() == r"$\Delta\mu_{Ni} (eV)$"
        assert ax.get_ylabel() == r"$\Delta\mu_{O} (eV)$"
        assert sorted(txt.get_text() for txt in ax.texts) == ["Ni_111_O", "Ni_111_clean"]

        # delu_Ni does not enter the surface energies, so the clean/adsorbed
        # boundary is a horizontal line at the equilibrium delu_O.
        clean = next(iter(self.metals_O_entry_dict["Ni"][1, 1, 1]))
        ads = self.metals_O_entry_dict["Ni"][1, 1, 1][clean][0]
        u_eq = float(self.Ni_analyzer.get_surface_equilibrium([clean, ads])[delu_O])
        horizontal = [line for line in ax.get_lines() if line.get_ydata()[0] == line.get_ydata()[-1]]
        assert any(line.get_ydata()[0] == approx(u_eq) for line in horizontal)
        # The O-adsorbed surface becomes unphysical (gamma < 0) near delu_O = 0
        dashed = [line for line in ax.get_lines() if line.get_linestyle() == "--"]
        assert len(dashed) == 3
        plt.close("all")

        _, ax = plt.subplots()
        ax_out = self.Ni_analyzer.surface_chempot_range_map(
            ["Ni", "O"],
            (1, 1, 1),
            [[-1, 0], [-8, 0]],
            incr=3,
            ax=ax,
            annotate=False,
            show_unphysical_only=True,
        )
        assert ax_out is ax
        assert len(ax.texts) == 0
        # Only the dashed unphysical-region markers are drawn
        assert all(line.get_linestyle() == "--" for line in ax.get_lines())
        plt.close("all")

    def test_get_surface_equilibrium_no_solution(self):
        clean = next(iter(self.metals_O_entry_dict["Pt"][1, 1, 1]))
        ads = self.metals_O_entry_dict["Pt"][1, 1, 1][clean][0]
        # Fixing delu_O leaves 2 equations with 1 unknown (gamma)
        with pytest.warns(UserWarning, match="No solution"):
            soln = self.Pt_analyzer.get_surface_equilibrium([clean, ads], delu_dict={Symbol("delu_O"): -3})
        assert not soln

    def test_stable_u_range_dict_values(self):
        delu_O = Symbol("delu_O")
        entries = self.metals_O_entry_dict["Ni"]
        clean111 = next(iter(entries[1, 1, 1]))
        ads111 = entries[1, 1, 1][clean111][0]
        u_eq = float(self.Ni_analyzer.get_surface_equilibrium([clean111, ads111])[delu_O])
        assert u_eq == approx(-5.771551, abs=1e-5)

        u_range, se_dict = self.Ni_analyzer.stable_u_range_dict(
            [0, -8], delu_O, no_doped=False, miller_index=(1, 1, 1), return_se_dict=True
        )
        assert set(u_range) == {clean111, ads111}
        assert_allclose(np.array(u_range[clean111], dtype=float), [-8, u_eq])
        assert_allclose(np.array(u_range[ads111], dtype=float), [u_eq, 0])
        gamma_clean = self.Ni_analyzer.surfe_dict[clean111]
        assert_allclose(np.array(se_dict[clean111], dtype=float), [gamma_clean, gamma_clean])
        assert float(se_dict[ads111][0]) == approx(gamma_clean)

        # dmu_at_0 inserts the chempot at which the adsorbed surface energy vanishes
        u_range0, se_dict0 = self.Ni_analyzer.stable_u_range_dict(
            [-8, 0], delu_O, no_doped=False, miller_index=(1, 1, 1), dmu_at_0=True, return_se_dict=True
        )
        assert len(u_range0[ads111]) == 3
        assert float(se_dict0[ads111][1]) == 0
        u_zero = float(u_range0[ads111][1])
        assert u_zero == approx(-0.681264, abs=1e-5)
        coeffs = self.Ni_analyzer.as_coeffs_dict[ads111]
        assert float(coeffs[1] + coeffs[delu_O] * u_zero) == approx(0, abs=1e-10)
        assert len(u_range0[clean111]) == 2

    def test_stable_u_range_dict_single_entry(self):
        # With no_doped, each facet only has one (clean) entry, which is stable
        # over the entire chempot range with a constant surface energy.
        delu_O = Symbol("delu_O")
        u_range, se_dict = self.Ni_analyzer.stable_u_range_dict([0, -5], delu_O, return_se_dict=True)
        assert len(u_range) == 2
        for entry, urange in u_range.items():
            assert urange == [-5, 0]
            assert_allclose(se_dict[entry], [self.Ni_analyzer.surfe_dict[entry]] * 2)

    def test_set_all_variables(self):
        delu_O = Symbol("delu_O")
        assert self.Ni_analyzer.list_of_chempots == [1, delu_O]
        assert self.Ni_analyzer.set_all_variables({delu_O: -3}, 0) == {1: 1, delu_O: -3}
        assert self.Ni_analyzer.set_all_variables(None, -2) == {1: 1, delu_O: -2}
        assert self.Cu_analyzer.list_of_chempots == [1]

    def test_non_stoichiometric(self):
        with open(f"{TEST_DIR}/ucell_entries.txt", encoding="utf-8") as file:
            ucell_entries = orjson.loads(file.read())
        MgO_ucell = ComputedStructureEntry.from_dict(ucell_entries["MgO"])
        Mg_ucell = ComputedStructureEntry.from_dict(ucell_entries["Mg"])
        entry_dict = get_entry_dict(f"{TEST_DIR}/MgO_slab_entries.txt")
        analyzer = SurfaceEnergyPlotter(entry_dict, MgO_ucell, ref_entries=[Mg_ucell])
        delu_Mg = Symbol("delu_Mg")
        assert analyzer.list_of_chempots == [1, delu_Mg]

        # Mg-rich and O-rich terminations have opposite slopes w.r.t. delu_Mg
        Mg5O6, Mg6O5 = entry_dict[1, 1, 1]
        assert Mg5O6.composition.reduced_formula == "Mg5O6"
        assert analyzer.as_coeffs_dict[Mg5O6][delu_Mg] == approx(0.0642299605)
        assert analyzer.as_coeffs_dict[Mg6O5][delu_Mg] == approx(-0.0642299605)

        # Dashed line for slabs that are non-stoichiometric w.r.t. the bulk
        _, ax = plt.subplots()
        ax = analyzer.chempot_vs_gamma_plot_one(ax, Mg5O6, delu_Mg, [-1, 0])
        assert ax.get_lines()[0].get_linestyle() == "--"
        plt.close("all")

    def test_sub_chempots(self):
        delu_O = Symbol("delu_O")
        assert sub_chempots({1: 0.5, delu_O: -0.1}, {delu_O: -2}) == approx(0.7)
        assert sub_chempots({1: 0.5, delu_O: -0.1}, {1: 1, delu_O: 0}) == approx(0.5)


class TestWorkFunctionAnalyzer(MatSciTest):
    def setup_method(self):
        self.kwargs = {
            "poscar_filename": f"{TEST_DIR}/CONTCAR.relax1.gz",
            "locpot_filename": f"{TEST_DIR}/LOCPOT.gz",
            "outcar_filename": f"{TEST_DIR}/OUTCAR.relax1.gz",
        }
        self.wf_analyzer = WorkFunctionAnalyzer.from_files(**self.kwargs)

    def test_shift(self):
        wf_analyzer_shift = WorkFunctionAnalyzer.from_files(shift=-0.25, blength=3.7, **self.kwargs)
        assert self.wf_analyzer.ave_bulk_p == approx(wf_analyzer_shift.ave_bulk_p, abs=1e-1)

    def test_is_converged(self):
        assert self.wf_analyzer.is_converged()
        # A vanishing tolerance over a wide window can never be satisfied
        assert not self.wf_analyzer.is_converged(min_points_frac=0.2, tol=0)

    def test_attributes(self):
        wf = self.wf_analyzer
        assert wf.efermi == approx(0.0603)
        assert wf.vacuum_locpot == approx(3.24599516, abs=1e-6)
        assert wf.work_function == approx(wf.vacuum_locpot - wf.efermi)
        assert wf.work_function == approx(3.18569516, abs=1e-6)
        assert wf.ave_bulk_p == approx(-13.0406972, abs=1e-6)
        assert len(wf.locpot_along_c) == len(wf.along_c) == 672
        assert wf.ave_locpot == approx((wf.vacuum_locpot - min(wf.locpot_along_c)) / 2)
        assert len(wf.slab_regions) == 1
        assert wf.slab_regions[0] == approx((0.028937, 0.471063), abs=1e-6)
        frac_z = [site.frac_coords[2] for site in wf.sorted_sites]
        assert frac_z == sorted(frac_z)

    @pytest.mark.parametrize(("shift", "expected"), [(1.25, 0.25), (-0.25, 0.75), (0.5, 0.5)])
    def test_shift_normalization(self, shift, expected):
        wf = WorkFunctionAnalyzer.from_files(shift=shift, **self.kwargs)
        assert wf.shift == approx(expected)
        # The work function is independent of the shift
        assert wf.work_function == approx(self.wf_analyzer.work_function)

    @pytest.mark.parametrize("label_energies", [True, False])
    def test_get_locpot_along_slab_plot(self, label_energies):
        wf = self.wf_analyzer
        ax = wf.get_locpot_along_slab_plot(label_energies=label_energies)
        assert isinstance(ax, plt.Axes)
        assert ax.get_xlim() == approx((0, 1))
        assert ax.get_ylim() == approx((min(wf.locpot_along_c), wf.vacuum_locpot + wf.ave_locpot * 0.2))
        assert ax.get_xlabel() == r"Fractional coordinates ($\hat{c}$)"
        assert ax.get_ylabel() == "Potential (eV)"
        # Raw signal and the averaged signal
        assert_allclose(ax.get_lines()[0].get_ydata(), wf.locpot_along_c)
        averaged = ax.get_lines()[1].get_ydata()
        assert min(averaged) == approx(wf.ave_bulk_p)
        assert max(averaged) == approx(wf.vacuum_locpot)
        if label_energies:
            texts = [txt.get_text() for txt in ax.texts]
            assert texts == [
                f"$V_{{vac}}={wf.vacuum_locpot:.2f}$",
                f"$E_F={wf.efermi:.2f}$",
                f"$V^{{interior}}_{{slab}}={wf.ave_bulk_p:.2f}$",
                rf"$\Phi={wf.work_function:.2f}$",
            ]
            # Single slab region with more vacuum above it: label work function there
            label_in_vac = (1 + wf.slab_regions[0][1]) / 2
            assert ax.texts[3].xy[0] == approx(label_in_vac)
        else:
            assert len(ax.texts) == 0
        plt.close("all")

    @pytest.mark.parametrize(
        ("slab_regions", "label_in_bulk", "label_in_vac"),
        [
            ([(0.0, 0.2), (0.7, 1.0)], 0.85, 0.45),
            ([(0.0, 0.4), (0.7, 1.0)], 0.2, 0.55),
            ([(0.6, 0.95)], 0.775, 0.3),
        ],
    )
    def test_get_labels_positions(self, slab_regions, label_in_bulk, label_in_vac):
        wf = self.wf_analyzer
        wf.slab_regions = slab_regions
        _, ax = plt.subplots()
        ax = wf.get_labels(ax)
        v_vac, e_f, v_bulk, phi = ax.texts
        assert v_vac.xy[0] == approx(label_in_bulk)
        assert e_f.xy[0] == approx(label_in_bulk)
        assert v_bulk.xy[0] == approx(label_in_vac)
        assert phi.xy[0] == approx(label_in_vac)
        assert phi.xy[1] == approx(wf.efermi + wf.ave_locpot * 0.05)
        plt.close("all")

    def test_get_locpot_along_slab_plot_two_regions(self):
        wf = self.wf_analyzer
        wf.slab_regions = [(0.0, 0.2), (0.7, 1.0)]
        ax = wf.get_locpot_along_slab_plot(plt=plt)
        averaged = ax.get_lines()[1].get_ydata()
        along_c = ax.get_lines()[1].get_xdata()
        # Everything at or below the first region and at or above the
        # second is flagged as slab interior
        in_slab = (along_c <= 0.2) | (along_c >= 0.7)
        assert_allclose(averaged[in_slab], wf.ave_bulk_p)
        plt.close("all")


class TestNanoscaleStability(MatSciTest):
    def setup_method(self):
        # Load all entries
        La_hcp_entry_dict = get_entry_dict(f"{TEST_DIR}/La_hcp_entries.txt")
        La_fcc_entry_dict = get_entry_dict(f"{TEST_DIR}/La_fcc_entries.txt")
        with open(f"{TEST_DIR}/ucell_entries.txt", encoding="utf-8") as txt_file:
            ucell_entries = orjson.loads(txt_file.read())
        La_hcp_ucell_entry = ComputedStructureEntry.from_dict(ucell_entries["La_hcp"])
        La_fcc_ucell_entry = ComputedStructureEntry.from_dict(ucell_entries["La_fcc"])

        # Set up the NanoscaleStabilityClass
        self.La_hcp_analyzer = SurfaceEnergyPlotter(La_hcp_entry_dict, La_hcp_ucell_entry)
        self.La_fcc_analyzer = SurfaceEnergyPlotter(La_fcc_entry_dict, La_fcc_ucell_entry)
        self.nanoscale_stability = NanoscaleStability([self.La_fcc_analyzer, self.La_hcp_analyzer])

    def test_stability_at_r(self):
        # Check that we have a different polymorph that is
        # stable below or above the equilibrium particle size
        radius = self.nanoscale_stability.solve_equilibrium_point(self.La_hcp_analyzer, self.La_fcc_analyzer) * 10

        # hcp phase of La particle should be the stable
        # polymorph above the equilibrium radius
        hcp_wulff = self.La_hcp_analyzer.wulff_from_chempot()
        bulk = self.La_hcp_analyzer.ucell_entry
        ghcp, _rhcp = self.nanoscale_stability.wulff_gform_and_r(hcp_wulff, bulk, radius + 10, from_sphere_area=True)
        fcc_wulff = self.La_fcc_analyzer.wulff_from_chempot()
        bulk = self.La_fcc_analyzer.ucell_entry
        gfcc, _rfcc = self.nanoscale_stability.wulff_gform_and_r(fcc_wulff, bulk, radius + 10, from_sphere_area=True)
        assert gfcc > ghcp

        # fcc phase of La particle should be the stable
        # polymorph below the equilibrium radius
        hcp_wulff = self.La_hcp_analyzer.wulff_from_chempot()
        bulk = self.La_hcp_analyzer.ucell_entry
        ghcp, _rhcp = self.nanoscale_stability.wulff_gform_and_r(hcp_wulff, bulk, radius - 10, from_sphere_area=True)
        fcc_wulff = self.La_fcc_analyzer.wulff_from_chempot()
        bulk = self.La_fcc_analyzer.ucell_entry
        gfcc, _rfcc = self.nanoscale_stability.wulff_gform_and_r(fcc_wulff, bulk, radius - 10, from_sphere_area=True)
        assert gfcc < ghcp

    def test_scaled_wulff(self):
        # Ensure for a given radius, the effective radius
        # of the Wulff shape is the same (correctly scaled)
        hcp_wulff = self.La_hcp_analyzer.wulff_from_chempot()
        fcc_wulff = self.La_fcc_analyzer.wulff_from_chempot()
        w1 = self.nanoscale_stability.scaled_wulff(hcp_wulff, 10)
        w2 = self.nanoscale_stability.scaled_wulff(fcc_wulff, 10)
        assert w1.effective_radius == approx(w2.effective_radius)
        assert w1.effective_radius == approx(10)
        assert w2.effective_radius == approx(10)

    def test_plot(self):
        ax = self.nanoscale_stability.plot_all_stability_map(max_r=10.0, labels=["La_hcp", "La_fcc"])

        assert isinstance(ax, plt.Axes)
        assert "Particle radius" in ax.get_xlabel()
        assert "$G_{form}$" in ax.get_ylabel()
        assert [line.get_label() for line in ax.get_lines()] == ["La_hcp", "La_fcc"]
        assert all(len(line.get_xdata()) == 50 for line in ax.get_lines())
        plt.close("all")

    def test_plot_one_stability_map(self):
        ax = self.nanoscale_stability.plot_one_stability_map(
            self.La_hcp_analyzer, 20, label="hcp", increments=5, r_units="Angstrom", e_units="eV"
        )
        assert ax.get_xlabel() == r"Particle radius ($\AA$)"
        assert ax.get_ylabel() == r"$G_{form}$ ($eV/\AA^3$)"
        (line,) = ax.get_lines()
        assert line.get_label() == "hcp"
        wulff = self.La_hcp_analyzer.wulff_from_chempot()
        gform, radius = self.nanoscale_stability.wulff_gform_and_r(
            wulff, self.La_hcp_analyzer.ucell_entry, 20, r_units="Angstrom", e_units="eV"
        )
        assert line.get_xdata()[-1] == approx(radius)
        assert line.get_ydata()[-1] == approx(gform)
        plt.close("all")

    def test_bulk_gform(self):
        bulk = self.La_hcp_analyzer.ucell_entry
        assert NanoscaleStability.bulk_gform(bulk) == approx(bulk.energy / bulk.structure.volume)

    def test_solve_equilibrium_point_units(self):
        r_nm = self.nanoscale_stability.solve_equilibrium_point(self.La_hcp_analyzer, self.La_fcc_analyzer)
        r_ang = self.nanoscale_stability.solve_equilibrium_point(
            self.La_hcp_analyzer, self.La_fcc_analyzer, units="Angstrom"
        )
        assert r_ang == approx(10 * r_nm)
        # Swapping the polymorphs gives the same equilibrium radius
        assert self.nanoscale_stability.solve_equilibrium_point(self.La_fcc_analyzer, self.La_hcp_analyzer) == approx(
            r_nm
        )

    def test_wulff_gform_and_r(self):
        wulff = self.La_hcp_analyzer.wulff_from_chempot()
        bulk = self.La_hcp_analyzer.ucell_entry
        radius = 30
        calc = self.nanoscale_stability.wulff_gform_and_r

        # Wulff-shape based estimate
        e_kev, r_nm = calc(wulff, bulk, radius)
        scaled = self.nanoscale_stability.scaled_wulff(wulff, radius)
        surf = sum(wulff.miller_energy_dict[hkl] * area for hkl, area in scaled.miller_area_dict.items())
        expected = NanoscaleStability.bulk_gform(bulk) * scaled.volume + surf
        assert e_kev == approx(expected / 1000)
        assert r_nm == approx(scaled.effective_radius / 10)

        # Spherical particle estimate
        e_sphere, r_sphere = calc(wulff, bulk, radius, from_sphere_area=True, r_units="Angstrom", e_units="eV")
        expected = (
            NanoscaleStability.bulk_gform(bulk) * 4 / 3 * np.pi * radius**3
            + wulff.weighted_surface_energy * 4 * np.pi * radius**2
        )
        assert e_sphere == approx(expected)
        assert r_sphere == radius

        # Normalizations
        e_norm, _ = calc(wulff, bulk, radius, from_sphere_area=True, r_units="Angstrom", e_units="eV", normalize=True)
        assert e_norm == approx(expected / (4 / 3 * np.pi * radius**3))
        e_atom, _ = calc(
            wulff, bulk, radius, from_sphere_area=True, r_units="Angstrom", e_units="eV", scale_per_atom=True
        )
        density = len(bulk.structure) / bulk.structure.volume
        assert e_atom == approx(expected / (density * 4 / 3 * np.pi * radius**3))


def get_entry_dict(filename):
    # helper to generate an entry_dict

    entry_dict = {}
    with open(filename, encoding="utf-8") as file:
        entries = orjson.loads(file.read())
    for entry in entries:
        sub_str = entry[25:]
        miller_index = []
        for idx, char in enumerate(sub_str):
            if char == "_":
                break
            if char == "-":
                continue
            t = int(char)
            if sub_str[idx - 1] == "-":
                t *= -1
            miller_index.append(t)
        hkl = tuple(miller_index)
        if hkl not in entry_dict:
            entry_dict[hkl] = {}
        entry = ComputedStructureEntry.from_dict(entries[entry])
        entry_dict[hkl][SlabEntry(entry.structure, entry.energy, hkl, label=entry)] = []

    return entry_dict


def load_O_adsorption():
    # Loads the dictionary for clean and O adsorbed Rh, Pt, and Ni entries

    # Load the adsorbate as an entry
    with open(f"{TEST_DIR}/isolated_O_entry.txt", encoding="utf-8") as file:
        isolated_O_entry = orjson.loads(file.read())
    O_entry = ComputedStructureEntry.from_dict(isolated_O_entry)

    # entry_dict for the adsorption case, O adsorption on Ni, Rh and Pt
    metals_O_entry_dict = {
        "Ni": {(1, 1, 1): {}, (1, 0, 0): {}},
        "Pt": {(1, 1, 1): {}},
        "Rh": {(1, 0, 0): {}},
    }

    with open(f"{TEST_DIR}/cs_entries_slabs.json", encoding="utf-8") as file:
        entries = orjson.loads(file.read())
    for key in entries:
        entry = ComputedStructureEntry.from_dict(entries[key])
        for el in metals_O_entry_dict:
            if el in key:
                if "111" in key:
                    clean = SlabEntry(entry.structure, entry.energy, (1, 1, 1), label=f"{key}_clean")
                    metals_O_entry_dict[el][1, 1, 1][clean] = []
                if "110" in key:
                    clean = SlabEntry(entry.structure, entry.energy, (1, 1, 0), label=f"{key}_clean")
                    metals_O_entry_dict[el][1, 1, 0][clean] = []
                if "100" in key:
                    clean = SlabEntry(entry.structure, entry.energy, (1, 0, 0), label=f"{key}_clean")
                    metals_O_entry_dict[el][1, 0, 0][clean] = []

    with open(f"{TEST_DIR}/cs_entries_o_ads.json", encoding="utf-8") as file:
        entries = orjson.loads(file.read())
    for key in entries:
        entry = ComputedStructureEntry.from_dict(entries[key])
        for el, val in metals_O_entry_dict.items():
            if el in key:
                if "111" in key:
                    clean = next(iter(val[1, 1, 1]))
                    ads = SlabEntry(
                        entry.structure,
                        entry.energy,
                        (1, 1, 1),
                        label=f"{key}_O",
                        adsorbates=[O_entry],
                        clean_entry=clean,
                    )
                    metals_O_entry_dict[el][1, 1, 1][clean] = [ads]
                if "110" in key:
                    clean = next(iter(val[1, 1, 0]))
                    ads = SlabEntry(
                        entry.structure,
                        entry.energy,
                        (1, 1, 0),
                        label=f"{key}_O",
                        adsorbates=[O_entry],
                        clean_entry=clean,
                    )
                    metals_O_entry_dict[el][1, 1, 0][clean] = [ads]
                if "100" in key:
                    clean = next(iter(val[1, 0, 0]))
                    ads = SlabEntry(
                        entry.structure,
                        entry.energy,
                        (1, 0, 0),
                        label=f"{key}_O",
                        adsorbates=[O_entry],
                        clean_entry=clean,
                    )
                    metals_O_entry_dict[el][1, 0, 0][clean] = [ads]

    return metals_O_entry_dict
