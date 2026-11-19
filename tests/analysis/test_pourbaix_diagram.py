from __future__ import annotations

import multiprocessing
from copy import deepcopy

import matplotlib.pyplot as plt
import numpy as np
import pytest
from monty.serialization import dumpfn, loadfn
from pytest import approx

from pymatgen.analysis.pourbaix_diagram import (
    IonEntry,
    MultiEntry,
    PourbaixDiagram,
    PourbaixEntry,
    PourbaixPlotter,
    ion_or_solid_comp_object,
)
from pymatgen.core.composition import Composition
from pymatgen.core.ion import Ion
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.util.testing import MatSciTest
from tests.testing import TEST_FILES_DIR

TEST_DIR = f"{TEST_FILES_DIR}/analysis/pourbaix_diagram"


class TestPourbaixEntry(MatSciTest):
    """Test all functions using a fictitious entry"""

    def setup_method(self):
        # comp = Composition("Mn2O3")
        self.sol_entry = ComputedEntry("Mn2O3", 49)
        ion = Ion.from_formula("MnO4-")
        self.ion_entry = IonEntry(ion, 25)
        self.px_ion = PourbaixEntry(self.ion_entry)
        self.px_sol = PourbaixEntry(self.sol_entry)
        self.px_ion.concentration = 1e-4

    def test_pourbaix_entry(self):
        assert self.px_ion.entry.energy == 25, "Wrong Energy!"
        assert self.px_ion.entry.name == "MnO4[-1]", "Wrong Entry!"
        assert self.px_sol.entry.energy == 49, "Wrong Energy!"
        assert self.px_sol.entry.name == "Mn2O3", "Wrong Entry!"
        # assert self.PxIon.energy == 25, "Wrong Energy!"
        # assert self.PxSol.energy == 49, "Wrong Energy!"
        assert self.px_ion.concentration == approx(1e-4), "Wrong concentration!"

    def test_calc_coeff_terms(self):
        assert self.px_ion.npH == -8, "Wrong npH!"
        assert self.px_ion.nPhi == -7, "Wrong nPhi!"
        assert self.px_ion.nH2O == 4, "Wrong nH2O!"

        assert self.px_sol.npH == -6, "Wrong npH!"
        assert self.px_sol.nPhi == -6, "Wrong nPhi!"
        assert self.px_sol.nH2O == 3, "Wrong nH2O!"

    def test_as_from_dict(self):
        dct = self.px_ion.as_dict()
        ion_entry = self.px_ion.from_dict(dct)
        assert ion_entry.entry.name == "MnO4[-1]", "Wrong Entry!"

        dct = self.px_sol.as_dict()
        sol_entry = self.px_sol.from_dict(dct)
        assert sol_entry.name == "Mn2O3(s)", "Wrong Entry!"
        assert sol_entry.energy == self.px_sol.energy, "as_dict and from_dict energies unequal"

        # Ensure computed entry data persists
        entry = ComputedEntry("TiO2", energy=-20, data={"test": "test"})
        pbx_entry = PourbaixEntry(entry=entry)
        dumpfn(pbx_entry, "pbx_entry.json")
        reloaded = loadfn("pbx_entry.json")
        assert isinstance(reloaded.entry, ComputedEntry)
        assert reloaded.entry.data is not None

    def test_energy_functions(self):
        # TODO: test these for values
        self.px_sol.energy_at_conditions(10, 0)
        self.px_sol.energy_at_conditions(np.array([1, 2, 3]), 0)
        self.px_sol.energy_at_conditions(10, np.array([1, 2, 3]))
        self.px_sol.energy_at_conditions(np.array([1, 2, 3]), np.array([1, 2, 3]))

    def test_multi_entry(self):
        # TODO: More robust multi-entry test
        m_entry = MultiEntry([self.px_sol, self.px_ion])
        for attr in ["energy", "composition", "nPhi"]:
            assert getattr(m_entry, attr) == getattr(self.px_sol, attr) + getattr(self.px_ion, attr)

        # As dict, from dict
        m_entry_dict = m_entry.as_dict()
        m_entry_new = MultiEntry.from_dict(m_entry_dict)
        assert m_entry_new.energy == m_entry.energy

    def test_multi_entry_repr(self):
        m_entry = MultiEntry([self.px_sol, self.px_ion])
        assert (
            repr(m_entry) == "PourbaixMultiEntry(energy=90.9717, npH=-14.0, nPhi=-13.0, nH2O=7.0, "
            "entry_id=[None, None], species='Mn2O3(s) + MnO4[-1]')"
        )

    def test_get_elt_fraction(self):
        entry = ComputedEntry("Mn2Fe3O3", 49)
        pb_entry = PourbaixEntry(entry)
        assert pb_entry.get_element_fraction("Fe") == approx(0.6)
        assert pb_entry.get_element_fraction("Mn") == approx(0.4)

    def test_energy_terms(self):
        # Solids have unit activity, so no concentration term
        assert self.px_sol.concentration == 1
        assert self.px_sol.conc_term == 0
        # E = E0 + conc_term - MU_H2O * nH2O with MU_H2O = -2.4583 eV, PREFAC = 0.0591
        assert self.px_sol.energy == approx(49 + 3 * 2.4583)
        assert self.px_ion.conc_term == approx(0.0591 * -4)
        assert self.px_ion.energy == approx(25 - 0.0591 * 4 + 4 * 2.4583)
        # Normalized by the number of non-H/O atoms (2 Mn in Mn2O3, 1 Mn in MnO4-)
        assert self.px_sol.normalization_factor == approx(0.5)
        assert self.px_sol.normalized_energy == approx(self.px_sol.energy / 2)
        assert self.px_ion.normalized_energy == approx(self.px_ion.energy)
        assert self.px_sol.energy_per_atom == approx(self.px_sol.energy / 5)
        assert self.px_sol.num_atoms == 5

    @pytest.mark.parametrize(
        ("pH", "V"),
        [(10, 0.5), (np.array([0, 7, 14]), 0.5), (7, np.array([-1, 0, 1])), (np.arange(3), np.arange(3))],
    )
    def test_energy_at_conditions(self, pH, V):
        # Mn2O3: npH = nPhi = -6
        expected = self.px_sol.energy - 6 * 0.0591 * np.asarray(pH) - 6 * np.asarray(V)
        assert self.px_sol.energy_at_conditions(pH, V) == approx(expected)
        assert self.px_sol.normalized_energy_at_conditions(pH, V) == approx(expected / 2)

    def test_names_and_repr(self):
        assert self.px_sol.name == self.px_sol.to_pretty_string() == "Mn2O3(s)"
        assert self.px_ion.name == self.px_ion.to_pretty_string() == "MnO4[-1]"
        assert self.px_sol.phase_type == "Solid"
        assert self.px_ion.phase_type == "Ion"
        assert self.px_ion.charge == -1
        assert {el.symbol for el in self.px_ion.elements} == {"Mn", "O"}
        assert repr(self.px_sol) == (
            "PourbaixEntry(Mn2 O3 with energy=56.3749, npH=-6.0, nPhi=-6.0, nH2O=3.0, entry_id=None)"
        )
        assert repr(self.ion_entry) == "IonEntry : Mn1 O4 with energy = 25.0000"

    def test_entry_id(self):
        # entry_id falls back to that of the wrapped entry, unless explicitly given
        entry = ComputedEntry("Mn2O3", 49, entry_id="mp-1172")
        assert PourbaixEntry(entry).entry_id == "mp-1172"
        assert PourbaixEntry(entry, entry_id="custom").entry_id == "custom"
        assert self.px_sol.entry_id is None

    def test_ion_entry_round_trip(self):
        ion_entry = IonEntry(Ion.from_formula("Mn[2+]"), -2.4, name="manganous")
        new = IonEntry.from_dict(ion_entry.as_dict())
        assert new.name == "manganous"
        assert new.energy == approx(-2.4)
        assert new.ion == Ion.from_formula("Mn[2+]")
        assert IonEntry(Ion.from_formula("Mn[2+]"), -2.4).name == "Mn[+2]"

        # Ion PourbaixEntry keeps its concentration through serialization
        pbx_entry = PourbaixEntry(ion_entry, entry_id="ion-mn", concentration=1e-3)
        new_pbx = PourbaixEntry.from_dict(pbx_entry.as_dict())
        assert pbx_entry.as_dict()["entry_type"] == "Ion"
        assert new_pbx.concentration == approx(1e-3)
        assert new_pbx.entry_id == "ion-mn"
        assert new_pbx.energy == approx(pbx_entry.energy)

    def test_multi_entry_weights(self):
        m_entry = MultiEntry([self.px_sol, self.px_ion], weights=[0.5, 2])
        assert m_entry.energy == approx(0.5 * self.px_sol.energy + 2 * self.px_ion.energy)
        assert m_entry.npH == approx(0.5 * -6 + 2 * -8)
        assert m_entry.nH2O == approx(0.5 * 3 + 2 * 4)
        assert m_entry.composition == Composition("Mn3O9.5")
        assert m_entry.phase_type == ["Solid", "Ion"]
        assert m_entry.name == "Mn2O3(s) + MnO4[-1]"

        new = MultiEntry.from_dict(m_entry.as_dict())
        assert new.weights == [0.5, 2]
        assert new.energy == approx(m_entry.energy)
        assert [entry.name for entry in new.entry_list] == ["Mn2O3(s)", "MnO4[-1]"]


class TestPourbaixDiagram:
    @classmethod
    def setup_class(cls):
        cls.test_data = loadfn(f"{TEST_DIR}/pourbaix_test_data.json")
        cls.pbx = PourbaixDiagram(cls.test_data["Zn"], filter_solids=True)
        cls.pbx_no_filter = PourbaixDiagram(cls.test_data["Zn"], filter_solids=False)

    def test_pourbaix_diagram(self):
        assert {entry.name for entry in self.pbx.stable_entries} == {
            "ZnO(s)",
            "Zn[2+]",
            "ZnHO2[-]",
            "ZnO2[2-]",
            "Zn(s)",
        }, "List of stable entries does not match"

        assert {entry.name for entry in self.pbx_no_filter.stable_entries} == {
            "ZnO(s)",
            "Zn[2+]",
            "ZnHO2[-]",
            "ZnO2[2-]",
            "Zn(s)",
            "ZnO2(s)",
            "ZnH(s)",
        }, "List of stable entries for unfiltered pbx does not match"

        pbx_low_conc = PourbaixDiagram(self.test_data["Zn"], conc_dict={"Zn": 1e-8}, filter_solids=True)
        assert {entry.name for entry in pbx_low_conc.stable_entries} == {
            "Zn(HO)2(aq)",
            "Zn[2+]",
            "ZnHO2[-]",
            "ZnO2[2-]",
            "Zn(s)",
        }

    def test_properties(self):
        assert len(self.pbx.unstable_entries) == 2

    def test_multicomponent(self):
        # Assure no ions get filtered at high concentration
        ag_n = [entry for entry in self.test_data["Ag-Te-N"] if "Te" not in entry.composition]
        highconc = PourbaixDiagram(ag_n, filter_solids=True, conc_dict={"Ag": 1e-5, "N": 1})
        entry_sets = [set(entry.entry_id) for entry in highconc.stable_entries]
        assert {"mp-124", "ion-17"} in entry_sets

        # Binary system
        pd_binary = PourbaixDiagram(
            self.test_data["Ag-Te"],
            filter_solids=True,
            comp_dict={"Ag": 0.5, "Te": 0.5},
            conc_dict={"Ag": 1e-8, "Te": 1e-8},
        )
        assert len(pd_binary.stable_entries) == 30
        test_entry = pd_binary.find_stable_entry(8, 2)
        assert "mp-499" in test_entry.entry_id

        # Find a specific multi-entry to test
        assert pd_binary.get_decomposition_energy(test_entry, 8, 2) == 0

        pd_ternary = PourbaixDiagram(self.test_data["Ag-Te-N"], filter_solids=True)
        assert len(pd_ternary.stable_entries) == 49

        # Fetch a solid entry and a ground state entry mixture
        ag_te_n = self.test_data["Ag-Te-N"][-1]
        ground_state_ag_with_ions = MultiEntry(
            [self.test_data["Ag-Te-N"][i] for i in [4, 18, 30]],
            weights=[1 / 3, 1 / 3, 1 / 3],
        )
        assert pd_ternary.get_decomposition_energy(ag_te_n, 2, -1) == approx(2.767822855765)
        assert pd_ternary.get_decomposition_energy(ag_te_n, 10, -2) == approx(3.756840056890625)
        assert pd_ternary.get_decomposition_energy(ground_state_ag_with_ions, 2, -1) == approx(0)

        # Test invocation of Pourbaix diagram from ternary data
        new_ternary = PourbaixDiagram(pd_ternary.all_entries)
        assert len(new_ternary.stable_entries) == 49
        assert new_ternary.get_decomposition_energy(ag_te_n, 2, -1) == approx(2.767822855765)
        assert new_ternary.get_decomposition_energy(ag_te_n, 10, -2) == approx(3.756840056890625)
        assert new_ternary.get_decomposition_energy(ground_state_ag_with_ions, 2, -1) == approx(0)

        # Test processing of multi-entries with degenerate reaction, produced
        # a bug in a prior implementation
        entries = [
            PourbaixEntry(ComputedEntry("VFe2Si", -1.8542253150000008), entry_id="mp-4595"),
            PourbaixEntry(ComputedEntry("Fe", 0), entry_id="mp-13"),
            PourbaixEntry(ComputedEntry("V2Ir2", -2.141851640000006), entry_id="mp-569250"),
            PourbaixEntry(
                IonEntry(Ion.from_formula("Fe[2+]"), -0.7683100214319288),
                entry_id="ion-0",
            ),
            PourbaixEntry(
                IonEntry(Ion.from_formula("Li[1+]"), -3.0697590542787156),
                entry_id="ion-12",
            ),
        ]
        comp_dict = Composition({"Fe": 1, "Ir": 1, "Li": 2, "Si": 1, "V": 2}).fractional_composition

        multi_entry = PourbaixDiagram.process_multientry(entries, prod_comp=comp_dict)
        assert multi_entry is None

    def test_comp_dict_with_h_or_o_is_stripped(self):
        # Regression test: comp_dict containing H or O used to leave _elt_comp
        # with those species, which then mismatched the non-H/O fractional
        # composition checked inside get_decomposition_energy and raised
        # ValueError for entries that should have been valid.
        pbx = PourbaixDiagram(
            self.test_data["Ag-Te"],
            filter_solids=True,
            comp_dict={"Ag": 0.5, "Te": 0.5, "O": 1, "H": 1},
        )
        assert pbx._elt_comp == {"Ag": 0.5, "Te": 0.5}
        entry = pbx.find_stable_entry(8, 2)
        # Must not raise.
        pbx.get_decomposition_energy(entry, 8, 2)

        with pytest.raises(ValueError, match="non-H/O element"):
            PourbaixDiagram(self.test_data["Ag-Te"], comp_dict={"H": 1, "O": 1})

    def test_entry_with_comp_dict_composition(self):
        # Regression test for #4709: a single entry with the same composition as
        # comp_dict must become a MultiEntry, so it can be a stable domain
        cu3p = PourbaixEntry(ComputedEntry("Cu3P", -1.0))
        multi_entry = PourbaixDiagram.process_multientry([cu3p], Composition("Cu3P"))
        assert multi_entry is not None
        assert multi_entry.weights == approx([1])
        # with another entry, that entry would need a coefficient of 0, so no MultiEntry
        cu = PourbaixEntry(ComputedEntry("Cu", 0.0))
        assert PourbaixDiagram.process_multientry([cu3p, cu], Composition("Cu3P")) is None
        assert PourbaixDiagram.process_multientry([cu, cu3p], Composition("Cu3P")) is None

        entries = [
            cu3p,
            PourbaixEntry(ComputedEntry("Cu", 0.0)),
            PourbaixEntry(ComputedEntry("P", 0.0)),
            PourbaixEntry(IonEntry(Ion.from_formula("Cu[2+]"), 65.49 / 96.485)),
            PourbaixEntry(IonEntry(Ion.from_formula("H2PO4[-]"), -1130.28 / 96.485)),
        ]
        pbx = PourbaixDiagram(entries, comp_dict={"Cu": 3, "P": 1})
        assert "Cu3P(s)" in [entry.name for entry in pbx.stable_entries]
        assert pbx.get_decomposition_energy(cu3p, 7, -1) == approx(0)

    def test_get_pourbaix_domains(self):
        domains = PourbaixDiagram.get_pourbaix_domains(self.test_data["Zn"])
        assert len(domains[0]) == 7

    def test_get_decomposition(self):
        # Test a stable entry to ensure that it's zero in the stable region
        entry = self.test_data["Zn"][12]  # Should correspond to mp-2133
        assert self.pbx.get_decomposition_energy(entry, 10, 1) == approx(0.0, 5), (
            "Decomposition energy of ZnO is not 0."
        )

        # Test an unstable entry to ensure that it's never zero
        entry = self.test_data["Zn"][11]
        ph, v = np.meshgrid(np.linspace(0, 14), np.linspace(-2, 4))
        result = self.pbx_no_filter.get_decomposition_energy(entry, ph, v)
        assert (result >= 0).all(), "Unstable energy has hull energy of 0 or less"

        # Test an unstable hydride to ensure HER correction works
        assert self.pbx.get_decomposition_energy(entry, -3, -2) == approx(3.6979147983333)
        # Test a list of pHs
        self.pbx.get_decomposition_energy(entry, np.linspace(0, 2, 5), 2)

        # Test a list of Vs
        self.pbx.get_decomposition_energy(entry, 4, np.linspace(-3, 3, 10))

        # Test a set of matching arrays
        ph, v = np.meshgrid(np.linspace(0, 14), np.linspace(-3, 3))
        self.pbx.get_decomposition_energy(entry, ph, v)

        # Test custom ions
        entries = self.test_data["C-Na-Sn"]
        ion = IonEntry(Ion.from_formula("NaO28H80Sn12C24+"), -161.676)
        custom_ion_entry = PourbaixEntry(ion, entry_id="some_ion")
        pbx = PourbaixDiagram(
            [*entries, custom_ion_entry],
            filter_solids=True,
            comp_dict={"Na": 1, "Sn": 12, "C": 24},
        )
        assert pbx.get_decomposition_energy(custom_ion_entry, 5, 2) == approx(2.1209002582, abs=1e-1)

    def test_get_stable_entry(self):
        entry = self.pbx.get_stable_entry(0, 0)
        assert entry.entry_id == "ion-0"

    @pytest.mark.parametrize(
        ("pH", "V", "name", "hull_energy"),
        [
            (0, 0, "Zn[2+]", -1.835004511978826),
            (7, -1.5, "Zn(s)", 0),
            (14, 0, "ZnO2[2-]", -2.687444511978826),
        ],
    )
    def test_stable_entry_and_hull_energy(self, pH, V, name, hull_energy):
        assert self.pbx.get_stable_entry(pH, V).name == name
        assert self.pbx.find_stable_entry(pH, V).name == name
        assert self.pbx.get_hull_energy(pH, V) == approx(hull_energy)
        # The stable entry sits on the hull, i.e. has zero decomposition energy
        stable = self.pbx.get_stable_entry(pH, V)
        assert stable.normalized_energy_at_conditions(pH, V) == approx(hull_energy)
        assert self.pbx.get_decomposition_energy(stable, pH, V) == approx(0, abs=1e-10)

    def test_get_hull_energy_vectorized(self):
        expected = [-1.835004511978826, 0, -2.687444511978826]
        assert self.pbx.get_hull_energy(np.array([0, 7, 14]), np.array([0, -1.5, 0])) == approx(expected)
        # Plain lists are accepted as documented
        assert self.pbx.get_hull_energy([0, 7, 14], [0, -1.5, 0]) == approx(expected)
        # Scalar broadcasting against an array
        assert self.pbx.get_hull_energy(0, np.array([0, 0])) == approx([expected[0]] * 2)

    def test_get_decomposition_energy_arrays(self):
        entry = self.test_data["Zn"][11]
        assert entry.name == "ZnO2(s)"
        ph_grid, v_grid = np.meshgrid(np.linspace(0, 14, 4), np.linspace(-2, 2, 3))
        result = self.pbx_no_filter.get_decomposition_energy(entry, ph_grid, v_grid)
        assert result.shape == (3, 4)
        expected = [
            self.pbx_no_filter.get_decomposition_energy(entry, ph, v)
            for ph, v in zip(ph_grid.ravel(), v_grid.ravel(), strict=True)
        ]
        assert result.ravel() == approx(expected)
        # ZnO2(s) is on the unfiltered hull only at high pH and potential
        assert (result >= -1e-10).all()
        assert result[0, 0] == approx(expected[0])
        assert result.min() == approx(0, abs=1e-10)

        # Entry whose non-H/O composition differs from the diagram's
        ag_entry = PourbaixEntry(ComputedEntry("Ag", 0))
        with pytest.raises(ValueError, match="Composition of stability entry does not match"):
            self.pbx.get_decomposition_energy(ag_entry, 7, 0)

    def test_entry_lists(self):
        # filter_solids drops the unstable ZnO2 and ZnH solids (and all but the
        # lowest-energy polymorphs) but keeps all supplied entries as unprocessed
        assert len(self.pbx.unprocessed_entries) == len(self.test_data["Zn"]) == 24
        assert len(self.pbx_no_filter.unprocessed_entries) == 24
        assert {entry.name for entry in self.pbx.all_entries} == {
            "Zn(s)",
            "ZnO(s)",
            "Zn[2+]",
            "ZnO2[2-]",
            "ZnHO[+]",
            "ZnHO2[-]",
            "Zn(HO)2(aq)",
        }
        assert set(self.pbx.stable_entries).isdisjoint(self.pbx.unstable_entries)
        assert len(self.pbx.stable_entries) + len(self.pbx.unstable_entries) == len(self.pbx.all_entries)
        assert self.pbx._multi_element is False

    def test_invalid_entries(self):
        entries = deepcopy(self.test_data["Zn"])
        entries[5].phase_type = "Gas"
        with pytest.raises(ValueError, match='phase type of either "Solid" or "Ion"'):
            PourbaixDiagram(entries)

        # Multi-element ion without a concentration cannot be assigned one from conc_dict
        multi_ion = PourbaixEntry(IonEntry(Ion.from_formula("AgTeO3[-]"), -2), concentration=0)
        with pytest.raises(ValueError, match="not compatible with multi-element ions"):
            PourbaixDiagram([*self.test_data["Ag-Te"], multi_ion])

    def test_generate_multielement_entries(self):
        # Brute-force multi-entry generation must contain every multi-entry kept by
        # the hull-based preprocessing used in __init__
        pbx = PourbaixDiagram(self.test_data["Ag-Te"], comp_dict={"Ag": 0.5, "Te": 0.5})
        assert pbx._multi_element
        generated = pbx._generate_multielement_entries(pbx._filtered_entries)
        assert len(generated) == 153
        assert all(isinstance(entry, MultiEntry) for entry in generated)
        generated_ids = {frozenset(entry.entry_id) for entry in generated}
        assert {frozenset(entry.entry_id) for entry in pbx.all_entries} <= generated_ids
        # Each multi-entry has the target 1:1 Ag:Te composition
        for entry in generated[:10]:
            comp = entry.composition
            assert comp["Ag"] == approx(comp["Te"])

        generated_parallel = pbx._generate_multielement_entries(pbx._filtered_entries, nproc=2)
        assert {frozenset(entry.entry_id) for entry in generated_parallel} == generated_ids

    def test_multielement_parallel(self):
        # Simple test to ensure that multiprocessing is working
        test_entries = self.test_data["Ag-Te-N"]
        nproc = multiprocessing.cpu_count()
        pbx = PourbaixDiagram(test_entries, filter_solids=True, nproc=nproc)
        assert len(pbx.stable_entries) == 49

    def test_solid_filter(self):
        entries = self.test_data["Zn"]
        pbx = PourbaixDiagram(entries, filter_solids=False)
        oxidized_phase = pbx.find_stable_entry(10, 2)
        assert oxidized_phase.name == "ZnO2(s)"

        entries = self.test_data["Zn"]
        pbx = PourbaixDiagram(entries, filter_solids=True)
        oxidized_phase = pbx.find_stable_entry(10, 2)
        assert oxidized_phase.name == "ZnO(s)"

    def test_serialization(self):
        dct = self.pbx.as_dict()
        new = PourbaixDiagram.from_dict(dct)
        assert {entry.name for entry in new.stable_entries} == {
            "ZnO(s)",
            "Zn[2+]",
            "ZnHO2[-]",
            "ZnO2[2-]",
            "Zn(s)",
        }, "List of stable entries does not match"

        # Test with unstable solid entries included (filter_solids=False), this should result in the
        # previously filtered entries being included
        dct = self.pbx_no_filter.as_dict()
        new = PourbaixDiagram.from_dict(dct)
        assert {entry.name for entry in new.stable_entries} == {
            "ZnO(s)",
            "Zn[2+]",
            "ZnHO2[-]",
            "ZnO2[2-]",
            "Zn(s)",
            "ZnO2(s)",
            "ZnH(s)",
        }, "List of stable entries for unfiltered pbx does not match"

        pd_binary = PourbaixDiagram(
            self.test_data["Ag-Te"],
            filter_solids=True,
            comp_dict={"Ag": 0.5, "Te": 0.5},
            conc_dict={"Ag": 1e-8, "Te": 1e-8},
        )
        new_binary = PourbaixDiagram.from_dict(pd_binary.as_dict())
        assert len(pd_binary.stable_entries) == len(new_binary.stable_entries)


class TestPourbaixPlotter:
    def setup_method(self):
        self.test_data = loadfn(f"{TEST_DIR}/pourbaix_test_data.json")
        self.pd = PourbaixDiagram(self.test_data["Zn"])
        self.plotter = PourbaixPlotter(self.pd)

    def test_plot_pourbaix(self):
        plotter = PourbaixPlotter(self.pd)
        # Default limits
        plotter.get_pourbaix_plot()
        # Non-standard limits
        plotter.get_pourbaix_plot(limits=[[-5, 4], [-2, 2]])

    def test_plot_entry_stability(self):
        entry = self.pd.all_entries[0]
        self.plotter.plot_entry_stability(entry, limits=[[-2, 14], [-3, 3]])

        # binary system
        pd_binary = PourbaixDiagram(self.test_data["Ag-Te"], comp_dict={"Ag": 0.5, "Te": 0.5})
        binary_plotter = PourbaixPlotter(pd_binary)
        ax = binary_plotter.plot_entry_stability(self.test_data["Ag-Te"][53])
        assert isinstance(ax, plt.Axes)

    @pytest.fixture(autouse=True)
    def _close_figures(self):
        yield
        plt.close("all")

    def test_pourbaix_plot_content(self):
        ax = self.plotter.get_pourbaix_plot(title="Zn", label_fontsize=12)
        assert ax.get_title() == "Zn"
        assert ax.get_xlabel() == "pH"
        assert ax.get_ylabel() == "E (V)"
        assert ax.get_xlim() == approx((-2, 16))
        assert ax.get_ylim() == approx((-3, 3))
        # 2 water lines + 2 neutral axes + one closed outline per stable domain
        n_domains = len(self.pd.stable_entries)
        assert n_domains == 5
        assert len(ax.lines) == 4 + n_domains
        h_line = ax.lines[0]
        assert h_line.get_linestyle() == "--"
        assert h_line.get_ydata() == approx([2 * 0.0591, -16 * 0.0591])
        assert ax.lines[1].get_ydata() == approx([2 * 0.0591 + 1.23, -16 * 0.0591 + 1.23])
        labels = {text.get_text() for text in ax.texts}
        assert len(labels) == n_domains
        assert "Zn$^{2+}$" in labels
        assert all(text.get_fontsize() == 12 for text in ax.texts)

        # Domain outlines are closed polygons of the stored vertices
        entry = self.pd.get_stable_entry(0, 0)
        vertices = self.plotter.domain_vertices(entry)
        assert len(vertices) >= 3
        outline = ax.lines[4 + list(self.pd._stable_domain_vertices).index(entry)]
        assert outline.get_xdata() == approx([*np.asarray(vertices)[:, 0], vertices[0][0]])

    def test_pourbaix_plot_no_extras(self):
        _, ax = plt.subplots()
        out = self.plotter.get_pourbaix_plot(
            limits=((0, 14), (-1, 1)), label_domains=False, show_water_lines=False, show_neutral_axes=False, ax=ax
        )
        assert out is ax
        assert len(ax.lines) == len(self.pd.stable_entries)
        assert len(ax.texts) == 0
        assert ax.get_xlim() == approx((0, 14))

    def test_plot_entry_stability_content(self):
        entry = self.test_data["Zn"][11]  # ZnO2(s)
        ax = self.plotter.plot_entry_stability(entry, pH_resolution=5, V_resolution=4, e_hull_max=2)
        mesh = ax.collections[-1]
        stability = np.asarray(mesh.get_array())
        assert stability.size == 5 * 4
        assert mesh.get_clim() == approx((0, 2))
        pH, V = np.mgrid[-2:16:5j, -3:3:4j]
        assert np.sort(stability.ravel()) == approx(np.sort(self.pd.get_decomposition_energy(entry, pH, V).ravel()))
        cbar_label = ax.figure.axes[-1].get_ylabel()
        assert cbar_label == "Stability of ZnO$_{2}$(s) (eV/atom)"

    def test_show(self, monkeypatch):
        shown = []
        monkeypatch.setattr(plt, "show", lambda: shown.append(plt.gca()))
        self.plotter.show(title="Zn")
        assert len(shown) == 1
        assert shown[0].get_title() == "Zn"


class TestIonOrSolidCompObject:
    def test_ion(self):
        # Test cations
        assert ion_or_solid_comp_object("Li+").charge == 1
        assert ion_or_solid_comp_object("Li[+]").charge == 1
        assert ion_or_solid_comp_object("Ca[2+]").charge == 2
        assert ion_or_solid_comp_object("Ca[+2]").charge == 2
        assert ion_or_solid_comp_object("Ca++").charge == 2
        assert ion_or_solid_comp_object("Ca[++]").charge == 2
        assert ion_or_solid_comp_object("Ca2+").charge == 1
        assert ion_or_solid_comp_object("C2O4-2").charge == -2

        # Test anions
        assert ion_or_solid_comp_object("Cl-").charge == -1
        assert ion_or_solid_comp_object("Cl[-]").charge == -1
        assert ion_or_solid_comp_object("SO4[-2]").charge == -2
        assert ion_or_solid_comp_object("SO4-2").charge == -2
        assert ion_or_solid_comp_object("SO42-").charge == -1
        assert ion_or_solid_comp_object("SO4--").charge == -2
        assert ion_or_solid_comp_object("SO4[--]").charge == -2
        assert ion_or_solid_comp_object("N3-").charge == -1

    def test_solid(self):
        # Test end with "(s)"
        assert ion_or_solid_comp_object("Fe2O3(s)") == Composition("Fe2O3")
        assert ion_or_solid_comp_object("Fe(s)") == Composition("Fe1")

        # Test end without "(s)"
        assert type(ion_or_solid_comp_object("Na2O")) is Composition
