from __future__ import annotations

import logging
import multiprocessing
import unittest
import warnings

import numpy as np
from monty.serialization import dumpfn, loadfn
from pytest import approx

from pymatgen.analysis.pourbaix_diagram import IonEntry, MultiEntry, PourbaixDiagram, PourbaixEntry, PourbaixPlotter
from pymatgen.core.composition import Composition
from pymatgen.core.ion import Ion
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

logger = logging.getLogger(__name__)


class TestPourbaixEntry(PymatgenTest):
    _multiprocess_shared_ = True
    """
    Test all functions using a fictitious entry
    """

    def setUp(self):
        # comp = Composition("Mn2O3")
        self.sol_entry = ComputedEntry("Mn2O3", 49)
        ion = Ion.from_formula("MnO4-")
        self.ion_entry = IonEntry(ion, 25)
        self.PxIon = PourbaixEntry(self.ion_entry)
        self.PxSol = PourbaixEntry(self.sol_entry)
        self.PxIon.concentration = 1e-4

    def test_pourbaix_entry(self):
        assert self.PxIon.entry.energy == 25, "Wrong Energy!"
        assert self.PxIon.entry.name == "MnO4[-1]", "Wrong Entry!"
        assert self.PxSol.entry.energy == 49, "Wrong Energy!"
        assert self.PxSol.entry.name == "Mn2O3", "Wrong Entry!"
        # assert self.PxIon.energy == 25, "Wrong Energy!"
        # assert self.PxSol.energy == 49, "Wrong Energy!"
        assert self.PxIon.concentration == 1e-4, "Wrong concentration!"

    def test_calc_coeff_terms(self):
        assert self.PxIon.npH == -8, "Wrong npH!"
        assert self.PxIon.nPhi == -7, "Wrong nPhi!"
        assert self.PxIon.nH2O == 4, "Wrong nH2O!"

        assert self.PxSol.npH == -6, "Wrong npH!"
        assert self.PxSol.nPhi == -6, "Wrong nPhi!"
        assert self.PxSol.nH2O == 3, "Wrong nH2O!"

    def test_to_from_dict(self):
        d = self.PxIon.as_dict()
        ion_entry = self.PxIon.from_dict(d)
        assert ion_entry.entry.name == "MnO4[-1]", "Wrong Entry!"

        d = self.PxSol.as_dict()
        sol_entry = self.PxSol.from_dict(d)
        assert sol_entry.name == "Mn2O3(s)", "Wrong Entry!"
        assert sol_entry.energy == self.PxSol.energy, "as_dict and from_dict energies unequal"

        # Ensure computed entry data persists
        entry = ComputedEntry("TiO2", energy=-20, data={"test": "test"})
        pbx_entry = PourbaixEntry(entry=entry)
        dumpfn(pbx_entry, "pbx_entry.json")
        reloaded = loadfn("pbx_entry.json")
        assert isinstance(reloaded.entry, ComputedEntry)
        assert reloaded.entry.data is not None

    def test_energy_functions(self):
        # TODO: test these for values
        self.PxSol.energy_at_conditions(10, 0)
        self.PxSol.energy_at_conditions(np.array([1, 2, 3]), 0)
        self.PxSol.energy_at_conditions(10, np.array([1, 2, 3]))
        self.PxSol.energy_at_conditions(np.array([1, 2, 3]), np.array([1, 2, 3]))

    def test_multi_entry(self):
        # TODO: More robust multi-entry test
        m_entry = MultiEntry([self.PxSol, self.PxIon])
        for attr in ["energy", "composition", "nPhi"]:
            assert getattr(m_entry, attr) == getattr(self.PxSol, attr) + getattr(self.PxIon, attr)

        # As dict, from dict
        m_entry_dict = m_entry.as_dict()
        m_entry_new = MultiEntry.from_dict(m_entry_dict)
        assert m_entry_new.energy == m_entry.energy

    def test_get_elt_fraction(self):
        entry = ComputedEntry("Mn2Fe3O3", 49)
        pb_entry = PourbaixEntry(entry)
        assert pb_entry.get_element_fraction("Fe") == approx(0.6)
        assert pb_entry.get_element_fraction("Mn") == approx(0.4)


class TestPourbaixDiagram(unittest.TestCase):
    _multiprocess_shared_ = True

    @classmethod
    def setUpClass(cls):
        cls.test_data = loadfn(f"{TEST_FILES_DIR}/pourbaix_test_data.json")
        cls.pbx = PourbaixDiagram(cls.test_data["Zn"], filter_solids=True)
        cls.pbx_nofilter = PourbaixDiagram(cls.test_data["Zn"], filter_solids=False)

    def test_pourbaix_diagram(self):
        assert {e.name for e in self.pbx.stable_entries} == {
            "ZnO(s)",
            "Zn[2+]",
            "ZnHO2[-]",
            "ZnO2[2-]",
            "Zn(s)",
        }, "List of stable entries does not match"

        assert {e.name for e in self.pbx_nofilter.stable_entries} == {
            "ZnO(s)",
            "Zn[2+]",
            "ZnHO2[-]",
            "ZnO2[2-]",
            "Zn(s)",
            "ZnO2(s)",
            "ZnH(s)",
        }, "List of stable entries for unfiltered pbx does not match"

        pbx_lowconc = PourbaixDiagram(self.test_data["Zn"], conc_dict={"Zn": 1e-8}, filter_solids=True)
        assert {e.name for e in pbx_lowconc.stable_entries} == {
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
        ag_n = [e for e in self.test_data["Ag-Te-N"] if "Te" not in e.composition]
        highconc = PourbaixDiagram(ag_n, filter_solids=True, conc_dict={"Ag": 1e-5, "N": 1})
        entry_sets = [set(e.entry_id) for e in highconc.stable_entries]
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
            PourbaixEntry(IonEntry(Ion.from_formula("Fe[2+]"), -0.7683100214319288), entry_id="ion-0"),
            PourbaixEntry(IonEntry(Ion.from_formula("Li[1+]"), -3.0697590542787156), entry_id="ion-12"),
        ]
        comp_dict = Composition({"Fe": 1, "Ir": 1, "Li": 2, "Si": 1, "V": 2}).fractional_composition

        multi_entry = PourbaixDiagram.process_multientry(entries, prod_comp=comp_dict)
        assert multi_entry is None

    def test_get_pourbaix_domains(self):
        domains = PourbaixDiagram.get_pourbaix_domains(self.test_data["Zn"])
        assert len(domains[0]) == 7

    def test_get_decomposition(self):
        # Test a stable entry to ensure that it's zero in the stable region
        entry = self.test_data["Zn"][12]  # Should correspond to mp-2133
        assert self.pbx.get_decomposition_energy(entry, 10, 1) == approx(
            0.0, 5
        ), "Decomposition energy of ZnO is not 0."

        # Test an unstable entry to ensure that it's never zero
        entry = self.test_data["Zn"][11]
        ph, v = np.meshgrid(np.linspace(0, 14), np.linspace(-2, 4))
        result = self.pbx_nofilter.get_decomposition_energy(entry, ph, v)
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
        custom_ion_entry = PourbaixEntry(ion, entry_id="my_ion")
        pbx = PourbaixDiagram(
            [*entries, custom_ion_entry],
            filter_solids=True,
            comp_dict={"Na": 1, "Sn": 12, "C": 24},
        )
        assert pbx.get_decomposition_energy(custom_ion_entry, 5, 2) == approx(2.1209002582, abs=1e-1)

    def test_get_stable_entry(self):
        entry = self.pbx.get_stable_entry(0, 0)
        assert entry.entry_id == "ion-0"

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
        assert {e.name for e in new.stable_entries} == {
            "ZnO(s)",
            "Zn[2+]",
            "ZnHO2[-]",
            "ZnO2[2-]",
            "Zn(s)",
        }, "List of stable entries does not match"

        # Test with unstable solid entries included (filter_solids=False), this should result in the
        # previously filtered entries being included
        dct = self.pbx_nofilter.as_dict()
        new = PourbaixDiagram.from_dict(dct)
        assert {e.name for e in new.stable_entries} == {
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


class TestPourbaixPlotter(unittest.TestCase):
    def setUp(self):
        warnings.simplefilter("ignore")
        self.test_data = loadfn(f"{TEST_FILES_DIR}/pourbaix_test_data.json")
        self.pd = PourbaixDiagram(self.test_data["Zn"])
        self.plotter = PourbaixPlotter(self.pd)

    def tearDown(self):
        warnings.simplefilter("default")

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
        plt = binary_plotter.plot_entry_stability(self.test_data["Ag-Te"][53])
        plt.close()
