from __future__ import annotations

import logging

import numpy as np
import pandas as pd
import pytest
from pytest import approx

from pymatgen.analysis.graphs import StructureGraph
from pymatgen.analysis.local_env import MinimumDistanceNN
from pymatgen.analysis.magnetism.heisenberg import HeisenbergMapper, HeisenbergModel
from pymatgen.core import Lattice
from pymatgen.core.structure import Structure
from tests.testing import TEST_FILES_DIR

TEST_DIR = f"{TEST_FILES_DIR}/analysis/magnetic_orderings"

# Same Boltzmann constant HeisenbergMapper.get_mft_temperature uses (meV/K).
K_BOLTZMANN = 0.0861733


class TestHeisenbergMapper:
    """End-to-end run on real Mn3Al DFT orderings.

    Mn3Al has no independently-known exchange constants -- the mapper is the only
    estimator -- so asserting specific J or Tc values here would be circular.
    These tests assert physical *invariants* instead (geometry, sublattice
    structure, finiteness, MSON round-trips), and that the mapper flags the fit
    as ill-conditioned rather than silently returning unphysical numbers. Exact
    numeric recovery is validated on synthetic known-Hamiltonian models below.
    """

    @classmethod
    def setup_class(cls):
        cls.Mn3Al = pd.read_json(f"{TEST_DIR}/Mn3Al.json")
        cls.compounds = [cls.Mn3Al]

        cls.hms = []
        for c in cls.compounds:
            ordered_structures = [Structure.from_dict(d) for d in c["structure"]]
            epa = list(c["energy_per_atom"])
            energies = [e * len(s) for (e, s) in zip(epa, ordered_structures, strict=True)]

            hm = HeisenbergMapper(ordered_structures, energies, cutoff=5.0, tol=0.02)
            cls.hms.append(hm)

    def test_graphs(self):
        for hm in self.hms:
            assert len(hm.sgraphs) == 7

    def test_sites(self):
        for hm in self.hms:
            # One sublattice label per site, one label list per ordering
            assert len(hm.site_labels) == len(hm.sgraphs)
            # Site 0 of the first ordering belongs to sublattice 0
            assert hm.site_labels[0][0] == 0

            # unique_site_ids groups each ordering's sites by sublattice id
            unique_site_ids = hm.unique_site_ids
            assert len(unique_site_ids) == len(hm.sgraphs)
            assert unique_site_ids[0][0, 1] == 0
            # Sublattice ids are shared across orderings and match the wyckoff labels
            used_ids = {sub for site_map in unique_site_ids for sub in site_map.values()}
            assert used_ids == set(hm.wyckoff_ids)

    def test_nn_interactions(self):
        for hm in self.hms:
            assert len(hm.nn_interactions) == 3
            assert hm.dists["nn"] == approx(2.51)

    def test_exchange_params_are_physical_invariants(self):
        # No ground-truth J for real DFT data: assert structure and finiteness.
        for hm in self.hms:
            ex_params = hm.get_exchange()
            assert set(hm.wyckoff_ids) == {0, 1}  # two Mn sublattices (1a, 2c)
            assert "E0" in ex_params
            assert all(np.isfinite(v) for v in ex_params.values())

    def test_ill_conditioned_fit_warns(self, caplog):
        # The Mn3Al system is over-parameterized: several orderings are
        # near-degenerate, so H is nearly singular (cond ~ 1e6) and the fitted
        # exchange constants are unphysical. The mapper must warn, not hand back
        # garbage silently.
        for hm in self.hms:
            with caplog.at_level(logging.WARNING):
                hm.get_exchange()
            assert "ill-conditioned" in caplog.text
            caplog.clear()

    def test_get_igraph(self):
        for hm in self.hms:
            igraph = hm.get_interaction_graph()
            assert len(igraph) == 6

    def test_heisenberg_model(self):
        for hm in self.hms:
            hmodel = hm.get_heisenberg_model()
            assert hmodel.formula == "Mn3Al"

    def test_as_from_dict_round_trip(self):
        # HeisenbergModel must survive repeated MSON round-trips. as_dict()
        # serializes ex_mat with jsanitize (a DataFrame becomes a nested dict),
        # so from_dict() must reconstruct the DataFrame from that dict.
        # https://github.com/materialsproject/pymatgen/issues/4664
        for hm in self.hms:
            model = hm.get_heisenberg_model()

            model_rt = HeisenbergModel.from_dict(model.as_dict())
            assert isinstance(model_rt.ex_mat, pd.DataFrame)
            assert model_rt.formula == model.formula

            # A second round-trip must be a no-op on the exchange matrix.
            model_rt2 = HeisenbergModel.from_dict(model_rt.as_dict())
            assert isinstance(model_rt2.ex_mat, pd.DataFrame)
            pd.testing.assert_frame_equal(model_rt.ex_mat, model_rt2.ex_mat)


class TestHeisenbergMapperMixedSupercells:
    """Round-trip on a J1-J2 square lattice with orderings living in different
    supercells of one parent cell.

    Total energies are assigned from a known Hamiltonian
    E = N * e0 - sum_<ij> J_ab s_i s_j, so the fit must recover e0 and the
    physical exchange constants independent of the supercell sizes.
    """

    E0 = -5.0  # eV per magnetic ion
    J_NN = 0.012  # eV
    J_NNN = -0.004  # eV
    NN_VECS = ((1, 0), (-1, 0), (0, 1), (0, -1))
    NNN_VECS = ((1, 1), (1, -1), (-1, 1), (-1, -1))

    @staticmethod
    def _structure(spins):
        spins = np.atleast_2d(np.asarray(spins, dtype=float))
        n_y, n_x = spins.shape
        lattice = Lattice.from_parameters(n_x, n_y, 10, 90, 90, 90)
        coords = [[x / n_x, y / n_y, 0.5] for y in range(n_y) for x in range(n_x)]
        magmoms = [spins[y, x] for y in range(n_y) for x in range(n_x)]
        return Structure(lattice, ["Fe"] * (n_x * n_y), coords, site_properties={"magmom": magmoms})

    def _energy(self, spins):
        spins = np.atleast_2d(np.asarray(spins, dtype=float))
        n_y, n_x = spins.shape
        e_ex = 0.0
        for y in range(n_y):
            for x in range(n_x):
                s_i = spins[y, x]
                nn = sum(s_i * spins[(y + dy) % n_y, (x + dx) % n_x] for dx, dy in self.NN_VECS)
                nnn = sum(s_i * spins[(y + dy) % n_y, (x + dx) % n_x] for dx, dy in self.NNN_VECS)
                e_ex -= 0.5 * (self.J_NN * nn + self.J_NNN * nnn)  # 1/2: bonds counted from both ends
        return n_x * n_y * self.E0 + e_ex

    def test_physical_exchange_recovery(self):
        # FM in the 1x1 primitive cell, stripe AFM in 1x2, Neel AFM in 2x2
        all_spins = [[[1]], [[1], [-1]], [[1, -1], [-1, 1]]]
        structures = [self._structure(spins) for spins in all_spins]
        energies = [self._energy(spins) for spins in all_spins]
        hm = HeisenbergMapper(structures, energies, cutoff=1.5, tol=0.02)

        ex_params = hm.get_exchange()
        assert ex_params["E0"] == approx(self.E0, abs=1e-8)
        assert ex_params["0-0-nn"] == approx(self.J_NN * 1000, abs=1e-6)
        assert ex_params["0-0-nnn"] == approx(self.J_NNN * 1000, abs=1e-6)

    def test_interaction_graph_consistent_across_orderings(self):
        # The fitted J_ij must be recoverable from the interaction graph of any
        # ordering, not just the first, even though the orderings live in
        # different-sized supercells.
        all_spins = [[[1]], [[1], [-1]], [[1, -1], [-1, 1]]]
        structures = [self._structure(spins) for spins in all_spins]
        energies = [self._energy(spins) for spins in all_spins]
        hm = HeisenbergMapper(structures, energies, cutoff=1.5, tol=0.02)
        hm.get_exchange()

        for ordering_index in range(len(structures)):
            igraph = hm.get_interaction_graph(ordering_index=ordering_index)
            weights = {round(d["weight"], 6) for *_, d in igraph.graph.edges(data=True)}
            assert weights == {self.J_NN * 1000, self.J_NNN * 1000}

    def test_incompatible_ordering_raises(self):
        fm, stripe = [[1]], [[1], [-1]]
        triangular = Structure(
            Lattice.from_parameters(1, 1, 10, 90, 90, 120),  # not a supercell of the square parent
            ["Fe"],
            [[0, 0, 0.5]],
            site_properties={"magmom": [1.0]},
        )
        structures = [self._structure(fm), self._structure(stripe), triangular]
        energies = [self._energy(fm), self._energy(stripe), -5.0]
        with pytest.raises(ValueError, match="parent cell"):
            HeisenbergMapper(structures, energies, cutoff=1.5, tol=0.02)


class TestHeisenbergMeanFieldTemperature:
    """Single magnetic sublattice (square lattice, nearest-neighbor coupling
    only) so both the average exchange and the mean-field critical temperature
    have closed forms to check against:

        <J> = z * J           (z = 4 nearest neighbors)
        Tc  = 2 |<J>| / (3 k_B)
    """

    E0 = -5.0  # eV per magnetic ion
    Z = 4  # nearest neighbors on the square lattice
    NN_VECS = ((1, 0), (-1, 0), (0, 1), (0, -1))

    @staticmethod
    def _structure(spins):
        spins = np.atleast_2d(np.asarray(spins, dtype=float))
        n_y, n_x = spins.shape
        lattice = Lattice.from_parameters(n_x, n_y, 10, 90, 90, 90)
        coords = [[x / n_x, y / n_y, 0.5] for y in range(n_y) for x in range(n_x)]
        magmoms = [spins[y, x] for y in range(n_y) for x in range(n_x)]
        return Structure(lattice, ["Fe"] * (n_x * n_y), coords, site_properties={"magmom": magmoms})

    def _energy(self, spins, j_nn):
        spins = np.atleast_2d(np.asarray(spins, dtype=float))
        n_y, n_x = spins.shape
        e_ex = 0.0
        for y in range(n_y):
            for x in range(n_x):
                s_i = spins[y, x]
                nn = sum(s_i * spins[(y + dy) % n_y, (x + dx) % n_x] for dx, dy in self.NN_VECS)
                e_ex -= 0.5 * j_nn * nn
        return n_x * n_y * self.E0 + e_ex

    def _mapper(self, j_nn):
        fm, afm = [[1, 1], [1, 1]], [[1, -1], [-1, 1]]
        structures = [self._structure(fm), self._structure(afm)]
        energies = [self._energy(fm, j_nn), self._energy(afm, j_nn)]
        # cutoff between NN (1.0) and NNN (sqrt(2)) keeps a single interaction
        return HeisenbergMapper(structures, energies, cutoff=1.2, tol=0.02)

    def test_single_sublattice_mft_matches_analytic(self):
        j_nn = 0.010  # eV
        hm = self._mapper(j_nn)
        assert len(hm.wyckoff_ids) == 1  # one magnetic sublattice

        j_avg = hm.estimate_exchange()
        assert j_avg == approx(self.Z * j_nn * 1000, abs=1e-6)  # <J> = z * J (meV)

        mft_t = hm.get_mft_temperature(j_avg)
        tc_expected = 2 * abs(self.Z * j_nn * 1000) / 3 / K_BOLTZMANN
        assert mft_t == approx(tc_expected, abs=1e-3)  # ~309.5 K

    @pytest.mark.parametrize(("j_nn", "fm_ground_state"), [(0.010, True), (-0.010, False)])
    def test_exchange_sign_follows_ground_state(self, j_nn, fm_ground_state):
        # J > 0 stabilizes FM (<J> > 0); J < 0 stabilizes AFM (<J> < 0).
        hm = self._mapper(j_nn)
        j_avg = hm.estimate_exchange()
        assert bool(j_avg > 0) == fm_ground_state
        assert j_avg == approx(self.Z * j_nn * 1000, abs=1e-6)


class TestHeisenbergTwoSublatticeCoupling:
    """Two inequivalent magnetic sublattices (Mn, Fe) arranged on a checkerboard
    so the parent cell has two Wyckoff orbits. The mapper must resolve a separate
    exchange constant for *each* sublattice pair (A-B nearest neighbor, A-A and
    B-B next-nearest neighbor) rather than averaging them, and the coupled
    mean-field temperature must follow from those constants.
    """

    E0 = -5.0  # eV per magnetic ion
    J_AB, J_AA, J_BB = 0.011, 0.004, -0.006  # eV, all distinct
    A, B = "Mn", "Fe"
    NN, NNN = 1.0, np.sqrt(2)  # a = sqrt(2): A-B nn = 1, A-A / B-B nnn = sqrt(2)

    # Orderings as (A-sublattice spins, B-sublattice spins) grids, in 1x1 and 2x2
    # supercells of the parent, chosen to constrain all three couplings + E0.
    ORDERINGS = (
        ([[1]], [[1]]),  # FM
        ([[1]], [[-1]]),  # A up, B down
        ([[1, -1], [-1, 1]], [[1, 1], [1, 1]]),  # A Neel, B FM
        ([[1, 1], [1, 1]], [[1, -1], [-1, 1]]),  # A FM, B Neel
        ([[1, -1], [-1, 1]], [[1, -1], [-1, 1]]),  # both Neel
    )

    @classmethod
    def _structure(cls, spins_a, spins_b):
        spins_a = np.atleast_2d(np.asarray(spins_a, dtype=float))
        spins_b = np.atleast_2d(np.asarray(spins_b, dtype=float))
        n_y, n_x = spins_a.shape
        a = np.sqrt(2)
        lattice = Lattice.from_parameters(n_x * a, n_y * a, 10, 90, 90, 90)
        species, coords, magmoms = [], [], []
        for y in range(n_y):
            for x in range(n_x):
                species += [cls.A, cls.B]
                coords += [[x / n_x, y / n_y, 0.5], [(x + 0.5) / n_x, (y + 0.5) / n_y, 0.5]]
                magmoms += [spins_a[y, x], spins_b[y, x]]
        return Structure(lattice, species, coords, site_properties={"magmom": magmoms})

    def _energy(self, struct):
        # Forward Heisenberg energy on the mapper's own neighbor graph: the
        # neighbor topology is a geometric fact, while the coupling assigned to
        # each (sublattice pair, shell) is the physics the mapper must invert.
        graph = StructureGraph.from_local_env_strategy(struct, MinimumDistanceNN(cutoff=1.5, get_all_sites=True))
        j_by_pair = {
            (frozenset((self.A, self.B)), "nn"): self.J_AB,
            (frozenset((self.A,)), "nnn"): self.J_AA,
            (frozenset((self.B,)), "nnn"): self.J_BB,
        }
        mags = struct.site_properties["magmom"]
        e_ex = 0.0
        for i in range(len(struct)):
            for cs in graph.get_connected_sites(i):
                j, dist = cs[2], cs[-1]
                if abs(dist - self.NN) < 0.05:
                    shell = "nn"
                elif abs(dist - self.NNN) < 0.05:
                    shell = "nnn"
                else:
                    continue
                pair = frozenset((struct[i].specie.symbol, struct[j].specie.symbol))
                e_ex -= 0.5 * j_by_pair[pair, shell] * mags[i] * mags[j]
        return len(struct) * self.E0 + e_ex

    def _mapper(self):
        structures = [self._structure(sa, sb) for sa, sb in self.ORDERINGS]
        energies = [self._energy(s) for s in structures]
        return HeisenbergMapper(structures, energies, cutoff=1.5, tol=0.02)

    def test_distinct_sublattice_couplings_resolved(self):
        hm = self._mapper()
        assert len(hm.wyckoff_ids) == 2  # Mn and Fe are distinct sublattices

        ex_params = hm.get_exchange()
        # One constant per sublattice pair, each recovered independently. Compare
        # the (order-independent) multiset of J values against the inputs.
        j_values = sorted(v for k, v in ex_params.items() if k != "E0")
        expected = sorted([self.J_AB * 1000, self.J_AA * 1000, self.J_BB * 1000])
        assert len(j_values) == 3
        assert j_values == approx(expected, abs=1e-6)
        assert ex_params["E0"] == approx(self.E0, abs=1e-8)

    def test_multi_sublattice_mft_matches_analytic(self):
        hm = self._mapper()
        hm.get_exchange()  # populates ex_params used by the multi-sublattice branch

        mft_t = hm.get_mft_temperature(hm.estimate_exchange())

        # get_mft_temperature builds omega = (2/3 k_B) * [[2 Jaa, Jab], [Jab, 2 Jbb]]
        # and takes its largest eigenvalue. That eigenvalue is symmetric in
        # (Jaa, Jbb), so it does not depend on which sublattice is labelled 0/1.
        jab, jaa, jbb = self.J_AB * 1000, self.J_AA * 1000, self.J_BB * 1000
        max_eig = (jaa + jbb) + np.sqrt((jaa - jbb) ** 2 + jab**2)
        tc_expected = 2 / 3 / K_BOLTZMANN * max_eig
        assert mft_t == approx(tc_expected, rel=1e-9)


class TestHeisenbergMapperFullCellSymmetry:
    """Site equivalence must be read from the full structure. Removing the
    nonmagnetic ions first can raise the apparent site symmetry and wrongly
    merge magnetic sublattices that are actually distinct.
    """

    lattice = Lattice.from_parameters(6, 4, 4, 90, 90, 90)

    def _ordering(self, spins, with_x=True):
        # Two Fe that look equivalent on their own; an off-center nonmagnetic Zn
        # breaks the symmetry between them.
        species, coords, magmoms = ["Fe", "Fe"], [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]], list(spins)
        if with_x:
            species, coords, magmoms = [*species, "Zn"], [*coords, [0.18, 0.0, 0.0]], [*magmoms, 0.0]
        return Structure(self.lattice, species, coords, site_properties={"magmom": magmoms})

    def test_nonmagnetic_ions_split_sublattices(self):
        structures = [self._ordering([3, 3]), self._ordering([3, -3])]
        hm = HeisenbergMapper(structures, [-10.0, -9.0], cutoff=3.5, tol=0.02)

        # The two Fe are distinct sublattices because of the off-center Zn.
        assert len(hm.wyckoff_ids) == 2
        assert hm.site_labels == [[0, 1], [0, 1]]

    def test_stripping_nonmagnetic_ions_would_merge_sublattices(self):
        # Without the nonmagnetic Zn the two Fe look equivalent (one sublattice).
        # This is the skewed result the full-cell symmetry analysis avoids.
        structures = [self._ordering([3, 3], with_x=False), self._ordering([3, -3], with_x=False)]
        hm = HeisenbergMapper(structures, [-10.0, -9.0], cutoff=3.5, tol=0.02)

        assert len(hm.wyckoff_ids) == 1
        assert hm.site_labels == [[0, 0], [0, 0]]
