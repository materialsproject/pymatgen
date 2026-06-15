from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from pytest import approx

from pymatgen.analysis.magnetism.heisenberg import HeisenbergMapper
from pymatgen.core import Lattice
from pymatgen.core.structure import Structure
from tests.testing import TEST_FILES_DIR

TEST_DIR = f"{TEST_FILES_DIR}/analysis/magnetic_orderings"


class TestHeisenbergMapper:
    @classmethod
    def setup_class(cls):
        cls.df = pd.read_json(f"{TEST_DIR}/mag_orderings_test_cases.json")

        # Good tests
        cls.Mn3Al = pd.read_json(f"{TEST_DIR}/Mn3Al.json")

        cls.compounds = [cls.Mn3Al]

        cls.hms = []
        for c in cls.compounds:
            ordered_structures = list(c["structure"])
            ordered_structures = [Structure.from_dict(d) for d in ordered_structures]
            epa = list(c["energy_per_atom"])
            energies = [e * len(s) for (e, s) in zip(epa, ordered_structures, strict=True)]

            hm = HeisenbergMapper(ordered_structures, energies, cutoff=5.0, tol=0.02)
            cls.hms.append(hm)

    def test_graphs(self):
        for hm in self.hms:
            struct_graphs = hm.sgraphs
            assert len(struct_graphs) == 7

    def test_sites(self):
        for hm in self.hms:
            unique_site_ids = hm.unique_site_ids
            # One site map per ordering
            assert len(unique_site_ids) == len(hm.sgraphs)
            # Site 0 of the first ordering belongs to global sublattice 0
            assert hm._global_orbit_id(unique_site_ids[0], 0) == 0
            # Global ids are shared across orderings and match the wyckoff labels
            used_ids = {gid for site_map in unique_site_ids for gid in site_map.values()}
            assert used_ids == set(hm.wyckoff_ids)

    def test_nn_interactions(self):
        for hm in self.hms:
            n_interacts = len(hm.nn_interactions)
            assert n_interacts == 3

            dists = hm.dists
            assert dists["nn"] == approx(2.51)

    def test_exchange_params(self):
        for hm in self.hms:
            ex_params = hm.get_exchange()
            assert ex_params["0-1-nn"] == approx(30.813729, abs=1e-2)

    def test_mean_field(self):
        for hm in self.hms:
            j_avg = hm.estimate_exchange()
            value = round(52.54997149705518, 3)
            assert round(j_avg, 3) == value

            mft_t = hm.get_mft_temperature(j_avg)
            assert mft_t == approx(21596.7, abs=1)

    def test_get_igraph(self):
        for hm in self.hms:
            igraph = hm.get_interaction_graph()
            assert len(igraph) == 6

    def test_heisenberg_model(self):
        for hm in self.hms:
            hmodel = hm.get_heisenberg_model()
            assert hmodel.formula == "Mn3Al"


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
