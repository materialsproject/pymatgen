from __future__ import annotations

import warnings
from collections import Counter

import numpy as np
import pandas as pd
import pytest
from pytest import approx
from scipy.constants import physical_constants

from pymatgen.analysis.magnetism.heisenberg import HeisenbergMapper, HeisenbergModel
from pymatgen.core import Lattice
from pymatgen.core.structure import Structure

# Taken from scipy rather than copied out of HeisenbergMapper.get_mft_temperature, so a
# wrong constant in the implementation shows up here instead of cancelling out.
K_BOLTZMANN = physical_constants["Boltzmann constant in eV/K"][0] * 1000  # meV/K


def _wyckoff_multiplicity(symbol: str) -> int:
    """Leading integer of a wyckoff symbol, e.g. '2c' -> 2."""
    return int("".join(char for char in symbol if char.isdigit()))


class TestHeisenbergMapperKnownHamiltonian:
    """Round-trip against a Hamiltonian whose exchange constants are known.

    Two magnetic species sit on alternating columns of a rectangular lattice::

        A  B  A  B      columns alternate along x, spaced D_AB
        |  |  |  |      chains run along y, spaced D_AA
        A  B  A  B

    giving three distinct interactions (A-A, B-B, A-B) without ever passing a cutoff.
    The chains are deliberately spaced much wider than the columns, so the two
    intra-sublattice interactions are nowhere near a site's overall shortest bond:
    they only survive because the neighbor search takes the nearest shell of each
    sublattice pair separately. The orderings below live in 1x1, 1x2 and 2x2
    supercells of the parent, so the fit must recover the same constants from cells
    of 2, 4 and 8 sites.

    Total energies are assigned from E = N * e0 - sum_<ij> J_ab s_i s_j, evaluated
    over explicit lattice vectors rather than over the mapper's own coupling graph,
    so the topology is ground truth here and not something the test borrows back
    from the code under test.
    """

    E0 = -5.0  # eV per magnetic ion
    J_AB, J_AA, J_BB = 0.011, 0.004, -0.006  # eV, all distinct
    A, B = "Mn", "Fe"

    # The chain spacing is 60% wider than the column spacing, far outside the 10% window
    # a single global nearest-neighbor distance would allow. It must stay below 2 * D_AB,
    # though, or the chain bond stops being the shortest A-A one.
    D_AB = 1.0  # A-B spacing along x
    D_AA = 1.6  # A-A and B-B spacing along y

    # (A spins, B spins) per ordering, indexed [y][x] in units of the parent cell.
    ORDERINGS = (
        ([[1]], [[1]]),  # FM, 1x1
        ([[1]], [[-1]]),  # A up / B down, 1x1
        ([[1], [-1]], [[1], [1]]),  # A chains antialigned, B FM, 1x2
        ([[1, 1], [1, 1]], [[1, -1], [-1, 1]]),  # A FM, B checkerboard, 2x2
    )

    @classmethod
    def _structure(cls, spins_a, spins_b, b_x=0.5):
        # b_x: fractional x of B within each column pair; off 0.5 it splits the A-B bond.
        spins_a = np.atleast_2d(np.asarray(spins_a, dtype=float))
        spins_b = np.atleast_2d(np.asarray(spins_b, dtype=float))
        n_y, n_x = spins_a.shape
        lattice = Lattice.from_parameters(2 * cls.D_AB * n_x, cls.D_AA * n_y, 10, 90, 90, 90)
        species, coords, magmoms = [], [], []
        for y in range(n_y):
            for x in range(n_x):
                species += [cls.A, cls.B]
                coords += [[x / n_x, y / n_y, 0.5], [(x + b_x) / n_x, y / n_y, 0.5]]
                magmoms += [spins_a[y, x], spins_b[y, x]]
        return Structure(lattice, species, coords, site_properties={"magmom": magmoms})

    def _energy(self, spins_a, spins_b):
        spins_a = np.atleast_2d(np.asarray(spins_a, dtype=float))
        spins_b = np.atleast_2d(np.asarray(spins_b, dtype=float))
        n_y, n_x = spins_a.shape
        e_ex = 0.0
        for y in range(n_y):
            for x in range(n_x):
                up, down = (y + 1) % n_y, (y - 1) % n_y
                # Chains along y couple a site to its own species.
                e_ex -= 0.5 * self.J_AA * spins_a[y, x] * (spins_a[up, x] + spins_a[down, x])
                e_ex -= 0.5 * self.J_BB * spins_b[y, x] * (spins_b[up, x] + spins_b[down, x])
                # Along x each A sits between the B of its own cell and the B of the cell
                # to its left; each B likewise between two A. 1/2: bonds counted twice.
                e_ex -= 0.5 * self.J_AB * spins_a[y, x] * (spins_b[y, x] + spins_b[y, (x - 1) % n_x])
                e_ex -= 0.5 * self.J_AB * spins_b[y, x] * (spins_a[y, x] + spins_a[y, (x + 1) % n_x])
        return 2 * n_x * n_y * self.E0 + e_ex

    def _mapper(self):
        structures = [self._structure(*spins) for spins in self.ORDERINGS]
        energies = [self._energy(*spins) for spins in self.ORDERINGS]
        return HeisenbergMapper(structures, energies)

    @staticmethod
    def _labels(hm):
        """J labels of the A-A, B-B and A-B interactions, keyed off species.

        Which orbit is labelled 0 and which 1 is an artifact of the symmetry analysis,
        so the assertions below never assume an order.
        """
        sub_ids = {
            site.specie.symbol: sub_id
            for site, sub_id in zip(hm.parent.magnetic_structure, hm.parent.sublattice_ids, strict=True)
        }
        a, b = sub_ids[TestHeisenbergMapperKnownHamiltonian.A], sub_ids[TestHeisenbergMapperKnownHamiltonian.B]
        return f"{a}-{a}-nn", f"{b}-{b}-nn", f"{min(a, b)}-{max(a, b)}-nn"

    def test_physical_exchange_recovery(self):
        hm = self._mapper()
        assert len(hm.sublattice_wyckoff_symbols) == 2  # Mn and Fe are distinct sublattices
        # The orderings really do live in different-sized supercells of the parent.
        assert sorted(len(struct) for struct in hm.magnetic_structures) == [2, 2, 4, 8]

        aa, bb, ab = self._labels(hm)
        ex_params, residual = hm.get_exchange()

        # One constant per sublattice pair, each recovered independently rather than
        # averaged together, and independent of the supercell each ordering lives in.
        assert set(ex_params) == {"E0", aa, bb, ab}
        assert ex_params[ab] == approx(self.J_AB * 1000, abs=1e-6)
        assert ex_params[aa] == approx(self.J_AA * 1000, abs=1e-6)
        assert ex_params[bb] == approx(self.J_BB * 1000, abs=1e-6)
        assert ex_params["E0"] == approx(self.E0, abs=1e-8)

        # The energies are exactly Heisenberg here, so the fit reproduces them.
        assert residual == approx(0, abs=1e-9)  # meV per ion; float noise only

    def test_residual_is_rms_energy_error_per_ion(self):
        # Four parameters, so the four orderings above make a square system the fit solves
        # exactly whatever the energies are. A fifth ordering plus an energy pushed off the
        # Heisenberg surface is what leaves a residual to look at.
        orderings = [*self.ORDERINGS, ([[1, 1], [1, -1]], [[1, 1], [1, 1]])]  # A stripe/B FM, 2x2
        structures = [self._structure(*spins) for spins in orderings]
        energies = [self._energy(*spins) for spins in orderings]
        energies[0] += 0.05  # eV, on the 2-ion FM cell
        hm = HeisenbergMapper(structures, energies)

        ex_params, residual = hm.get_exchange()
        assert residual > 0

        # The residual is a root mean square over the orderings, in meV per magnetic ion:
        # averaging rather than summing is what keeps it from growing with the number of
        # orderings, and makes it comparable between materials.
        j_columns = [col for col in hm.ex_mat.columns if col != "E"]
        H = hm.ex_mat[j_columns].to_numpy(dtype=float)
        E = hm.ex_mat["E"].to_numpy(dtype=float)
        # ex_params reports E0 in eV and the J_ij in meV; ex_mat is all eV per ion.
        params = np.array([ex_params["E0"], *(ex_params[col] / 1000 for col in j_columns[1:])])
        assert residual == approx(np.sqrt(np.mean((H @ params - E) ** 2)) * 1000, rel=1e-9)

    def test_repeated_hamiltonian_rows_are_kept(self):
        # Reversing every moment leaves each product s_i.s_j untouched, so this ordering
        # repeats an earlier row of the Hamiltonian while its energy sits somewhere else.
        # That gap is the Heisenberg model failing to tell the two orderings apart, and
        # the residual is what reports it - dropping the repeated row would hide it.
        flipped = ([[-1]], [[1]])  # ORDERINGS[1] with every moment reversed
        orderings = [*self.ORDERINGS, flipped]
        structures = [self._structure(*spins) for spins in orderings]
        energies = [self._energy(*spins) for spins in orderings]
        energies[-1] += 0.05  # eV of non-Heisenberg energy, on the 2-ion cell

        hm = HeisenbergMapper(structures, energies)
        assert len(hm.ex_mat) == len(orderings)  # every ordering has a row
        assert hm.ex_mat.drop(columns="E").round(10).duplicated().any()  # and one is a repeat

        assert hm.get_exchange()[1] > 0

    def test_unconstrained_interactions_are_dropped_with_warning(self):
        # Three orderings constrain E0 and two J_ij, one short of the three sublattice pairs
        # the parent has. Even without a cutoff, the nearest-neighbor coupling of one whole
        # pair has to go, and the caller must hear of it.
        orderings = self.ORDERINGS[:3]
        structures = [self._structure(*spins) for spins in orderings]
        energies = [self._energy(*spins) for spins in orderings]
        with pytest.warns(UserWarning, match="orderings constrain only"):
            hm = HeisenbergMapper(structures, energies)

        aa, bb, ab = self._labels(hm)
        j_columns = [col for col in hm.ex_mat.columns if col not in ("E", "E0")]
        assert len(j_columns) == 2
        assert ab in j_columns  # the short A-B bond survives; one of the 1.6 A chains is cut
        assert len({aa, bb} & set(j_columns)) == 1

    def test_degenerate_orderings_are_dropped(self):
        # The same state in a 1x1 and in a 2x2 cell has the same energy per magnetic ion,
        # so screening keeps only one of the two. Per-ion normalization is what makes
        # cells of different size comparable, and hence what makes them degenerate here.
        orderings = [*self.ORDERINGS, ([[1, 1], [1, 1]], [[1, 1], [1, 1]])]  # the FM ordering again, 2x2
        structures = [self._structure(*spins) for spins in orderings]
        energies = [self._energy(*spins) for spins in orderings]
        hm = HeisenbergMapper(structures, energies)

        assert len(hm.orderings) == len(self.ORDERINGS)  # the duplicate was dropped
        # sorted by energy per magnetic ion, ground state first; the total energies of cells
        # of different size need not follow that order
        assert hm.energies_per_magnetic_ion == sorted(hm.energies_per_magnetic_ion)

    def test_sublattices_follow_parent_wyckoff_orbits(self):
        hm = self._mapper()
        # A and B occupy two inequivalent orbits of the parent cell, and every ordering
        # draws its labels from those orbits.
        assert len(hm.sublattice_wyckoff_symbols) == 2
        assert {sub_id for sub_ids in hm.sublattice_ids for sub_id in sub_ids} == set(hm.sublattice_wyckoff_symbols)

        for ordering in hm.orderings:
            # Labels are indexed against the magnetic-only cell the graph is built from,
            # so no None placeholders survive from the parent's full-cell labelling.
            assert len(ordering.sublattice_ids) == len(ordering.magnetic_structure)

            # Each sublattice is filled in proportion to its wyckoff multiplicity, i.e.
            # every ordering is a whole number of parent cells. Independent of which
            # orbit happens to be labelled 0.
            counts = Counter(ordering.sublattice_ids)
            n_cells = {
                counts[sub_id] / _wyckoff_multiplicity(symbol)
                for sub_id, symbol in hm.sublattice_wyckoff_symbols.items()
            }
            assert len(n_cells) == 1

    def test_heisenberg_model(self):
        hm = self._mapper()
        hmodel = hm.get_heisenberg_model()

        assert {sp.symbol for struct in hmodel.magnetic_structures for sp in struct.composition} == {self.A, self.B}

        # The fit travels with the model, residual included.
        assert set(hmodel.ex_params) == {"E0", *self._labels(hm)}
        assert hmodel.residual == approx(0, abs=1e-9)

        # Total energies and the per-magnetic-ion energies the fit used both travel with it.
        assert hmodel.energies == hm.energies
        assert hmodel.energies_per_magnetic_ion == approx(
            [e / len(struct) for e, struct in zip(hm.energies, hm.magnetic_structures, strict=True)]
        )

    def test_as_from_dict_round_trip(self):
        # HeisenbergModel must survive repeated MSON round-trips. as_dict()
        # serializes ex_mat with jsanitize (a DataFrame becomes a nested dict),
        # so from_dict() must reconstruct the DataFrame from that dict.
        # https://github.com/materialsproject/pymatgen/issues/4664
        model = self._mapper().get_heisenberg_model()

        model_rt = HeisenbergModel.from_dict(model.as_dict())
        assert isinstance(model_rt.ex_mat, pd.DataFrame)
        assert model_rt.formula == model.formula
        assert model_rt.structures == model.structures
        assert model_rt.magnetic_structures == model.magnetic_structures
        assert model_rt.residual == model.residual
        assert model_rt.energies == model.energies
        assert model_rt.energies_per_magnetic_ion == model.energies_per_magnetic_ion

        # A second round-trip must be a no-op on the exchange matrix.
        model_rt2 = HeisenbergModel.from_dict(model_rt.as_dict())
        assert isinstance(model_rt2.ex_mat, pd.DataFrame)
        pd.testing.assert_frame_equal(model_rt.ex_mat, model_rt2.ex_mat)

    def test_from_dict_rejects_legacy_serialization(self):
        # A dict serialized with HeisenbergModel <0.2 (pre parent-cell-sublattice
        # refactor) has `unique_site_ids`/`wyckoff_ids`/`sgraphs` instead of
        # `sublattice_ids`/`magnetic_structures`, which don't map onto the new
        # parent-cell sublattices. Loading it must fail with a clear explanation,
        # not a bare KeyError.
        dct = self._mapper().get_heisenberg_model().as_dict()
        del dct["sublattice_ids"]
        with pytest.raises(ValueError, match="HeisenbergModel"):
            HeisenbergModel.from_dict(dct)

    def test_shells_are_counted_within_a_sublattice_pair(self):
        hm = self._mapper()
        aa, bb, ab = self._labels(hm)

        # The A-A and B-B chain bonds are the 'nn' of their own pair even though they
        # are longer than the A-B bond, which is the 'nn' of its pair.
        assert hm.dists[ab] == approx(self.D_AB, abs=0.01)
        assert hm.dists[aa] == hm.dists[bb] == approx(self.D_AA, abs=0.01)

    def test_without_cutoff_tol_still_splits_shells(self):
        # Moving B off-center splits the A-B bond into 0.96 and 1.04 Angstrom. Both lie in
        # the nearest shell of the pair, but they differ by more than tol, so they are two
        # interactions even without a cutoff.
        structures = [self._structure(*spins, b_x=0.48) for spins in self.ORDERINGS]
        energies = [self._energy(*spins) for spins in self.ORDERINGS]
        hm = HeisenbergMapper(structures, energies, tol=0.05)

        _aa, _bb, ab = self._labels(hm)
        assert hm.dists[ab] == approx(0.96, abs=1e-6)
        assert hm.dists[ab.replace("-nn", "-nnn")] == approx(1.04, abs=1e-6)

    def test_bond_is_labelled_with_the_shell_it_was_grouped_into(self):
        # B at x=0.45 gives A-B bonds of 0.90, 1.10, 1.84 and 1.94 Angstrom. tol=0.95 groups
        # the first three into one shell and starts the next at 1.94, so 1.84 is nn although
        # it lies closer to the nnn distance.
        parent = self._structure([[1]], [[1]], b_x=0.45)
        structures = [self._structure(*spins, b_x=0.45) for spins in self.ORDERINGS]
        energies = [self._energy(*spins) for spins in self.ORDERINGS]
        hm = HeisenbergMapper(structures, energies, parent=parent, cutoff=1.95, tol=0.95)

        _aa, _bb, ab = self._labels(hm)
        assert hm.dists[ab.replace("-nn", "-nnn")] == approx(np.hypot(1.1, self.D_AA), abs=1e-6)
        i_id, j_id = map(int, ab.split("-")[:2])
        assert hm._interaction_label(i_id, j_id, np.hypot(0.9, self.D_AA)) == ab

    def test_relaxation_does_not_change_the_couplings(self):
        # Every ordering relaxes B off-center, splitting the A-B coupling into 0.9 and
        # 1.1 Angstrom. On the relaxed cells 1.1 falls outside the 10% nearest-shell
        # window and half the A-B couplings would drop out of every row. Built on the
        # parent geometry, the couplings stay those of the parent and the fit is exact.
        parent = self._structure([[1]], [[1]])
        structures = [self._structure(*spins, b_x=0.45) for spins in self.ORDERINGS]
        energies = [self._energy(*spins) for spins in self.ORDERINGS]
        hm = HeisenbergMapper(structures, energies, parent=parent)

        aa, bb, ab = self._labels(hm)
        assert hm.dists[ab] == approx(self.D_AB, abs=1e-6)
        ex_params, residual = hm.get_exchange()
        assert ex_params[ab] == approx(self.J_AB * 1000, abs=1e-6)
        assert ex_params[aa] == approx(self.J_AA * 1000, abs=1e-6)
        assert ex_params[bb] == approx(self.J_BB * 1000, abs=1e-6)
        assert residual == approx(0, abs=1e-9)

    def test_cutoff_groups_bonds_into_shells_by_tol(self):
        # Up to 2.1 Angstrom each pair has two bond lengths: A-B at 1.0 and 1.89, A-A and
        # B-B at 1.6 (chains) and 2.0 (along x). tol decides whether those are two shells.
        structures = [self._structure(*spins) for spins in self.ORDERINGS]
        energies = [self._energy(*spins) for spins in self.ORDERINGS]

        hm = HeisenbergMapper(structures, energies, cutoff=2.1, tol=0.05)
        aa, _bb, ab = self._labels(hm)
        aa_nnn, ab_nnn = aa.replace("-nn", "-nnn"), ab.replace("-nn", "-nnn")
        assert hm.dists[ab_nnn] == approx(np.hypot(self.D_AB, self.D_AA), abs=1e-6)
        assert hm.dists[aa_nnn] == approx(2 * self.D_AB, abs=1e-6)

        # 2.0 - 1.6 = 0.4 is within tol=0.5, so the A-A bonds merge into one shell; the
        # A-B bonds are 0.89 apart and stay separate.
        hm = HeisenbergMapper(structures, energies, cutoff=2.1, tol=0.5)
        assert aa_nnn not in hm.dists
        assert ab_nnn in hm.dists

    def test_interaction_graph_consistent_across_orderings(self):
        # The fitted J_ij must be recoverable from the interaction graph of any
        # ordering, not just the first, even though the orderings live in
        # different-sized supercells.
        hm = self._mapper()
        hm.get_exchange()
        expected = {self.J_AA * 1000, self.J_BB * 1000, self.J_AB * 1000}

        for ordering_index in range(len(self.ORDERINGS)):
            igraph = hm.get_interaction_graph(ordering_index=ordering_index)
            weights = {round(data["weight"], 6) for *_, data in igraph.graph.edges(data=True)}
            assert weights == expected

    def test_incompatible_ordering_raises(self):
        triangular = Structure(
            Lattice.from_parameters(1, 1, 10, 90, 90, 120),  # not a supercell of the parent
            [self.A, self.B],
            [[0, 0, 0.5], [0.5, 0.5, 0.5]],
            site_properties={"magmom": [1.0, 1.0]},
        )
        structures = [self._structure(*self.ORDERINGS[0]), self._structure(*self.ORDERINGS[1]), triangular]
        energies = [self._energy(*self.ORDERINGS[0]), self._energy(*self.ORDERINGS[1]), -5.0]
        with pytest.raises(ValueError, match="parent cell"):
            HeisenbergMapper(structures, energies)

    def test_too_few_orderings_raises_value_error(self):
        # ValueError, not SystemExit: an unusable input must not tear down the
        # interpreter session the mapper is being called from.
        structures = [self._structure(*self.ORDERINGS[0])]
        with pytest.raises(ValueError, match="at least 2 unique orderings"):
            HeisenbergMapper(structures, [self._energy(*self.ORDERINGS[0])])

    def test_non_structure_parent_raises_clear_error(self):
        # parent moved to third position when this class was refactored (was
        # HeisenbergMapper(ordered_structures, energies, cutoff, tol)); a caller still
        # passing cutoff positionally now hands a number to `parent`. That must fail
        # with a pointed error, not deep inside unrelated Structure/symmetry code.
        structures = [self._structure(*spins) for spins in self.ORDERINGS]
        energies = [self._energy(*spins) for spins in self.ORDERINGS]
        with pytest.raises(TypeError, match="parent must be a Structure"):
            HeisenbergMapper(structures, energies, 5.0)

    def test_coupling_off_the_parent_geometry_raises(self):
        # A coupling shorter than any parent shell of its pair, or of a pair the parent does
        # not have, cannot be labelled; that must say so rather than leak a StopIteration.
        hm = self._mapper()
        _aa, _bb, ab = self._labels(hm)
        i_id, j_id = map(int, ab.split("-")[:2])
        with pytest.raises(ValueError, match="No interaction of sublattices"):
            hm._interaction_label(i_id, j_id, hm.dists[ab] / 2)
        with pytest.raises(ValueError, match="No interaction of sublattices"):
            hm._interaction_label(i_id, 99, hm.dists[ab])

        # The coupling graphs carry each ordering's own moments.
        for ordering, graph in zip(hm.orderings, hm.coupling_graphs, strict=True):
            assert graph.structure.site_properties["magmom"] == ordering.magnetic_structure.site_properties["magmom"]

    def test_inferred_parent_warns(self):
        structures = [self._structure(*spins) for spins in self.ORDERINGS]
        energies = [self._energy(*spins) for spins in self.ORDERINGS]

        with pytest.warns(UserWarning, match="No `parent` cell supplied"):
            HeisenbergMapper(structures, energies)

        # An explicit parent is the documented way out, and must stay silent.
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            HeisenbergMapper(structures, energies, parent=self._structure(*self.ORDERINGS[0]))
        assert not [warning for warning in caught if "No `parent` cell supplied" in str(warning.message)]

    def test_ill_conditioned_fit_warns(self):
        # Near-degenerate orderings make H nearly singular, so tiny energy differences
        # blow up into unphysical exchange constants. The mapper must warn rather than
        # hand the numbers back silently. Perturbing one row of the exchange matrix is
        # the only way to reach that state here: the rows are sums of spin products, so
        # no choice of collinear ordering puts two of them ~1e-7 apart.
        hm = self._mapper()

        # The honest fit does not warn.
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            hm.get_exchange()
        assert not [warning for warning in caught if "ill-conditioned" in str(warning.message)]

        ex_mat = hm.ex_mat.copy()
        ex_mat.iloc[1] = ex_mat.iloc[0] + 1e-7
        hm.ex_mat = ex_mat

        with pytest.warns(UserWarning, match="ill-conditioned"):
            hm.get_exchange()


class TestHeisenbergMeanFieldTemperature:
    """Single magnetic sublattice (square lattice, nearest-neighbor coupling only)
    so both the average exchange and the mean-field critical temperature have
    closed forms to check against:

        <J> = z * J           (z = 4 nearest neighbors)
        Tc  = 2 |<J>| / (3 k_B)

    The next-nearest neighbors sit at sqrt(2) a, outside the default
    nearest-neighbor search, so a single interaction survives without a cutoff.
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
                e_ex -= 0.5 * j_nn * nn  # 1/2: bonds counted from both ends
        return n_x * n_y * self.E0 + e_ex

    def _mapper(self, j_nn):
        fm, afm = [[1, 1], [1, 1]], [[1, -1], [-1, 1]]
        structures = [self._structure(fm), self._structure(afm)]
        energies = [self._energy(fm, j_nn), self._energy(afm, j_nn)]
        return HeisenbergMapper(structures, energies)

    def test_single_sublattice_mft_matches_analytic(self):
        j_nn = 0.010  # eV
        hm = self._mapper(j_nn)
        assert len(hm.sublattice_wyckoff_symbols) == 1  # one magnetic sublattice

        j_avg = hm.estimate_exchange()
        assert j_avg == approx(self.Z * j_nn * 1000, abs=1e-6)  # <J> = z * J (meV)

        mft_t = hm.get_mft_temperature(j_avg)
        tc_expected = 2 * abs(self.Z * j_nn * 1000) / 3 / K_BOLTZMANN
        assert mft_t == approx(tc_expected, abs=1e-3)  # ~309.5 K

    def test_single_interaction_cannot_be_fitted(self):
        # One interaction cannot constrain E0 and a J at once, so there is nothing for
        # the least-squares fit to solve and the mapper must say so rather than guess.
        hm = self._mapper(0.010)
        assert len(hm.interactions) == 1

        with pytest.raises(ValueError, match="needs at least 2"):
            hm.get_exchange()

    @pytest.mark.parametrize(("j_nn", "fm_ground_state"), [(0.010, True), (-0.010, False)])
    def test_exchange_sign_follows_ground_state(self, j_nn, fm_ground_state):
        # J > 0 stabilizes FM (<J> > 0); J < 0 stabilizes AFM (<J> < 0).
        hm = self._mapper(j_nn)
        j_avg = hm.estimate_exchange()
        assert bool(j_avg > 0) == fm_ground_state
        assert j_avg == approx(self.Z * j_nn * 1000, abs=1e-6)


class TestHeisenbergMapperFullCellSymmetry:
    """Site equivalence must be read from the full structure. Removing the
    nonmagnetic ions first can raise the apparent site symmetry and wrongly
    merge magnetic sublattices that are actually distinct.
    """

    lattice = Lattice.from_parameters(6, 4, 4, 90, 90, 90)

    def _ordering(self, spins, with_zn=True):
        # Two Fe that look equivalent on their own; an off-center nonmagnetic Zn
        # breaks the symmetry between them.
        species, coords, magmoms = ["Fe", "Fe"], [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]], list(spins)
        if with_zn:
            species, coords, magmoms = [*species, "Zn"], [*coords, [0.18, 0.0, 0.0]], [*magmoms, 0.0]
        return Structure(self.lattice, species, coords, site_properties={"magmom": magmoms})

    def test_nonmagnetic_ions_split_sublattices(self):
        structures = [self._ordering([3, 3]), self._ordering([3, -3])]
        hm = HeisenbergMapper(structures, [-10.0, -9.0])

        # The two Fe are distinct sublattices because of the off-center Zn. They share a
        # wyckoff symbol, so the sublattice ids -- not the symbols -- are what separate
        # them. Which Fe is labelled 0 is arbitrary, only the split is asserted.
        assert len(hm.sublattice_wyckoff_symbols) == 2
        assert all(len(set(sub_ids)) == 2 for sub_ids in hm.sublattice_ids)

    def test_equivalent_sites_share_a_sublattice(self):
        # Control for the test above: with the Zn gone the two Fe genuinely are
        # equivalent and belong to one sublattice. Stripping the nonmagnetic ions
        # before the symmetry analysis would make the case above look like this one.
        structures = [self._ordering([3, 3], with_zn=False), self._ordering([3, -3], with_zn=False)]
        hm = HeisenbergMapper(structures, [-10.0, -9.0])

        assert len(hm.sublattice_wyckoff_symbols) == 1
        assert all(set(sub_ids) == {0} for sub_ids in hm.sublattice_ids)


class TestHeisenbergMapperZeroMomentIon:
    """A magnetic species that relaxes to zero moment throughout one ordering must
    stay on the magnetic lattice there. Magnetic sites are selected by species
    pooled over *all* orderings, so those ions contribute zero terms instead of
    vanishing and leaving that ordering with a different site count, graph
    topology and energy per magnetic ion than its siblings.
    """

    lattice = Lattice.from_parameters(2, 2, 10, 90, 90, 90)
    species = ("Fe", "Mn", "Fe", "Mn")
    coords = ((0, 0, 0.5), (0.5, 0, 0.5), (0, 0.5, 0.5), (0.5, 0.5, 0.5))

    def _ordering(self, magmoms):
        return Structure(
            self.lattice,
            list(self.species),
            [list(coord) for coord in self.coords],
            site_properties={"magmom": list(magmoms)},
        )

    def test_quenched_species_stays_on_the_magnetic_lattice(self):
        # Both Mn relaxed to zero moment in the second ordering, so on its own that
        # ordering looks like an Fe-only compound.
        structures = [self._ordering([3, 2, 3, 2]), self._ordering([3, 0, -3, 0])]
        hm = HeisenbergMapper(structures, [-20.0, -19.0])

        assert hm.parent.magn_species == {"Fe", "Mn"}
        # Without pooling, the second cell would carry two magnetic sites instead of
        # four and its energy per magnetic ion would be off by a factor of two.
        assert [len(struct) for struct in hm.magnetic_structures] == [4, 4]
        assert [len(sub_ids) for sub_ids in hm.sublattice_ids] == [4, 4]
        assert hm.energies == [-20.0, -19.0]
        assert hm.energies_per_magnetic_ion == [-5.0, -4.75]
