from __future__ import annotations

import logging
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
    over explicit lattice vectors rather than over the mapper's own neighbor graph,
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
    def _structure(cls, spins_a, spins_b):
        spins_a = np.atleast_2d(np.asarray(spins_a, dtype=float))
        spins_b = np.atleast_2d(np.asarray(spins_b, dtype=float))
        n_y, n_x = spins_a.shape
        lattice = Lattice.from_parameters(2 * cls.D_AB * n_x, cls.D_AA * n_y, 10, 90, 90, 90)
        species, coords, magmoms = [], [], []
        for y in range(n_y):
            for x in range(n_x):
                species += [cls.A, cls.B]
                coords += [[x / n_x, y / n_y, 0.5], [(x + 0.5) / n_x, y / n_y, 0.5]]
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
        return HeisenbergMapper(structures, energies, tol=0.02)

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
        assert residual == approx(0, abs=1e-12)

    def test_residual_is_rms_energy_error_per_ion(self):
        # Four parameters, so the four orderings above make a square system the fit solves
        # exactly whatever the energies are. A fifth ordering plus an energy pushed off the
        # Heisenberg surface is what leaves a residual to look at.
        orderings = [*self.ORDERINGS, ([[1, 1], [1, -1]], [[1, 1], [1, 1]])]  # A stripe/B FM, 2x2
        structures = [self._structure(*spins) for spins in orderings]
        energies = [self._energy(*spins) for spins in orderings]
        energies[0] += 0.05  # eV, on the 2-ion FM cell
        hm = HeisenbergMapper(structures, energies, tol=0.02)

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

    def test_degenerate_orderings_are_dropped(self):
        # The same state in a 1x1 and in a 2x2 cell has the same energy per magnetic ion,
        # so screening keeps only one of the two. Per-ion normalization is what makes
        # cells of different size comparable, and hence what makes them degenerate here.
        orderings = [*self.ORDERINGS, ([[1, 1], [1, 1]], [[1, 1], [1, 1]])]  # the FM ordering again, 2x2
        structures = [self._structure(*spins) for spins in orderings]
        energies = [self._energy(*spins) for spins in orderings]
        hm = HeisenbergMapper(structures, energies, tol=0.02)

        assert len(hm.orderings) == len(self.ORDERINGS)  # the duplicate was dropped
        assert hm.energies == sorted(hm.energies)  # sorted, ground state first

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
            assert all(isinstance(sub_id, int) for sub_id in ordering.sublattice_ids)

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

        assert hmodel.formula == "MnFe"
        assert {sp.symbol for struct in hmodel.magnetic_structures for sp in struct.composition} == {self.A, self.B}

        # The fit travels with the model, residual included.
        assert set(hmodel.ex_params) == {"E0", *self._labels(hm)}
        assert hmodel.residual == approx(0, abs=1e-12)

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

        # A second round-trip must be a no-op on the exchange matrix.
        model_rt2 = HeisenbergModel.from_dict(model_rt.as_dict())
        assert isinstance(model_rt2.ex_mat, pd.DataFrame)
        pd.testing.assert_frame_equal(model_rt.ex_mat, model_rt2.ex_mat)

    def test_shells_are_counted_within_a_sublattice_pair(self):
        hm = self._mapper()
        aa, bb, ab = self._labels(hm)

        # The A-A and B-B chain bonds are the 'nn' of their own pair even though they
        # are longer than the A-B bond, which is the 'nn' of its pair.
        assert hm.dists[ab] == approx(self.D_AB, abs=0.01)
        assert hm.dists[aa] == hm.dists[bb] == approx(self.D_AA, abs=0.01)

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
            HeisenbergMapper(structures, energies, tol=0.02)

    def test_ill_conditioned_fit_warns(self, caplog):
        # Near-degenerate orderings make H nearly singular, so tiny energy differences
        # blow up into unphysical exchange constants. The mapper must warn rather than
        # hand the numbers back silently. Perturbing one row of the exchange matrix is
        # the only way to reach that state here: the rows are sums of spin products, so
        # no choice of collinear ordering puts two of them ~1e-7 apart.
        hm = self._mapper()

        with caplog.at_level(logging.WARNING):
            hm.get_exchange()
        assert "ill-conditioned" not in caplog.text  # the honest fit does not warn
        caplog.clear()

        ex_mat = hm.ex_mat.copy()
        ex_mat.iloc[1] = ex_mat.iloc[0] + 1e-7
        hm.ex_mat = ex_mat

        with caplog.at_level(logging.WARNING):
            hm.get_exchange()
        assert "ill-conditioned" in caplog.text

    def test_multi_sublattice_mft_pins_current_formula(self):
        # NOT an independent physical check. get_mft_temperature's multi-sublattice
        # branch uses bare J_ij with no coordination numbers, and adds each diagonal
        # term twice, so it does not reduce to the single-sublattice formula
        # Tc = 2<J>/(3 k_B) that TestHeisenbergMeanFieldTemperature validates. This
        # test pins the matrix the branch currently builds,
        #     omega = (2 / 3 k_B) * [[2 Jaa, Jab], [Jab, 2 Jbb]],
        # so any change to it is deliberate rather than accidental. Its largest
        # eigenvalue is symmetric in (Jaa, Jbb), hence independent of the labelling.
        hm = self._mapper()
        hm.get_exchange()  # populates ex_params used by the multi-sublattice branch

        mft_t = hm.get_mft_temperature(hm.estimate_exchange())

        jab, jaa, jbb = self.J_AB * 1000, self.J_AA * 1000, self.J_BB * 1000
        max_eig = (jaa + jbb) + np.sqrt((jaa - jbb) ** 2 + jab**2)
        assert mft_t == approx(2 / 3 / K_BOLTZMANN * max_eig, rel=1e-5)


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
        return HeisenbergMapper(structures, energies, tol=0.02)

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
        assert len(hm.nn_interactions) == 1

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
        hm = HeisenbergMapper(structures, [-10.0, -9.0], tol=0.02)

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
        hm = HeisenbergMapper(structures, [-10.0, -9.0], tol=0.02)

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
        hm = HeisenbergMapper(structures, [-20.0, -19.0], tol=0.02)

        assert hm.parent.magn_species == {"Fe", "Mn"}
        # Without pooling, the second cell would carry two magnetic sites instead of
        # four and its energy per magnetic ion would be off by a factor of two.
        assert [len(struct) for struct in hm.magnetic_structures] == [4, 4]
        assert [len(sub_ids) for sub_ids in hm.sublattice_ids] == [4, 4]
        assert hm.energies == [-5.0, -4.75]
