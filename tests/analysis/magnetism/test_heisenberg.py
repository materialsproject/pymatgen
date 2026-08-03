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
from tests.testing import TEST_FILES_DIR

TEST_DIR = f"{TEST_FILES_DIR}/analysis/magnetic_orderings"

# Taken from scipy rather than copied out of HeisenbergMapper.get_mft_temperature, so a
# wrong constant in the implementation shows up here instead of cancelling out.
K_BOLTZMANN = physical_constants["Boltzmann constant in eV/K"][0] * 1000  # meV/K


def _wyckoff_multiplicity(symbol: str) -> int:
    """Leading integer of a wyckoff symbol, e.g. '2c' -> 2."""
    return int("".join(char for char in symbol if char.isdigit()))


class TestHeisenbergMapper:
    """End-to-end run on real Mn3Al DFT orderings.

    Mn3Al has no independently-known exchange constants -- the mapper is the only
    estimator -- so asserting specific J or Tc values here would be circular. These
    tests assert physical invariants instead (sublattice structure, geometry, the
    <J> fallback, MSON round-trips). Exact numeric recovery is validated against a
    known Hamiltonian in TestHeisenbergMapperKnownHamiltonian below.

    No cutoff is passed, so each site keeps only its nearest-neighbor shell. In
    Mn3Al that shell holds a single sublattice pair, which leaves the mapper with
    too few interactions to invert and sends it down the <J> path.
    """

    @classmethod
    def setup_class(cls):
        cls.Mn3Al = pd.read_json(f"{TEST_DIR}/Mn3Al.json")
        cls.structures = [Structure.from_dict(dct) for dct in cls.Mn3Al["structure"]]
        cls.energies = [
            epa * len(struct) for epa, struct in zip(cls.Mn3Al["energy_per_atom"], cls.structures, strict=True)
        ]
        cls.hm = HeisenbergMapper(cls.structures, cls.energies, tol=0.02)

    def test_degenerate_orderings_are_dropped(self):
        # Different initial configs relax to the same state and show up as orderings
        # sharing an energy per magnetic ion; only one of each survives screening.
        n_magnetic = [sum(site.specie.symbol == "Mn" for site in struct) for struct in self.structures]
        energies_per_ion = [energy / n for energy, n in zip(self.energies, n_magnetic, strict=True)]

        assert len(self.hm.orderings) < len(self.structures)  # something was dropped
        assert len(self.hm.orderings) == len({round(energy, 6) for energy in energies_per_ion})
        assert self.hm.energies == sorted(self.hm.energies)  # sorted, ground state first

    def test_sublattices_follow_parent_wyckoff_orbits(self):
        hm = self.hm
        # Mn occupies two inequivalent orbits of the parent cell.
        assert sorted(hm.sublattice_wyckoff_symbols.values()) == ["1a", "2c"]

        # Every ordering draws its labels from those orbits, and every orbit is used.
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

    def test_single_interaction_falls_back_to_javg(self):
        hm = self.hm
        # The nearest-neighbor shell of Mn3Al holds exactly one sublattice pair.
        [(pair, labels)] = hm.nn_interactions.items()
        assert set(pair) == set(hm.sublattice_wyckoff_symbols)  # the two orbits are coupled
        assert labels == [f"{pair[0]}-{pair[1]}-nn"]
        assert hm.dists[labels[0]] == approx(2.51, abs=0.01)

        # One interaction cannot constrain E0 and a J at once, so the mapper reports the
        # average exchange from E_AFM - E_FM rather than inverting a degenerate system.
        ex_params = hm.get_exchange()
        assert set(ex_params) == {"<J>"}
        assert np.isfinite(ex_params["<J>"])

    def test_interaction_graph_carries_javg_on_every_bond(self):
        hm = self.hm
        igraph = hm.get_interaction_graph()

        # In the <J> fallback every nearest-neighbor bond carries the same constant, and
        # none are dropped: the interaction graph covers the whole neighbor graph.
        assert igraph.graph.number_of_edges() == hm.nn_graphs[0].graph.number_of_edges()
        weights = {round(data["weight"], 6) for *_, data in igraph.graph.edges(data=True)}
        assert weights == {round(hm.get_exchange()["<J>"], 6)}

    def test_heisenberg_model(self):
        hmodel = self.hm.get_heisenberg_model()
        assert hmodel.formula == "Mn3Al"

        # structures keep every ion; magnetic_structures keep only the magnetic species,
        # and those are what the graphs and sublattice labels are indexed against.
        assert {sp.symbol for struct in hmodel.magnetic_structures for sp in struct.composition} == {"Mn"}
        assert all(
            len(mag_struct) < len(struct)
            for struct, mag_struct in zip(hmodel.structures, hmodel.magnetic_structures, strict=True)
        )

    def test_as_from_dict_round_trip(self):
        # HeisenbergModel must survive repeated MSON round-trips. as_dict()
        # serializes ex_mat with jsanitize (a DataFrame becomes a nested dict),
        # so from_dict() must reconstruct the DataFrame from that dict.
        # https://github.com/materialsproject/pymatgen/issues/4664
        model = self.hm.get_heisenberg_model()

        model_rt = HeisenbergModel.from_dict(model.as_dict())
        assert isinstance(model_rt.ex_mat, pd.DataFrame)
        assert model_rt.formula == model.formula
        assert model_rt.structures == model.structures
        assert model_rt.magnetic_structures == model.magnetic_structures

        # A second round-trip must be a no-op on the exchange matrix.
        model_rt2 = HeisenbergModel.from_dict(model_rt.as_dict())
        assert isinstance(model_rt2.ex_mat, pd.DataFrame)
        pd.testing.assert_frame_equal(model_rt.ex_mat, model_rt2.ex_mat)


class TestHeisenbergMapperKnownHamiltonian:
    """Round-trip against a Hamiltonian whose exchange constants are known.

    Two magnetic species sit on alternating columns of a rectangular lattice::

        A  B  A  B      columns alternate along x, spaced D_AB
        |  |  |  |      chains run along y, spaced D_AA
        A  B  A  B

    so a site's nearest-neighbor shell holds both its own-species chain neighbors
    and its cross-species column neighbors. That gives three distinct interactions
    (A-A, B-B, A-B) without ever passing a cutoff, and the orderings below live in
    1x1, 1x2 and 2x2 supercells of the parent, so the fit must recover the same
    constants from cells of 2, 4 and 8 sites.

    Total energies are assigned from E = N * e0 - sum_<ij> J_ab s_i s_j, evaluated
    over explicit lattice vectors rather than over the mapper's own neighbor graph,
    so the topology is ground truth here and not something the test borrows back
    from the code under test.
    """

    E0 = -5.0  # eV per magnetic ion
    J_AB, J_AA, J_BB = 0.011, 0.004, -0.006  # eV, all distinct
    A, B = "Mn", "Fe"

    # D_AA must stay inside the default nearest-neighbor search's 10% window around
    # D_AB, or the chain bonds drop out of the graph; and the two must differ by more
    # than the mapper's tol, or they collapse into one shell.
    D_AB = 1.0  # A-B spacing along x
    D_AA = 1.05  # A-A and B-B spacing along y

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
        ex_params = hm.get_exchange()

        # One constant per sublattice pair, each recovered independently rather than
        # averaged together, and independent of the supercell each ordering lives in.
        assert set(ex_params) == {"E0", aa, bb, ab}
        assert ex_params[ab] == approx(self.J_AB * 1000, abs=1e-6)
        assert ex_params[aa] == approx(self.J_AA * 1000, abs=1e-6)
        assert ex_params[bb] == approx(self.J_BB * 1000, abs=1e-6)
        assert ex_params["E0"] == approx(self.E0, abs=1e-8)

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
