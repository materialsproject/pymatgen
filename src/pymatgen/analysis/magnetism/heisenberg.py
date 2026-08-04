"""
This module implements a simple algorithm for extracting nearest neighbor
exchange parameters by mapping low energy magnetic orderings to a Heisenberg
model.
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from ast import literal_eval
from functools import cached_property
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from monty.dev import deprecated
from monty.json import MSONable, jsanitize
from monty.serialization import dumpfn

from pymatgen.analysis.graphs import StructureGraph
from pymatgen.analysis.local_env import MinimumDistanceNN
from pymatgen.analysis.magnetism import CollinearMagneticStructureAnalyzer, Ordering
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

if TYPE_CHECKING:
    from typing import Self

# __author__ = "ncfrey"
# __version__ = "0.1"
# __maintainer__ = "Nathan C. Frey"
# __email__ = "ncfrey@lbl.gov"
# __status__ = "Development"
# __date__ = "June 2019"

__author__ = "Luca Frey"
__version__ = "0.2"
__maintainer__ = "Luca Frey"
__email__ = "luca.frey@student.kit.edu"
__status__ = "Development"
__date__ = "July 2026"

logger = logging.getLogger(__name__)

DEFAULT_TOL = 0.05  # Angstrom tolerance
DISTANCE_ROUND_DECIMALS = 2  # Round distances to this many decimals when matching interactions


def _analyzer(structure: Structure, **kwargs) -> CollinearMagneticStructureAnalyzer:
    """CollinearMagneticStructureAnalyzer factory with the settings used throughout this module.

    ``make_primitive=False`` keeps the cell as given, ``threshold=0.0`` retains any nonzero
    moment, and ``threshold_nonmag=100.0`` always zeroes out nonmagnetic ions, so induced
    magnetic moments don't change the number of magnetic sites per unit cell between orderings.

    Extra keyword arguments are passed through to the analyzer.
    """
    # Last value wins on conflicting keys, so the defaults above are overridden by kwargs.
    return CollinearMagneticStructureAnalyzer(
        structure, **{"make_primitive": False, "threshold": 0.0, "threshold_nonmag": 100.0, **kwargs}
    )


class MagneticOrdering(ABC):
    """A single collinear magnetic configuration on a common parent lattice.

    Base class holding one ordering's structure, its magnetic-only reduction,
    its NN graph, and the parent-sublattice labels of its magnetic sites.
    """

    def __init__(
        self,
        structure: Structure,
        magn_species: set[str],
        cutoff: float = 0,
        tol: float = DEFAULT_TOL,
    ):
        self.analyzer = _analyzer(structure)
        # The analyzer works on a copy, so the caller's structure is never mutated. On
        # this copy every site carries a float 'magmom' site property, whether the input
        # supplied moments as site properties or as species spins, and induced moments on
        # nonmagnetic ions have been zeroed out. Moment magnitudes are left untouched.
        self.structure = self.analyzer.structure
        self.magn_species = magn_species
        self.cutoff = cutoff
        self.tol = tol

        self.sublattice_ids: list[int] = []
        self.sublattice_wyckoff_symbols: dict[int, str] = {}

    @cached_property
    def magnetic_structure(self) -> Structure:
        """Magnetic-only cell. nn_graph and sublattice_ids are both indexed against this.

        Membership is decided by *species*, not by moment, using the magnetic species
        pooled over every ordering (see ``HeisenbergMapper._initialize_orderings``). A
        magnetic ion whose moment happens to relax to zero therefore stays on the lattice
        and contributes a zero term, rather than dropping out and changing this ordering's
        site count, graph topology and per-ion energy relative to its siblings.
        """
        return Structure.from_sites([site for site in self.structure if site.specie.symbol in self.magn_species])

    @staticmethod
    def _nonmagnetic(structure) -> Structure:
        """Nonmagnetic copy of a structure (moments zeroed, all ions kept)."""
        s0 = _analyzer(structure).get_nonmagnetic_structure(make_primitive=False)
        if "wyckoff" in s0.site_properties:
            s0.remove_site_property("wyckoff")
        return s0

    @staticmethod
    def _magnetic_species(structure) -> set[str]:
        """Symbols of the species carrying a moment in this structure.

        Used to pool the magnetic species across all orderings; the pooled set then
        defines the magnetic sublattice for every ordering and for the parent.

        A site counts as magnetic if it ends up with a nonzero moment after ``_analyzer``
        has processed it, i.e. it is of a species listed in the analyzer's default magmoms
        *and* was given a nonzero moment (``threshold=0.0`` retains any nonzero value).
        Induced moments on nonmagnetic ions are zeroed out first (``threshold_nonmag=100.0``)
        and so never contribute a spurious species here.
        """
        magnetic = _analyzer(structure).get_structure_with_only_magnetic_atoms(make_primitive=False)
        return {site.specie.symbol for site in magnetic}

    @cached_property
    def nn_graph(self) -> StructureGraph:
        """NN StructureGraph of the magnetic-only structure.

        If self.cutoff is set, the graph includes all neighbors within that distance.
        Otherwise, only the nearest neighbors of each site are included.
        """
        strategy = MinimumDistanceNN(cutoff=self.cutoff, get_all_sites=True) if self.cutoff else MinimumDistanceNN()
        return StructureGraph.from_local_env_strategy(self.magnetic_structure, strategy=strategy)

    @abstractmethod
    def set_sublattice_ids(self) -> None:
        """
        Sets the following attributes for this ordering:
            - `self.sublattice_ids` = [sublattice id for each *magnetic* site], aligned with
              `self.magnetic_structure` and hence with `self.nn_graph`. Nonmagnetic sites are
              not represented at all, so these are always plain ints - never None.
            - `self.sublattice_wyckoff_symbols` = {sublattice id: wyckoff symbol}

        For `ParentOrdering`, this defines the sublattice ids by grouping symmetrically
        equivalent sites into sublattices (nonmagnetic ions present for symmetry analysis)
        and keeping only the orbits of the magnetic species.

        For `RelaxedOrdering`, this maps the ordering's sites onto the parent's to find
        which sublattice each one belongs to.

        `ParentOrdering` additionally stores the *full-cell* ids (None on nonmagnetic sites)
        as a `sublattice_id` site property, which is what carries the labels across to each
        `RelaxedOrdering` during structure matching.
        """


class ParentOrdering(MagneticOrdering):
    """The nonmagnetic parent cell that defines the sublattices.

    Sublattices are its symmetry orbits (computed with nonmagnetic ions present),
    restricted to the magnetic species. Because the parent carries no moments it cannot
    tell which of its species are magnetic; that is pooled from the orderings and passed
    in as ``magn_species``.
    """

    def __init__(
        self,
        structure: Structure,
        magn_species: set[str],
        cutoff: float = 0,
        tol: float = DEFAULT_TOL,
    ):
        # Strip the moments *before* super().__init__, so the structure this ordering keeps
        # is the one its substructure and graph get derived from.
        super().__init__(self._nonmagnetic(structure), magn_species, cutoff, tol)
        self.set_sublattice_ids()

    def set_sublattice_ids(self):
        symmetrized_parent = SpacegroupAnalyzer(self.structure).get_symmetrized_structure()

        # Full-cell ids, None on the nonmagnetic sites. These exist to be stamped as a site
        # property, which is how the labels reach each RelaxedOrdering via structure matching.
        full_ids: list[int | None] = [None] * len(self.structure)
        for indices, wyckoff_symbol in zip(
            symmetrized_parent.equivalent_indices, symmetrized_parent.wyckoff_symbols, strict=True
        ):
            if symmetrized_parent[indices[0]].specie.symbol not in self.magn_species:
                continue
            sub_id = len(self.sublattice_wyckoff_symbols)
            # A sublattice always has one wyckoff symbol, but several sublattices can share
            # one (e.g. two species sitting on the same wyckoff position).
            self.sublattice_wyckoff_symbols[sub_id] = wyckoff_symbol
            for index in indices:
                full_ids[index] = sub_id

        self.structure.add_site_property("sublattice_id", full_ids)

        # Drop the nonmagnetic sites to line the ids up with magnetic_structure. Both select
        # on species, so the surviving ids are exactly the non-None ones, in order.
        self.sublattice_ids = [sub_id for sub_id in full_ids if sub_id is not None]


class RelaxedOrdering(MagneticOrdering):
    """One DFT-relaxed magnetic ordering with its energy and parent back-reference.

    The parent is optional at construction: the orderings have to exist, and be screened
    and sorted, before a parent can be inferred from them, so ``HeisenbergMapper`` attaches
    it afterwards via ``set_parent``.
    """

    def __init__(
        self,
        structure: Structure,
        energy: float,
        magn_species: set[str],
        parent_ordering: ParentOrdering | None = None,
        cutoff: float = 0,
        tol: float = DEFAULT_TOL,
    ):
        super().__init__(structure, magn_species, cutoff, tol)
        self.energy = energy  # total energy, as supplied by the caller
        self.parent_ordering = parent_ordering
        if parent_ordering is not None:
            self.set_sublattice_ids()

    @property
    def energy_per_magnetic_ion(self) -> float:
        """Total energy divided by the number of magnetic ions (eV).

        This is what makes orderings living in different-sized supercells comparable.
        """
        return self.energy / len(self.magnetic_structure)

    def set_parent(self, parent_ordering: ParentOrdering) -> None:
        """Attach the parent that defines the sublattices, and label this ordering."""
        self.parent_ordering = parent_ordering
        self.set_sublattice_ids()

    def set_sublattice_ids(self):
        matcher = StructureMatcher(primitive_cell=False, attempt_supercell=True)

        # Fold this ordering's full geometry onto the tagged parent cell (the nonmagnetic
        # ions help pin down a unique mapping). get_s2_like_s1 returns the parent's sites,
        # carrying their 'sublattice_id', in this ordering's site order - so matched_parent[i]
        # is the parent site that site i of this ordering sits on.
        matched_parent = matcher.get_s2_like_s1(self._nonmagnetic(self.structure), self.parent_ordering.structure)
        if matched_parent is None:
            raise ValueError(
                "This ordering is not a supercell of the parent cell; it cannot be mapped "
                "onto the parent sublattices. Pass an explicit `parent` cell that all "
                "orderings share."
            )

        # Keep only the magnetic sites, in order, so the ids line up with magnetic_structure
        # and the graph built from it. A parent site is labelled iff its species is magnetic,
        # which is the same test magnetic_structure applies.
        self.sublattice_ids = [
            int(site.properties["sublattice_id"])
            for site in matched_parent
            if site.properties["sublattice_id"] is not None
        ]
        self.sublattice_wyckoff_symbols = self.parent_ordering.sublattice_wyckoff_symbols


class HeisenbergMapper:
    """Compute exchange parameters from low energy magnetic orderings.

    Sublattices are defined once, on a paramagnetic *parent cell* (its symmetry
    orbits / Wyckoff positions). Every magnetic ordering is treated as a spin
    sample drawn on that parent lattice, so each magnetic site is labelled with
    the parent sublattice it belongs to. Because the labels live in the parent
    cell - not in any one ordering's supercell - orderings that occupy
    different-sized supercells share a single, consistent set of exchange
    parameters.

    Attributes:
        orderings (list[RelaxedOrdering]): The screened orderings, sorted by energy per
            magnetic ion. Each owns its magnetic-only structure, its NN graph and its
            parent-sublattice labels.
        parent (ParentOrdering): Nonmagnetic parent cell that defines the sublattices.
        nn_interactions (dict): {(i, j): [J label, ...]} - the distinct interactions of
            each sublattice pair, ordered from the nearest shell outwards.
        dists (dict): {J label: interaction distance in Angstrom}.
        ex_mat (DataFrame): Heisenberg Hamiltonian (per magnetic ion) for each ordering.
        ex_params (dict): Exchange parameter values. The J_ij are in meV/muB^2 (they
            multiply the raw moments, see get_exchange); the included 'E0' offset is in
            eV per magnetic ion. get_interaction_graph carries the per-bond J_ij in meV.
    """

    def __init__(self, ordered_structures, energies, parent=None, cutoff=0, tol: float = DEFAULT_TOL):
        """Exchange parameters are computed by mapping to a classical Heisenberg
        model. n+1 unique orderings are required to compute n exchange parameters.

        First run a MagneticOrderings Flow to obtain low energy collinear magnetic
        orderings and find the magnetic ground state, enumerate magnetic
        states with the ground state as the input structure and do static
        calculations for these orderings. The orderings may live in different
        supercells - they only need to be commensurate with a common parent cell.

        ***IMPORTANT NOTE***

        If the parent is not supplied, it is inferred as the primitive cell of the lowest-energy ordering.
        Not supplying a parent cell is only safe if the lowest-energy ordering still preserves the symmetry of the paramagnetic parent.
        In most cases, relaxation will lower the symmetry and the parent cell must be supplied explicitly.

        Args:
            ordered_structures (list): Structure objects with magmoms.
            energies (list): Total energies of each relaxed magnetic structure.
            parent (Structure): Paramagnetic parent cell whose symmetry defines the
                magnetic sublattices. If None, it is inferred as the primitive cell of
                the lowest-energy ordering. Defaults to None.
            cutoff (float): Cutoff in Angstrom for nearest neighbor search. Defaults to 0,
                which keeps only each site's nearest neighbors, i.e. one shell per
                sublattice pair; with a cutoff, every interaction up to it is kept.
            tol (float): Tolerance (in Angstrom) on nearest neighbor distances
                being equal.
        """
        # Save original copies of inputs
        self.ordered_structures_ = ordered_structures
        self.energies_ = energies
        self.parent_ = parent

        self.cutoff = cutoff
        self.tol = tol

        # These attributes are set by internal methods, listed here for clarity.
        self.orderings = self.parent = None  # set by _initialize_orderings
        self.nn_interactions = self.dists = None  # set by _set_interactions
        self.ex_mat = self.ex_params = None  # set by _build_exchange_mat and get_exchange

        self._initialize_orderings(ordered_structures, energies, parent)
        self._set_interactions()
        self._build_exchange_mat()

    def _initialize_orderings(self, ordered_structures, energies, parent):
        """Build the RelaxedOrdering objects and the ParentOrdering that labels them.

        Sets self.orderings and self.parent.

        The magnetic species are pooled over *all* orderings and then used to define the
        magnetic sublattice of every ordering and of the parent. Pooling matters: if an ion
        happens to relax to a zero moment in one ordering, it must stay on the magnetic
        lattice there (contributing a zero term) rather than vanishing and leaving that
        ordering with a different site count, graph topology and per-ion energy than its
        siblings.

        This function does:
         - Build a set of magnetic species pooled over all orderings.
         - Build a RelaxedOrdering for each ordering, with its magnetic-only structure and NN graph. Since the parent is not yet known, the sublattice ids are not set yet.
         - Drop duplicate/degenerate orderings and sort by energy per magnetic ion using HeisenbergScreener.
         - Build the ParentOrdering from the lowest-energy ordering (or the explicit parent if supplied), and set its sublattice ids.

        Args:
            ordered_structures (list): Structure objects with magmoms.
            energies (list): Total energies of each relaxed magnetic structure.
            parent (Structure | None): Explicit paramagnetic parent cell, or None to infer
                it from the lowest-energy ordering.

        Raises:
            SystemExit: If fewer than 2 unique orderings remain after screening.
        """
        # Pool the magnetic species over all orderings, to make sure a species that relaxes to zero moment in one ordering still counts as magnetic there.
        # Since the magmom of this species is zero, it contributes a zero term to the Heisenberg Hamiltonian, but it stays on the lattice and keeps the site count and graph topology consistent with the other orderings.
        magn_species = set().union(*(MagneticOrdering._magnetic_species(struct) for struct in ordered_structures))

        orderings = [
            RelaxedOrdering(struct, energy, magn_species, cutoff=self.cutoff, tol=self.tol)
            for struct, energy in zip(ordered_structures, energies, strict=True)
        ]

        # Drop duplicate/degenerate orderings and sort by energy per magnetic ion.
        self.orderings = HeisenbergScreener(orderings, screen=False).screened_orderings

        if len(self.orderings) < 2:
            raise SystemExit("HeisenbergMapper needs at least 2 unique orderings.")

        # The nonmagnetic ions are kept in the parent: site equivalence is read from its
        # symmetry, and removing them first can raise the apparent site symmetry and merge
        # sublattices that are actually distinct.
        reference = parent if parent is not None else self.orderings[0].structure
        self.parent = ParentOrdering(
            MagneticOrdering._nonmagnetic(reference).get_primitive_structure(),
            magn_species,
            cutoff=self.cutoff,
            tol=self.tol,
        )

        # The parent cell defines the sublattices; label every magnetic site in every
        # ordering with the parent sublattice it belongs to.
        for ordering in self.orderings:
            ordering.set_parent(self.parent)

    @property
    def structures(self):
        """list[Structure]: Each ordering with all ions retained."""
        return [ordering.structure for ordering in self.orderings]

    @property
    def magnetic_structures(self):
        """list[Structure]: Magnetic-only structure of each ordering."""
        return [ordering.magnetic_structure for ordering in self.orderings]

    @property
    def energies(self):
        """list[float]: Energy per magnetic ion (eV) of each ordering."""
        return [ordering.energy_per_magnetic_ion for ordering in self.orderings]

    @property
    def nn_graphs(self):
        """list[StructureGraph]: NN graph of each ordering's magnetic-only structure."""
        return [ordering.nn_graph for ordering in self.orderings]

    @property
    def sublattice_ids(self):
        """list[list[int]]: sublattice_ids[k][i] is the sublattice id of magnetic site i
        in ordering k, aligned with that ordering's graph.
        """
        return [ordering.sublattice_ids for ordering in self.orderings]

    @property
    def parent_nn_graph(self):
        """StructureGraph: NN graph of the magnetic-only parent structure."""
        return self.parent.nn_graph

    @property
    def parent_sublattice_ids(self):
        """list[int]: Sublattice id of each magnetic site in the parent."""
        return self.parent.sublattice_ids

    @property
    def sublattice_wyckoff_symbols(self):
        """dict[int, str]: Maps each sublattice id to its wyckoff symbol."""
        return self.parent.sublattice_wyckoff_symbols

    @staticmethod
    def _order_sublattice_ids(i_id, j_id):
        """Returns the sublattice_ids in the order (i, j) with i <= j. This is the key used to look up the interactions of a sublattice pair."""
        return tuple(sorted((i_id, j_id)))

    def _interaction_label(self, i_id, j_id, dist):
        """Return the J label of an interaction, or None if it matches no interaction.

        Args:
            i_id (int): sublattice id of the ith site
            j_id (int): sublattice id of the jth site
            dist (float): distance (Angstrom) between the sites

        Returns:
            str | None: '<i>-<j>-<shell>' label, e.g. '0-1-nn'.
        """
        for label in self.nn_interactions.get(self._order_sublattice_ids(i_id, j_id), ()):
            if abs(dist - self.dists[label]) <= self.tol:
                return label
        logger.warning(
            f"Interaction between sublattices {i_id} and {j_id} at {dist:.2f} Angstrom does not match any interaction in nn_interactions built from the parent:\n{self.nn_interactions}\n"
        )
        return None

    def _set_interactions(self):
        """Set self.dists and self.nn_interactions describing the distinct interactions.

        An interaction is a (sublattice pair, neighbor shell) combination, labelled
        '<i>-<j>-<shell>' with i <= j and shell one of 'nn', 'nnn', 'nnnn', ... Shells are
        counted *within* a sublattice pair, not globally: '0-0-nn' and '0-1-nn' are the
        nearest 0-0 and 0-1 interaction respectively, even when one of them is the longer interaction.

        Every distinct interaction the parent graph contains is kept, so which interactions are
        parameterized follows entirely from how that graph is built: with cutoff=0 each
        site keeps only its nearest neighbors (one nn-shell per sublattice pair, but every
        pair the lattice realizes), with cutoff > 0 it keeps every interaction up to the cutoff.

        Distances and connectivity are read from the parent ordering; see
        _initialize_orderings() for how the parent is defined and how the sublattices are
        labelled. The connectivity comes from the parent's StructureGraph, which is built
        from the magnetic-only parent structure - the nonmagnetic ions are ignored for the
        graph, but kept in the parent structure to preserve the true site symmetry.
        """
        nn_graph = self.parent.nn_graph
        sub_ids = self.parent.sublattice_ids

        # Distinct interaction distances of each sublattice pair. Every site is visited, since
        # with cutoff=0 the two ends of an interaction need not both count it as a nearest neighbor.
        pair_dists: dict[tuple[int, int], set[float]] = {}
        for i in range(len(nn_graph)):
            for cs in nn_graph.get_connected_sites(i):
                pair = self._order_sublattice_ids(sub_ids[i], sub_ids[cs[2]])
                pair_dists.setdefault(pair, set()).add(round(cs[-1], DISTANCE_ROUND_DECIMALS))

        self.dists = {}
        self.nn_interactions = {}
        for pair, dists in sorted(pair_dists.items()):
            labels = []
            for dist in sorted(dists):
                # Distances equal within tol are the same shell; the smallest one names it.
                if labels and dist - self.dists[labels[-1]] <= self.tol:
                    continue
                label = f"{pair[0]}-{pair[1]}-{'n' * (len(labels) + 2)}"
                self.dists[label] = dist
                labels.append(label)
            self.nn_interactions[pair] = labels

    def _exchange_columns(self):
        """Column labels for the exchange matrix: 'E', 'E0' then one per J_ij.

        The J columns are ordered by increasing interaction length, so when the system has more
        interactions than the orderings can constrain, the longest-ranged ones are dropped.
        """
        j_columns = sorted(self.dists, key=self.dists.get)
        # Keep at most n interactions for n+1 orderings
        return ["E", "E0", *j_columns][: len(self.orderings) + 1]

    def _build_exchange_mat(self):
        """Build the Heisenberg Hamiltonian, one row per ordering, by summing the
        signed products S_i . S_j over each graph. Sets self.ex_mat.

        Each row is normalised per magnetic ion so that orderings living in
        different-sized supercells share a single linear system (the energies are
        per magnetic ion too, see RelaxedOrdering.energy_per_magnetic_ion).
        """
        columns = self._exchange_columns()
        j_columns = [c for c in columns if c not in ("E", "E0")]

        if len(j_columns) < 2:  # too few interactions to fit; get_exchange raises on this
            self.ex_mat = pd.DataFrame(columns=columns)
            return

        rows = []
        seen = set()  # de-duplicate identical interaction rows so no ordering is weighted twice
        for ordering in self.orderings:
            nn_graph = ordering.nn_graph
            sub_ids = ordering.sublattice_ids
            magmoms = nn_graph.structure.site_properties["magmom"]
            n_sites = len(nn_graph.structure)

            row = dict.fromkeys(columns, 0.0)
            for i in range(len(nn_graph.graph.nodes)):
                s_i = magmoms[i]
                for cs in nn_graph.get_connected_sites(i):
                    col = self._interaction_label(sub_ids[i], sub_ids[cs[2]], round(cs[-1], DISTANCE_ROUND_DECIMALS))
                    # Two ways an interaction has no column in row:
                    # 1. col is None, meaning this
                    # ordering's geometry has drifted off the parent lattice (a real
                    # problem, warned about in _interaction_label), or
                    # 2. col is just not in row, but is a genuine
                    # interaction that _exchange_columns truncated away because the
                    # orderings cannot constrain that many parameters (by design).
                    if col in row:
                        row[col] -= s_i * magmoms[cs[2]]

            # Extensive sums -> per magnetic ion, with the 1/2 Heisenberg factor for double counting.
            for c in j_columns:
                row[c] /= 2 * n_sites

            key = tuple(round(row[c], 10) for c in j_columns)
            if key in seen:
                continue
            seen.add(key)
            row["E0"] = 1.0  # nonmagnetic contribution (per ion)
            row["E"] = ordering.energy_per_magnetic_ion
            rows.append(row)

        ex_mat = pd.DataFrame(rows, columns=columns)

        # Drop interaction columns that never appear (all zero) to keep H full rank
        j_columns = [c for c in j_columns if not (ex_mat[c] == 0).all()]

        # Every ordering is kept: get_exchange fits the parameters by least squares, so
        # surplus orderings average out the noise on the energies instead of being
        # discarded to square the system.
        self.ex_mat = ex_mat[["E", "E0", *j_columns]].reset_index(drop=True)

    def get_exchange(self):
        """
        Take Heisenberg Hamiltonian and corresponding energy for each row and
        solve for the exchange parameters in the least-squares sense.

        With exactly as many orderings as parameters this reproduces the solution of the
        square linear system of equations; with more it is a fit over all of them compensating
        for the noise in the DFT-energies.

        The rows of ex_mat multiply the raw magnetic moments in muB rather than normalized
        spins, so the fitted J_ij come out in meV/muB^2: it takes J_ij * m_i * m_j to get
        the energy of a bond. Leaving the moments unnormalized is deliberate: orderings
        relax to different moment magnitudes, and it is the bond energy that should follow
        the product m_i * m_j, while the interaction strength itself stays the same.
        Normalizing every moment to 1 would fold that ordering-to-ordering variation into
        the fitted J_ij as noise. Codes and papers working with
        normalized spins (VAMPIRE, UppASD, TB2J) instead report J in meV for
        E = -sum_<ij> J_ij e_i.e_j with |e| = 1; get_interaction_graph does that conversion.

        Returns:
            ex_params (dict[str, float]): Exchange parameters. The J_ij are in meV/muB^2;
                the included 'E0' offset is in eV per magnetic ion (its column in ex_mat
                is 1.0 and the fitted energies are per magnetic ion, so it comes out
                intensive, which is what lets orderings of different size share one fit).
            residual (float): Sum of squared residuals of the least-squares fit, in
                (meV per magnetic ion)^2. It is a residual on the fitted energies, so
                unlike the J_ij it keeps the per-magnetic-ion normalisation.

        Raises:
            ValueError: If the orderings constrain fewer than two exchange interactions,
                leaving nothing for the fit to solve.
        """
        ex_mat = self.ex_mat
        E = ex_mat[["E"]]
        col_names = [c for c in ex_mat.columns if c != "E"]

        # E0 and a single J cannot be fitted from the energies alone, so there is nothing
        # to solve here.
        if len(col_names) < 3:
            raise ValueError(
                f"Exchange matrix holds {len(col_names) - 1} interaction(s) besides E0; a least-squares "
                "fit needs at least 2. Supply orderings whose nearest-neighbor shell holds more than "
                "one sublattice pair, or set a cutoff so that further shells are included."
            )

        # Fit E0 and the J_ij to the energies of every ordering
        H = np.array(ex_mat.loc[:, ex_mat.columns != "E"].values).astype(float)

        # Warn when the fit is ill-conditioned: near-degenerate orderings or an
        # over-parameterized model make H nearly singular, so tiny energy
        # differences blow up into unphysical exchange parameters.
        cond = np.linalg.cond(H)
        if cond > 1e5:
            logger.warning(
                f"Exchange matrix is ill-conditioned (cond={cond:.1e}); the fitted exchange "
                "parameters are unreliable. The input orderings are near-degenerate or the "
                "model has more parameters than the orderings can constrain. Supply more, "
                "more-distinct orderings."
            )

        j_ij, residuals, rank, _singular = np.linalg.lstsq(H, E, rcond=None)

        # lstsq only fills residuals for an overdetermined full-rank fit; for a square,
        # underdetermined or rank-deficient one it returns an empty array, so compute the
        # sum of squared residuals here. It sits on the energy side of the system, not the
        # J_ij side, so it is a squared energy in (eV per magnetic ion)^2.
        residual = float(residuals[0]) if residuals.size else float(np.sum((H @ j_ij - np.asarray(E)) ** 2))

        # A least-squares fit does not fail on a system it cannot determine,
        # it silently returns the minimum-norm solution. cond above cannot see this case:
        # with fewer rows than parameters it stays finite.
        if rank < H.shape[1]:
            logger.warning(
                f"Exchange matrix is rank deficient (rank {rank} for {H.shape[1]} parameters); "
                "the orderings do not constrain every exchange parameter, and the values "
                "returned are one of infinitely many fits. Supply more distinct orderings."
            )

        j_ij[1:] *= 1000  # J_ij in meV/muB^2 (E0 stays in eV per magnetic ion)
        residual *= 1e6  # convert to (meV per magnetic ion)^2
        ex_params = {j_name: j[0] for j_name, j in zip(col_names, j_ij.tolist(), strict=True)}

        self.ex_params = ex_params
        self.residual = residual
        return self.ex_params, self.residual

    def get_low_energy_orderings(self):
        """Find lowest energy FM and AFM orderings to compute E_AFM - E_FM.

        Returns:
            fm_struct (Structure): fm structure with 'magmom' site property
            afm_struct (Structure): afm structure with 'magmom' site property
            fm_e (float): fm energy
            afm_e (float): afm energy
        """
        fm_struct, afm_struct = None, None
        mag_min = np.inf
        mag_max = 0.001
        fm_e = afm_e = fm_e_min = afm_e_min = 0

        for magnetic_ordering in self.orderings:
            s = magnetic_ordering.magnetic_structure
            e = magnetic_ordering.energy_per_magnetic_ion
            ordering = _analyzer(s).ordering
            magmoms = s.site_properties["magmom"]

            # Try to find matching orderings first
            if ordering == Ordering.FM and e < fm_e_min:
                fm_struct = s
                mag_max = abs(sum(magmoms))
                fm_e = e
                fm_e_min = e

            if ordering == Ordering.AFM and e < afm_e_min:
                afm_struct = s
                afm_e = e
                mag_min = abs(sum(magmoms))
                afm_e_min = e

        # Brute force search for closest thing to FM and AFM
        if not fm_struct or not afm_struct:
            for magnetic_ordering in self.orderings:
                s = magnetic_ordering.magnetic_structure
                e = magnetic_ordering.energy_per_magnetic_ion
                magmoms = s.site_properties["magmom"]

                if abs(sum(magmoms)) > mag_max:  # FM ground state
                    fm_struct = s
                    fm_e = e
                    mag_max = abs(sum(magmoms))

                # AFM ground state
                if abs(sum(magmoms)) < mag_min:
                    afm_struct = s
                    afm_e = e
                    mag_min = abs(sum(magmoms))
                    afm_e_min = e
                elif abs(sum(magmoms)) == 0 and mag_min == 0 and e < afm_e_min:
                    afm_struct = s
                    afm_e = e
                    afm_e_min = e

        return fm_struct, afm_struct, fm_e, afm_e

    @deprecated(
        get_exchange,
        message=(
            "<J> is a single average in meV/magnetic ion rather than a set of shell-resolved "
            "J_ij, and it is only defined for a pair of FM/AFM orderings."
        ),
        category=DeprecationWarning,
        deadline=(2027, 8, 1),
    )
    def estimate_exchange(self, fm_struct=None, afm_struct=None, fm_e=None, afm_e=None):
        """Estimate <J> for a structure based on low energy FM and AFM orderings.

        .. deprecated::
            Use :meth:`get_exchange` instead, which fits shell-resolved J_ij over all
            supplied orderings.

        Args:
            fm_struct (Structure): fm structure with 'magmom' site property
            afm_struct (Structure): afm structure with 'magmom' site property
            fm_e (float): fm energy per magnetic ion
            afm_e (float): afm energy per magnetic ion

        Returns:
            float: Average J exchange parameter (meV / magnetic ion)
        """
        # Get low energy orderings if not supplied
        if any(arg is None for arg in [fm_struct, afm_struct, fm_e, afm_e]):
            fm_struct, afm_struct, fm_e, afm_e = self.get_low_energy_orderings()

        magmoms = fm_struct.site_properties["magmom"]
        m_avg = np.mean([np.sqrt(m**2) for m in magmoms])

        # If m_avg for FM config is < 1 we won't get sensible results.
        if m_avg < 1:
            logger.warning(
                "Local magnetic moments are small (< 1 muB / atom). The exchange parameters may "
                "be wrong, but <J> and the mean field critical temperature estimate may be OK."
            )

        delta_e = afm_e - fm_e  # J > 0 -> FM
        j_avg = delta_e / (m_avg**2)  # eV / magnetic ion
        j_avg *= 1000  # meV / ion

        return j_avg

    @deprecated(
        message=(
            "<J> is in units of meV/magnetic ion, the multi-sublattice branch double counts the "
            "diagonal entries of omega, and the result is only a crude estimate of the true "
            "critical temperature."
        ),
        category=DeprecationWarning,
        deadline=(2027, 8, 1),
    )
    def get_mft_temperature(self, j_avg):
        """
        Crude mean field estimate of critical temperature based on <J> for
        one sublattice, or solving the coupled equations for a multi-sublattice
        material.

        .. deprecated::
            No direct replacement; use a Monte Carlo solver (e.g. VAMPIRE) on the
            exchange parameters from :meth:`get_exchange` for a reliable T_c.

        Args:
            j_avg (float): Average exchange parameter (meV / magnetic ion)

        Returns:
            float: Critical temperature mft_t (K)
        """
        # Number of magnetic sublattices = number of parent orbits
        n_sub_lattices = len(self.sublattice_wyckoff_symbols)
        k_boltzmann = 0.0861733  # meV/K

        # Only 1 magnetic sublattice
        if n_sub_lattices == 1:
            mft_t = 2 * abs(j_avg) / 3 / k_boltzmann

        else:  # multiple magnetic sublattices
            omega = np.zeros((n_sub_lattices, n_sub_lattices))
            ex_params = {k: v for (k, v) in self.ex_params.items() if k != "E0"}  # ignore E0
            for k, j_val in ex_params.items():
                # split into i, j sublattice ids (cut the shell identifier)
                i, j = (int(num) for num in k.split("-")[:2])
                omega[i, j] += j_val
                omega[j, i] += j_val

            omega = omega * 2 / 3 / k_boltzmann
            # omega is symmetric by construction, so use eigvalsh to guarantee
            # real eigenvalues (np.linalg.eig can return complex128 with zero
            # imaginary part depending on the LAPACK build).
            eigen_vals = np.linalg.eigvalsh(omega)
            mft_t = max(eigen_vals)

        if mft_t > 1500:  # Not sensible!
            logger.warning(
                "This mean field estimate is too high! Probably the true low energy orderings were not given as inputs."
            )

        return mft_t

    def get_interaction_graph(self, filename=None, ordering_index=0):
        """Get a StructureGraph with edges and weights that correspond to exchange
        interactions and J_ij values, respectively.

        Edge weights are in meV, for the normalized-spin Hamiltonian
        E = -sum_<ij> J_ij e_i.e_j, i.e. this ordering's moments are folded into the
        fitted meV/muB^2 parameters (see get_exchange). That is the convention VAMPIRE,
        UppASD and TB2J use, so the weights can go straight into a ucf file.

        Args:
            filename (str): if not None, save interaction graph to filename.
            ordering_index (int): Which ordering (and its supercell) to build the
                graph for. Site indices and the J_ij lookup use this ordering's
                sublattice labels. Defaults to 0 (the lowest-energy ordering).

        Returns:
            StructureGraph: Exchange interaction graph.
        """
        if self.ex_params is None:
            self.get_exchange()

        ordering = self.orderings[ordering_index]
        structure = ordering.magnetic_structure
        nn_graph = ordering.nn_graph
        magmoms = structure.site_properties["magmom"]
        sub_ids = ordering.sublattice_ids

        igraph = StructureGraph.from_empty_graph(
            structure, edge_weight_name="exchange_constant", edge_weight_units="meV"
        )

        # J_ij exchange interaction matrix
        for i in range(len(nn_graph.graph.nodes)):
            for c in nn_graph.get_connected_sites(i):
                jimage = c[1]  # relative integer coordinates of atom j
                j = c[2]  # index of neighbor
                j_exc = self._get_j_exc(sub_ids[i], sub_ids[j], c[-1])
                # Weights are in meV, so fold this ordering's moments into the fitted
                # meV/muB^2 parameters. Only the magnitudes enter: the relative orientation
                # of the two sites belongs to the spin vectors, not to J_ij.
                j_exc *= abs(magmoms[i] * magmoms[j])
                # Only add interactions the fit actually parameterized. Unparameterized
                # interactions (no matching sublattice-pair/shell in ex_params, exactly
                # the interactions _build_exchange_mat also ignores) get j_exc == 0, which
                # StructureGraph.add_edge silently drops via its falsy-weight
                # guard, leaving an edge with weight=None that breaks downstream
                # per-interaction consumers.
                if not j_exc:
                    continue
                igraph.add_edge(i, j, to_jimage=jimage, weight=j_exc, warn_duplicates=False)

        if filename:
            if not filename.endswith(".json"):
                filename += ".json"
            dumpfn(igraph, filename)

        return igraph

    def _get_j_exc(self, i_id, j_id, dist):
        """Look up the exchange parameter between two sublattices at a distance.

        Args:
            i_id (int): sublattice id of the ith site
            j_id (int): sublattice id of the jth site
            dist (float): distance (Angstrom) between sites (10E-2 precision)

        Returns:
            float: Exchange parameter J_exc in meV/muB^2 (0 if no matching interaction).
                Multiply by the two moments to get meV, as get_interaction_graph does.
        """
        label = self._interaction_label(i_id, j_id, round(dist, DISTANCE_ROUND_DECIMALS))

        return self.ex_params.get(label, 0)

    def get_heisenberg_model(self):
        """Save results of mapping to a HeisenbergModel object.

        Returns:
            HeisenbergModel: MSONable object.
        """
        ex_params, residual = self.get_exchange()
        return HeisenbergModel(
            formula=str(self.ordered_structures_[0].reduced_formula),
            structures=self.structures,
            magnetic_structures=self.magnetic_structures,
            energies=self.energies,
            cutoff=self.cutoff,
            tol=self.tol,
            nn_graphs=self.nn_graphs,
            sublattice_ids=self.sublattice_ids,
            sublattice_wyckoff_symbols=self.sublattice_wyckoff_symbols,
            nn_interactions=self.nn_interactions,
            dists=self.dists,
            ex_mat=self.ex_mat,
            ex_params=ex_params,
            residual=residual,
            igraph=self.get_interaction_graph(),
        )


class HeisenbergScreener:
    """Clean and screen magnetic orderings."""

    def __init__(self, orderings: list[RelaxedOrdering], screen=False):
        """Pre-processes magnetic orderings for HeisenbergMapper.
        It prioritizes low-energy orderings with large and localized magnetic moments.

        Args:
            orderings (list[RelaxedOrdering]): The orderings to screen. Each one already
                owns its magnetic-only substructure and its per-magnetic-ion energy.
            screen (bool): Try to screen out high energy and low-spin configurations.

        Attributes:
            screened_orderings (list[RelaxedOrdering]): Deduplicated orderings, sorted by
                energy per magnetic ion.
        """
        orderings = self._do_cleanup(orderings)

        # If there are more than 2 orderings, we want to perform a
        # screening to prioritize well-behaved ones
        if screen and len(orderings) > 2:
            orderings = self._do_screen(orderings)

        self.screened_orderings = orderings

    @staticmethod
    def _do_cleanup(orderings: list[RelaxedOrdering]):
        """Drop duplicate/degenerate orderings and sort by energy per magnetic ion.

        Sometimes different initial configs relax to the same state; those show up as
        orderings with (near-)identical energies per magnetic ion.

        Args:
            orderings (list[RelaxedOrdering]): The orderings to clean up.

        Returns:
            list[RelaxedOrdering]: Deduplicated, sorted by energy per magnetic ion.
        """
        e_tol = 6  # 10^-6 eV/atom tol on energies
        energies = [round(ordering.energy_per_magnetic_ion, e_tol) for ordering in orderings]

        remove_list = []
        for idx, energy in enumerate(energies):
            if idx not in remove_list:
                for i_check, e_check in enumerate(energies):
                    if idx != i_check and i_check not in remove_list and energy == e_check:
                        remove_list.append(i_check)

        keep = [idx for idx in range(len(energies)) if idx not in remove_list]
        keep.sort(key=lambda idx: energies[idx])

        return [orderings[idx] for idx in keep]

    @staticmethod
    def _do_screen(orderings: list[RelaxedOrdering]):
        """Screen and sort magnetic orderings based on some criteria.

        Prioritize low energy orderings and large, localized magmoms. _do_cleanup should be
        run first to deduplicate and sort the orderings by energy.

        Args:
            orderings (list[RelaxedOrdering]): Cleaned up orderings, sorted by energy.

        Returns:
            list[RelaxedOrdering]: The ground and first excited state, followed by the
                remaining orderings sorted by how few small moments they carry.
        """

        def n_below_1ub(ordering):
            magmoms = ordering.magnetic_structure.site_properties["magmom"]
            return sum(abs(magmom) < 1 for magmom in magmoms)

        # Keep the ground and first excited state fixed to capture the low-energy spectrum,
        # and prioritize the rest by having fewer magmoms < 1 uB.
        return [*orderings[:2], *sorted(orderings[2:], key=n_below_1ub)]


class HeisenbergModel(MSONable):
    """
    Store a Heisenberg model fit to low-energy magnetic orderings.
    Intended to be generated by HeisenbergMapper.get_heisenberg_model().
    """

    def __init__(
        self,
        formula=None,
        structures=None,
        magnetic_structures=None,
        energies=None,
        cutoff=None,
        tol=None,
        nn_graphs=None,
        sublattice_ids=None,
        sublattice_wyckoff_symbols=None,
        nn_interactions=None,
        dists=None,
        ex_mat=None,
        ex_params=None,
        residual=None,
        igraph=None,
    ):
        """
        Args:
            formula (str): Reduced formula of compound.
            structures (list): Each ordering with all ions retained, with magmoms.
            magnetic_structures (list): Magnetic-only cell of each ordering. nn_graphs and
                sublattice_ids are indexed against these, not against structures.
            energies (list): Energies of each relaxed magnetic structure.
            cutoff (float): Cutoff in Angstrom for nearest neighbor search.
            tol (float): Tolerance (in Angstrom) on nearest neighbor distances being equal.
            nn_graphs (list): StructureGraph objects.
            sublattice_ids (list[list[int]]): sublattice_ids[k][i] is the sublattice id of
                site i in ordering k.
            sublattice_wyckoff_symbols (dict): Maps each sublattice id to its wyckoff symbol.
            nn_interactions (dict): {(i, j): [J label, ...]} - the distinct interactions
                of each sublattice pair, ordered from the nearest shell outwards.
            dists (dict): {J label: interaction distance in Angstrom}.
            ex_mat (DataFrame): Heisenberg Hamiltonian (per magnetic ion) for each ordering.
            ex_params (dict): Exchange parameter values. The J_ij are in meV/muB^2 (they
                multiply the raw moments, see HeisenbergMapper.get_exchange); the included
                'E0' offset is in eV per magnetic ion.
            residual (float): Sum of squared residuals of the fit that produced ex_params,
                in (meV per magnetic ion)^2.
            igraph (StructureGraph): Exchange interaction graph, edge weights in meV.
        """
        self.formula = formula
        self.structures = structures
        self.magnetic_structures = magnetic_structures
        self.energies = energies
        self.cutoff = cutoff
        self.tol = tol
        self.nn_graphs = nn_graphs
        self.sublattice_ids = sublattice_ids
        self.sublattice_wyckoff_symbols = sublattice_wyckoff_symbols
        self.nn_interactions = nn_interactions
        self.dists = dists
        self.ex_mat = ex_mat
        self.ex_params = ex_params
        self.residual = residual
        self.igraph = igraph

    def as_dict(self):
        """Because some dicts have int keys, some sanitization is required for JSON compatibility."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "@version": __version__,
            "formula": self.formula,
            "structures": [struct.as_dict() for struct in self.structures],
            "magnetic_structures": [struct.as_dict() for struct in self.magnetic_structures],
            "energies": self.energies,
            "cutoff": self.cutoff,
            "tol": self.tol,
            "nn_graphs": [nn_graph.as_dict() for nn_graph in self.nn_graphs],
            "sublattice_ids": self.sublattice_ids,
            "dists": self.dists,
            "ex_params": self.ex_params,
            "residual": self.residual,
            "igraph": self.igraph.as_dict(),
            # Sanitize int keys / DataFrame
            "ex_mat": jsanitize(self.ex_mat.to_dict()),
            "nn_interactions": jsanitize(self.nn_interactions),
            "sublattice_wyckoff_symbols": jsanitize(self.sublattice_wyckoff_symbols),
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Create a HeisenbergModel from a dict."""
        # Reconstitute the tuple/int-keyed dicts that jsanitize stringified
        nn_interactions = {literal_eval(pair): labels for pair, labels in dct["nn_interactions"].items()}
        sublattice_wyckoff_symbols = {literal_eval(k): v for k, v in dct["sublattice_wyckoff_symbols"].items()}

        structures = [Structure.from_dict(v) for v in dct["structures"]]
        magnetic_structures = [Structure.from_dict(v) for v in dct["magnetic_structures"]]
        nn_graphs = [StructureGraph.from_dict(v) for v in dct["nn_graphs"]]
        igraph = StructureGraph.from_dict(dct["igraph"])

        # Reconstitute the exchange matrix DataFrame
        ex_mat = pd.DataFrame.from_dict(dct["ex_mat"]) if dct["ex_mat"] else pd.DataFrame(columns=["E", "E0"])

        return cls(
            formula=dct["formula"],
            structures=structures,
            magnetic_structures=magnetic_structures,
            energies=dct["energies"],
            cutoff=dct["cutoff"],
            tol=dct["tol"],
            nn_graphs=nn_graphs,
            sublattice_ids=dct["sublattice_ids"],
            sublattice_wyckoff_symbols=sublattice_wyckoff_symbols,
            nn_interactions=nn_interactions,
            dists=dct["dists"],
            ex_mat=ex_mat,
            ex_params=dct["ex_params"],
            residual=dct["residual"],
            igraph=igraph,
        )
