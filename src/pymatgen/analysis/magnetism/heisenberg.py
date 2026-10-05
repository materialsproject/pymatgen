"""
This module implements a simple algorithm for extracting nearest neighbor
exchange parameters by mapping low energy magnetic orderings to a Heisenberg
model.

Migrating from HeisenbergMapper 0.1
-----------------------------------
Sublattices used to be derived from the orderings themselves, which restricted
all orderings to a single common supercell. They are now defined once on a
paramagnetic *parent cell* and every ordering is mapped onto it, so orderings in
different supercells can be fitted together. This changed the public surface:

* ``HeisenbergMapper(ordered_structures, energies, cutoff, tol)`` gained a
  ``parent`` argument in third position, so positional ``cutoff``/``tol`` move
  along: ``HeisenbergMapper(ordered_structures, energies, parent, cutoff, tol)``.
  Supply the paramagnetic parent cell whenever you have it; leaving it as None
  infers it from the lowest-energy ordering and warns.
* ``sgraphs`` is gone. Each ordering is now a :class:`RelaxedOrdering` in
  ``mapper.orderings``, carrying its own ``structure``, ``magnetic_structure``,
  ``energy`` and ``coupling_graph``; ``mapper.coupling_graphs`` gives the list of graphs.
  Note these are built on the magnetic-only structure, not the full one, at the
  parent's positions (``ordering.ideal_magnetic_structure``), not the relaxed ones.
* ``nn_interactions`` is now ``interactions`` and maps each sublattice pair
  ``(i, j)`` to its J labels (``'<i>-<j>-nn'``, ``'<i>-<j>-nnn'``, ...), instead
  of mapping ``'nn'``/``'nnn'``/``'nnnn'`` to site pairs.
* ``unique_site_ids`` and ``wyckoff_ids`` are gone. Sublattices are now the
  symmetry orbits of the parent: ``mapper.sublattice_ids`` labels the magnetic
  sites of each ordering and ``mapper.sublattice_wyckoff_symbols`` maps a
  sublattice id to its Wyckoff symbol.
* ``ordered_structures`` (the screened, energy-sorted list) is now
  ``mapper.structures``. ``mapper.energies`` now holds the *total* energies of
  the screened orderings, as passed to the constructor, rather than energies
  per magnetic ion; those moved to ``mapper.energies_per_magnetic_ion``. The
  same split applies to :class:`HeisenbergModel`. The unmodified constructor
  inputs are still ``ordered_structures_``/``energies_``.
* ``get_exchange`` now returns ``(ex_params, residual)`` rather than
  ``ex_params`` alone, and is a least-squares fit over all orderings instead of
  an exactly-determined solve - so supplying more orderings than parameters is
  now useful, and the ``{"<J>": ...}`` fallback for under-determined systems is
  gone. ``residual`` is the RMS fit residual in meV per magnetic ion and is also
  stored on ``mapper.residual`` and on :class:`HeisenbergModel`. It reports an
  ill-conditioned or rank-deficient fit as a ``UserWarning`` rather than through
  the module logger, so the two signals that the returned parameters are
  untrustworthy can be filtered or turned into errors with ``warnings``.
* ``estimate_exchange``, ``get_low_energy_orderings`` and ``get_mft_temperature`` are deprecated. Use
  ``get_exchange`` for shell-resolved ``J_ij``, and a Monte Carlo solver (e.g.
  VAMPIRE, via :class:`HeisenbergModel`) rather than the mean field estimate for
  a critical temperature.
* ``HeisenbergScreener(ordered_structures, energies)`` now takes the list of
  :class:`RelaxedOrdering` objects and exposes ``screened_orderings`` in place
  of ``screened_structures``/``screened_energies``.
* Couplings are found on the parent geometry, mapped into each ordering's
  supercell, not on the relaxed orderings: relaxation of the orderings no longer
  decides which pairs couple or which shell they fall in. An inferred parent is
  itself a relaxed cell, so pass the unrelaxed one to get the full benefit.
"""

from __future__ import annotations

import logging
import warnings
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

__author__ = "Luguza, ncfrey"
__version__ = "0.2"
__maintainer__ = "Luca Frey, Nathan C. Frey"
__email__ = "luca.frey@student.kit.edu, ncfrey@lbl.gov"
__status__ = "Development"
__date__ = "July 2026"

logger = logging.getLogger(__name__)

# Gap (Angstrom) between consecutive coupling lengths that starts a new shell in the same sublattice pair.
DEFAULT_TOL = 0.02

# Symmetry tolerance (Angstrom) for finding the parent's sublattices, as in SpacegroupAnalyzer.
DEFAULT_SYMPREC = 0.01


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


class SublatticeMinimumDistanceNN(MinimumDistanceNN):
    """Nearest-neighbor strategy taking the nearest shell of every sublattice pair.

    Like ``MinimumDistanceNN``, but the (1 + tol) window is applied per neighbor sublattice
    instead of to the site's overall nearest-neighbor distance, so a long intra-sublattice
    bond is not lost to a shorter inter-sublattice one.
    """

    def __init__(self, sublattice_ids: list[int], tol: float = 0.1, cutoff: float = 10) -> None:
        """
        Args:
            sublattice_ids (list[int]): Sublattice id of each site, indexed against the
                structure this strategy will be applied to.
            tol (float): Relative tolerance on neighbor distances being in the same shell,
                applied within a sublattice pair. Defaults to 0.1.
            cutoff (float): Cutoff radius in Angstrom to look for trial neighbors in.
                Defaults to 10.
        """
        super().__init__(tol=tol, cutoff=cutoff, get_all_sites=False)
        self.sublattice_ids = sublattice_ids

    def get_nn_info(self, structure: Structure, n: int) -> list[dict]:
        """Nearest neighbors of site n, one shell per sublattice the neighbors belong to.

        Args:
            structure (Structure): input structure.
            n (int): index of the site to find neighbors of.

        Returns:
            list[dict]: dicts with the neighbor site, its image, weight and site index.
        """
        by_sublattice: dict[int, list] = {}
        for nn in structure.get_neighbors(structure[n], self.cutoff):
            site_index = self._get_original_site(structure, nn)
            by_sublattice.setdefault(self.sublattice_ids[site_index], []).append((nn, site_index))

        siw = []
        for neighbors in by_sublattice.values():
            min_dist = min(nn.nn_distance for nn, _ in neighbors)
            for nn, site_index in neighbors:
                if nn.nn_distance < (1 + self.tol) * min_dist:
                    siw.append(
                        {
                            "site": nn,
                            "image": self._get_image(structure, nn),
                            "weight": min_dist / nn.nn_distance,
                            "site_index": site_index,
                        }
                    )
        return siw 


class MagneticOrdering(ABC):
    """A single collinear magnetic configuration on a common parent lattice.

    Base class holding one ordering's structure, its magnetic-only reduction,
    its coupling graph, and the parent-sublattice labels of its magnetic sites.
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
        # Magnetic-only cell at the parent's positions, site for site like magnetic_structure
        # and with its moments. The coupling graph is built on it.
        self.ideal_magnetic_structure: Structure | None = None

    @cached_property
    def magnetic_structure(self) -> Structure:
        """Magnetic-only cell with sites of the magnetic species, including any whose moment relaxed to zero.

        coupling_graph and sublattice_ids are indexed against this.
        """
        return Structure.from_sites([site for site in self.structure if site.specie.symbol in self.magn_species])

    @staticmethod
    def _nonmagnetic(structure) -> Structure:
        """Nonmagnetic copy of a structure (moments zeroed, all ions kept)."""
        s0 = _analyzer(structure).get_nonmagnetic_structure(make_primitive=False)
        if "wyckoff" in s0.site_properties:
            s0.remove_site_property("wyckoff")
        return s0

    @abstractmethod
    def set_sublattice_ids(self) -> None:
        """Set sublattice_ids (one per site of magnetic_structure), sublattice_wyckoff_symbols
        and ideal_magnetic_structure.
        """

    @cached_property
    def coupling_graph(self) -> StructureGraph:
        """Graph of the coupled site pairs, built on ideal_magnetic_structure.

        With a cutoff, every pair within it; otherwise SublatticeMinimumDistanceNN.
        """
        if self.cutoff:
            strategy = MinimumDistanceNN(cutoff=self.cutoff, get_all_sites=True)
        else:
            strategy = SublatticeMinimumDistanceNN(self.sublattice_ids)
        return StructureGraph.from_local_env_strategy(self.ideal_magnetic_structure, strategy=strategy)


class ParentOrdering(MagneticOrdering):
    """The nonmagnetic parent cell; its symmetry orbits of the magnetic species are the sublattices."""

    def __init__(
        self,
        structure: Structure,
        magn_species: set[str],
        cutoff: float = 0,
        tol: float = DEFAULT_TOL,
        symprec: float = DEFAULT_SYMPREC,
    ):
        super().__init__(self._nonmagnetic(structure), magn_species, cutoff, tol)
        self.symprec = symprec
        self.set_sublattice_ids()

    def set_sublattice_ids(self):
        symmetrized_parent = SpacegroupAnalyzer(self.structure, symprec=self.symprec).get_symmetrized_structure()

        # Full-cell ids, None on nonmagnetic sites; RelaxedOrdering reads them via structure matching.
        full_ids: list[int | None] = [None] * len(self.structure)
        for indices, wyckoff_symbol in zip(
            symmetrized_parent.equivalent_indices, symmetrized_parent.wyckoff_symbols, strict=True
        ):
            if symmetrized_parent[indices[0]].specie.symbol not in self.magn_species:
                continue
            sub_id = len(self.sublattice_wyckoff_symbols)
            self.sublattice_wyckoff_symbols[sub_id] = wyckoff_symbol
            for index in indices:
                full_ids[index] = sub_id

        self.structure.add_site_property("sublattice_id", full_ids)

        self.sublattice_ids = [sub_id for sub_id in full_ids if sub_id is not None]
        self.ideal_magnetic_structure = self.magnetic_structure  # the parent is its own ideal cell


class RelaxedOrdering(MagneticOrdering):
    """One DFT-relaxed magnetic ordering with its energy, labelled by the parent's sublattices."""

    def __init__(
        self,
        structure: Structure,
        energy: float,
        parent_ordering: ParentOrdering,
        cutoff: float = 0,
        tol: float = DEFAULT_TOL,
    ):
        super().__init__(structure, parent_ordering.magn_species, cutoff, tol)
        self.energy = energy  # total energy, as supplied by the caller
        self.parent_ordering = parent_ordering
        self.set_sublattice_ids()

    @property
    def energy_per_magnetic_ion(self) -> float:
        """Total energy divided by the number of magnetic ions (eV)."""
        return self.energy / len(self.magnetic_structure)

    def set_sublattice_ids(self):
        matcher = StructureMatcher(primitive_cell=False, attempt_supercell=True)

        # matched_parent[i] is the parent site, with its 'sublattice_id', that site i sits on.
        matched_parent = matcher.get_s2_like_s1(self._nonmagnetic(self.structure), self.parent_ordering.structure)
        if matched_parent is None:
            raise ValueError(
                "This ordering is not a supercell of the parent cell; it cannot be mapped "
                "onto the parent sublattices. Pass an explicit `parent` cell that all "
                "orderings share."
            )

        # A partial or wrong match would silently shift every label.
        aligned = len(matched_parent) == len(self.structure) and all(
            parent_site.specie.symbol == site.specie.symbol
            for parent_site, site in zip(matched_parent, self.structure, strict=True)
        )
        if not aligned:
            raise ValueError(
                "The parent cell was matched onto this ordering, but the matched sites do "
                "not line up with the ordering's species site by site, so the sublattice "
                "labels would land on the wrong sites. Pass an explicit `parent` cell that "
                "all orderings share."
            )

        ideal = Structure.from_sites([site for site in matched_parent if site.specie.symbol in self.magn_species])
        ideal.add_site_property("magmom", self.magnetic_structure.site_properties["magmom"])
        self.ideal_magnetic_structure = ideal
        self.sublattice_ids = [int(sub_id) for sub_id in ideal.site_properties["sublattice_id"]]
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
            magnetic ion. Each owns its magnetic-only structure, its coupling graph and its
            parent-sublattice labels.
        parent (ParentOrdering): Nonmagnetic parent cell that defines the sublattices.
        interactions (dict): {(i, j): [J label, ...]}, nearest shell first.
        dists (dict): {J label: length (Angstrom) at which its shell starts}.
        ex_mat (DataFrame): Heisenberg Hamiltonian (per magnetic ion) for each ordering.
        ex_params (dict): Exchange parameter values. The J_ij are in meV/muB^2 (they
            multiply the raw moments, see get_exchange); the included 'E0' offset is in
            eV per magnetic ion. get_interaction_graph carries the per-bond J_ij in meV.
    """

    def __init__(
        self,
        ordered_structures,
        energies,
        parent=None,
        cutoff=0,
        tol: float = DEFAULT_TOL,
        symprec: float = DEFAULT_SYMPREC,
    ):
        """Map collinear orderings, possibly in different supercells of a common parent
        cell, onto a classical Heisenberg model.

        Args:
            ordered_structures (list): Structure objects with magmoms.
            energies (list): Total energies of each relaxed magnetic structure.
            parent (Structure): Paramagnetic parent cell whose symmetry defines the
                magnetic sublattices; it is reduced to its primitive cell. If None, the
                lowest-energy ordering is used, which is only correct if relaxation kept
                the parent's symmetry. Defaults to None.
            cutoff (float): Cutoff in Angstrom for the bond search. Defaults to 0,
                which keeps the nearest shell of every sublattice pair.
            tol (float): Gap (Angstrom) between consecutive bond lengths of a sublattice
                pair that starts a new shell. Defaults to 0.02.
            symprec (float): Symmetry tolerance (Angstrom) for finding the parent's
                sublattices. Raise it, or pass a symmetrized parent, if sublattices that
                should be equivalent come out split. Defaults to 0.01.
        """
        if parent is not None and not isinstance(parent, Structure):
            raise TypeError(
                f"parent must be a Structure or None, got {type(parent).__name__}. "
                "Note the constructor signature changed: parent now comes third, "
                "before cutoff/tol - see the module docstring's migration guide."
            )

        # Save original copies of inputs
        self.ordered_structures_ = ordered_structures
        self.energies_ = energies
        self.parent_ = parent

        self.cutoff = cutoff
        self.tol = tol
        self.symprec = symprec

        # These attributes are set by internal methods, listed here for clarity.
        # Set by _initialize_orderings.
        self.orderings = None
        self.parent = None
        # Set by _set_interactions.
        self.interactions: dict[tuple[int, int], list[str]] | None = None
        """J labels of each sublattice pair (i, j), i <= j, nearest shell first.
        E.g. {(0, 1): ['0-1-nn', '0-1-nnn']}."""
        self.dists: dict[str, float] | None = None
        """Length (Angstrom) at which each J label's shell starts. E.g. {'0-1-nn': 3.0, '0-1-nnn': 4.2}."""
        # Set by _build_exchange_mat.
        self.ex_mat = None
        # Set by get_exchange.
        self.ex_params = None
        self.residual = None

        self._initialize_orderings(ordered_structures, energies, parent)
        self._set_interactions()
        self._build_exchange_mat()

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
        """list[float]: Total energy (eV) of each ordering, as supplied to the constructor."""
        return [ordering.energy for ordering in self.orderings]

    @property
    def energies_per_magnetic_ion(self):
        """list[float]: Energy per magnetic ion (eV) of each ordering - the energies the fit uses."""
        return [ordering.energy_per_magnetic_ion for ordering in self.orderings]

    @property
    def coupling_graphs(self):
        """list[StructureGraph]: Coupling graph of each ordering, on its ideal magnetic-only cell."""
        return [ordering.coupling_graph for ordering in self.orderings]

    @property
    def sublattice_ids(self):
        """list[list[int]]: sublattice_ids[k][i] is the sublattice id of magnetic site i
        in ordering k, aligned with that ordering's graph.
        """
        return [ordering.sublattice_ids for ordering in self.orderings]

    @property
    def sublattice_wyckoff_symbols(self):
        """dict[int, str]: Maps each sublattice id to its wyckoff symbol."""
        return self.parent.sublattice_wyckoff_symbols

    def _initialize_orderings(self, ordered_structures, energies, parent):
        """Set self.parent and self.orderings (screened and sorted by energy per magnetic ion).

        Raises:
            ValueError: If fewer than 2 unique orderings remain after screening.
        """
        # Pooled, so an ion whose moment relaxed to zero in one ordering stays on its lattice.
        magn_species = set()
        for struct in ordered_structures:
            for site in _analyzer(struct).structure:
                if site.properties["magmom"] != 0:
                    magn_species.add(site.specie.symbol)

        if parent is None:
            warnings.warn(
                "No `parent` cell supplied; the magnetic sublattices are inferred from the "
                "primitive cell of the lowest-energy ordering. This is only correct if that "
                "ordering still has the symmetry of the paramagnetic parent - relaxation "
                "usually lowers it, which silently splits or merges sublattices. Pass an "
                "explicit `parent` unless you have checked that it does not.",
                UserWarning,
                # _initialize_orderings <- __init__ <- caller
                stacklevel=3,
            )

            # Lowest energy per magnetic ion, rounded and tie-broken like HeisenbergScreener.
            def energy_per_magnetic_ion(idx):
                n_magnetic = sum(site.specie.symbol in magn_species for site in ordered_structures[idx])
                return round(energies[idx] / n_magnetic, 6)

            parent = ordered_structures[min(range(len(ordered_structures)), key=energy_per_magnetic_ion)]

        self.parent = ParentOrdering(
            MagneticOrdering._nonmagnetic(parent).get_primitive_structure(),
            magn_species,
            cutoff=self.cutoff,
            tol=self.tol,
            symprec=self.symprec,
        )

        orderings = [
            RelaxedOrdering(struct, energy, self.parent, cutoff=self.cutoff, tol=self.tol)
            for struct, energy in zip(ordered_structures, energies, strict=True)
        ]

        # Drop duplicate/degenerate orderings and sort by energy per magnetic ion.
        self.orderings = HeisenbergScreener(orderings, screen=False).screened_orderings

        if len(self.orderings) < 2:
            raise ValueError("HeisenbergMapper needs at least 2 unique orderings.")

    def _set_interactions(self):
        """Set self.dists and self.interactions from the parent's coupling graph.

        An interaction is a shell of a sublattice pair, labelled '<i>-<j>-<shell>' with
        i <= j and shell 'nn', 'nnn', ... counted within that pair.
        """
        coupling_graph = self.parent.coupling_graph
        sub_ids = self.parent.sublattice_ids

        pair_dists: dict[tuple[int, int], list[float]] = {}
        for i in range(len(coupling_graph)):
            for neighbor in coupling_graph.get_connected_sites(i):
                sub_id_pair = tuple(sorted((sub_ids[i], sub_ids[neighbor.index])))
                pair_dists.setdefault(sub_id_pair, []).append(neighbor.dist)

        self.dists = {}
        self.interactions = {}
        for sub_id_pair, dists in sorted(pair_dists.items()):
            labels = []
            previous = None
            for dist in sorted(dists):
                # A gap of more than tol to the previous bond starts a shell, named by its shortest bond.
                if previous is None or dist - previous > self.tol:
                    label = f"{sub_id_pair[0]}-{sub_id_pair[1]}-{'n' * (len(labels) + 2)}"
                    self.dists[label] = dist
                    labels.append(label)
                previous = dist
            self.interactions[sub_id_pair] = labels

    def _interaction_label(self, i_id, j_id, dist):
        """Look up the J label of a bond from the sublattices at its ends and its length.

        Each shell starts at its length in self.dists, and the bond gets the last shell
        starting at or below its length.

        Example, with self.dists = {'0-1-nn': 3.0, '0-1-nnn': 4.2}:
            _interaction_label(0, 1, 3.05) -> '0-1-nn'
            _interaction_label(1, 0, 4.2)  -> '0-1-nnn'
            _interaction_label(0, 1, 2.5)  -> ValueError, no shell that short

        Args:
            i_id (int): sublattice id of the ith site
            j_id (int): sublattice id of the jth site
            dist (float): distance (Angstrom) between the sites

        Returns:
            str: '<i>-<j>-<shell>' label, e.g. '0-1-nn'.

        Raises:
            ValueError: If the pair has no shell at or below dist.
        """
        pair = tuple(sorted((i_id, j_id)))
        label = None
        for shell in self.interactions.get(pair, []):  # nearest shell first
            if self.dists[shell] <= dist + 1e-6:  # 1e-6 absorbs float noise
                label = shell
        if label is None:
            raise ValueError(
                f"No interaction of sublattices {i_id} and {j_id} at {dist:.4f} Angstrom in the parent; "
                f"its interactions are {self.interactions}. Couplings must come from a graph built on "
                "the parent geometry (ideal_magnetic_structure)."
            )
        return label

    def _build_exchange_mat(self):
        """Set self.ex_mat: one row per ordering, -sum m_i m_j per J label, normalised to number of magnetic ions.

        Keeps at most (n - 1) J columns for n orderings, dropping those with the longest bond length
        (self.dists) with a UserWarning.
        """
        # J columns ordered by increasing interaction length, so truncation drops the longest.
        j_columns = sorted(self.dists, key=self.dists.get)
        columns = ["E", "E0", *j_columns]

        rows = []
        for ordering in self.orderings:
            coupling_graph = ordering.coupling_graph
            sub_ids = ordering.sublattice_ids
            magmoms = ordering.magnetic_structure.site_properties["magmom"]
            n_sites = len(ordering.magnetic_structure)

            row = dict.fromkeys(columns, 0.0)
            for i in range(len(coupling_graph.graph.nodes)):
                m_i = magmoms[i]
                for neighbor in coupling_graph.get_connected_sites(i):
                    col = self._interaction_label(sub_ids[i], sub_ids[neighbor.index], neighbor.dist)
                    row[col] -= m_i * magmoms[neighbor.index]

            # Extensive sums -> normalised per magnetic ion, with the 1/2 Heisenberg factor for double counting.
            for c in j_columns:
                row[c] /= 2 * n_sites

            row["E0"] = 1.0  # nonmagnetic contribution (per ion)
            row["E"] = ordering.energy_per_magnetic_ion
            rows.append(row)

        ex_mat = pd.DataFrame(rows, columns=columns)

        # Drop interaction columns that never appear (all zero) to keep H full rank, before
        # truncating so they do not use up a slot.
        j_columns = [c for c in j_columns if not (ex_mat[c] == 0).all()]

        # Keep at most (n - 1) J_ij for n orderings (E0 is the nth parameter).
        n_j_max = len(self.orderings) - 1
        if len(j_columns) > n_j_max:
            dropped = j_columns[n_j_max:]
            j_columns = j_columns[:n_j_max]
            remedy = "Supply more orderings or lower the cutoff." if self.cutoff else "Supply more orderings."
            warnings.warn(
                f"{len(self.orderings)} orderings constrain only {n_j_max} exchange interactions; "
                f"left out of the fit and the interaction graph: {dropped}. {remedy}",
                UserWarning,
                # _build_exchange_mat <- __init__ <- caller
                stacklevel=3,
            )

        self.ex_mat = ex_mat[["E", "E0", *j_columns]].reset_index(drop=True)

    def get_exchange(self):
        """Fit E0 and the J_ij to the energies of all orderings by least squares.

        The J_ij multiply the raw moments, so they are in meV/muB^2 (bond energy
        J_ij * m_i * m_j); this keeps moment magnitudes that differ between orderings out
        of the fitted J_ij. get_interaction_graph converts to meV for normalized spins.

        Returns:
            ex_params (dict[str, float]): J_ij in meV/muB^2 and 'E0' in eV per magnetic ion.
            residual (float): RMS fit residual in meV per magnetic ion.

        Raises:
            ValueError: If fewer than two exchange interactions are left to fit.
        """
        ex_mat = self.ex_mat
        E = ex_mat[["E"]]
        col_names = [c for c in ex_mat.columns if c != "E"]

        if len(col_names) < 3:
            raise ValueError(
                f"Exchange matrix holds {len(col_names) - 1} interaction(s) besides E0; a least-squares "
                "fit needs at least 2. Supply more orderings if interactions were left out of the "
                "fit, otherwise set a cutoff so that further shells are included."
            )

        H = np.array(ex_mat.loc[:, ex_mat.columns != "E"].values).astype(float)

        # Column-normalized, so the threshold does not depend on the moment magnitudes.
        cond = np.linalg.cond(H / np.linalg.norm(H, axis=0))
        if cond > 1e5:
            warnings.warn(
                f"Exchange matrix is ill-conditioned (cond={cond:.1e}); the fitted exchange "
                "parameters are unreliable. The input orderings are near-degenerate or the "
                "model has more parameters than the orderings can constrain. Supply more, "
                "more-distinct orderings.",
                UserWarning,
                stacklevel=2,
            )

        j_ij, residuals, rank, _singular = np.linalg.lstsq(H, E, rcond=None)

        # lstsq leaves residuals empty unless the fit is overdetermined and full rank.
        ssr = float(residuals[0]) if residuals.size else float(np.sum((H @ j_ij - np.asarray(E)) ** 2))
        # Divided by the row count, not the degrees of freedom, so a square system works.
        residual = float(np.sqrt(ssr / H.shape[0]))

        # lstsq returns a minimum-norm solution here rather than failing.
        if rank < H.shape[1]:
            warnings.warn(
                f"Exchange matrix is rank deficient (rank {rank} for {H.shape[1]} parameters); "
                "the orderings do not constrain every exchange parameter, and the values "
                "returned are one of infinitely many fits. Supply more distinct orderings.",
                UserWarning,
                stacklevel=2,
            )

        residual *= 1000  # convert to meV per magnetic ion
        ex_params = {
            name: value[0] if name == "E0" else value[0] * 1000  # J_ij in meV/muB^2, E0 in eV per ion
            for name, value in zip(col_names, j_ij.tolist(), strict=True)
        }

        self.ex_params = ex_params
        self.residual = residual
        return self.ex_params, self.residual

    def get_interaction_graph(self, filename=None, ordering_index=0):
        """Get a StructureGraph with edges and weights that correspond to exchange
        interactions and J_ij values, respectively.

        Edge weights are in meV for normalized spins, E = -sum_<ij> J_ij e_i.e_j, as used
        by VAMPIRE, UppASD and TB2J.

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
        structure = ordering.ideal_magnetic_structure.copy()  # the returned graph must not share it
        magmoms = structure.site_properties["magmom"]
        coupling_graph = ordering.coupling_graph
        sub_ids = ordering.sublattice_ids

        igraph = StructureGraph.from_empty_graph(
            structure, edge_weight_name="exchange_constant", edge_weight_units="meV"
        )

        # J_ij exchange interaction matrix
        for i in range(len(coupling_graph.graph.nodes)):
            for neighbor in coupling_graph.get_connected_sites(i):
                j = neighbor.index
                j_exc = self.ex_params.get(self._interaction_label(sub_ids[i], sub_ids[j], neighbor.dist), 0)
                j_exc *= abs(magmoms[i] * magmoms[j])  # meV/muB^2 -> meV
                # Interactions left out of the fit; add_edge would store them with weight None.
                if not j_exc:
                    continue
                igraph.add_edge(i, j, to_jimage=neighbor.jimage, weight=j_exc, warn_duplicates=False)

        if filename:
            if not filename.endswith(".json"):
                filename += ".json"
            dumpfn(igraph, filename)

        return igraph

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
            energies_per_magnetic_ion=self.energies_per_magnetic_ion,
            cutoff=self.cutoff,
            tol=self.tol,
            symprec=self.symprec,
            coupling_graphs=self.coupling_graphs,
            sublattice_ids=self.sublattice_ids,
            sublattice_wyckoff_symbols=self.sublattice_wyckoff_symbols,
            interactions=self.interactions,
            dists=self.dists,
            ex_mat=self.ex_mat,
            ex_params=ex_params,
            residual=residual,
            igraph=self.get_interaction_graph(),
        )

    @deprecated(
        get_exchange,
        message="It only serves the deprecated estimate_exchange.",
        category=DeprecationWarning,
        deadline=(2027, 8, 1),
    )
    def get_low_energy_orderings(self):
        """Find lowest energy FM and AFM orderings to compute E_AFM - E_FM.

        .. deprecated::
            Only used by the deprecated :meth:`estimate_exchange`.

        Returns:
            fm_struct (Structure): fm structure with 'magmom' site property
            afm_struct (Structure): afm structure with 'magmom' site property
            fm_e (float): fm energy
            afm_e (float): afm energy
        """
        fm_struct = None
        afm_struct = None
        mag_min = np.inf
        mag_max = 0.001
        fm_e = 0
        afm_e = 0
        fm_e_min = 0
        afm_e_min = 0

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
        energies_per_magnetic_ion=None,
        cutoff=None,
        tol=None,
        symprec=None,
        coupling_graphs=None,
        sublattice_ids=None,
        sublattice_wyckoff_symbols=None,
        interactions=None,
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
            magnetic_structures (list): Magnetic-only cell of each ordering. coupling_graphs and
                sublattice_ids are indexed against these, not against structures. The
                coupling graphs hold the same sites at the parent's positions.
            energies (list): Total energy (eV) of each relaxed magnetic structure.
            energies_per_magnetic_ion (list): Energy per magnetic ion (eV) of each relaxed
                magnetic structure, the energies the exchange parameters were fitted to.
            cutoff (float): Cutoff in Angstrom for nearest neighbor search.
            tol (float): Gap (Angstrom) between consecutive bond lengths that separates shells.
            symprec (float): Symmetry tolerance (Angstrom) the parent's sublattices were found with.
            coupling_graphs (list): Coupling graph of each ordering, built on its magnetic
                sites at the parent's positions.
            sublattice_ids (list[list[int]]): sublattice_ids[k][i] is the sublattice id of
                site i in ordering k.
            sublattice_wyckoff_symbols (dict): Maps each sublattice id to its wyckoff symbol.
            interactions (dict): {(i, j): [J label, ...]} - the distinct interactions
                of each sublattice pair, ordered from the nearest shell outwards.
            dists (dict): {J label: interaction distance in Angstrom}.
            ex_mat (DataFrame): Heisenberg Hamiltonian (per magnetic ion) for each ordering.
            ex_params (dict): Exchange parameter values. The J_ij are in meV/muB^2 (they
                multiply the raw moments, see HeisenbergMapper.get_exchange); the included
                'E0' offset is in eV per magnetic ion.
            residual (float): RMS residual of the fit that produced ex_params, in meV per
                magnetic ion.
            igraph (StructureGraph): Exchange interaction graph, edge weights in meV.
        """
        self.formula = formula
        self.structures = structures
        self.magnetic_structures = magnetic_structures
        self.energies = energies
        self.energies_per_magnetic_ion = energies_per_magnetic_ion
        self.cutoff = cutoff
        self.tol = tol
        self.symprec = symprec
        self.coupling_graphs = coupling_graphs
        self.sublattice_ids = sublattice_ids
        self.sublattice_wyckoff_symbols = sublattice_wyckoff_symbols
        self.interactions = interactions
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
            "energies_per_magnetic_ion": self.energies_per_magnetic_ion,
            "cutoff": self.cutoff,
            "tol": self.tol,
            "symprec": self.symprec,
            "coupling_graphs": [coupling_graph.as_dict() for coupling_graph in self.coupling_graphs],
            "sublattice_ids": self.sublattice_ids,
            "dists": self.dists,
            "ex_params": self.ex_params,
            "residual": self.residual,
            "igraph": self.igraph.as_dict(),
            # Sanitize int keys / DataFrame
            "ex_mat": jsanitize(self.ex_mat.to_dict()),
            "interactions": jsanitize(self.interactions),
            "sublattice_wyckoff_symbols": jsanitize(self.sublattice_wyckoff_symbols),
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Create a HeisenbergModel from a dict."""
        if "sublattice_ids" not in dct or "magnetic_structures" not in dct:
            raise ValueError(
                f"This dict was serialized with HeisenbergModel {dct.get('@version', '<0.2')}, which "
                "predates the parent-cell sublattice refactor (see this module's migration guide) and "
                "cannot be loaded by this version - it is missing `sublattice_ids`/`magnetic_structures` "
                "(pre-0.2 dicts have `unique_site_ids`/`wyckoff_ids`/`sgraphs` instead, which don't map "
                "onto the new parent-cell sublattices). Recompute the HeisenbergModel from the original "
                "orderings with the current HeisenbergMapper."
            )

        # Reconstitute the tuple/int-keyed dicts that jsanitize stringified
        interactions = {literal_eval(pair): labels for pair, labels in dct["interactions"].items()}
        sublattice_wyckoff_symbols = {literal_eval(k): v for k, v in dct["sublattice_wyckoff_symbols"].items()}

        structures = [Structure.from_dict(v) for v in dct["structures"]]
        magnetic_structures = [Structure.from_dict(v) for v in dct["magnetic_structures"]]
        coupling_graphs = [StructureGraph.from_dict(v) for v in dct["coupling_graphs"]]
        igraph = StructureGraph.from_dict(dct["igraph"])

        # Older serializations store ex_mat as a string.
        ex_mat = dct["ex_mat"]
        if isinstance(ex_mat, str):
            try:
                ex_mat = literal_eval(ex_mat)
            except (SyntaxError, ValueError):  # empty or unparsable string
                ex_mat = None
        ex_mat = pd.DataFrame.from_dict(ex_mat) if ex_mat else pd.DataFrame(columns=["E", "E0"])

        return cls(
            formula=dct["formula"],
            structures=structures,
            magnetic_structures=magnetic_structures,
            energies=dct["energies"],
            energies_per_magnetic_ion=dct["energies_per_magnetic_ion"],
            cutoff=dct["cutoff"],
            tol=dct["tol"],
            symprec=dct["symprec"],
            coupling_graphs=coupling_graphs,
            sublattice_ids=dct["sublattice_ids"],
            sublattice_wyckoff_symbols=sublattice_wyckoff_symbols,
            interactions=interactions,
            dists=dct["dists"],
            ex_mat=ex_mat,
            ex_params=dct["ex_params"],
            residual=dct["residual"],
            igraph=igraph,
        )
