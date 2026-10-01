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
* ``estimate_exchange`` and ``get_mft_temperature`` are deprecated. Use
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

# Tolerance (Angstrom) within which two coupling distances count as the same shell. Shells
# are read off the parent, where symmetry-equivalent couplings have identical lengths, so
# tol only has to separate genuinely distinct distances.
DEFAULT_TOL = 0.02


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

    ``MinimumDistanceNN`` is sublattice-blind: it keeps only the neighbors closer than
    (1 + tol) times a site's *overall* nearest-neighbor distance. A sublattice pair whose
    bond is longer than that then never enters the graph, so a structure whose sublattices
    sit close together but are internally widely spaced yields only the inter-sublattice
    bonds and no intra-sublattice ones.

    Here the trial neighbors of a site are grouped by the sublattice they belong to and the
    same window is applied within each group, so every sublattice pair the lattice realizes
    within ``cutoff`` contributes its nearest shell, whatever the other pairs do.
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
        # The site's own sublattice is just another group here, so intra-sublattice bonds
        # survive even when they are much longer than the shortest inter-sublattice one.
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
        """Magnetic-only cell. coupling_graph and sublattice_ids are both indexed against this.

        Membership is decided by *species*, not by moment, using the magnetic species
        pooled over every ordering (see ``HeisenbergMapper._initialize_orderings``). A
        magnetic ion whose moment happens to relax to zero therefore stays on the lattice
        and contributes a zero term, rather than dropping out and changing this ordering's
        site count, graph topology and per-ion energy relative to its siblings.
        """
        return self._magnetic_only(self.structure)

    def _magnetic_only(self, structure: Structure) -> Structure:
        """Copy of structure keeping only the sites of the pooled magnetic species."""
        return Structure.from_sites([site for site in structure if site.specie.symbol in self.magn_species])

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
    def coupling_graph(self) -> StructureGraph:
        """Coupling graph of the magnetic-only structure: an edge for every pair of sites
        whose exchange coupling enters the Heisenberg model.

        If self.cutoff is set, the graph includes all neighbors within that distance.
        Otherwise each site keeps the nearest shell of every sublattice pair it takes part
        in, so no pair drops out of the graph for being longer-ranged than another
        (see SublatticeMinimumDistanceNN).
        """
        # Cached, so a graph built before the labels exist would be silently wrong.
        if self.ideal_magnetic_structure is None:
            raise ValueError("Call set_sublattice_ids() before the coupling graph is built.")
        if self.cutoff:
            strategy = MinimumDistanceNN(cutoff=self.cutoff, get_all_sites=True)
        else:
            strategy = SublatticeMinimumDistanceNN(self.sublattice_ids)
        return StructureGraph.from_local_env_strategy(self.ideal_magnetic_structure, strategy=strategy)

    @abstractmethod
    def set_sublattice_ids(self) -> None:
        """
        Sets the following attributes for this ordering:
            - `self.sublattice_ids` = [sublattice id for each *magnetic* site], aligned with
              `self.magnetic_structure` and hence with `self.coupling_graph`. Nonmagnetic sites are
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
        self.ideal_magnetic_structure = self.magnetic_structure  # the parent is its own ideal cell


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

        # get_s2_like_s1 returns the parent's sites in this ordering's site order, but only
        # those it could match, and it matches on geometry within StructureMatcher's
        # tolerances. On relaxed cells a partial or wrong match would shift every label -
        # and every J_ij derived from it - without erroring, so require the two to line up
        # species by species.
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

        # matched_parent is the parent supercell in this ordering's site order, each site in
        # the periodic image of the ordering site it matches. Its magnetic sites are thus this
        # ordering's magnetic sites at their parent positions, with their sublattice ids.
        ideal = self._magnetic_only(matched_parent)
        # Replace the parent's zero moments with this ordering's, so the graph built on it
        # carries the spins it describes.
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
        interactions (dict): {(i, j): [J label, ...]} - the distinct interactions of
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
        Not supplying a parent cell is only safe if the lowest-energy ordering still
        preserves the symmetry of the paramagnetic parent.
        In most cases, relaxation will lower the symmetry and the parent cell must be supplied explicitly.

        Args:
            ordered_structures (list): Structure objects with magmoms.
            energies (list): Total energies of each relaxed magnetic structure.
            parent (Structure): Paramagnetic parent cell whose symmetry defines the
                magnetic sublattices. Reduced to its primitive cell either way, so a
                deliberately-supplied supercell parent is analyzed on its primitive
                cell too, not as given. If None, it is inferred as the primitive cell
                of the lowest-energy ordering. Defaults to None.
            cutoff (float): Cutoff in Angstrom for the bond search. Defaults to 0,
                which keeps the nearest shell of every sublattice pair; with a cutoff,
                every bond up to it is kept.
            tol (float): Bond lengths of a sublattice pair within tol (Angstrom) of a
                shell's shortest bond belong to that shell.
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

        # These attributes are set by internal methods, listed here for clarity.
        self.orderings = self.parent = None  # set by _initialize_orderings
        self.interactions = self.dists = None  # set by _set_interactions
        self.ex_mat = self.ex_params = self.residual = None  # set by _build_exchange_mat and get_exchange

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
         - Build a RelaxedOrdering for each ordering, with its magnetic-only structure and
           coupling graph. Since the parent is not yet known, the sublattice ids are not set yet.
         - Drop duplicate/degenerate orderings and sort by energy per magnetic ion using
           HeisenbergScreener.
         - Build the ParentOrdering from the lowest-energy ordering (or the explicit parent
           if supplied), and set its sublattice ids.

        Args:
            ordered_structures (list): Structure objects with magmoms.
            energies (list): Total energies of each relaxed magnetic structure.
            parent (Structure | None): Explicit paramagnetic parent cell, or None to infer
                it from the lowest-energy ordering.

        Raises:
            ValueError: If fewer than 2 unique orderings remain after screening.
        """
        # Pool the magnetic species over all orderings, to make sure a species that relaxes to
        # zero moment in one ordering still counts as magnetic there. Since the magmom of this
        # species is zero, it contributes a zero term to the Heisenberg Hamiltonian, but it stays
        # on the lattice and keeps the site count and graph topology consistent with the others.
        magn_species = set().union(*(MagneticOrdering._magnetic_species(struct) for struct in ordered_structures))

        orderings = [
            RelaxedOrdering(struct, energy, magn_species, cutoff=self.cutoff, tol=self.tol)
            for struct, energy in zip(ordered_structures, energies, strict=True)
        ]

        # Drop duplicate/degenerate orderings and sort by energy per magnetic ion.
        self.orderings = HeisenbergScreener(orderings, screen=False).screened_orderings

        if len(self.orderings) < 2:
            raise ValueError("HeisenbergMapper needs at least 2 unique orderings.")

        # The nonmagnetic ions are kept in the parent: site equivalence is read from its
        # symmetry, and removing them first can raise the apparent site symmetry and merge
        # sublattices that are actually distinct.
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
        """list[float]: Total energy (eV) of each ordering, as supplied to the constructor."""
        return [ordering.energy for ordering in self.orderings]

    @property
    def energies_per_magnetic_ion(self):
        """list[float]: Energy per magnetic ion (eV) of each ordering - the energies the fit uses.

        The magnetic ions are counted over the magnetic species pooled across all orderings,
        so an ion that relaxed to zero moment in one ordering still counts there.
        """
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
    def parent_coupling_graph(self):
        """StructureGraph: Coupling graph of the magnetic-only parent structure."""
        return self.parent.coupling_graph

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
        """Returns the sublattice_ids in the order (i, j) with i <= j.

        This is the key used to look up the interactions of a sublattice pair.
        """
        return tuple(sorted((i_id, j_id)))

    def _interaction_label(self, i_id, j_id, dist):
        """Return the J label of a coupling: the shell of its sublattice pair that dist falls in.

        Every coupling graph is built on the parent geometry (see ideal_magnetic_structure),
        so dist is one of the parent's own distances up to float noise. It belongs to the last
        shell starting at or below it, the rule _set_interactions grouped it by.

        Args:
            i_id (int): sublattice id of the ith site
            j_id (int): sublattice id of the jth site
            dist (float): distance (Angstrom) between the sites

        Returns:
            str: '<i>-<j>-<shell>' label, e.g. '0-1-nn'.

        Raises:
            ValueError: If the parent has no coupling of this sublattice pair at or below
                dist, i.e. the coupling did not come from the parent geometry.
        """
        labels = self.interactions.get(self._order_sublattice_ids(i_id, j_id), [])
        label = next((label for label in reversed(labels) if self.dists[label] <= dist + 1e-6), None)
        if label is None:
            raise ValueError(
                f"No interaction of sublattices {i_id} and {j_id} at {dist:.4f} Angstrom in the parent; "
                f"its interactions are {self.interactions}. Couplings must come from a graph built on "
                "the parent geometry (ideal_magnetic_structure)."
            )
        return label

    def _set_interactions(self):
        """Set self.dists and self.interactions describing the distinct interactions.

        An interaction is a (sublattice pair, neighbor shell) combination, labelled
        '<i>-<j>-<shell>' with i <= j and shell one of 'nn', 'nnn', 'nnnn', ... Shells are
        counted *within* a sublattice pair, not globally: '0-0-nn' and '0-1-nn' are the
        nearest 0-0 and 0-1 interaction respectively, even when one of them is the longer interaction.

        With cutoff=0 the parent graph holds the nearest shell of every sublattice pair (see
        SublatticeMinimumDistanceNN); with cutoff > 0 every bond up to the cutoff. Either way a
        pair's bond lengths are grouped into shells: a bond within tol of the current shell's
        shortest bond joins it, a longer one starts the next shell.

        Distances and connectivity are read from the parent ordering; see
        _initialize_orderings() for how the parent is defined and how the sublattices are
        labelled. The connectivity comes from the parent's StructureGraph, which is built
        from the magnetic-only parent structure - the nonmagnetic ions are ignored for the
        graph, but kept in the parent structure to preserve the true site symmetry.
        """
        coupling_graph = self.parent.coupling_graph
        sub_ids = self.parent.sublattice_ids

        # Bond lengths of each sublattice pair. Every site is visited, since with cutoff=0 the
        # two ends of an interaction need not both count it as a nearest neighbor.
        pair_dists: dict[tuple[int, int], list[float]] = {}
        for i in range(len(coupling_graph)):
            for conn_site in coupling_graph.get_connected_sites(i):
                pair = self._order_sublattice_ids(sub_ids[i], sub_ids[conn_site[2]])
                pair_dists.setdefault(pair, []).append(conn_site[-1])

        self.dists = {}
        self.interactions = {}
        for pair, dists in sorted(pair_dists.items()):
            labels = []
            for dist in sorted(dists):
                # The shortest bond names the shell.
                if labels and dist - self.dists[labels[-1]] <= self.tol:
                    continue
                label = f"{pair[0]}-{pair[1]}-{'n' * (len(labels) + 2)}"
                self.dists[label] = dist
                labels.append(label)
            self.interactions[pair] = labels

    def _build_exchange_mat(self):
        """Build the Heisenberg Hamiltonian, one row per ordering, by summing the
        signed products S_i . S_j over each graph. Sets self.ex_mat.

        Each row is normalised per magnetic ion so that orderings living in
        different-sized supercells share a single linear system (the energies are
        per magnetic ion too, see RelaxedOrdering.energy_per_magnetic_ion).

        n orderings constrain at most E0 and n - 1 J_ij; any longer-ranged interactions are
        left out of the fit with a UserWarning.
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
                s_i = magmoms[i]
                for conn_site in coupling_graph.get_connected_sites(i):
                    col = self._interaction_label(sub_ids[i], sub_ids[conn_site[2]], conn_site[-1])
                    row[col] -= s_i * magmoms[conn_site[2]]

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

        # Keep at most n - 1 J_ij for n orderings (E0 is the nth parameter).
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

        # Every ordering is kept: get_exchange fits the parameters by least squares, so
        # surplus orderings average out the noise on the energies instead of being
        # discarded to square the system. That includes orderings whose rows coincide while
        # their energies do not - HeisenbergScreener has already dropped the ones that
        # agree on both, so what is left is the Heisenberg model failing to tell two
        # orderings apart, which is precisely what the residual is there to report.
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
            residual (float): Root-mean-square residual of the least-squares fit, in meV
                per magnetic ion. It is a residual on the fitted energies, so unlike the
                J_ij it keeps the per-magnetic-ion normalisation, and averaging over the
                orderings instead of summing makes it independent of how many orderings
                went into the fit. Both normalisations together make it an intensive
                measure of how well the Heisenberg model describes the energies, i.e.
                one that is comparable between materials.

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
                "fit needs at least 2. Supply more orderings if interactions were left out of the "
                "fit, otherwise set a cutoff so that further shells are included."
            )

        # Fit E0 and the J_ij to the energies of every ordering
        H = np.array(ex_mat.loc[:, ex_mat.columns != "E"].values).astype(float)

        # Warn when the fit is ill-conditioned: near-degenerate orderings or an
        # over-parameterized model make H nearly singular, so tiny energy
        # differences blow up into unphysical exchange parameters. Judge that on
        # column-normalized H - the E0 column is exactly 1 while the J columns are sums of
        # m_i * m_j, so the condition number of H itself tracks the moment magnitudes as
        # much as the degeneracy of the orderings, and a fixed threshold on it would fire
        # or not depending on the species involved.
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

        # lstsq only fills residuals for an overdetermined full-rank fit; for a square,
        # underdetermined or rank-deficient one it returns an empty array, so compute the
        # sum of squared residuals here. It sits on the energy side of the system, not the
        # J_ij side, so it is a squared energy in (eV per magnetic ion)^2.
        ssr = float(residuals[0]) if residuals.size else float(np.sum((H @ j_ij - np.asarray(E)) ** 2))

        # Take the root mean square over the orderings rather than the raw sum: the sum
        # grows with the number of orderings, while the RMS stays a typical per-ordering
        # energy error and can be compared between materials fitted from different numbers
        # of orderings. Divided by the row count, not by the degrees of freedom, so that a
        # square system (as many orderings as parameters) does not divide by zero.
        residual = float(np.sqrt(ssr / H.shape[0]))

        # A least-squares fit does not fail on a system it cannot determine,
        # it silently returns the minimum-norm solution. cond above cannot see this case:
        # with fewer rows than parameters it stays finite.
        if rank < H.shape[1]:
            warnings.warn(
                f"Exchange matrix is rank deficient (rank {rank} for {H.shape[1]} parameters); "
                "the orderings do not constrain every exchange parameter, and the values "
                "returned are one of infinitely many fits. Supply more distinct orderings.",
                UserWarning,
                stacklevel=2,
            )

        residual *= 1000  # convert to meV per magnetic ion
        # Keyed by column name rather than by position, so that reordering ex_mat's columns
        # cannot silently convert the offset and leave a J_ij in eV.
        ex_params = {
            name: value[0] if name == "E0" else value[0] * 1000  # J_ij in meV/muB^2, E0 in eV per ion
            for name, value in zip(col_names, j_ij.tolist(), strict=True)
        }

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
        # The ideal cell: the edges come from it, and it keeps the ucf geometry consistent
        # with the fit rather than with one ordering's relaxation.
        structure = ordering.ideal_magnetic_structure.copy()  # the returned graph must not share it
        magmoms = structure.site_properties["magmom"]
        coupling_graph = ordering.coupling_graph
        sub_ids = ordering.sublattice_ids

        igraph = StructureGraph.from_empty_graph(
            structure, edge_weight_name="exchange_constant", edge_weight_units="meV"
        )

        # J_ij exchange interaction matrix
        for i in range(len(coupling_graph.graph.nodes)):
            for c in coupling_graph.get_connected_sites(i):
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
            dist (float): distance (Angstrom) between sites

        Returns:
            float: Exchange parameter J_exc in meV/muB^2 (0 if the interaction was left
                out of the fit).
                Multiply by the two moments to get meV, as get_interaction_graph does.
        """
        label = self._interaction_label(i_id, j_id, dist)

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
            energies_per_magnetic_ion=self.energies_per_magnetic_ion,
            cutoff=self.cutoff,
            tol=self.tol,
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
            tol (float): Tolerance (in Angstrom) on bond lengths being in the same shell.
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
            residual (float): Root-mean-square residual of the fit that produced ex_params,
                in meV per magnetic ion. Intensive in both the cell size and the number of
                orderings, so it is comparable between materials.
            igraph (StructureGraph): Exchange interaction graph, edge weights in meV.
        """
        self.formula = formula
        self.structures = structures
        self.magnetic_structures = magnetic_structures
        self.energies = energies
        self.energies_per_magnetic_ion = energies_per_magnetic_ion
        self.cutoff = cutoff
        self.tol = tol
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

        # Reconstitute the exchange matrix DataFrame. as_dict() calls .to_dict() on it,
        # but serializations written before that did (and by #4664, which jsanitized the
        # DataFrame directly) may store a (JSON/repr) string instead. Accept both forms
        # and fall back to an empty matrix when ex_mat is empty.
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
