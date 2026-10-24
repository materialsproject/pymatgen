"""
Fit the exchange constants of a classical Heisenberg model to the energies of
collinear magnetic orderings.

The model lives on a paramagnetic *parent cell* (:class:`ParentCell`): its
symmetry orbits of sites are the magnetic sublattices, and its symmetry orbits
of bonds are the :class:`Interaction` s, one exchange constant J each. Every
ordering (:class:`ComputedOrdering`) is mapped onto the parent and contributes
one row of spin products per interaction; :class:`HeisenbergMapper` fits the
J's to all rows by least squares.

Migrating from HeisenbergMapper 0.1
-----------------------------------
Sublattices used to be derived from the orderings themselves, which restricted
all orderings to a single common supercell. They are now defined once on a
paramagnetic *parent cell* and every ordering is mapped onto it, so orderings in
different supercells can be fitted together. This changed the public surface:

* ``HeisenbergMapper(ordered_structures, energies, cutoff, tol)`` gained a
  ``parent`` argument in third position, so positional ``cutoff`` moves along:
  ``HeisenbergMapper(ordered_structures, energies, parent, cutoff)``. Supply the
  paramagnetic parent cell whenever you have it; leaving it as None infers it
  from the lowest-energy ordering and warns. ``tol`` is deprecated and ignored.
* ``sgraphs`` is gone. Each ordering is now a :class:`ComputedOrdering` in
  ``mapper.orderings``, carrying its own ``structure``, ``magnetic_structure``,
  ``energy`` and ``coupling_graph``; ``mapper.coupling_graphs`` gives the list of graphs.
  Note these are built on the magnetic-only structure, not the full one, at the
  parent's positions (``ordering.ideal_magnetic_structure``), not the relaxed ones.
* J labels are now ``'J1(Fe1-Mn1)'``, ``'J2(Fe1-Fe1)'``, ...: every
  interaction is numbered by length across all sublattice pairs, and the pair
  is named per element as in a CIF file. They replace ``'<i>-<j>-nn'``,
  ``'<i>-<j>-nnn'``, ... ``nn_interactions`` is now ``interactions`` and maps
  each sublattice pair ``(i, j)`` to its J labels instead of mapping
  ``'nn'``/``'nnn'``/``'nnnn'`` to site pairs.
* ``unique_site_ids`` and ``wyckoff_ids`` are gone. Sublattices are now the
  symmetry orbits of the parent: ``mapper.sublattice_ids`` labels the magnetic
  sites of each ordering, and ``mapper.sublattice_names`` and
  ``mapper.sublattice_wyckoff_symbols`` map a sublattice id to its name
  (``'Fe1'``) and its Wyckoff symbol.
* ``ordered_structures`` (the screened, energy-sorted list) is now
  ``mapper.structures``. ``mapper.energies`` now holds the *total* energies of
  the screened orderings, as passed to the constructor, rather than energies
  per magnetic ion; those moved to ``mapper.energies_per_magnetic_ion``. The
  same split applies to :class:`HeisenbergModel`. The unmodified constructor
  inputs are still ``ordered_structures_``/``energies_``.
* ``get_exchange`` now returns ``(ex_params, residual)`` rather than
  ``ex_params`` alone. It is a least-squares fit over all orderings instead of
  an exactly-determined solve. The ``{"<J>": ...}`` fallback for a single
  interaction is gone: that J is now fitted together with E0. ``residual`` is
  the RMS fit residual in meV per magnetic ion and is also stored on
  ``mapper.residual`` and on :class:`HeisenbergModel`. It is only meaningful
  with substantially more orderings than fitted parameters; otherwise the extra
  parameters overfit, lowering the residual without improving the fit.
  ``get_exchange`` reports an ill-conditioned fit with a ``UserWarning``
  instead of the module logger.
* Interactions are no longer cut to n - 1 for n orderings, longest first.
  Instead, each interaction the orderings cannot separate from E0 and the
  shorter interactions is left out, with a ``UserWarning``.
* ``estimate_exchange``, ``get_low_energy_orderings`` and ``get_mft_temperature`` are deprecated. Use
  ``get_exchange`` for bond-resolved ``J_ij``, and a Monte Carlo solver (e.g.
  VAMPIRE, via :class:`HeisenbergModel`) rather than the mean field estimate for
  a critical temperature.
* ``HeisenbergScreener(ordered_structures, energies)`` now takes the list of
  :class:`ComputedOrdering` objects and exposes ``screened_orderings`` in place
  of ``screened_structures``/``screened_energies``.
* Couplings are found on the parent geometry, mapped into each ordering's
  supercell, not on the relaxed orderings: relaxation of the orderings no longer
  decides which pairs couple or which J they belong to. An inferred parent is
  itself a relaxed cell, so pass the unrelaxed one to get the full benefit.
* Every sublattice pair couples through at least its nearest bonds, also with a
  ``cutoff``. Without a cutoff, the nearest bonds are taken per pair rather than
  once per site, and are only the shortest ones, not all within 10% of them.
* Bonds share an exchange constant only if a symmetry operation of the parent
  maps one onto the other, not because they have the same length. Bonds of
  equal length that are not symmetry-equivalent get separate J's.
"""

from __future__ import annotations

import logging
import warnings
from ast import literal_eval
from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import cached_property
from itertools import combinations_with_replacement
from typing import TYPE_CHECKING, NamedTuple

import numpy as np
import pandas as pd
from monty.dev import deprecated
from monty.json import MSONable, jsanitize
from monty.serialization import dumpfn

from pymatgen.analysis.graphs import StructureGraph
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

# Distance (Angstrom) up to which the parent looks for a sublattice pair's shortest bonds.
MAX_SEARCH_DIST = 10

# Float noise (Angstrom) absorbed when comparing bond lengths, so the parent and every ordering agree.
DIST_EPSILON = 1e-6

# Symmetry tolerance (Angstrom) for finding the parent's sublattices, as in SpacegroupAnalyzer.
DEFAULT_SYMPREC = 0.01

# 'sublattice_id' of nonmagnetic parent sites. Not None, which Structure.from_sites warns about.
_NONMAGNETIC_SUBLATTICE_ID = -1

# Moment (muB) above which zeroing a species outside default_magmoms warns; induced moments on
# anions typically stay below it.
_ZEROED_MOMENT_WARNING_THRESHOLD = 0.5

# Condition number of the column-normalized exchange matrix above which the fit warns.
_MAX_CONDITION_NUMBER = 1e5

# Converts the fitted energies from eV to meV.
_MEV_PER_EV = 1000

# Decimals to which two orderings' energies and spin products must agree to count as duplicates.
_DUPLICATE_DECIMALS = 6


def _analyzer(structure: Structure, **kwargs) -> CollinearMagneticStructureAnalyzer:
    """CollinearMagneticStructureAnalyzer factory with the settings used throughout this module.

    Oxidation states are removed first: the analyzer matches each site's species string
    against default_magmoms, so 'Fe2+' would not match the key 'Fe' and its moment would be
    zeroed. ``make_primitive=False`` keeps the cell as given, ``threshold=0.0`` retains any
    nonzero moment, and ``threshold_nonmag=100.0`` always zeroes out nonmagnetic ions, so
    induced magnetic moments don't change the number of magnetic sites per unit cell between
    orderings.

    Extra keyword arguments are passed through to the analyzer.
    """
    structure = structure.copy()
    # Species spins would be dropped with the oxidation states, so move them to 'magmom' first.
    spins = [getattr(specie, "spin", None) for specie in structure.species]
    if any(spins) and "magmom" not in structure.site_properties:
        structure.add_site_property("magmom", [spin or 0 for spin in spins])

    # Last value wins on conflicting keys, so the defaults above are overridden by kwargs.
    return CollinearMagneticStructureAnalyzer(
        structure.remove_oxidation_states(),
        **{"make_primitive": False, "threshold": 0.0, "threshold_nonmag": 100.0, **kwargs},
    )


def _nonmagnetic(structure: Structure) -> Structure:
    """Nonmagnetic copy of a structure (moments zeroed, all ions kept)."""
    nonmagnetic = _analyzer(structure).get_nonmagnetic_structure(make_primitive=False)
    if "wyckoff" in nonmagnetic.site_properties:
        nonmagnetic.remove_site_property("wyckoff")
    return nonmagnetic


def _magnetic_sites(sites, magnetic_species: set[str]) -> list:
    """The sites whose element is in magnetic_species, in their original order."""
    return [site for site in sites if site.specie.symbol in magnetic_species]


def _as_image(cells) -> tuple[int, int, int]:
    """A whole-cell lattice translation as a tuple of ints."""
    return tuple(int(cell) for cell in cells)


class Bond(NamedTuple):
    """Bond from parent site start to parent site end, image parent cells away."""

    start: int
    end: int
    image: tuple[int, int, int]

    def reversed(self) -> Bond:
        """The same bond, seen from its end site."""
        return Bond(self.end, self.start, (-self.image[0], -self.image[1], -self.image[2]))


@dataclass(frozen=True)
class Interaction:
    """One exchange constant J: a symmetry orbit of parent bonds, in both directions."""

    label: str  # e.g. 'J1(Fe1-Mn1)': numbered by length, then the names of the two sublattices
    sublattices: tuple[int, int]  # i <= j
    length: float  # Angstrom, mean over the orbit; symmetry makes its bonds equally long, up to symprec
    bonds: frozenset[Bond]


class ParentCell:
    """Paramagnetic parent cell: its site orbits are the sublattices, its bond orbits the interactions."""

    def __init__(
        self,
        structure: Structure,
        magnetic_species: set[str],
        cutoff: float = 0,
        symprec: float = DEFAULT_SYMPREC,
    ):
        """
        Args:
            structure (Structure): Parent cell. Moments are dropped and it is reduced to its
                primitive cell.
            magnetic_species (set[str]): Elements whose sites form the magnetic sublattices.
            cutoff (float): Every bond within this distance (Angstrom) couples. A sublattice pair
                without such a bond couples through its shortest bonds. Defaults to 0, which
                keeps only the shortest bonds of each pair.
            symprec (float): Symmetry tolerance (Angstrom) for the sublattices and interactions.
                Defaults to 0.01.
        """
        self.structure = _nonmagnetic(structure).get_primitive_structure()
        self.magnetic_species = magnetic_species
        self.cutoff = cutoff
        self.symprec = symprec

        analyzer = SpacegroupAnalyzer(self.structure, symprec=symprec)
        self._symmetry_ops = analyzer.get_symmetry_operations()

        # Sublattices: the symmetry orbits of the magnetic sites. StructureMatcher carries the
        # 'sublattice_id' and 'parent_index' site properties into every ordering.
        sublattice_ids, self.sublattice_names, self.sublattice_wyckoff_symbols = self._sublattices(
            analyzer.get_symmetrized_structure()
        )
        self.structure.add_site_property("sublattice_id", sublattice_ids)
        self.structure.add_site_property("parent_index", list(range(len(self.structure))))

        # Interactions: the symmetry orbits of the coupled bonds, one exchange constant each.
        self.interactions = self._interactions()
        self._interaction_of = {bond: interaction for interaction in self.interactions for bond in interaction.bonds}
        self._max_bond_length = max((self._length(bond) for bond in self._interaction_of), default=0.0)

        # Coverage: every sublattice pair should have at least one interaction.
        coupled_pairs = {interaction.sublattices for interaction in self.interactions}
        all_pairs = combinations_with_replacement(sorted(self.sublattice_names.keys()), 2)
        uncoupled_pairs = [pair for pair in all_pairs if pair not in coupled_pairs]
        if uncoupled_pairs:
            warnings.warn(
                f"The sublattice pairs {uncoupled_pairs} have no neighbor within "
                f"{max(cutoff, MAX_SEARCH_DIST)} Angstrom, so the model has no exchange "
                "constant between them.",
                UserWarning,
                # ParentCell.__init__ <- HeisenbergMapper.__init__ <- caller
                stacklevel=3,
            )

    @cached_property
    def magnetic_structure(self) -> Structure:
        """Magnetic sites of the parent; sublattice_ids is indexed against it."""
        return Structure.from_sites(_magnetic_sites(self.structure, self.magnetic_species))

    @property
    def sublattice_ids(self) -> list[int]:
        """Sublattice id of each site of magnetic_structure."""
        return self.magnetic_structure.site_properties["sublattice_id"]

    def coupled_bonds(self, structure: Structure):
        """Yield (i, j, image, interaction) for every coupled bond between sites of structure,
        once from each end.

        Args:
            structure (Structure): Sites at parent positions that carry 'parent_index', in any
                supercell of the parent.
        """
        frac_coords = self.structure.frac_coords
        for i, j, image, _length in zip(
            *structure.get_neighbor_list(self._max_bond_length + DIST_EPSILON), strict=True
        ):
            start, end = structure[i].properties["parent_index"], structure[j].properties["parent_index"]
            vector = structure.lattice.get_cartesian_coords(structure[j].frac_coords + image - structure[i].frac_coords)
            # The bond vector, not the site positions, so a shift of the whole structure drops out.
            cells = self.structure.lattice.get_fractional_coords(vector) - (frac_coords[end] - frac_coords[start])
            interaction = self._interaction_of.get(Bond(start, end, _as_image(np.round(cells))))
            if interaction is not None:
                yield int(i), int(j), _as_image(image), interaction

    def _sublattices(self, symmetrized) -> tuple[list[int], dict[int, str], dict[int, str]]:
        """The sublattices: the symmetry orbits of the magnetic sites of symmetrized.

        Args:
            symmetrized (SymmetrizedStructure): The parent cell with its symmetry orbits.

        Returns:
            sublattice_ids (list[int]): Sublattice id of each site, _NONMAGNETIC_SUBLATTICE_ID
                for nonmagnetic ones.
            sublattice_names (dict[int, str]): Name of each sublattice id, numbered per element
                as in a CIF file: 'Fe1', 'Fe2', 'Mn1'.
            sublattice_wyckoff_symbols (dict[int, str]): Wyckoff symbol of each sublattice id.
        """
        sublattice_ids = [_NONMAGNETIC_SUBLATTICE_ID] * len(symmetrized)
        sublattice_names: dict[int, str] = {}
        sublattice_wyckoff_symbols: dict[int, str] = {}
        magnetic_orbits = [
            (indices, wyckoff_symbol)
            for indices, wyckoff_symbol in zip(symmetrized.equivalent_indices, symmetrized.wyckoff_symbols, strict=True)
            if symmetrized[indices[0]].specie.symbol in self.magnetic_species
        ]
        n_sublattices_of = Counter()
        for sublattice_id, (indices, wyckoff_symbol) in enumerate(magnetic_orbits):
            element = symmetrized[indices[0]].specie.symbol
            n_sublattices_of[element] += 1
            sublattice_names[sublattice_id] = f"{element}{n_sublattices_of[element]}"
            sublattice_wyckoff_symbols[sublattice_id] = wyckoff_symbol
            for index in indices:
                sublattice_ids[index] = sublattice_id
        return sublattice_ids, sublattice_names, sublattice_wyckoff_symbols

    def _interactions(self) -> list[Interaction]:
        """One Interaction per symmetry orbit of the coupled bonds.

        For each magnetic site and neighbor sublattice, every bond within cutoff couples, or, if
        there is none, the shortest bonds. Searched up to cutoff, but at least MAX_SEARCH_DIST.
        Each orbit holds every image of its bonds under the parent's symmetry operations, in both
        directions. Orbits are labelled 'J1(Fe1-Mn1)', 'J2(Fe1-Fe1)', ... by length; orbits of
        equal length are ordered by sublattice names, then by number of bonds, then by smallest
        bond.
        """
        sublattice_of = self.structure.site_properties["sublattice_id"]
        frac_coords = self.structure.frac_coords

        # Select the coupled bonds.
        candidates = defaultdict(list)  # (site, neighbor sublattice) -> [(length, bond)]
        for start, end, image, length in zip(
            *self.structure.get_neighbor_list(max(self.cutoff, MAX_SEARCH_DIST)), strict=True
        ):
            if _NONMAGNETIC_SUBLATTICE_ID not in (sublattice_of[start], sublattice_of[end]):
                bond = Bond(int(start), int(end), _as_image(image))
                candidates[start, sublattice_of[end]].append((float(length), bond))

        selected = []
        for length_and_bond in candidates.values():
            reach = max(self.cutoff, min(length for length, _bond in length_and_bond))
            selected += [bond for length, bond in length_and_bond if length <= reach + DIST_EPSILON]

        # Group them into symmetry orbits.
        orbits: list[frozenset[Bond]] = []
        for bond in selected:
            if any(bond in orbit for orbit in orbits):
                continue
            orbit = set()
            for op in self._symmetry_ops:
                start, start_cell = self._locate(op.operate(frac_coords[bond.start]))
                end, end_cell = self._locate(op.operate(frac_coords[bond.end] + bond.image))
                image = Bond(start, end, _as_image(end_cell - start_cell))
                orbit |= {image, image.reversed()}
            orbits.append(frozenset(orbit))

        # Label them by length.
        def sublattices(orbit):
            bond = next(iter(orbit))
            return tuple(sorted((sublattice_of[bond.start], sublattice_of[bond.end])))

        def names(orbit):
            return sorted(self.sublattice_names[sub_id] for sub_id in sublattices(orbit))

        def mean_length(orbit):
            return float(np.mean([self._length(bond) for bond in orbit]))

        # Lengths equal up to float noise count as equal, so their order does not hinge on it.
        orbits.sort(key=lambda orbit: (round(mean_length(orbit) / DIST_EPSILON), names(orbit), len(orbit), min(orbit)))
        return [
            Interaction(f"J{number}({'-'.join(names(orbit))})", sublattices(orbit), mean_length(orbit), orbit)
            for number, orbit in enumerate(orbits, start=1)
        ]

    def _locate(self, frac_coords) -> tuple[int, np.ndarray]:
        """The parent site at frac_coords, within symprec, and the parent cell it lies in."""
        offsets = frac_coords - self.structure.frac_coords
        cells = np.round(offsets)
        misfits = np.linalg.norm(self.structure.lattice.get_cartesian_coords(offsets - cells), axis=1)
        index = int(np.argmin(misfits))
        if misfits[index] > self.symprec:
            raise ValueError(
                f"A symmetry operation of the parent maps a site {misfits[index]:.3f} Angstrom away from "
                f"any site, beyond symprec={self.symprec}. Pass a symmetrized parent or raise symprec."
            )
        return index, cells[index].astype(int)

    def _length(self, bond: Bond) -> float:
        """Length (Angstrom) of bond in the parent geometry."""
        frac_coords = self.structure.frac_coords
        return float(
            np.linalg.norm(
                self.structure.lattice.get_cartesian_coords(
                    frac_coords[bond.end] + bond.image - frac_coords[bond.start]
                )
            )
        )


class ComputedOrdering:
    """One collinear magnetic ordering with its computed energy, mapped onto the parent cell."""

    def __init__(
        self,
        structure: Structure,
        energy: float,
        parent: ParentCell,
        matcher: StructureMatcher,
        default_magmoms: dict | None = None,
    ):
        """
        Args:
            structure (Structure): The ordering, with moments, in any supercell of the parent.
            energy (float): Its total energy (eV).
            parent (ParentCell): The parent cell it is mapped onto.
            matcher (StructureMatcher): Maps the parent onto the ordering. It must have
                primitive_cell=False and attempt_supercell=True.
            default_magmoms (dict): Elements that can carry a moment, as in
                CollinearMagneticStructureAnalyzer. Defaults to None, which uses
                default_magmoms.yaml.
        """
        settings = matcher.as_dict()
        if settings["primitive_cell"] or not settings["attempt_supercell"]:
            raise ValueError(
                "The `matcher` must be a StructureMatcher with primitive_cell=False and "
                "attempt_supercell=True, so it can match each ordering as a supercell of the "
                "parent cell."
            )

        # The analyzer works on a copy, so the caller's structure is never mutated. On this
        # copy every site carries a float 'magmom' site property, whether the input supplied
        # moments as site properties or as species spins, and induced moments on nonmagnetic
        # ions are zeroed out.
        self.structure = _analyzer(structure, default_magmoms=default_magmoms).structure
        self.energy = energy
        self.parent = parent
        self.ideal_magnetic_structure = self._magnetic_sites_on_parent(matcher)

    @cached_property
    def magnetic_structure(self) -> Structure:
        """Sites of the magnetic species, including any whose moment relaxed to zero."""
        return Structure.from_sites(_magnetic_sites(self.structure, self.parent.magnetic_species))

    @property
    def energy_per_magnetic_ion(self) -> float:
        """Total energy divided by the number of magnetic ions (eV)."""
        return self.energy / len(self.magnetic_structure)

    @property
    def sublattice_ids(self) -> list[int]:
        """Parent sublattice of each site of magnetic_structure."""
        return self.ideal_magnetic_structure.site_properties["sublattice_id"]

    def coupled_bonds(self):
        """Yield (i, j, image, interaction) for every coupled bond between sites of
        ideal_magnetic_structure, once from each end.
        """
        return self.parent.coupled_bonds(self.ideal_magnetic_structure)

    def spin_products(self, signs_only: bool = False) -> dict[str, float]:
        """-sum m_i m_j over the bonds of each interaction, per magnetic ion, keyed by label.

        Args:
            signs_only (bool): Use the sign of each moment instead of its value.
        """
        moments = np.array(self.ideal_magnetic_structure.site_properties["magmom"])
        if signs_only:
            moments = np.sign(moments)
        products = dict.fromkeys((interaction.label for interaction in self.parent.interactions), 0.0)
        for i, j, _image, interaction in self.coupled_bonds():
            products[interaction.label] -= moments[i] * moments[j]
        # Each bond is met once from each end.
        return {label: total / (2 * len(moments)) for label, total in products.items()}

    @cached_property
    def coupling_graph(self) -> StructureGraph:
        """Graph of the coupled bonds, on ideal_magnetic_structure."""
        graph = StructureGraph.from_empty_graph(self.ideal_magnetic_structure)
        for i, j, image, _interaction in self.coupled_bonds():
            graph.add_edge(i, j, to_jimage=image, warn_duplicates=False)
        return graph

    def _magnetic_sites_on_parent(self, matcher: StructureMatcher) -> Structure:
        """Magnetic sites at the parent's positions, with this ordering's moments and the
        parent's 'sublattice_id' and 'parent_index'.
        """
        # matched_parent[i] is the parent site, with its labels, that site i sits on.
        matched_parent = matcher.get_s2_like_s1(_nonmagnetic(self.structure), self.parent.structure)
        if matched_parent is None:
            raise ValueError(
                "This ordering is not a supercell of the parent cell; it cannot be mapped "
                "onto the parent sublattices. Pass an explicit `parent` cell that all "
                "orderings share, or, if relaxation distorted the ordering, a `matcher` "
                "with looser tolerances."
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
                "all orderings share, or, if you loosened the `matcher` tolerances, tighten "
                "them."
            )

        ideal = Structure.from_sites(_magnetic_sites(matched_parent, self.parent.magnetic_species))
        ideal.add_site_property("magmom", self.magnetic_structure.site_properties["magmom"])
        return ideal


class HeisenbergMapper:
    """Compute exchange parameters from the energies of collinear magnetic orderings.

    The model is defined once, on a paramagnetic *parent cell*: its symmetry orbits of
    sites are the magnetic sublattices, and its symmetry orbits of bonds are the
    interactions, one exchange constant J each. Every ordering is mapped onto the parent,
    so orderings in different supercells share one set of exchange constants, and each
    ordering contributes one row of spin products to a least-squares fit.

    Compute all orderings as static (SCF) calculations in the parent's geometry. Relaxing
    each ordering on its own mixes lattice energy into the energy differences, and the fit
    puts it into the J's.

    Attributes:
        parent (ParentCell): Defines the sublattices and the interactions.
        orderings (list[ComputedOrdering]): The screened orderings, sorted by energy per
            magnetic ion.
        ex_mat (DataFrame): Heisenberg Hamiltonian (per magnetic ion) for each ordering.
        ex_params (dict): Exchange parameter values. The J_ij are in meV/muB^2 (they
            multiply the raw moments, see get_exchange); the included 'E0' offset is in
            eV per magnetic ion.
    """

    def __init__(
        self,
        ordered_structures,
        energies,
        parent=None,
        cutoff=0,
        tol: float | None = None,
        symprec: float = DEFAULT_SYMPREC,
        matcher: StructureMatcher | None = None,
        default_magmoms: dict | None = None,
    ):
        """Map collinear orderings, possibly in different supercells of a common parent
        cell, onto a classical Heisenberg model.

        Args:
            ordered_structures (list): Magnetically ordered structures, with magmoms.
            energies (list): Total energy (eV) of each ordering.
            parent (Structure): Paramagnetic parent cell whose symmetry defines the
                magnetic sublattices and the interactions; it is reduced to its primitive
                cell. If None, the lowest-energy ordering is used, which is only correct if
                that ordering kept the parent's symmetry. Defaults to None.
            cutoff (float): Every bond within this distance (Angstrom) couples. A sublattice
                pair without such a bond still couples through its shortest bonds, so the
                default of 0 keeps only the shortest bonds of each pair.
            tol (float): Deprecated and ignored. Symmetry decides which bonds share a J, so
                there are no shells of near-equal length any more.
            symprec (float): Symmetry tolerance (Angstrom) for finding the parent's
                sublattices and interactions. Raise it, or pass a symmetrized parent, if
                sublattices that should be equivalent come out split. Defaults to 0.01.
            matcher (StructureMatcher): Maps each ordering onto the parent cell. If
                relaxation distorted an ordering beyond the defaults, loosen its ltol, stol
                or angle_tol as far as needed and keep its other settings. It must have
                primitive_cell=False and attempt_supercell=True, so each ordering is matched
                as a supercell of the parent. Defaults to
                StructureMatcher(primitive_cell=False, attempt_supercell=True).
            default_magmoms (dict): Its keys are the elements that can carry a moment, as
                in CollinearMagneticStructureAnalyzer; oxidation states are ignored, so key
                by element ("Fe", not "Fe2+"). Moments on any other element are zeroed, and
                those sites drop out of the fit. The values are only used as moments for a
                structure that carries no moments at all. It replaces default_magmoms.yaml
                rather than extending it, so list every magnetic element. Defaults to None,
                which uses default_magmoms.yaml.
        """
        if parent is not None and not isinstance(parent, Structure):
            raise TypeError(
                f"parent must be a Structure or None, got {type(parent).__name__}. "
                "Note the constructor signature changed: parent now comes third, "
                "before cutoff - see the module docstring's migration guide."
            )
        if tol is not None:
            warnings.warn(
                "`tol` is deprecated and ignored: symmetry decides which bonds share a J.",
                DeprecationWarning,
                stacklevel=2,
            )

        if matcher is None:
            matcher = StructureMatcher(primitive_cell=False, attempt_supercell=True)

        # Save original copies of inputs
        self.ordered_structures_ = ordered_structures
        self.energies_ = energies
        self.parent_ = parent

        self.cutoff = cutoff
        self.symprec = symprec
        self.matcher = matcher
        self.default_magmoms = default_magmoms
        # Set by get_exchange.
        self.ex_params = None
        self.residual = None

        # The magnetic species: the elements that carry a moment in any ordering.
        magnetic_species = self._magnetic_species(ordered_structures, default_magmoms)

        # The model: sublattices and interactions, defined on the parent cell.
        if parent is None:
            parent = self._inferred_parent(ordered_structures, energies, magnetic_species)
        self.parent = ParentCell(parent, magnetic_species, cutoff=cutoff, symprec=symprec)

        # The data: every ordering mapped onto the parent, deduplicated and sorted by energy per magnetic ion.
        orderings = [
            ComputedOrdering(struct, energy, self.parent, matcher, default_magmoms=default_magmoms)
            for struct, energy in zip(ordered_structures, energies, strict=True)
        ]
        self.orderings = HeisenbergScreener(orderings, screen=False).screened_orderings
        if len(self.orderings) < 2:
            raise ValueError("HeisenbergMapper needs at least 2 unique orderings.")

        # The fit's linear system: one row per ordering, one column per interaction it can resolve.
        self.ex_mat = self._exchange_matrix()

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

    @property
    def sublattice_names(self):
        """dict[int, str]: Maps each sublattice id to its name, e.g. 'Fe1'."""
        return self.parent.sublattice_names

    @property
    def interactions(self):
        """dict[tuple[int, int], list[str]]: J labels of each sublattice pair (i, j), i <= j,
        nearest first. E.g. {(0, 1): ['J1(Fe1-Mn1)', 'J4(Fe1-Mn1)']}.
        """
        labels = defaultdict(list)
        for interaction in self.parent.interactions:
            labels[interaction.sublattices].append(interaction.label)
        return dict(labels)

    @property
    def dists(self):
        """dict[str, float]: Bond length (Angstrom) of each J label."""
        return {interaction.label: interaction.length for interaction in self.parent.interactions}

    @staticmethod
    def _magnetic_species(ordered_structures, default_magmoms):
        """Elements of default_magmoms with a moment in any ordering. Pooled over all orderings,
        so an ion whose moment relaxed to zero in one ordering stays on its sublattice there.

        Warns about moments above _ZEROED_MOMENT_WARNING_THRESHOLD on any other species, since
        the fit zeroes them.
        """
        magnetic_species = set()
        largest_zeroed_moment: dict[str, float] = {}
        for struct in ordered_structures:
            # threshold_nonmag=0.0 keeps the moments on other species, so they can be reported.
            analyzer = _analyzer(struct, default_magmoms=default_magmoms, threshold_nonmag=0.0)
            for site in analyzer.structure:
                species = site.species_string
                moment = abs(site.properties["magmom"])
                if species in analyzer.default_magmoms:
                    if moment:
                        magnetic_species.add(site.specie.symbol)
                elif moment > _ZEROED_MOMENT_WARNING_THRESHOLD:
                    largest_zeroed_moment[species] = max(largest_zeroed_moment.get(species, 0.0), moment)

        if largest_zeroed_moment:
            listed = ", ".join(
                f"{species} (up to {moment:.2f} muB)" for species, moment in sorted(largest_zeroed_moment.items())
            )
            warnings.warn(
                f"Moments on {listed} are zeroed because these species are not in `default_magmoms`, "
                "so their sites drop out of the fit. Pass a `default_magmoms` that includes them if "
                "they are magnetic.",
                UserWarning,
                # _magnetic_species <- __init__ <- caller
                stacklevel=3,
            )
        return magnetic_species

    @staticmethod
    def _inferred_parent(ordered_structures, energies, magnetic_species):
        """The ordering with the lowest energy per magnetic ion, as a stand-in parent, with a UserWarning."""
        warnings.warn(
            "No `parent` cell supplied; the magnetic sublattices are inferred from the "
            "primitive cell of the lowest-energy ordering. This is only correct if that "
            "ordering still has the symmetry of the paramagnetic parent - relaxation "
            "usually lowers it, which silently splits or merges sublattices. Pass an "
            "explicit `parent` unless you have checked that it does not.",
            UserWarning,
            # _inferred_parent <- __init__ <- caller
            stacklevel=3,
        )

        energies_per_magnetic_ion = [
            energy / len(_magnetic_sites(struct, magnetic_species))
            for struct, energy in zip(ordered_structures, energies, strict=True)
        ]
        return ordered_structures[int(np.argmin(energies_per_magnetic_ion))]

    def _exchange_matrix(self):
        """One row per ordering: its energy per magnetic ion 'E', 'E0' = 1 and its spin products.

        Drops each J whose column depends linearly on E0 and the shorter J's, with a UserWarning.
        """
        # Rows: per ordering, the energy to fit and the spin products each J multiplies.
        rows = [
            {"E": ordering.energy_per_magnetic_ion, "E0": 1.0, **ordering.spin_products()}
            for ordering in self.orderings
        ]
        ex_mat = pd.DataFrame(rows)

        # Columns: shortest first, each J the orderings can separate from E0 and the shorter J's.
        # Signs, not moments: small moment differences between orderings would hide a dependence.
        sign_mat = pd.DataFrame([{"E0": 1.0, **ordering.spin_products(signs_only=True)} for ordering in self.orderings])
        kept = ["E0"]
        dropped = []
        for label in sorted(self.dists, key=self.dists.get):
            if np.linalg.matrix_rank(sign_mat[[*kept, label]]) > len(kept):
                kept.append(label)
            else:
                dropped.append(label)

        if dropped:
            lower_cutoff = " or lower the cutoff" if self.cutoff else ""
            warnings.warn(
                f"The orderings cannot separate {dropped} from E0 and the shorter interactions; "
                "they are left out of the fit and the interaction graph. Add orderings that tell "
                f"them apart, possibly in larger supercells{lower_cutoff}.",
                UserWarning,
                # _exchange_matrix <- __init__ <- caller
                stacklevel=3,
            )

        return ex_mat[["E", *kept]]

    def get_exchange(self):
        """Fit E0 and the J_ij to the energies of all orderings by least squares.

        The J_ij multiply the raw moments, so they are in meV/muB^2 (bond energy
        J_ij * m_i * m_j); this keeps moment magnitudes that differ between orderings out
        of the fitted J_ij.

        The residual is only meaningful with substantially more orderings than fitted
        parameters; otherwise the extra parameters overfit, lowering the residual without
        improving the fit.

        Returns:
            ex_params (dict[str, float]): J_ij in meV/muB^2 and 'E0' in eV per magnetic ion.
            residual (float): RMS fit residual in meV per magnetic ion.

        Raises:
            ValueError: If no exchange interaction is left to fit.
        """
        col_names = [c for c in self.ex_mat.columns if c != "E"]
        if col_names == ["E0"]:
            raise ValueError(
                "Exchange matrix holds no interaction besides E0: the orderings cannot separate "
                "any interaction from E0, so the energies constrain no J."
            )

        design_matrix = self.ex_mat[col_names].to_numpy(dtype=float)
        energies = self.ex_mat["E"].to_numpy(dtype=float)

        # Column-normalized, so the threshold does not depend on the moment magnitudes.
        condition_number = np.linalg.cond(design_matrix / np.linalg.norm(design_matrix, axis=0))
        if condition_number > _MAX_CONDITION_NUMBER:
            warnings.warn(
                f"Exchange matrix is ill-conditioned (cond={condition_number:.1e}); the fitted exchange "
                "parameters are unreliable. The input orderings are near-degenerate or the "
                "model has more parameters than the orderings can constrain. Supply more, "
                "more-distinct orderings.",
                UserWarning,
                stacklevel=2,
            )

        params = np.linalg.lstsq(design_matrix, energies, rcond=None)[0]

        # Averaged over the orderings, not the degrees of freedom, so a square system works.
        residual = _MEV_PER_EV * float(np.sqrt(np.mean((design_matrix @ params - energies) ** 2)))  # meV per ion
        ex_params = {
            name: value if name == "E0" else _MEV_PER_EV * value  # J_ij in meV/muB^2, E0 in eV per ion
            for name, value in zip(col_names, params.tolist(), strict=True)
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
                graph for. Defaults to 0 (the lowest-energy ordering).

        Returns:
            StructureGraph: Exchange interaction graph.
        """
        if self.ex_params is None:
            self.get_exchange()

        ordering = self.orderings[ordering_index]
        structure = ordering.ideal_magnetic_structure.copy()  # the returned graph must not share it
        magmoms = structure.site_properties["magmom"]
        igraph = StructureGraph.from_empty_graph(
            structure, edge_weight_name="exchange_constant", edge_weight_units="meV"
        )

        for i, j, image, interaction in ordering.coupled_bonds():
            left_out_of_fit = interaction.label not in self.ex_params
            if left_out_of_fit:
                continue
            exchange_constant = self.ex_params[interaction.label] * abs(magmoms[i] * magmoms[j])  # meV
            # add_edge stores a zero weight as no weight at all.
            if exchange_constant == 0:
                continue
            igraph.add_edge(i, j, to_jimage=image, weight=exchange_constant, warn_duplicates=False)

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
            symprec=self.symprec,
            coupling_graphs=self.coupling_graphs,
            sublattice_ids=self.sublattice_ids,
            sublattice_names=self.sublattice_names,
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
            ordering = _analyzer(s, default_magmoms=self.default_magmoms).ordering
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
            "<J> is a single average in meV/magnetic ion rather than a set of bond-resolved "
            "J_ij, and it is only defined for a pair of FM/AFM orderings."
        ),
        category=DeprecationWarning,
        deadline=(2027, 8, 1),
    )
    def estimate_exchange(self, fm_struct=None, afm_struct=None, fm_e=None, afm_e=None):
        """Estimate <J> for a structure based on low energy FM and AFM orderings.

        .. deprecated::
            Use :meth:`get_exchange` instead, which fits bond-resolved J_ij over all
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
        j_avg *= _MEV_PER_EV  # meV / ion

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
            sublattices_of = {interaction.label: interaction.sublattices for interaction in self.parent.interactions}
            for k, j_val in ex_params.items():
                i, j = sublattices_of[k]
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

    def __init__(self, orderings: list[ComputedOrdering], screen=False):
        """Pre-processes magnetic orderings for HeisenbergMapper.
        It prioritizes low-energy orderings with large and localized magnetic moments.

        Args:
            orderings (list[ComputedOrdering]): The orderings to screen. Each one already
                owns its magnetic-only substructure and its per-magnetic-ion energy.
            screen (bool): Try to screen out high energy and low-spin configurations.

        Attributes:
            screened_orderings (list[ComputedOrdering]): Deduplicated orderings, sorted by
                energy per magnetic ion.
        """
        orderings = self._do_cleanup(orderings)

        # If there are more than 2 orderings, we want to perform a
        # screening to prioritize well-behaved ones
        if screen and len(orderings) > 2:
            orderings = self._do_screen(orderings)

        self.screened_orderings = orderings

    @staticmethod
    def _do_cleanup(orderings: list[ComputedOrdering]):
        """Drop duplicate orderings and sort by energy per magnetic ion.

        Two orderings are duplicates if both their energies per magnetic ion and their spin
        products agree to _DUPLICATE_DECIMALS decimals, e.g. one state submitted in two
        supercells. The first is kept.

        Args:
            orderings (list[ComputedOrdering]): The orderings to clean up.

        Returns:
            list[ComputedOrdering]: Deduplicated, sorted by energy per magnetic ion.
        """
        unique: dict[tuple, ComputedOrdering] = {}
        for ordering in orderings:
            spin_products = (round(product, _DUPLICATE_DECIMALS) for product in ordering.spin_products().values())
            energy = round(ordering.energy_per_magnetic_ion, _DUPLICATE_DECIMALS)
            unique.setdefault((energy, *spin_products), ordering)
        return sorted(unique.values(), key=lambda ordering: ordering.energy_per_magnetic_ion)

    @staticmethod
    def _do_screen(orderings: list[ComputedOrdering]):
        """Screen and sort magnetic orderings based on some criteria.

        Prioritize low energy orderings and large, localized magmoms.

        Args:
            orderings (list[ComputedOrdering]): Cleaned up orderings, sorted by energy.

        Returns:
            list[ComputedOrdering]: The ground and first excited state, followed by the
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
        symprec=None,
        coupling_graphs=None,
        sublattice_ids=None,
        sublattice_names=None,
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
            energies (list): Total energy (eV) of each ordering.
            energies_per_magnetic_ion (list): Energy per magnetic ion (eV) of each ordering,
                the energies the exchange parameters were fitted to.
            cutoff (float): Every bond within this distance (Angstrom) couples.
            symprec (float): Symmetry tolerance (Angstrom) the parent's sublattices and
                interactions were found with.
            coupling_graphs (list): Coupling graph of each ordering, built on its magnetic
                sites at the parent's positions.
            sublattice_ids (list[list[int]]): sublattice_ids[k][i] is the sublattice id of
                site i in ordering k.
            sublattice_names (dict): Maps each sublattice id to its name, e.g. 'Fe1'.
            sublattice_wyckoff_symbols (dict): Maps each sublattice id to its wyckoff symbol.
            interactions (dict): {(i, j): [J label, ...]} - the interactions of each
                sublattice pair, nearest first.
            dists (dict): {J label: bond length in Angstrom}.
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
        self.symprec = symprec
        self.coupling_graphs = coupling_graphs
        self.sublattice_ids = sublattice_ids
        self.sublattice_names = sublattice_names
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
            "sublattice_names": jsanitize(self.sublattice_names),
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
        sublattice_names = {literal_eval(k): v for k, v in dct["sublattice_names"].items()}
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
            symprec=dct["symprec"],
            coupling_graphs=coupling_graphs,
            sublattice_ids=dct["sublattice_ids"],
            sublattice_names=sublattice_names,
            sublattice_wyckoff_symbols=sublattice_wyckoff_symbols,
            interactions=interactions,
            dists=dct["dists"],
            ex_mat=ex_mat,
            ex_params=dct["ex_params"],
            residual=dct["residual"],
            igraph=igraph,
        )
