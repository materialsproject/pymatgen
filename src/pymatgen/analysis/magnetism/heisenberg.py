"""
This module implements a simple algorithm for extracting nearest neighbor
exchange parameters by mapping low energy magnetic orderings to a Heisenberg
model.
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from ast import literal_eval
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
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


def _analyzer(structure: Structure, **kwargs) -> CollinearMagneticStructureAnalyzer:
    """CollinearMagneticStructureAnalyzer factory with the settings used throughout this module.

    ``make_primitive=False`` keeps the cell as given, ``threshold=0.0`` retains any nonzero
    moment, and ``threshold_nonmag=100.0`` always zeroes out nonmagnetic ions, so induced
    magnetic moments don't change the number of magnetic sites per unit cell between orderings.

    Extra keyword arguments are passed through to the analyzer.
    """
    return CollinearMagneticStructureAnalyzer(
        structure, **{"make_primitive": False, "threshold": 0.0, "threshold_nonmag": 100.0, **kwargs} # last value wins if there are multiple conflicting values for the same key -> defaults overridden by kwargs
    )



class MagneticOrdering(ABC):
    """A single collinear magnetic configuration on a common parent lattice.

    Base class holding one ordering's structure, its magnetic-only reduction,
    its NN graph, and the parent-sublattice labels of its magnetic sites.
    """

    def __init__(self, structure: Structure, cutoff: float = 0, tol: float = DEFAULT_TOL):
        self.structure = structure
        self.analyzer = _analyzer(structure)
        self.cutoff = cutoff
        self.tol = tol
        self.nn_graph = self._get_nn_graph()
        self.magn_species: set[str] = {site.specie.symbol for site in self.magn_substructure}

        self.sublattice_ids_dict: dict[tuple[int, ...], int] = {}
        self.sublattice_labels: list[int] | None = None
        self.sublattice_wyckoff_symbols: dict[int, str] = {}

    @property
    def magn_substructure(self) -> Structure:
        return self._magnetic_only(self.structure)

    @staticmethod
    def _nonmagnetic(structure) -> Structure:
        """Nonmagnetic copy of a structure (moments zeroed, all ions kept)."""
        s0 = _analyzer(structure).get_nonmagnetic_structure(make_primitive=False)
        if "wyckoff" in s0.site_properties:
            s0.remove_site_property("wyckoff")
        return s0

    @staticmethod
    def _magnetic_only(structure) -> Structure:
        """Magnetic-only copy (nonmagnetic ions removed).

        Whether a site is magnetic is decided by its *species*: only species listed in
        CollinearMagneticStructureAnalyzer's default magmoms are treated as magnetic.

        - Magnetic-species sites keep their existing magmom (any nonzero value is retained,
          since ``threshold=0.0``). If the input has no magmom site property, the default magmoms are assigned (``replace_all_if_undefined``).
        - Nonmagnetic-species sites always get a magmom of zero. Because ``threshold_nonmag=100.0``,
          even induced moments on nonmagnetic ions are zeroed out, so the number of magnetic sites
          per cell stays constant across orderings.

        Finally, all sites with zero magmom are dropped (even if magnetic-species sites with magmom=0. or magmom=None) and only sites with nonzero magmom are kept.
        """
        # overwrite_magmom_mode: if the structure has no magmom site property, assign default magmoms to all magnetic species
        return _analyzer(structure).get_structure_with_only_magnetic_atoms(make_primitive=False)

    @staticmethod
    def _sublattice_ids_dict(labels: list[int]) -> dict[tuple[int, ...], int]:
        """Group site indices by sublattice id: {tuple(site indices): sublattice id}."""
        groups: dict[int, list[int]] = {}
        for idx, sub_id in enumerate(labels):
            groups.setdefault(sub_id, []).append(idx)
        return {tuple(indices): sub_id for sub_id, indices in groups.items()}
    
    def _get_nn_graph(self) -> StructureGraph:
        """
        NN StructureGraph of the magnetic-only structure.
        
        If self.cutoff is set, the graph includes all neighbors within that distance. 
        Otherwise, only the nearest neighbors of each site are included.
        """

        strategy = MinimumDistanceNN(cutoff=self.cutoff, get_all_sites=True) if self.cutoff else MinimumDistanceNN()  # only NN
        return StructureGraph.from_local_env_strategy(self.magn_substructure, strategy=strategy)
    
    @abstractmethod
    def set_sublattice_ids(self) -> dict[tuple[int, ...], int]:
        """
        Sets the following properties for this ordering:
            - `self.sublattice_ids_dict` = {tuple(site indices): sublattice id}
            - `self.sublattice_wyckoff_symbols` = {sublattice id: wyckoff symbol}
            - `self.sublattice_ids` = [sublattice id for each magnetic site in the ordering or None for nonmagnetic sites]

        All nonmagnetic (`CollinearMagneticStructureAnalyzer` default_magmom=0) sites are grouped into one sublattice with label `None`.
        
        For `ParentOrdering`, this defines the sublattice ids by grouping symmetrically equivalent sites of ions into sublattices (nonmagnetic ions present for symmetry analysis). 
        
        For `RelaxedOrdering`, this is a mapping from the ordering's site indices to the sublattice ids via the finding the corresponding site in the parent for each site in the relaxed ordering.
        
        Finally, the sublattice ids are stored as a site property in the ordering's structure for easy access.
        """


class ParentOrdering(MagneticOrdering):
    """The nonmagnetic parent cell that defines the sublattices.

    Sublattices are its symmetry orbits (computed with nonmagnetic ions present),
    restricted to the magnetic species.

    """

    def __init__(self, structure: Structure, cutoff: float = 0, tol: float = DEFAULT_TOL):
        super().__init__(structure, cutoff, tol)
        self.structure = self._nonmagnetic(structure)  # parent is always nonmagnetic
        self.set_sublattice_ids()

    def set_sublattice_ids(self):
        symmetrized_parent = SpacegroupAnalyzer(self.structure).get_symmetrized_structure()
        self.sublattice_ids = [None] * len(self.structure)
        for indices, wyckoff_symbol in zip(symmetrized_parent.equivalent_indices, symmetrized_parent.wyckoff_symbols, strict=True):
            assert all(symmetrized_parent[i].specie == symmetrized_parent[indices[0]].specie for i in indices)
            if symmetrized_parent[indices[0]].specie.symbol not in self.magn_species:
                continue
            sub_id = len(self.sublattice_wyckoff_symbols)
            self.sublattice_wyckoff_symbols[sub_id] = wyckoff_symbol # one sublattice always has one wyckoff symbol, but multiple sublattices can share the same wyckoff symbol (e.g. two species on the same wyckoff position)
            for index in indices:
                self.sublattice_ids[index] = sub_id

        self.structure.add_site_property("sublattice_id", self.sublattice_ids)

        self.sublattice_ids_dict = self._sublattice_ids_dict(self.sublattice_ids)


class RelaxedOrdering(MagneticOrdering):
    """One DFT-relaxed magnetic ordering with its energy and parent back-reference."""

    def __init__(self, structure: Structure, energy: float, parent_ordering: ParentOrdering, cutoff: float = 0, tol: float = DEFAULT_TOL):
        super().__init__(structure, cutoff, tol)
        self.parent_ordering = parent_ordering
        self.energy = energy
        self.set_sublattice_ids()

    def set_sublattice_ids(self):
        matcher = StructureMatcher(primitive_cell=False, attempt_supercell=True)

        matched_parent = matcher.get_s2_like_s1(self._nonmagnetic(self.structure), self.parent_ordering.structure)
        if matched_parent is None:
            raise ValueError(
                "The RelaxedOrdering is not a supercell of the ParentOrdering cell; it "
                "cannot be mapped onto the parent sublattices."
            )
        self.sublattice_ids = [int(matched_parent[i].properties["sublattice_id"]) for i in range(len(self.structure))]




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
        sgraphs (list): StructureGraph objects, one per ordering.
        parent (Structure): Nonmagnetic parent cell that defines the sublattices.
        site_labels (list[list[int]]): site_labels[k][i] is the parent
            sublattice id of site i in ordering k.
        sublattice_wyckoff_symbols (dict): Maps each sublattice id to its wyckoff symbol.
        nn_interactions (dict): {i: j} pairs of NN interactions between sublattices.
        dists (dict): NN, NNN, and NNNN interaction distances.
        ex_mat (DataFrame): Heisenberg Hamiltonian (per magnetic ion) for each ordering.
        ex_params (dict): Exchange parameter values. The J_ij are in meV; the
            included 'E0' offset stays in eV.
    """

    def __init__(self, ordered_structures, energies, parent, cutoff=0, tol: float = DEFAULT_TOL):
        """Exchange parameters are computed by mapping to a classical Heisenberg
        model. n+1 unique orderings are required to compute n exchange parameters.

        First run a MagneticOrderingsWF to obtain low energy collinear magnetic
        orderings and find the magnetic ground state, enumerate magnetic
        states with the ground state as the input structure and do static
        calculations for these orderings. The orderings may live in different
        supercells - they only need to be commensurate with a common parent cell.

        Args:
            ordered_structures (list): Structure objects with magmoms.
            energies (list): Total energies of each relaxed magnetic structure.
            parent (Structure): Paramagnetic parent cell whose symmetry defines the
                magnetic sublattices. 
            cutoff (float): Cutoff in Angstrom for nearest neighbor search.
                Defaults to 0 (only NN, no NNN, etc.).
            tol (float): Tolerance (in Angstrom) on nearest neighbor distances
                being equal and structure matching.
        """
        # Save original copies of inputs
        self.ordered_structures_ = ordered_structures
        self.energies_ = energies
        self.parent_ = parent

        # Sanitize inputs: standardize moments, convert energies to per magnetic
        # ion, drop duplicates and sort by energy. The screener keeps both a
        # magnetic-only copy (for the graphs) and a full copy (all ions, needed
        # to determine site symmetry) of every ordering.
        hs = HeisenbergScreener(ordered_structures, energies, screen=False)
        self.ordered_structures = hs.screened_structures
        self.full_structures = hs.screened_full_structures
        self.energies = hs.screened_energies
        self.cutoff = cutoff
        self.tol = tol

        # Graph representation of each (magnetic-only) ordering
        self.sgraphs = [self._get_graph(cutoff, s) for s in self.ordered_structures]

        # The parent cell defines the sublattices; label every magnetic site in
        # every ordering with the parent sublattice it belongs to.
        self.parent = self._get_parent(parent, self.full_structures)
        self.parent_sgraph = self._get_graph(cutoff, self._magnetic_only(self.parent))
        self.sublattice_wyckoff_symbols, self.site_labels, self.parent_site_labels = self._label_orderings(self.parent, self.full_structures)

        # These attributes are set by internal methods
        self.nn_interactions = self.dists = self.ex_mat = self.ex_params = None

        if len(self.sgraphs) < 2:
            raise SystemExit("We need at least 2 unique orderings.")

        self._get_nn_dict()
        self._build_exchange_matmat()

    @staticmethod
    def _get_graph(cutoff, structure):
        """Generate graph representations of a magnetic structure with nearest
        neighbor bonds. Right now this only works for MinimumDistanceNN.

        Args:
            cutoff (float): Cutoff in Angstrom for nearest neighbor search.
            structure (Structure): A single magnetic structure.

        Returns:
            StructureGraph: The graph representation of the structure.
        """
        strategy = MinimumDistanceNN(cutoff=cutoff, get_all_sites=True) if cutoff else MinimumDistanceNN()  # only NN
        return StructureGraph.from_local_env_strategy(structure, strategy=strategy)


    # @classmethod
    # def _get_parent(cls, parent, full_structures):
    #     """Return the full primitive parent cell that defines the sublattices.

    #     The nonmagnetic ions are kept: site equivalence is read from this parent's
    #     symmetry, and removing the nonmagnetic ions first can raise the apparent
    #     site symmetry and merge sublattices that are actually distinct.

    #     Args:
    #         parent (Structure | None): Explicit paramagnetic parent (all ions). If
    #             None, the primitive cell of the lowest-energy ordering is used.
    #         full_structures (list[Structure]): Orderings with all ions retained.

    #     Returns:
    #         Structure: Nonmagnetic primitive parent cell, nonmagnetic ions included.
    #     """
    #     ref = parent if parent is not None else full_structures[0]
    #     return cls._nonmagnetic(ref).get_primitive_structure()

    # @classmethod
    # def _label_orderings(cls, parent, full_structures):
    #     """Label every magnetic site in every ordering with its parent sublattice.

    #     The magnetic sublattices are the parent's symmetry orbits, computed on the
    #     *full* cell (nonmagnetic ions present) and then restricted to the magnetic
    #     species. Each ordering is folded back onto the full parent so its sites
    #     inherit the parent labels; the labels of the nonmagnetic sites are
    #     discarded. This keeps the labels consistent across orderings that live in
    #     different supercells, while preserving the true site symmetry.

    #     Args:
    #         parent (Structure): Nonmagnetic parent cell (all ions).
    #         full_structures (list[Structure]): Orderings with all ions retained and
    #             a 'magmom' site property on every site.

    #     Returns:
    #         tuple[dict, list[list[int]], list[int]]:
    #             - sublattice_wyckoff_symbols maps each magnetic sublattice id to its wyckoff symbol.
    #             - site_labels[k][i] is the sublattice id of magnetic site i in
    #               ordering k (aligned with the magnetic-only graphs).
    #             - parent_site_labels is the sublattice id of each site in the parent structure.

    #     Raises:
    #         ValueError: If an ordering is not a supercell of the parent cell.
    #     """
    #     # Species that carry a moment in any ordering are the magnetic species.
    #     mag_species = {site.specie.symbol for s in full_structures for site in s if abs(site.properties["magmom"]) > 0}

    #     # Sublattices = the parent's symmetry orbits restricted to magnetic
    #     # species. Nonmagnetic sites get the sentinel -1 (never used).
    #     symm_parent = SpacegroupAnalyzer(parent).get_symmetrized_structure()
    #     parent_labels = [-1] * len(parent)
    #     sublattice_wyckoff_symbols = {}
    #     for indices, symbol in zip(symm_parent.equivalent_indices, symm_parent.wyckoff_symbols, strict=True):
    #         if symm_parent[indices[0]].specie.symbol not in mag_species:
    #             continue
    #         sub_id = len(sublattice_wyckoff_symbols)
    #         sublattice_wyckoff_symbols[sub_id] = symbol
    #         for index in indices:
    #             parent_labels[index] = sub_id

    #     # Carry the parent labels across to each ordering by folding the ordering's
    #     # full geometry onto the tagged parent cell (the nonmagnetic ions help pin
    #     # down a unique mapping). This also validates that every ordering is a
    #     # genuine supercell of the parent.
    #     tagged_parent = parent.copy()
    #     tagged_parent.add_site_property("sublattice_id", parent_labels)
    #     matcher = StructureMatcher(primitive_cell=False, attempt_supercell=True)

    #     site_labels = []
    #     for idx, full in enumerate(full_structures):
    #         matched = matcher.get_s2_like_s1(cls._nonmagnetic(full), tagged_parent)
    #         if matched is None:
    #             raise ValueError(
    #                 f"Ordering {idx} is not a supercell of the parent cell; it "
    #                 "cannot be mapped onto the parent sublattices. Pass an explicit "
    #                 "`parent` cell that all orderings share."
    #             )
    #         # matched[i] is the parent site that ordering site i folds onto. Keep
    #         # only the magnetic sites, in order, so the labels align with the
    #         # magnetic-only structure used to build the graphs.
    #         magmoms = full.site_properties["magmom"]
    #         site_labels.append(
    #             [int(matched[i].properties["sublattice_id"]) for i in range(len(full)) if abs(magmoms[i]) > 0]
    #         )

    #     parent_labels = [int(label) for label in parent_labels if label >= 0]  # drop nonmagnetic sites

    #     return wyckoff_ids, site_labels, parent_labels


    def _shell_of(self, dist):
        """Return 'nn'/'nnn'/'nnnn' for a distance, or None if it matches no shell."""
        for shell in ("nn", "nnn", "nnnn"):
            if self.dists[shell] and abs(dist - self.dists[shell]) <= self.tol:
                return shell
        return None

    def _get_nn_dict(self):
        """Set self.dists and self.nn_interactions describing the unique NN, NNN and
        NNNN interactions between sublattices.

        Distances and connectivity are read from the parent ordering. See _get_parent() for how the parent is defined and how the sublattices are labelled. The connectivity is
    read from the parent StructureGraph, which is built from the magnetic-only parent structure (nonmagnetic ions are ignored for the graph, but they are kept in the parent structure to preserve the true site symmetry).
        """
        tol = self.tol
        sgraph = self.parent_sgraph
        labels = self.parent_site_labels

        # One representative site index per sublattice
        reps = {}
        for idx, sub_id in enumerate(labels):
            reps.setdefault(sub_id, idx)

        # Collect neighbor distances (up to NNNN) across the sublattices
        all_dists = []
        for i in reps.values():
            dists = sorted({round(cs[-1], 2) for cs in sgraph.get_connected_sites(i)})
            all_dists += dists[:3]

        # Collapse distances that are equal within tol, keep up to NNNN, pad with 0
        all_dists = sorted(set(all_dists))
        rm_list = [idx for idx, d in enumerate(all_dists[:-1], start=1) if abs(d - all_dists[idx]) < tol]
        all_dists = [d for idx, d in enumerate(all_dists) if idx not in rm_list]
        all_dists = [*all_dists, 0, 0, 0][:3]
        self.dists = dict(zip(("nn", "nnn", "nnnn"), all_dists, strict=True))

        # Determine which sublattice pairs interact at each shell
        nn_interactions = {"nn": {}, "nnn": {}, "nnnn": {}}
        for i_id, i in reps.items():
            for cs in sgraph.get_connected_sites(i):
                shell = self._shell_of(round(cs[-1], 2))
                if shell is not None:
                    nn_interactions[shell][i_id] = labels[cs[2]]

        self.nn_interactions = nn_interactions

    def _exchange_columns(self):
        """Column labels for the exchange matrix: 'E', 'E0' then one per J_ij."""
        columns = ["E", "E0"]
        for shell, pairs in self.nn_interactions.items():
            for i_id, j_id in pairs.items():
                col = f"{i_id}-{j_id}-{shell}"
                rev = f"{j_id}-{i_id}-{shell}"
                if col not in columns and rev not in columns:
                    columns.append(col)
        # Keep at most n interactions for n+1 orderings
        return columns[: len(self.sgraphs) + 1]

    def _build_exchange_matmat(self):
        """Build the Heisenberg Hamiltonian, one row per ordering, by summing the
        signed products S_i . S_j over each graph. Sets self.ex_mat.

        Each row is normalised per magnetic ion so that orderings living in
        different-sized supercells share a single linear system (the energies are
        already per magnetic ion, see HeisenbergScreener._do_cleanup).
        """
        columns = self._exchange_columns()
        j_columns = [c for c in columns if c not in ("E", "E0")]

        if len(j_columns) < 2:  # Only <J> can be calculated
            self.ex_mat = pd.DataFrame(columns=columns)
            return

        rows = []
        seen = set()  # de-duplicate identical interaction rows to avoid singularity
        for k, sgraph in enumerate(self.sgraphs):
            labels = self.site_labels[k]
            magmoms = sgraph.structure.site_properties["magmom"]
            n_sites = len(sgraph.structure)

            row = dict.fromkeys(columns, 0.0)
            for i in range(len(sgraph.graph.nodes)):
                s_i = magmoms[i]
                i_id = labels[i]
                for cs in sgraph.get_connected_sites(i):
                    shell = self._shell_of(round(cs[-1], 2))
                    if shell is None:
                        continue
                    j_id = labels[cs[2]]
                    col = f"{i_id}-{j_id}-{shell}"
                    rev = f"{j_id}-{i_id}-{shell}"
                    if col in row:
                        row[col] -= s_i * magmoms[cs[2]]
                    elif rev in row:
                        row[rev] -= s_i * magmoms[cs[2]]

            # Extensive sums -> per magnetic ion, with the 1/2 Heisenberg factor.
            # Each bond is visited from both ends, so 2 * n_sites also removes the
            # double counting.
            for c in j_columns:
                row[c] /= 2 * n_sites

            key = tuple(round(row[c], 10) for c in j_columns)
            if key in seen:
                continue
            seen.add(key)
            row["E0"] = 1.0  # nonmagnetic contribution (per ion)
            row["E"] = self.energies[k]
            rows.append(row)

        ex_mat = pd.DataFrame(rows, columns=columns)

        # Drop interaction columns that never appear (all zero) to keep H invertible
        j_columns = [c for c in j_columns if (ex_mat[c] != 0).any()]
        ex_mat = ex_mat[["E", "E0", *j_columns]]

        # Square system: one row per parameter (E0 + the J_ij)
        n_params = ex_mat.shape[1] - 1
        self.ex_mat = ex_mat.iloc[:n_params].reset_index(drop=True)

    def get_exchange(self):
        """
        Take Heisenberg Hamiltonian and corresponding energy for each row and
        solve for the exchange parameters.

        Returns:
            dict[str, float]: Exchange parameters. The J_ij are in meV; the
                included 'E0' offset stays in eV.
        """
        ex_mat = self.ex_mat
        E = ex_mat[["E"]]
        j_names = [j for j in ex_mat.columns if j != "E"]

        # Only 1 NN interaction -> estimate <J> from E_AFM - E_FM
        if len(j_names) < 3:
            ex_params = {"<J>": self.estimate_exchange()}
            self.ex_params = ex_params
            return ex_params

        # Solve the linear system for E0 and the J_ij
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

        j_ij = np.dot(np.linalg.inv(H), E)

        j_ij[1:] *= 1000  # J_ij in meV (E0 stays in eV)
        ex_params = {j_name: j[0] for j_name, j in zip(j_names, j_ij.tolist(), strict=True)}

        self.ex_params = ex_params
        return ex_params

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

        for s, e in zip(self.ordered_structures, self.energies, strict=True):
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
            for s, e in zip(self.ordered_structures, self.energies, strict=True):
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

        # Convert to magnetic structures with 'magmom' site property
        fm_struct = _analyzer(fm_struct).get_structure_with_only_magnetic_atoms(make_primitive=False)

        afm_struct = _analyzer(afm_struct).get_structure_with_only_magnetic_atoms(make_primitive=False)

        return fm_struct, afm_struct, fm_e, afm_e

    def estimate_exchange(self, fm_struct=None, afm_struct=None, fm_e=None, afm_e=None):
        """Estimate <J> for a structure based on low energy FM and AFM orderings.

        Args:
            fm_struct (Structure): fm structure with 'magmom' site property
            afm_struct (Structure): afm structure with 'magmom' site property
            fm_e (float): fm energy per magnetic ion
            afm_e (float): afm energy per magnetic ion

        Returns:
            float: Average J exchange parameter (meV)
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

    def get_mft_temperature(self, j_avg):
        """
        Crude mean field estimate of critical temperature based on <J> for
        one sublattice, or solving the coupled equations for a multi-sublattice
        material.

        Args:
            j_avg (float): Average exchange parameter (meV)

        Returns:
            float: Critical temperature mft_t (K)
        """
        # Number of magnetic sublattices = number of parent orbits
        n_sub_lattices = len(self.wyckoff_ids)
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

        structure = self.ordered_structures[ordering_index]
        sgraph = self.sgraphs[ordering_index]
        labels = self.site_labels[ordering_index]

        igraph = StructureGraph.from_empty_graph(
            structure, edge_weight_name="exchange_constant", edge_weight_units="meV"
        )

        if "<J>" in self.ex_params:  # Only <J> is available
            logger.warning("Only <J> is available. The interaction graph will not tell you much.")

        # J_ij exchange interaction matrix
        for i in range(len(sgraph.graph.nodes)):
            for c in sgraph.get_connected_sites(i):
                jimage = c[1]  # relative integer coordinates of atom j
                j = c[2]  # index of neighbor
                j_exc = self._get_j_exc(labels[i], labels[j], c[-1])
                # Only add bonds the fit actually parameterized. Unparameterized
                # bonds (no matching sublattice-pair/shell in ex_params, exactly
                # the bonds _build_exchange_mat also ignores) get j_exc == 0, which
                # StructureGraph.add_edge silently drops via its falsy-weight
                # guard, leaving an edge with weight=None that breaks downstream
                # per-bond consumers.
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
            float: Exchange parameter J_exc in meV (0 if no matching interaction).
        """
        shell = self._shell_of(round(dist, 2))
        order = f"-{shell}" if shell else ""
        j_ij = f"{i_id}-{j_id}{order}"
        j_ji = f"{j_id}-{i_id}{order}"

        if j_ij in self.ex_params:
            j_exc = self.ex_params[j_ij]
        elif j_ji in self.ex_params:
            j_exc = self.ex_params[j_ji]
        else:
            j_exc = 0

        # Check if only averaged NN <J> values are available
        if "<J>" in self.ex_params and order == "-nn":
            j_exc = self.ex_params["<J>"]

        return j_exc

    def get_heisenberg_model(self):
        """Save results of mapping to a HeisenbergModel object.

        Returns:
            HeisenbergModel: MSONable object.
        """
        return HeisenbergModel(
            formula=str(self.ordered_structures_[0].reduced_formula),
            structures=self.ordered_structures,
            energies=self.energies,
            cutoff=self.cutoff,
            tol=self.tol,
            sgraphs=self.sgraphs,
            site_labels=self.site_labels,
            wyckoff_ids=self.wyckoff_ids,
            nn_interactions=self.nn_interactions,
            dists=self.dists,
            ex_mat=self.ex_mat,
            ex_params=self.get_exchange(),
            javg=self.estimate_exchange(),
            igraph=self.get_interaction_graph(),
        )


class HeisenbergScreener:
    """Clean and screen magnetic orderings."""

    def __init__(self, structures, energies, screen=False):
        """Pre-processes magnetic orderings and energies for HeisenbergMapper.
        It prioritizes low-energy orderings with large and localized magnetic moments.

        Args:
            structures (list): Structure objects with magnetic moments.
            energies (list): Total energies of magnetic orderings.
            screen (bool): Try to screen out high energy and low-spin configurations.

        Attributes:
            screened_structures (list): Sorted magnetic-only structures.
            screened_full_structures (list): The same orderings with nonmagnetic
                ions retained, aligned index-for-index with screened_structures.
            screened_energies (list): Sorted energies.
        """
        # Cleanup
        full_structures, structures, energies = self._do_cleanup(structures, energies)

        n_structures = len(structures)

        # If there are more than 2 structures, we want to perform a
        # screening to prioritize well-behaved orderings
        if screen and n_structures > 2:
            full_structures, structures, energies = self._do_screen(full_structures, structures, energies)

        self.screened_structures = structures
        self.screened_full_structures = full_structures
        self.screened_energies = energies

    @staticmethod
    def _do_cleanup(structures, energies):
        """Sanitize input structures and energies.

        Takes magnetic structures and performs the following operations
        - Standardizes magmoms (every site gets a ['magmom'] site prop)
        - Builds a magnetic-only copy of each ordering, keeping a full copy too
        - Converts total energies -> energy / magnetic ion
        - Checks for duplicate/degenerate orderings
        - Sorts by energy

        Args:
            structures (list): Structure objects with magmoms.
            energies (list): Corresponding energies.

        Returns:
            full_structures (list): Sanitized orderings with all ions retained.
            ordered_structures (list): The same orderings with only magnetic ions.
            ordered_energies (list): Sorted energies (per magnetic ion).
        """
        # Standardize moments (zero threshold so small moments are preserved).
        # Keep both a full copy (all ions, needed for correct site symmetry) and
        # a magnetic-only copy (used to build the interaction graphs).
        # Make sure no induced moments are present (threshold_nonmag=100).
        analyzers = [_analyzer(s) for s in structures]
        full_structures = [analyzer.structure for analyzer in analyzers]
        ordered_structures = [
            analyzer.get_structure_with_only_magnetic_atoms(make_primitive=False) for analyzer in analyzers
        ]

        # Convert to energies / magnetic ion
        energies = [e / len(s) for (e, s) in zip(energies, ordered_structures, strict=True)]

        # Check for duplicate / degenerate states (sometimes different initial
        # configs relax to the same state)
        remove_list = []
        e_tol = 6  # 10^-6 eV/atom tol on energies

        for idx, energy in enumerate(energies):
            energy = round(energy, e_tol)
            if idx not in remove_list:
                for i_check, e_check in enumerate(energies):
                    e_check = round(e_check, e_tol)
                    if idx != i_check and i_check not in remove_list and energy == e_check:
                        remove_list.append(i_check)

        # Drop duplicates, then sort by energy, keeping the three lists aligned.
        keep = [idx for idx in range(len(energies)) if idx not in remove_list]
        keep.sort(key=lambda idx: energies[idx])

        full_structures = [full_structures[idx] for idx in keep]
        ordered_structures = [ordered_structures[idx] for idx in keep]
        ordered_energies = [energies[idx] for idx in keep]

        return full_structures, ordered_structures, ordered_energies

    @staticmethod
    def _do_screen(full_structures, structures, energies):
        """Screen and sort magnetic orderings based on some criteria.

        Prioritize low energy orderings and large, localized magmoms. do_clean should be run first to sanitize inputs.

        Args:
            full_structures (list): Orderings with all ions retained.
            structures (list): The same orderings with only magnetic ions.
            energies (list): Energies.

        Returns:
            screened_full_structures (list): Sorted full structures.
            screened_structures (list): Sorted magnetic-only structures.
            screened_energies (list): Sorted energies.
        """
        magmoms = [struct.site_properties["magmom"] for struct in structures]
        n_below_1ub = [sum(abs(m) < 1 for m in ms) for ms in magmoms]

        df_mag = pd.DataFrame(
            {
                "full": full_structures,
                "structure": structures,
                "energy": energies,
                "magmoms": magmoms,
                "n_below_1ub": n_below_1ub,
            }
        )

        # keep the ground and first excited state fixed to capture the
        # low-energy spectrum
        df_high_energy = df_mag.iloc[2:]

        # Prioritize structures with fewer magmoms < 1 uB
        df_high_energy = df_high_energy.sort_values(by="n_below_1ub")

        index = [0, 1, *df_high_energy.index]

        # sort
        df_mag = df_mag.reindex(index)
        return (
            list(df_mag["full"].values),
            list(df_mag["structure"].values),
            list(df_mag["energy"].values),
        )


class HeisenbergModel(MSONable):
    """
    Store a Heisenberg model fit to low-energy magnetic orderings.
    Intended to be generated by HeisenbergMapper.get_heisenberg_model().
    """

    def __init__(
        self,
        formula=None,
        structures=None,
        energies=None,
        cutoff=None,
        tol=None,
        sgraphs=None,
        site_labels=None,
        wyckoff_ids=None,
        nn_interactions=None,
        dists=None,
        ex_mat=None,
        ex_params=None,
        javg=None,
        igraph=None,
    ):
        """
        Args:
            formula (str): Reduced formula of compound.
            structures (list): Structure objects with magmoms.
            energies (list): Energies of each relaxed magnetic structure.
            cutoff (float): Cutoff in Angstrom for nearest neighbor search.
            tol (float): Tolerance (in Angstrom) on nearest neighbor distances being equal.
            sgraphs (list): StructureGraph objects.
            site_labels (list[list[int]]): site_labels[k][i] is the sublattice id of
                site i in ordering k.
            wyckoff_ids (dict): Maps each sublattice id to its wyckoff position.
            nn_interactions (dict): {i: j} pairs of NN interactions between sublattices.
            dists (dict): NN, NNN, and NNNN interaction distances.
            ex_mat (DataFrame): Heisenberg Hamiltonian (per magnetic ion) for each ordering.
            ex_params (dict): Exchange parameter values. The J_ij are in meV; the
                included 'E0' offset stays in eV.
            javg (float): <J> exchange param (meV).
            igraph (StructureGraph): Exchange interaction graph.
        """
        self.formula = formula
        self.structures = structures
        self.energies = energies
        self.cutoff = cutoff
        self.tol = tol
        self.sgraphs = sgraphs
        self.site_labels = site_labels
        self.wyckoff_ids = wyckoff_ids
        self.nn_interactions = nn_interactions
        self.dists = dists
        self.ex_mat = ex_mat
        self.ex_params = ex_params
        self.javg = javg
        self.igraph = igraph

    @property
    def unique_site_ids(self):
        """list[dict]: One {tuple(site indices): sublattice id} map per ordering."""
        return [HeisenbergMapper._site_ids_dict(labels) for labels in self.site_labels]

    def as_dict(self):
        """Because some dicts have int keys, some sanitization is required for JSON compatibility."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "@version": __version__,
            "formula": self.formula,
            "structures": [struct.as_dict() for struct in self.structures],
            "energies": self.energies,
            "cutoff": self.cutoff,
            "tol": self.tol,
            "sgraphs": [sgraph.as_dict() for sgraph in self.sgraphs],
            "site_labels": self.site_labels,
            "dists": self.dists,
            "ex_params": self.ex_params,
            "javg": self.javg,
            "igraph": self.igraph.as_dict(),
            # Sanitize int keys / DataFrame
            "ex_mat": jsanitize(self.ex_mat.to_dict()),
            "nn_interactions": jsanitize(self.nn_interactions),
            "wyckoff_ids": jsanitize(self.wyckoff_ids),
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Create a HeisenbergModel from a dict."""
        # Reconstitute the int-keyed dicts that jsanitize stringified
        nn_interactions = {
            shell: {literal_eval(i): j for i, j in pairs.items()} for shell, pairs in dct["nn_interactions"].items()
        }
        wyckoff_ids = {literal_eval(k): v for k, v in dct["wyckoff_ids"].items()}

        structures = [Structure.from_dict(v) for v in dct["structures"]]
        sgraphs = [StructureGraph.from_dict(v) for v in dct["sgraphs"]]
        igraph = StructureGraph.from_dict(dct["igraph"])

        # Reconstitute the exchange matrix DataFrame
        ex_mat = pd.DataFrame.from_dict(dct["ex_mat"]) if dct["ex_mat"] else pd.DataFrame(columns=["E", "E0"])

        return cls(
            formula=dct["formula"],
            structures=structures,
            energies=dct["energies"],
            cutoff=dct["cutoff"],
            tol=dct["tol"],
            sgraphs=sgraphs,
            site_labels=dct["site_labels"],
            wyckoff_ids=wyckoff_ids,
            nn_interactions=nn_interactions,
            dists=dct["dists"],
            ex_mat=ex_mat,
            ex_params=dct["ex_params"],
            javg=dct["javg"],
            igraph=igraph,
        )
