"""This module defines tools to generate and analyze phase diagrams."""

from __future__ import annotations

import itertools
import json
import logging
import math
import os
import re
import warnings
from collections import defaultdict
from functools import lru_cache
from typing import TYPE_CHECKING, no_type_check

import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
from matplotlib import cm
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.font_manager import FontProperties
from monty.json import MontyDecoder, MSONable
from scipy import interpolate
from scipy.optimize import minimize
from scipy.spatial import ConvexHull
from tqdm import tqdm

from pymatgen.analysis.reaction_calculator import Reaction, ReactionError
from pymatgen.core import DummySpecies, Element, get_el_sp
from pymatgen.core.composition import Composition
from pymatgen.entries import Entry
from pymatgen.util.coord import Simplex, in_coord_list
from pymatgen.util.due import Doi, due
from pymatgen.util.plotting import pretty_plot
from pymatgen.util.string import htmlify, latexify

if TYPE_CHECKING:
    from collections.abc import Collection, Iterator, Sequence
    from io import StringIO
    from typing import Any, Literal

    from numpy.typing import ArrayLike
    from typing_extensions import Self

logger = logging.getLogger(__name__)

with open(os.path.join(os.path.dirname(__file__), "..", "util", "plotly_pd_layouts.json"), encoding="utf-8") as file:
    plotly_layouts = json.load(file)


class PDEntry(Entry):
    """
    An object encompassing all relevant data for phase diagrams.

    Attributes:
        composition (Composition): The composition associated with the PDEntry.
        energy (float): The energy associated with the entry.
        name (str): A name for the entry. This is the string shown in the phase diagrams.
            By default, this is the reduced formula for the composition, but can be
            set to some other string for display purposes.
        attribute (MSONable): A arbitrary attribute. Can be used to specify that the
            entry is a newly found compound, or to specify a particular label for
            the entry, etc. An attribute can be anything but must be MSONable.
    """

    def __init__(
        self,
        composition: Composition,
        energy: float,
        name: str | None = None,
        attribute: object = None,
    ):
        """
        Args:
            composition (Composition): Composition
            energy (float): Energy for composition.
            name (str): Optional parameter to name the entry. Defaults
                to the reduced chemical formula.
            attribute: Optional attribute of the entry. Must be MSONable.
        """
        super().__init__(composition, energy)
        self.name = name or self.reduced_formula
        self.attribute = attribute

    def __repr__(self):
        name = ""
        if self.name != self.reduced_formula:
            name = f" ({self.name})"
        return f"{type(self).__name__} : {self.composition}{name} with energy = {self.energy:.4f}"

    @property
    def energy(self) -> float:
        """The entry's energy."""
        return self._energy

    def as_dict(self):
        """Get MSONable dict representation of PDEntry."""
        return super().as_dict() | {"name": self.name, "attribute": self.attribute}

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """
        Args:
            dct (dict): dictionary representation of PDEntry.

        Returns:
            PDEntry
        """
        return cls(
            composition=Composition(dct["composition"]),
            energy=dct["energy"],
            name=dct.get("name"),
            attribute=dct.get("attribute"),
        )


class GrandPotPDEntry(PDEntry):
    """
    A grand potential pd entry object encompassing all relevant data for phase
    diagrams. Chemical potentials are given as a element-chemical potential
    dict.
    """

    def __init__(self, entry, chempots, name=None):
        """
        Args:
            entry: A PDEntry-like object.
            chempots: Chemical potential specification as {Element: float}.
            name: Optional parameter to name the entry. Defaults to the reduced
                chemical formula of the original entry.
        """
        super().__init__(
            entry.composition,
            entry.energy,
            name or entry.name,
            getattr(entry, "attribute", None),
        )
        # NOTE if we init GrandPotPDEntry from ComputedEntry _energy is the
        # corrected energy of the ComputedEntry hence the need to keep
        # the original entry to not lose data.
        self.original_entry = entry
        self.original_comp = self._composition
        self.chempots = chempots

    @property
    def composition(self) -> Composition:
        """The composition after removing free species.

        Returns:
            Composition
        """
        return Composition({el: self._composition[el] for el in self._composition.elements if el not in self.chempots})

    @property
    def chemical_energy(self):
        """The chemical energy term mu*N in the grand potential.

        Returns:
            The chemical energy term mu*N in the grand potential
        """
        return sum(self._composition[el] * pot for el, pot in self.chempots.items())

    @property
    def energy(self) -> float:
        """Grand potential energy."""
        return self._energy - self.chemical_energy

    def __repr__(self):
        output = [
            (
                f"GrandPotPDEntry with original composition {self.original_entry.composition}, "
                f"energy = {self.original_entry.energy:.4f}, "
            ),
            "chempots = " + ", ".join(f"mu_{el} = {mu:.4f}" for el, mu in self.chempots.items()),
        ]
        return "".join(output)

    def as_dict(self):
        """Get MSONable dict representation of GrandPotPDEntry."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "entry": self.original_entry.as_dict(),
            "chempots": {el.symbol: u for el, u in self.chempots.items()},
            "name": self.name,
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """
        Args:
            dct (dict): dictionary representation of GrandPotPDEntry.

        Returns:
            GrandPotPDEntry
        """
        chempots = {Element(symbol): u for symbol, u in dct["chempots"].items()}
        entry = MontyDecoder().process_decoded(dct["entry"])
        return cls(entry, chempots, dct["name"])


class TransformedPDEntry(PDEntry):
    """
    This class represents a TransformedPDEntry, which allows for a PDEntry to be
    transformed to a different composition coordinate space. It is used in the
    construction of phase diagrams that do not have elements as the terminal
    compositions.
    """

    # Tolerance for determining if amount of a composition is positive.
    amount_tol = 1e-5

    def __init__(self, entry, sp_mapping, name=None):
        """
        Args:
            entry (PDEntry): Original entry to be transformed.
            sp_mapping ({Composition: DummySpecies}): dictionary mapping Terminal Compositions to Dummy Species.
        """
        super().__init__(
            entry.composition,
            entry.energy,
            name or entry.name,
            getattr(entry, "attribute", None),
        )
        self.original_entry = entry
        self.sp_mapping = sp_mapping

        self.rxn = Reaction(list(self.sp_mapping), [self._composition])
        self.rxn.normalize_to(self.original_entry.composition)

        # NOTE We only allow reactions that have positive amounts of reactants.
        if not all(self.rxn.get_coeff(comp) <= TransformedPDEntry.amount_tol for comp in self.sp_mapping):
            raise TransformedPDEntryError("Only reactions with positive amounts of reactants allowed")

    @property
    def composition(self) -> Composition:
        """The composition in the dummy species space.

        Returns:
            Composition
        """
        # NOTE this is not infallible as the original entry is mutable and an
        # end user could choose to normalize or change the original entry.
        # However, the risk of this seems low.
        factor = self._composition.num_atoms / self.original_entry.composition.num_atoms

        trans_comp = {self.sp_mapping[comp]: -self.rxn.get_coeff(comp) for comp in self.sp_mapping}

        trans_comp = {k: v * factor for k, v in trans_comp.items() if v > TransformedPDEntry.amount_tol}

        return Composition(trans_comp)

    def __repr__(self):
        output = [
            f"TransformedPDEntry {self.composition}",
            f" with original composition {self.original_entry.composition}",
            f", energy = {self.original_entry.energy:.4f}",
        ]
        return "".join(output)

    def as_dict(self):
        """Get MSONable dict representation of TransformedPDEntry."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "sp_mapping": self.sp_mapping,
            **self.original_entry.as_dict(),
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """
        Args:
            dct (dict): dictionary representation of TransformedPDEntry.

        Returns:
            TransformedPDEntry
        """
        sp_mapping = dct["sp_mapping"]
        del dct["sp_mapping"]
        entry = MontyDecoder().process_decoded(dct)
        return cls(entry, sp_mapping)


class TransformedPDEntryError(Exception):
    """An exception class for TransformedPDEntry."""


@due.dcite(Doi("10.1021/cm702327g"), description="Phase Diagram from First Principles Calculations")
@due.dcite(
    Doi("10.1016/j.elecom.2010.01.010"),
    description="Thermal stabilities of delithiated olivine MPO4 (M=Fe, Mn) cathodes "
    "investigated using first principles calculations",
)
class PhaseDiagram(MSONable):
    """
    Simple phase diagram class taking in elements and entries as inputs.
    The algorithm is based on the work in the following papers:

    1. S. P. Ong, L. Wang, B. Kang, and G. Ceder, Li-Fe-P-O2 Phase Diagram from
        First Principles Calculations. Chem. Mater., 2008, 20(5), 1798-1807.
        doi:10.1021/cm702327g

    2. S. P. Ong, A. Jain, G. Hautier, B. Kang, G. Ceder, Thermal stabilities
        of delithiated olivine MPO4 (M=Fe, Mn) cathodes investigated using first
        principles calculations. Electrochem. Comm., 2010, 12(3), 427-430.
        doi:10.1016/j.elecom.2010.01.010

    Attributes:
        dim (int): The dimensionality of the phase diagram.
        elements: Elements in the phase diagram.
        el_refs: List of elemental references for the phase diagrams. These are
            entries corresponding to the lowest energy element entries for simple
            compositional phase diagrams.
        all_entries: All entries provided for Phase Diagram construction. Note that this
            does not mean that all these entries are actually used in the phase
            diagram. For example, this includes the positive formation energy
            entries that are filtered out before Phase Diagram construction.
        qhull_entries: Actual entries used in convex hull. Excludes all positive formation
            energy entries.
        qhull_data: Data used in the convex hull operation. This is essentially a matrix of
            composition data and energy per atom values created from qhull_entries.
        facets: Facets of the phase diagram in the form of  [[1,2,3],[4,5,6]...].
            For a ternary, it is the indices (references to qhull_entries and
            qhull_data) for the vertices of the phase triangles. Similarly
            extended to higher D simplices for higher dimensions.
        simplices: The simplices of the phase diagram as a list of np.ndarray, i.e.,
            the list of stable compositional coordinates in the phase diagram.
    """

    # Tolerance for determining if formation energy is positive.
    formation_energy_tol = 1e-11
    numerical_tol = 1e-8

    def __init__(
        self,
        entries: Sequence[PDEntry] | set[PDEntry],
        elements: Sequence[Element] = (),
        *,
        computed_data: dict[str, Any] | None = None,
    ) -> None:
        """
        Args:
            entries (list[PDEntry]): A list of PDEntry-like objects having an
                energy, energy_per_atom and composition.
            elements (list[Element]): Optional list of elements in the phase
                diagram. If set to None, the elements are determined from
                the entries themselves and are sorted alphabetically.
                If specified, element ordering (e.g. for pd coordinates)
                is preserved.
            computed_data (dict): A dict containing pre-computed data. This allows
                PhaseDiagram object to be reconstituted without performing the
                expensive convex hull computation. The dict is the output from the
                PhaseDiagram._compute() method and is stored in PhaseDiagram.computed_data
                when generated for the first time.
        """
        if not entries:
            raise ValueError("Unable to build phase diagram without entries.")

        self.elements = elements
        self.entries = entries
        if computed_data is None:
            computed_data = self._compute()
        else:
            computed_data = MontyDecoder().process_decoded(computed_data)
            if not isinstance(computed_data, dict):
                raise TypeError(f"computed_data should be dict, got {type(computed_data).__name__}")

            # Update keys to be Element objects in case they are strings in pre-computed data
            computed_data["el_refs"] = [(Element(el_str), entry) for el_str, entry in computed_data["el_refs"]]
        self.computed_data = computed_data
        self.facets = computed_data["facets"]
        self.simplexes = computed_data["simplexes"]
        self.all_entries = computed_data["all_entries"]
        self.qhull_data = computed_data["qhull_data"]
        self.dim = computed_data["dim"]
        self.el_refs = dict(computed_data["el_refs"])
        self.qhull_entries = tuple(computed_data["qhull_entries"])
        self._qhull_spaces = tuple(frozenset(e.elements) for e in self.qhull_entries)
        self._stable_entries = tuple({self.qhull_entries[idx] for idx in set(itertools.chain(*self.facets))})
        self._stable_spaces = tuple(frozenset(e.elements) for e in self._stable_entries)

    def as_dict(self):
        """Get MSONable dict representation of PhaseDiagram."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "all_entries": [e.as_dict() for e in self.all_entries],
            "elements": [e.as_dict() for e in self.elements],
            "computed_data": self.computed_data,
        }

    @classmethod
    def from_dict(cls, dct: dict[str, Any]) -> Self:
        """
        Args:
            dct (dict): dictionary representation of PhaseDiagram.

        Returns:
            PhaseDiagram
        """
        entries = [MontyDecoder().process_decoded(entry) for entry in dct["all_entries"]]
        elements = [Element.from_dict(elem) for elem in dct["elements"]]
        computed_data = dct.get("computed_data")
        return cls(entries, elements, computed_data=computed_data)

    def _compute(self) -> dict[str, Any]:
        if self.elements == ():
            self.elements = sorted({els for e in self.entries for els in e.elements})

        elements = list(self.elements)
        dim = len(elements)

        entries = sorted(self.entries, key=lambda e: e.composition.reduced_composition)

        el_refs: dict[Element, PDEntry] = {}
        min_entries: list[PDEntry] = []
        all_entries: list[PDEntry] = []
        for composition, group_iter in itertools.groupby(entries, key=lambda e: e.composition.reduced_composition):
            group = list(group_iter)
            min_entry = min(group, key=lambda e: e.energy_per_atom)
            if composition.is_element:
                el_refs[composition.elements[0]] = min_entry
            min_entries.append(min_entry)
            all_entries.extend(group)

        if missing := set(elements) - set(el_refs):
            raise ValueError(f"Missing terminal entries for elements {sorted(map(str, missing))}")
        if extra := set(el_refs) - set(elements):
            raise ValueError(f"There are more terminal elements than dimensions: {sorted(map(str, extra))}")

        data = np.array(
            [[e.composition.get_atomic_fraction(el) for el in elements] + [e.energy_per_atom] for e in min_entries]
        )

        # Use only entries with negative formation energy
        vec = [el_refs[el].energy_per_atom for el in elements] + [-1]
        form_e = -np.dot(data, vec)
        idx = np.where(form_e < -PhaseDiagram.formation_energy_tol)[0].tolist()

        # Add the elemental references
        idx.extend([min_entries.index(el) for el in el_refs.values()])

        qhull_entries = [min_entries[idx] for idx in idx]
        qhull_data = data[idx][:, 1:]

        # Add an extra point to enforce full dimensionality.
        # This point will be present in all upper hull facets.
        extra_point = np.zeros(dim) + 1 / dim
        extra_point[-1] = np.max(qhull_data) + 1
        qhull_data = np.concatenate([qhull_data, [extra_point]], axis=0)

        if dim == 1:
            facets = [qhull_data.argmin(axis=0)]
        else:
            facets = get_facets(qhull_data)
            final_facets = []
            for facet in facets:
                # Skip facets that include the extra point
                if max(facet) == len(qhull_data) - 1:
                    continue
                mat = qhull_data[facet]
                mat[:, -1] = 1
                if abs(np.linalg.det(mat)) > 1e-14:
                    final_facets.append(facet)
            facets = final_facets

        simplexes = [Simplex(qhull_data[facet, :-1]) for facet in facets]
        self.elements = elements
        return {
            "facets": facets,
            "simplexes": simplexes,
            "all_entries": all_entries,
            "qhull_data": qhull_data,
            "dim": dim,
            # Dictionary with Element keys is not JSON-serializable
            "el_refs": list(el_refs.items()),
            "qhull_entries": qhull_entries,
        }

    def pd_coords(self, comp: Composition) -> np.ndarray:
        """
        The phase diagram is generated in a reduced dimensional space
        (n_elements - 1). This function returns the coordinates in that space.
        These coordinates are compatible with the stored simplex objects.

        Args:
            comp (Composition): A composition

        Returns:
            The coordinates for a given composition in the PhaseDiagram's basis
        """
        if set(comp.elements) - set(self.elements):
            raise ValueError(f"{comp} has elements not in the phase diagram {', '.join(map(str, self.elements))}")
        return np.array([comp.get_atomic_fraction(el) for el in self.elements[1:]])

    @property
    def all_entries_hulldata(self):
        """The ndarray used to construct the convex hull."""
        data = [
            [e.composition.get_atomic_fraction(el) for el in self.elements] + [e.energy_per_atom]
            for e in self.all_entries
        ]
        return np.array(data)[:, 1:]

    @property
    def unstable_entries(self) -> set[Entry]:
        """
        Returns:
            set[Entry]: unstable entries in the phase diagram. Includes positive formation energy entries.
        """
        return {e for e in self.all_entries if e not in self.stable_entries}

    @property
    def stable_entries(self) -> set[Entry]:
        """
        Returns:
            set[Entry]: of stable entries in the phase diagram.
        """
        return set(self._stable_entries)

    @lru_cache(1)  # noqa: B019
    def _get_stable_entries_in_space(self, space) -> list[Entry]:
        """
        Args:
            space (set[Element]): set of Element objects.

        Returns:
            list[Entry]: stable entries in the space.
        """
        return [e for e, s in zip(self._stable_entries, self._stable_spaces, strict=True) if space.issuperset(s)]

    def get_reference_energy(self, comp: Composition) -> float:
        """Sum of elemental reference energies over all elements in a composition.

        Args:
            comp (Composition): Input composition.

        Returns:
            float: Reference energy
        """
        return sum(comp[el] * self.el_refs[el].energy_per_atom for el in comp.elements)

    def get_reference_energy_per_atom(self, comp: Composition) -> float:
        """Sum of elemental reference energies over all elements in a composition.

        Args:
            comp (Composition): Input composition.

        Returns:
            float: Reference energy per atom
        """
        return self.get_reference_energy(comp) / comp.num_atoms

    def get_form_energy(self, entry: PDEntry) -> float:
        """Get the formation energy for an entry (NOT normalized) from the
        elemental references.

        Args:
            entry (PDEntry): A PDEntry-like object.

        Returns:
            float: Formation energy from the elemental references.
        """
        comp = entry.composition
        return entry.energy - self.get_reference_energy(comp)

    def get_form_energy_per_atom(self, entry: PDEntry) -> float:
        """Get the formation energy per atom for an entry from the
        elemental references.

        Args:
            entry (PDEntry): An PDEntry-like object

        Returns:
            Formation energy **per atom** from the elemental references.
        """
        return self.get_form_energy(entry) / entry.composition.num_atoms

    def __repr__(self) -> str:
        symbols = [el.symbol for el in self.elements]
        output = [
            f"{'-'.join(symbols)} phase diagram",
            f"{len(self.stable_entries)} stable phases: ",
            ", ".join(entry.name for entry in sorted(self.stable_entries, key=str)),
        ]
        return "\n".join(output)

    @lru_cache(1)  # noqa: B019
    def _get_facet_and_simplex(self, comp: Composition) -> tuple[Simplex, Simplex]:
        """Get any facet that a composition falls into. Cached so successive
        calls at same composition are fast.

        Args:
            comp (Composition): A composition
        """
        coord = self.pd_coords(comp)
        for facet, simplex in zip(self.facets, self.simplexes, strict=True):
            if simplex.in_simplex(coord, PhaseDiagram.numerical_tol / 10):
                return facet, simplex

        raise RuntimeError(f"No facet found for {comp = }")

    def _get_all_facets_and_simplexes(self, comp):
        """Get all facets that a composition falls into.

        Args:
            comp (Composition): A composition
        """
        coords = self.pd_coords(comp)

        all_facets = [
            facet
            for facet, simplex in zip(self.facets, self.simplexes, strict=True)
            if simplex.in_simplex(coords, PhaseDiagram.numerical_tol / 10)
        ]

        if not all_facets:
            raise RuntimeError(f"No facets found for {comp = }")

        return all_facets

    def _get_facet_chempots(self, facet: list[int]) -> dict[Element, float]:
        """
        Calculates the chemical potentials for each element within a facet.

        Args:
            facet (list): Indices of the entries in the facet.

        Returns:
            dict[Element, float]: Chemical potentials for each element in the facet.
        """
        comp_list = [self.qhull_entries[idx].composition for idx in facet]
        energy_list = [self.qhull_entries[idx].energy_per_atom for idx in facet]
        atom_frac_mat = [[c.get_atomic_fraction(e) for e in self.elements] for c in comp_list]
        chempots = np.linalg.solve(atom_frac_mat, energy_list)

        return dict(zip(self.elements, chempots, strict=True))

    def _get_simplex_intersections(self, c1, c2):
        """Get coordinates of the intersection of the tie line between two compositions
        and the simplexes of the PhaseDiagram.

        Args:
            c1: Reduced dimension coordinates of first composition
            c2: Reduced dimension coordinates of second composition

        Returns:
            Array of the intersections between the tie line and the simplexes of
            the PhaseDiagram
        """
        intersections = [c1, c2]
        for sc in self.simplexes:
            intersections.extend(sc.line_intersection(c1, c2))

        return np.array(intersections)

    def get_decomposition(self, comp: Composition) -> dict[PDEntry, float]:
        """
        Provides the decomposition at a particular composition.

        Args:
            comp (Composition): A composition

        Returns:
            Decomposition as a dict of {PDEntry: amount} where amount
            is the amount of the fractional composition.
        """
        facet, simplex = self._get_facet_and_simplex(comp)
        decomp_amts = simplex.bary_coords(self.pd_coords(comp))
        return {
            self.qhull_entries[f]: amt
            for f, amt in zip(facet, decomp_amts, strict=True)
            if abs(amt) > PhaseDiagram.numerical_tol
        }

    def get_decomp_and_hull_energy_per_atom(self, comp: Composition) -> tuple[dict[PDEntry, float], float]:
        """
        Args:
            comp (Composition): Input composition.

        Returns:
            Energy of lowest energy equilibrium at desired composition per atom
        """
        decomp = self.get_decomposition(comp)
        return decomp, sum(e.energy_per_atom * n for e, n in decomp.items())

    def get_hull_energy_per_atom(self, comp: Composition, **kwargs) -> float:
        """
        Args:
            comp (Composition): Input composition.

        Returns:
            Energy of lowest energy equilibrium at desired composition.
        """
        return self.get_decomp_and_hull_energy_per_atom(comp, **kwargs)[1]

    def get_hull_energy(self, comp: Composition) -> float:
        """
        Args:
            comp (Composition): Input composition.

        Returns:
            Energy of lowest energy equilibrium at desired composition. Not
                normalized by atoms, i.e. E(Li4O2) = 2 * E(Li2O)
        """
        return comp.num_atoms * self.get_hull_energy_per_atom(comp)

    def get_decomp_and_e_above_hull(
        self,
        entry: PDEntry,
        allow_negative: bool = False,
        check_stable: bool = True,
        on_error: Literal["raise", "warn", "ignore"] = "raise",
    ) -> tuple[dict[PDEntry, float], float] | tuple[None, None]:
        """
        Provides the decomposition and energy above convex hull for an entry.
        Due to caching, can be much faster if entries with the same composition
        are processed together.

        Args:
            entry (PDEntry): A PDEntry like object
            allow_negative (bool): Whether to allow negative e_above_hulls. Used to
                calculate equilibrium reaction energies. Defaults to False.
            check_stable (bool): Whether to first check whether an entry is stable.
                In normal circumstances, this is the faster option since checking for
                stable entries is relatively fast. However, if you have a huge proportion
                of unstable entries, then this check can slow things down. You should then
                set this to False.
            on_error ('raise' | 'warn' | 'ignore'): What to do if no valid decomposition was
                found. 'raise' will throw ValueError. 'warn' will print return (None, None).
                'ignore' just returns (None, None). Defaults to 'raise'.

        Raises:
            ValueError: If on_error is 'raise' and no valid decomposition exists in this
                phase diagram for given entry.

        Returns:
            tuple[decomp, energy_above_hull]: The decomposition is provided
                as a dict of {PDEntry: amount} where amount is the amount of the
                fractional composition. Stable entries should have energy above
                convex hull of 0. The energy is given per atom.
        """
        # Avoid computation for stable_entries.
        # NOTE scaled duplicates of stable_entries will not be caught.
        if check_stable and entry in self.stable_entries:
            return {entry: 1.0}, 0.0

        try:
            decomp, hull_energy = self.get_decomp_and_hull_energy_per_atom(entry.composition)
        except Exception as exc:
            if on_error == "raise":
                raise ValueError(f"Unable to get decomposition for {entry}") from exc
            if on_error == "warn":
                warnings.warn(f"Unable to get decomposition for {entry}, encountered {exc}")
            return None, None
        e_above_hull = entry.energy_per_atom - hull_energy

        if allow_negative or e_above_hull >= -PhaseDiagram.numerical_tol:
            return decomp, e_above_hull

        msg = f"No valid decomposition found for {entry}! (e_h: {e_above_hull})"
        if on_error == "raise":
            raise ValueError(msg)
        if on_error == "warn":
            warnings.warn(msg)
        return None, None  # 'ignore' and 'warn' case

    def get_e_above_hull(self, entry: PDEntry, **kwargs: Any) -> float | None:
        """
        Provides the energy above convex hull for an entry.

        Args:
            entry (PDEntry): A PDEntry like object.
            **kwargs: Passed to get_decomp_and_e_above_hull().

        Returns:
            float | None: Energy above convex hull of entry. Stable entries should have
                energy above hull of 0. The energy is given per atom.
        """
        return self.get_decomp_and_e_above_hull(entry, **kwargs)[1]

    def get_equilibrium_reaction_energy(self, entry: PDEntry) -> float | None:
        """
        Provides the reaction energy of a stable entry from the neighboring
        equilibrium stable entries (also known as the inverse distance to
        hull).

        Args:
            entry (PDEntry): A PDEntry like object

        Returns:
            float | None: Equilibrium reaction energy of entry. Stable entries should have
                equilibrium reaction energy <= 0. The energy is given per atom.
        """
        elem_space = entry.elements

        # NOTE scaled duplicates of stable_entries will not be caught.
        if entry not in self._get_stable_entries_in_space(frozenset(elem_space)):
            raise ValueError(
                f"{entry} is unstable, the equilibrium reaction energy is available only for stable entries."
            )

        if entry.is_element:
            return 0

        entries = [e for e in self._get_stable_entries_in_space(frozenset(elem_space)) if e != entry]
        mod_pd = PhaseDiagram(entries, elements=elem_space)

        return mod_pd.get_decomp_and_e_above_hull(entry, allow_negative=True)[1]

    def get_decomp_and_phase_separation_energy(
        self,
        entry: PDEntry,
        space_limit: int = 200,
        stable_only: bool = False,
        tols: Sequence[float] = (1e-8,),
        maxiter: int = 1000,
        **kwargs: Any,
    ) -> tuple[dict[PDEntry, float], float] | tuple[None, None]:
        """
        Provides the combination of entries in the PhaseDiagram that gives the
        lowest formation enthalpy with the same composition as the given entry
        excluding entries with the same composition and the energy difference
        per atom between the given entry and the energy of the combination found.

        For unstable entries that are not polymorphs of stable entries (or completely
        novel entries) this is simply the energy above (or below) the convex hull.

        For entries with the same composition as one of the stable entries in the
        phase diagram setting `stable_only` to `False` (Default) allows for entries
        not previously on the convex hull to be considered in the combination.
        In this case the energy returned is what is referred to as the decomposition
        enthalpy in:

        1. Bartel, C., Trewartha, A., Wang, Q., Dunn, A., Jain, A., Ceder, G.,
            A critical examination of compound stability predictions from
            machine-learned formation energies, npj Computational Materials 6, 97 (2020)

        For stable entries setting `stable_only` to `True` returns the same energy
        as `get_equilibrium_reaction_energy`. This function is based on a constrained
        optimization rather than recalculation of the convex hull making it
        algorithmically cheaper. However, if `tol` is too loose there is potential
        for this algorithm to converge to a different solution.

        Args:
            entry (PDEntry): A PDEntry like object.
            space_limit (int): The maximum number of competing entries to consider
                before calculating a second convex hull to reducing the complexity
                of the optimization.
            stable_only (bool): Only use stable materials as competing entries.
            tols (list[float]): Tolerances for convergence of the SLSQP optimization
                when finding the equilibrium reaction. Tighter tolerances tested first.
            maxiter (int): The maximum number of iterations of the SLSQP optimizer
                when finding the equilibrium reaction.
            **kwargs: Passed to get_decomp_and_e_above_hull.

        Returns:
            tuple[decomp, energy]: The decomposition is given as a dict of {PDEntry, amount}
                for all entries in the decomp reaction where amount is the amount of the
                fractional composition. The phase separation energy is given per atom.
        """
        entry_frac = entry.composition.fractional_composition
        entry_elems = frozenset(entry_frac.elements)

        # Handle elemental materials
        if entry.is_element:
            return self.get_decomp_and_e_above_hull(entry, allow_negative=True, **kwargs)

        # Select space to compare against
        if stable_only:
            compare_entries = self._get_stable_entries_in_space(entry_elems)
        else:
            compare_entries = [
                e for e, s in zip(self.qhull_entries, self._qhull_spaces, strict=True) if entry_elems.issuperset(s)
            ]

        # get memory ids of entries with the same composition.
        same_comp_mem_ids = [
            id(c)
            for c in compare_entries
            # NOTE use this construction to avoid calls to fractional_composition
            if (
                len(entry_frac) == len(c.composition)
                and all(
                    abs(v - c.composition.get_atomic_fraction(el)) <= Composition.amount_tolerance
                    for el, v in entry_frac.items()
                )
            )
        ]

        if not any(id(e) in same_comp_mem_ids for e in self._get_stable_entries_in_space(entry_elems)):
            return self.get_decomp_and_e_above_hull(entry, allow_negative=True, **kwargs)

        # take entries with negative e_form and different compositions as competing entries
        competing_entries = {c for c in compare_entries if id(c) not in same_comp_mem_ids}

        # NOTE SLSQP optimizer doesn't scale well for > 300 competing entries.
        if len(competing_entries) > space_limit and not stable_only:
            warnings.warn(
                f"There are {len(competing_entries)} competing entries "
                f"for {entry.composition} - Calculating inner hull to discard additional unstable entries"
            )

            reduced_space = competing_entries - {*self._get_stable_entries_in_space(entry_elems)} | {
                *self.el_refs.values()
            }

            # NOTE calling PhaseDiagram is only reasonable if the composition has fewer than 5 elements
            inner_hull = PhaseDiagram(reduced_space)

            competing_entries = inner_hull.stable_entries | {*self._get_stable_entries_in_space(entry_elems)}
            competing_entries = {c for c in compare_entries if id(c) not in same_comp_mem_ids}

        if len(competing_entries) > space_limit:
            warnings.warn(
                f"There are {len(competing_entries)} competing entries "
                f"for {entry.composition} - Using SLSQP to find decomposition likely to be slow"
            )

        decomp = _get_slsqp_decomp(entry.composition, competing_entries, tols, maxiter)

        # find the minimum alternative formation energy for the decomposition
        decomp_enthalpy = np.sum([c.energy_per_atom * amt for c, amt in decomp.items()])

        decomp_enthalpy = entry.energy_per_atom - decomp_enthalpy

        return decomp, decomp_enthalpy

    def get_phase_separation_energy(self, entry, **kwargs):
        """
        Provides the energy to the convex hull for the given entry. For stable entries
        already in the phase diagram the algorithm provides the phase separation energy
        which is referred to as the decomposition enthalpy in:

        1. Bartel, C., Trewartha, A., Wang, Q., Dunn, A., Jain, A., Ceder, G.,
            A critical examination of compound stability predictions from
            machine-learned formation energies, npj Computational Materials 6, 97 (2020)

        Args:
            entry (PDEntry): A PDEntry like object
            **kwargs: Keyword args passed to `get_decomp_and_decomp_energy`
                space_limit (int): The maximum number of competing entries to consider.
                stable_only (bool): Only use stable materials as competing entries
                tol (float): The tolerance for convergence of the SLSQP optimization
                    when finding the equilibrium reaction.
                maxiter (int): The maximum number of iterations of the SLSQP optimizer
                    when finding the equilibrium reaction.

        Returns:
            phase separation energy per atom of entry. Stable entries should have
            energies <= 0, Stable elemental entries should have energies = 0 and
            unstable entries should have energies > 0. Entries that have the same
            composition as a stable energy may have positive or negative phase
            separation energies depending on their own energy.
        """
        return self.get_decomp_and_phase_separation_energy(entry, **kwargs)[1]

    def get_composition_chempots(self, comp):
        """Get the chemical potentials for all elements at a given composition.

        Args:
            comp (Composition): Composition

        Returns:
            Dictionary of chemical potentials.
        """
        facet = self._get_facet_and_simplex(comp)[0]
        return self._get_facet_chempots(facet)

    def get_all_chempots(self, comp):
        """Get chemical potentials at a given composition.

        Args:
            comp (Composition): Composition

        Returns:
            Chemical potentials.
        """
        all_facets = self._get_all_facets_and_simplexes(comp)

        chempots = {}
        for facet in all_facets:
            facet_name = "-".join(self.qhull_entries[j].name for j in facet)
            chempots[facet_name] = self._get_facet_chempots(facet)

        return chempots

    def get_transition_chempots(self, element):
        """Get the critical chemical potentials for an element in the Phase
        Diagram.

        Args:
            element: An element. Has to be in the PD in the first place.

        Returns:
            A sorted sequence of critical chemical potentials, from less
            negative to more negative.
        """
        if element not in self.elements:
            raise ValueError("get_transition_chempots can only be called with elements in the phase diagram.")

        critical_chempots = []
        for facet in self.facets:
            chempots = self._get_facet_chempots(facet)
            critical_chempots.append(chempots[element])

        clean_pots = []
        for c in sorted(critical_chempots):
            if len(clean_pots) == 0 or abs(c - clean_pots[-1]) > PhaseDiagram.numerical_tol:
                clean_pots.append(c)
        clean_pots.reverse()
        return tuple(clean_pots)

    def get_critical_compositions(self, comp1, comp2):
        """Get the critical compositions along the tieline between two
        compositions. I.e. where the decomposition products change.
        The endpoints are also returned.

        Args:
            comp1 (Composition): First composition to define the tieline
            comp2 (Composition): Second composition to define the tieline

        Returns:
            [(Composition)]: list of critical compositions. All are of
                the form x * comp1 + (1-x) * comp2
        """
        n1 = comp1.num_atoms
        n2 = comp2.num_atoms
        pd_els = self.elements

        # NOTE the reduced dimensionality Simplexes don't use the
        # first element in the PD
        c1 = self.pd_coords(comp1)
        c2 = self.pd_coords(comp2)

        # NOTE none of the projections work if c1 == c2, so just
        # return *copies* of the inputs
        if np.all(c1 == c2):
            return [comp1.copy(), comp2.copy()]

        intersections = self._get_simplex_intersections(c1, c2)

        # find position along line
        line = c2 - c1
        line /= np.sum(line**2) ** 0.5
        proj = np.dot(intersections - c1, line)

        # only take compositions between endpoints
        proj = proj[
            np.logical_and(proj > -self.numerical_tol, proj < proj[1] + self.numerical_tol)  # proj[1] is |c2-c1|
        ]
        proj.sort()

        # only unique compositions
        valid = np.ones(len(proj), dtype=bool)
        valid[1:] = proj[1:] > proj[:-1] + self.numerical_tol
        proj = proj[valid]

        ints = c1 + line * proj[:, None]

        # reconstruct full-dimensional composition array
        cs = np.concatenate([np.array([1 - np.sum(ints, axis=-1)]).T, ints], axis=-1)

        # mixing fraction when compositions are normalized
        x = proj / np.dot(c2 - c1, line)

        # mixing fraction when compositions are not normalized
        x_unnormalized = x * n1 / (n2 + x * (n1 - n2))
        num_atoms = n1 + (n2 - n1) * x_unnormalized
        cs *= num_atoms[:, None]

        return [Composition((elem, val) for elem, val in zip(pd_els, m, strict=True)) for m in cs]

    def get_element_profile(self, element, comp, comp_tol=1e-5):
        """
        Provides the element evolution data for a composition. For example, can be used
        to analyze Li conversion voltages by varying mu_Li and looking at the phases
        formed. Also can be used to analyze O2 evolution by varying mu_O2.

        Args:
            element: An element. Must be in the phase diagram.
            comp: A Composition
            comp_tol: The tolerance to use when calculating decompositions.
                Phases with amounts less than this tolerance are excluded.
                Defaults to 1e-5.

        Returns:
            Evolution data as a list of dictionaries of the following format:
            [ {'chempot': -10.487582, 'evolution': -2.0,
            'reaction': Reaction Object], ...]
        """
        element = get_el_sp(element)

        if element not in self.elements:
            raise ValueError("get_transition_chempots can only be called with elements in the phase diagram.")

        gc_comp = Composition({el: amt for el, amt in comp.items() if el != element})
        el_ref = self.el_refs[element]
        el_comp = Composition(element.symbol)
        evolution = []

        for cc in self.get_critical_compositions(el_comp, gc_comp)[1:]:
            decomp_entries = list(self.get_decomposition(cc))
            decomp = [k.composition for k in decomp_entries]
            rxn = Reaction([comp], [*decomp, el_comp])
            rxn.normalize_to(comp)
            c = self.get_composition_chempots(cc + el_comp * 1e-5)[element]
            amt = -rxn.coeffs[rxn.all_comp.index(el_comp)]
            evolution.append(
                {
                    "chempot": c,
                    "evolution": amt,
                    "element_reference": el_ref,
                    "reaction": rxn,
                    "entries": decomp_entries,
                    "critical_composition": cc,
                }
            )
        return evolution

    def get_chempot_range_map(
        self, elements: Sequence[Element], referenced: bool = True, joggle: bool = True
    ) -> dict[Element, list[Simplex]]:
        """Get a chemical potential range map for each stable entry.

        Args:
            elements: Sequence of elements to be considered as independent variables.
                e.g. if you want to show the stability ranges
                of all Li-Co-O phases with respect to mu_Li and mu_O, you will supply
                [Element("Li"), Element("O")]
            referenced: If True, gives the results with a reference being the
                energy of the elemental phase. If False, gives absolute values.
            joggle (bool): Whether to joggle the input to avoid precision
                errors.

        Returns:
            Returns a dict of the form {entry: [simplices]}. The list of
            simplices are the sides of the N-1 dim polytope bounding the
            allowable chemical potential range of each entry.
        """
        all_chempots = []
        for facet in self.facets:
            chempots = self._get_facet_chempots(facet)
            all_chempots.append([chempots[el] for el in self.elements])

        inds = [self.elements.index(el) for el in elements]

        if referenced:
            el_energies = {el: self.el_refs[el].energy_per_atom for el in elements}
        else:
            el_energies = dict.fromkeys(elements, 0)

        chempot_ranges = defaultdict(list)
        vertices = [list(range(len(self.elements)))]

        if len(all_chempots) > len(self.elements):
            vertices = get_facets(all_chempots, joggle=joggle)

        for ufacet in vertices:
            for combi in itertools.combinations(ufacet, 2):
                data1 = self.facets[combi[0]]
                data2 = self.facets[combi[1]]
                common_ent_ind = set(data1).intersection(set(data2))
                if len(common_ent_ind) == len(elements):
                    common_entries = [self.qhull_entries[idx] for idx in common_ent_ind]
                    data = np.array(
                        [[all_chempots[ii][jj] - el_energies[self.elements[jj]] for jj in inds] for ii in combi]
                    )
                    sim = Simplex(data)
                    for entry in common_entries:
                        chempot_ranges[entry].append(sim)

        return chempot_ranges

    def getmu_vertices_stability_phase(self, target_comp, dep_elt, tol_en=1e-2):
        """Get a set of chemical potentials corresponding to the vertices of
        the simplex in the chemical potential phase diagram.
        The simplex is built using all elements in the target_composition
        except dep_elt.
        The chemical potential of dep_elt is computed from the target
        composition energy.
        This method is useful to get the limiting conditions for
        defects computations for instance.

        Args:
            target_comp: A Composition object
            dep_elt: the element for which the chemical potential is computed
                from the energy of the stable phase at the target composition
            tol_en: a tolerance on the energy to set

        Returns:
            [{Element: mu}]: An array of conditions on simplex vertices for
            which each element has a chemical potential set to a given
            value. "absolute" values (i.e., not referenced to element energies)
        """
        mu_ref = np.array([self.el_refs[elem].energy_per_atom for elem in self.elements if elem != dep_elt])
        chempot_ranges = self.get_chempot_range_map([elem for elem in self.elements if elem != dep_elt])

        for elem in self.elements:
            if elem not in target_comp.elements:
                target_comp += Composition({elem: 0.0})

        coeff = [-target_comp[elem] for elem in self.elements if elem != dep_elt]

        for elem, chempots in chempot_ranges.items():
            if elem.composition.reduced_composition == target_comp.reduced_composition:
                multiplier = elem.composition[dep_elt] / target_comp[dep_elt]
                ef = elem.energy / multiplier
                all_coords = []
                for simplex in chempots:
                    for v in simplex._coords:
                        elements = [elem for elem in self.elements if elem != dep_elt]
                        res = {}
                        for idx, el in enumerate(elements):
                            res[el] = v[idx] + mu_ref[idx]
                        res[dep_elt] = (np.dot(v + mu_ref, coeff) + ef) / target_comp[dep_elt]
                        already_in = False
                        for di in all_coords:
                            dict_equals = True
                            for k in di:
                                if abs(di[k] - res[k]) > tol_en:
                                    dict_equals = False
                                    break
                            if dict_equals:
                                already_in = True
                                break
                        if not already_in:
                            all_coords.append(res)

                return all_coords
        return None

    def get_chempot_range_stability_phase(self, target_comp, open_elt):
        """Get a set of chemical potentials corresponding to the max and min
        chemical potential of the open element for a given composition. It is
        quite common to have for instance a ternary oxide (e.g., ABO3) for
        which you want to know what are the A and B chemical potential leading
        to the highest and lowest oxygen chemical potential (reducing and
        oxidizing conditions). This is useful for defect computations.

        Args:
            target_comp: A Composition object
            open_elt: Element that you want to constrain to be max or min

        Returns:
            dict[Element, (float, float)]: A dictionary of the form {Element: (min_mu, max_mu)}
            where min_mu and max_mu are the minimum and maximum chemical potentials
            for the given element (as "absolute" values, i.e. not referenced to 0).
        """
        mu_ref = np.array([self.el_refs[elem].energy_per_atom for elem in self.elements if elem != open_elt])
        chempot_ranges = self.get_chempot_range_map([elem for elem in self.elements if elem != open_elt])
        for elem in self.elements:
            if elem not in target_comp.elements:
                target_comp += Composition({elem: 0.0})

        coeff = [-target_comp[elem] for elem in self.elements if elem != open_elt]
        max_open = -float("inf")
        min_open = float("inf")
        max_mus = min_mus = None

        for elem, chempots in chempot_ranges.items():
            if elem.composition.reduced_composition == target_comp.reduced_composition:
                multiplier = elem.composition[open_elt] / target_comp[open_elt]
                ef = elem.energy / multiplier
                all_coords = []
                for s in chempots:
                    for v in s._coords:
                        all_coords.append(v)
                        test_open = (np.dot(v + mu_ref, coeff) + ef) / target_comp[open_elt]
                        if test_open > max_open:
                            max_open = test_open
                            max_mus = v
                        if test_open < min_open:
                            min_open = test_open
                            min_mus = v

        elems = [elem for elem in self.elements if elem != open_elt]
        res = {}

        for idx, el in enumerate(elems):
            res[el] = (min_mus[idx] + mu_ref[idx], max_mus[idx] + mu_ref[idx])

        res[open_elt] = (min_open, max_open)
        return res

    def get_plot(
        self,
        show_unstable: float = 0.2,
        backend: Literal["plotly", "matplotlib"] = "plotly",
        ternary_style: Literal["2d", "3d"] = "2d",
        label_stable: bool = True,
        label_unstable: bool = True,
        ordering: Sequence[str] | None = None,
        energy_colormap=None,
        process_attributes: bool = False,
        ax: plt.Axes = None,
        label_uncertainties: bool = False,
        fill: bool = True,
        **kwargs,
    ):
        """
        Convenient wrapper for PDPlotter. Initializes a PDPlotter object and calls
        get_plot() with provided combined arguments.

        Plotting is only supported for phase diagrams with <=4 elements (unary,
        binary, ternary, or quaternary systems).

        Args:
            show_unstable (float): Whether unstable (above the hull) phases will be
                plotted. If a number > 0 is entered, all phases with
                e_hull < show_unstable (eV/atom) will be shown.
            backend ("plotly" | "matplotlib"): Python package to use for plotting.
                Defaults to "plotly".
            ternary_style ("2d" | "3d"): Ternary phase diagrams are typically plotted in
                two-dimensions (2d), but can be plotted in three dimensions (3d) to visualize
                the depth of the hull. This argument only applies when backend="plotly".
                Defaults to "2d".
            label_stable: Whether to label stable compounds.
            label_unstable: Whether to label unstable compounds.
            ordering: Ordering of vertices (matplotlib backend only).
            energy_colormap: Colormap for coloring energy (matplotlib backend only).
            process_attributes: Whether to process the attributes (matplotlib
                backend only).
            ax: Existing Axes object if plotting multiple phase diagrams (matplotlib backend only).
            label_uncertainties: Whether to add error bars to the hull (plotly
                backend only). For binaries, this also shades the hull with the
                uncertainty window.
            fill: Whether to shade the hull. For ternary_2d and quaternary plots, this
                colors facets arbitrarily for visual clarity. For ternary_3d plots, this
                shades the hull by formation energy (plotly backend only).
            **kwargs (dict): Keyword args passed to PDPlotter.get_plot(). Can be used to customize markers
                etc. If not set, the default is { "markerfacecolor": "#4daf4a", "markersize": 10, "linewidth": 3 }
        """
        plotter = PDPlotter(self, show_unstable=show_unstable, backend=backend, ternary_style=ternary_style)
        return plotter.get_plot(
            label_stable=label_stable,
            label_unstable=label_unstable,
            ordering=ordering,
            energy_colormap=energy_colormap,
            process_attributes=process_attributes,
            ax=ax,
            label_uncertainties=label_uncertainties,
            fill=fill,
            **kwargs,
        )


@due.dcite(Doi("10.1021/cm702327g"), description="Phase Diagram from First Principles Calculations")
@due.dcite(
    Doi("10.1016/j.elecom.2010.01.010"),
    description="Thermal stabilities of delithiated olivine MPO4 (M=Fe, Mn) cathodes "
    "investigated using first principles calculations",
)
class GrandPotentialPhaseDiagram(PhaseDiagram):
    """
    A class representing a Grand potential phase diagram. Grand potential phase
    diagrams are essentially phase diagrams that are open to one or more
    components. To construct such phase diagrams, the relevant free energy is
    the grand potential, which can be written as the Legendre transform of the
    Gibbs free energy as follows.

    Grand potential = G - u_X N_X

    The algorithm is based on the work in the following papers:

    1. S. P. Ong, L. Wang, B. Kang, and G. Ceder, Li-Fe-P-O2 Phase Diagram from
       First Principles Calculations. Chem. Mater., 2008, 20(5), 1798-1807.
       doi:10.1021/cm702327g

    2. S. P. Ong, A. Jain, G. Hautier, B. Kang, G. Ceder, Thermal stabilities
       of delithiated olivine MPO4 (M=Fe, Mn) cathodes investigated using first
       principles calculations. Electrochem. Comm., 2010, 12(3), 427-430.
       doi:10.1016/j.elecom.2010.01.010
    """

    def __init__(self, entries, chempots, elements=None, *, computed_data=None):
        """Standard constructor for grand potential phase diagram.

        Args:
            entries ([PDEntry]): A list of PDEntry-like objects having an
                energy, energy_per_atom and composition.
            chempots ({Element: float}): Specify the chemical potentials
                of the open elements.
            elements ([Element]): Optional list of elements in the phase
                diagram. If set to None, the elements are determined from
                the entries themselves.
            computed_data (dict): A dict containing pre-computed data. This allows
                PhaseDiagram object to be reconstituted without performing the
                expensive convex hull computation. The dict is the output from the
                PhaseDiagram._compute() method and is stored in PhaseDiagram.computed_data
                when generated for the first time.
        """
        if elements is None:
            elements = {els for entry in entries for els in entry.elements}

        self.chempots = {get_el_sp(el): u for el, u in chempots.items()}
        elements = set(elements) - set(self.chempots)

        all_entries = [
            GrandPotPDEntry(entry, self.chempots) for entry in entries if len(elements.intersection(entry.elements)) > 0
        ]

        super().__init__(all_entries, elements, computed_data=None)

    def __repr__(self):
        chemsys = "-".join(el.symbol for el in self.elements)
        chempots = ", ".join(f"mu_{el} = {mu:.4f}" for el, mu in self.chempots.items())

        output = [
            f"{chemsys} GrandPotentialPhaseDiagram with {chempots = }",
            f"{len(self.stable_entries)} stable phases: ",
            ", ".join(entry.name for entry in self.stable_entries),
        ]
        return "".join(output)

    def as_dict(self):
        """Get MSONable dict representation of GrandPotentialPhaseDiagram."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "all_entries": [entry.as_dict() for entry in self.all_entries],
            "chempots": self.chempots,
            "elements": [entry.as_dict() for entry in self.elements],
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """
        Args:
            dct (dict): dictionary representation of GrandPotentialPhaseDiagram.

        Returns:
            GrandPotentialPhaseDiagram
        """
        entries = MontyDecoder().process_decoded(dct["all_entries"])
        elements = MontyDecoder().process_decoded(dct["elements"])
        return cls(entries, dct["chempots"], elements)


class CompoundPhaseDiagram(PhaseDiagram):
    """
    Generates phase diagrams from compounds as terminations instead of
    elements.
    """

    # Tolerance for determining if amount of a composition is positive.
    amount_tol = 1e-5

    def __init__(self, entries, terminal_compositions, normalize_terminal_compositions=True):
        """Initialize a CompoundPhaseDiagram.

        Args:
            entries ([PDEntry]): Sequence of input entries. For example,
               if you want a Li2O-P2O5 phase diagram, you might have all
               Li-P-O entries as an input.
            terminal_compositions (list[Composition]): Terminal compositions of
                phase space. In the Li2O-P2O5 example, these will be the
                Li2O and P2O5 compositions.
            normalize_terminal_compositions (bool): Whether to normalize the
                terminal compositions to a per atom basis. If normalized,
                the energy above hulls will be consistent
                for comparison across systems. Non-normalized terminals are
                more intuitive in terms of compositional breakdowns.
        """
        self.original_entries = entries
        self.terminal_compositions = terminal_compositions
        self.normalize_terminals = normalize_terminal_compositions
        p_entries, species_mapping = self.transform_entries(entries, terminal_compositions)
        self.species_mapping = species_mapping
        super().__init__(p_entries, elements=species_mapping.values())

    def transform_entries(self, entries, terminal_compositions):
        """
        Method to transform all entries to the composition coordinate in the
        terminal compositions. If the entry does not fall within the space
        defined by the terminal compositions, they are excluded. For example,
        Li3PO4 is mapped into a Li2O:1.5, P2O5:0.5 composition. The terminal
        compositions are represented by DummySpecies.

        Args:
            entries: Sequence of all input entries
            terminal_compositions: Terminal compositions of phase space.

        Returns:
            Sequence of TransformedPDEntries falling within the phase space.
        """
        new_entries = []
        if self.normalize_terminals:
            terminal_compositions = [c.fractional_composition for c in terminal_compositions]

        # Map terminal compositions to unique dummy species.
        sp_mapping = {}
        for idx, comp in enumerate(terminal_compositions):
            sp_mapping[comp] = DummySpecies("X" + chr(102 + idx))

        for entry in entries:
            if getattr(entry, "attribute", None) is None:
                entry.attribute = getattr(entry, "entry_id", None)

            try:
                transformed_entry = TransformedPDEntry(entry, sp_mapping)
                new_entries.append(transformed_entry)
            except ReactionError:
                # If the reaction can't be balanced, the entry does not fall
                # into the phase space. We ignore them.
                pass
            except TransformedPDEntryError:
                # If the reaction has negative amounts for reactants the
                # entry does not fall into the phase space.
                pass

        return new_entries, sp_mapping

    def as_dict(self):
        """Get MSONable dict representation of CompoundPhaseDiagram."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "original_entries": [entry.as_dict() for entry in self.original_entries],
            "terminal_compositions": [c.as_dict() for c in self.terminal_compositions],
            "normalize_terminal_compositions": self.normalize_terminals,
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """
        Args:
            dct (dict): dictionary representation of CompoundPhaseDiagram.

        Returns:
            CompoundPhaseDiagram
        """
        entries = MontyDecoder().process_decoded(dct["original_entries"])
        terminal_compositions = MontyDecoder().process_decoded(dct["terminal_compositions"])
        return cls(entries, terminal_compositions, dct["normalize_terminal_compositions"])


class PatchedPhaseDiagram(PhaseDiagram):
    """
    Computing the Convex Hull of a large set of data in multiple dimensions is
    highly expensive. This class acts to breakdown large chemical spaces into
    smaller chemical spaces which can be computed much more quickly due to having
    both reduced dimensionality and data set sizes.

    Attributes:
        subspaces ({str: {Element, }}): Dictionary of the sets of elements for each of the
            PhaseDiagrams within the PatchedPhaseDiagram.
        pds ({str: PhaseDiagram}): Dictionary of PhaseDiagrams within the
            PatchedPhaseDiagram.
        all_entries (list[PDEntry]): All entries provided for Phase Diagram construction.
            Note that this does not mean that all these entries are actually used in
            the phase diagram. For example, this includes the positive formation energy
            entries that are filtered out before Phase Diagram construction.
        min_entries (list[PDEntry]): List of the lowest energy entries for each composition
            in the data provided for Phase Diagram construction.
        el_refs (list[PDEntry]): List of elemental references for the phase diagrams.
            These are entries corresponding to the lowest energy element entries for
            simple compositional phase diagrams.
        elements (list[Element]): List of elements in the phase diagram.
    """

    def __init__(
        self,
        entries: Sequence[PDEntry] | set[PDEntry],
        elements: Sequence[Element] | None = None,
        keep_all_spaces: bool = False,
        verbose: bool = False,
    ) -> None:
        """
        Args:
            entries (list[PDEntry]): A list of PDEntry-like objects having an
                energy, energy_per_atom and composition.
            elements (list[Element], optional): Optional list of elements in the phase
                diagram. If set to None, the elements are determined from
                the entries themselves and are sorted alphabetically.
                If specified, element ordering (e.g. for pd coordinates)
                is preserved.
            keep_all_spaces (bool): Pass True to keep chemical spaces that are subspaces
                of other spaces.
            verbose (bool): Whether to show progress bar during convex hull construction.
        """
        if elements is None:
            elements = sorted({els for entry in entries for els in entry.elements})

        self.dim = len(elements)

        entries = sorted(entries, key=lambda e: e.composition.reduced_composition)

        el_refs: dict[Element, PDEntry] = {}
        min_entries = []
        all_entries: list[PDEntry] = []
        for composition, group_iter in itertools.groupby(entries, key=lambda e: e.composition.reduced_composition):
            group = list(group_iter)
            min_entry = min(group, key=lambda e: e.energy_per_atom)
            if composition.is_element:
                el_refs[composition.elements[0]] = min_entry
            min_entries.append(min_entry)
            all_entries.extend(group)

        if len(el_refs) < self.dim:
            missing = set(elements) - set(el_refs)
            raise ValueError(f"Missing terminal entries for elements {sorted(map(str, missing))}")
        if len(el_refs) > self.dim:
            extra = set(el_refs) - set(elements)
            raise ValueError(f"There are more terminal elements than dimensions: {extra}")

        data = np.array(
            [
                [*(entry.composition.get_atomic_fraction(el) for el in elements), entry.energy_per_atom]
                for entry in min_entries
            ]
        )

        # Use only entries with negative formation energy
        vec = [el_refs[el].energy_per_atom for el in elements] + [-1]
        form_e = -np.dot(data, vec)
        inds = np.where(form_e < -PhaseDiagram.formation_energy_tol)[0].tolist()

        # Add the elemental references
        inds.extend([min_entries.index(el) for el in el_refs.values()])

        qhull_entries = tuple(min_entries[idx] for idx in inds)
        # make qhull spaces frozensets since they become keys to self.pds dict and frozensets are hashable
        # prevent repeating elements in chemical space and avoid the ordering problem (i.e. Fe-O == O-Fe automatically)
        qhull_spaces = tuple(frozenset(entry.elements) for entry in qhull_entries)

        # Get all unique chemical spaces
        spaces = {s for s in qhull_spaces if len(s) > 1}

        # Remove redundant chemical spaces
        spaces = self.remove_redundant_spaces(spaces, keep_all_spaces)

        # TODO comprhys: refactor to have self._compute method to allow serialization
        self.spaces = sorted(spaces, key=len, reverse=True)  # Calculate pds for smaller dimension spaces last
        self.qhull_entries = qhull_entries
        self._qhull_spaces = qhull_spaces
        self.pds = dict(self._get_pd_patch_for_space(s) for s in tqdm(self.spaces, disable=not verbose))
        self.all_entries = all_entries
        self.el_refs = el_refs
        self.elements = elements

        # Add terminal elements as we may not have PD patches including them
        # NOTE add el_refs in case no multielement entries are present for el
        _stable_entries = {se for pd in self.pds.values() for se in pd._stable_entries}
        self._stable_entries = tuple(_stable_entries | {*self.el_refs.values()})
        self._stable_spaces = tuple(frozenset(entry.elements) for entry in self._stable_entries)

    def __repr__(self):
        return f"{type(self).__name__} covering {len(self.spaces)} sub-spaces"

    def __len__(self):
        return len(self.spaces)

    def __getitem__(self, item: frozenset[Element]) -> PhaseDiagram:
        return self.pds[item]

    def __setitem__(self, key: frozenset[Element], value: PhaseDiagram) -> None:
        self.pds[key] = value

    def __delitem__(self, key: frozenset[Element]) -> None:
        del self.pds[key]

    def __iter__(self) -> Iterator[PhaseDiagram]:
        return iter(self.pds.values())

    def __contains__(self, item: frozenset[Element]) -> bool:
        return item in self.pds

    def as_dict(self) -> dict[str, Any]:
        """Write the entries and elements used to construct the PatchedPhaseDiagram
        to a dictionary.

        NOTE unlike PhaseDiagram the computation involved in constructing the
        PatchedPhaseDiagram is not saved on serialization. This is done because
        hierarchically calling the `PhaseDiagram.as_dict()` method would break the
        link in memory between entries in overlapping patches leading to a
        ballooning of the amount of memory used.

        NOTE For memory efficiency the best way to store patched phase diagrams is
        via pickling. As this allows all the entries in overlapping patches to share
        the same id in memory when unpickling.

        Returns:
            dict[str, Any]: MSONable dictionary representation of PatchedPhaseDiagram.
        """
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "all_entries": [entry.as_dict() for entry in self.all_entries],
            "elements": [entry.as_dict() for entry in self.elements],
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Reconstruct PatchedPhaseDiagram from dictionary serialization.

        NOTE unlike PhaseDiagram the computation involved in constructing the
        PatchedPhaseDiagram is not saved on serialization. This is done because
        hierarchically calling the `PhaseDiagram.as_dict()` method would break the
        link in memory between entries in overlapping patches leading to a
        ballooning of the amount of memory used.

        NOTE For memory efficiency the best way to store patched phase diagrams is
        via pickling. As this allows all the entries in overlapping patches to share
        the same id in memory when unpickling.

        Args:
            dct (dict): dictionary representation of PatchedPhaseDiagram.

        Returns:
            PatchedPhaseDiagram
        """
        entries = [MontyDecoder().process_decoded(entry) for entry in dct["all_entries"]]
        elements = [Element.from_dict(elem) for elem in dct["elements"]]
        return cls(entries, elements)

    @staticmethod
    def remove_redundant_spaces(spaces, keep_all_spaces=False):
        if keep_all_spaces or len(spaces) <= 1:
            return spaces

        # Sort spaces by size in descending order and pre-compute lengths
        sorted_spaces = sorted(spaces, key=len, reverse=True)

        result = []
        for idx, space_i in enumerate(sorted_spaces):
            if not any(space_i.issubset(larger_space) for larger_space in sorted_spaces[:idx]):
                result.append(space_i)

        return result

    # NOTE following methods are inherited unchanged from PhaseDiagram:
    # __repr__,
    # all_entries_hulldata,
    # unstable_entries,
    # stable_entries,
    # get_form_energy(),
    # get_form_energy_per_atom(),
    # get_hull_energy(),
    # get_e_above_hull(),
    # get_decomp_and_e_above_hull(),
    # get_decomp_and_phase_separation_energy(),
    # get_phase_separation_energy()

    def get_pd_for_entry(self, entry: Entry | Composition) -> PhaseDiagram:
        """Get the possible phase diagrams for an entry.

        Args:
            entry (PDEntry | Composition): A PDEntry or Composition-like object

        Returns:
            PhaseDiagram: phase diagram that the entry is part of

        Raises:
            ValueError: If no suitable PhaseDiagram is found for the entry.
        """
        entry_space = frozenset(entry.elements) if isinstance(entry, Composition) else frozenset(entry.elements)

        try:
            return self.pds[entry_space]
        except KeyError:
            for space, pd in self.pds.items():
                if space.issuperset(entry_space):
                    return pd

        raise ValueError(f"No suitable PhaseDiagrams found for {entry}.")

    def get_decomposition(self, comp: Composition) -> dict[PDEntry, float]:
        """See PhaseDiagram.

        Args:
            comp (Composition): A composition

        Returns:
            Decomposition as a dict of {PDEntry: amount} where amount
            is the amount of the fractional composition.
        """
        try:
            pd = self.get_pd_for_entry(comp)
            return pd.get_decomposition(comp)
        except ValueError as exc:
            # NOTE warn when stitching across pds is being used
            warnings.warn(f"{exc} Using SLSQP to find decomposition")
            competing_entries = self._get_stable_entries_in_space(frozenset(comp.elements))
            return _get_slsqp_decomp(comp, competing_entries)

    def get_equilibrium_reaction_energy(self, entry: Entry) -> float:
        """See PhaseDiagram.

        NOTE this is only approximately the same as the what we would get
        from `PhaseDiagram` as we make use of the slsqp approach inside
        get_phase_separation_energy().

        Args:
            entry (PDEntry): A PDEntry like object

        Returns:
            Equilibrium reaction energy of entry. Stable entries should have
            equilibrium reaction energy <= 0. The energy is given per atom.
        """
        return self.get_phase_separation_energy(entry, stable_only=True)

    def get_decomp_and_e_above_hull(
        self,
        entry: PDEntry,
        allow_negative: bool = False,
        check_stable: bool = False,
        on_error: Literal["raise", "warn", "ignore"] = "raise",
    ) -> tuple[dict[PDEntry, float], float] | tuple[None, None]:
        """Same as method on parent class PhaseDiagram except check_stable defaults to False
        for speed. See https://github.com/materialsproject/pymatgen/issues/2840 for details.
        """
        return super().get_decomp_and_e_above_hull(
            entry=entry, allow_negative=allow_negative, check_stable=check_stable, on_error=on_error
        )

    def _get_pd_patch_for_space(self, space: frozenset[Element]) -> tuple[frozenset[Element], PhaseDiagram]:
        """
        Args:
            space (frozenset[Element]): chemical space of the form A-B-X.

        Returns:
            space, PhaseDiagram for the given chemical space
        """
        space_entries = [e for e, s in zip(self.qhull_entries, self._qhull_spaces, strict=True) if space.issuperset(s)]

        return space, PhaseDiagram(space_entries)

    # NOTE the following functions are not implemented for PatchedPhaseDiagram

    def _get_facet_and_simplex(self):
        """Not Implemented - See PhaseDiagram."""
        raise NotImplementedError("_get_facet_and_simplex() not implemented for PatchedPhaseDiagram")

    def _get_all_facets_and_simplexes(self):
        """Not Implemented - See PhaseDiagram."""
        raise NotImplementedError("_get_all_facets_and_simplexes() not implemented for PatchedPhaseDiagram")

    def _get_facet_chempots(self):
        """Not Implemented - See PhaseDiagram."""
        raise NotImplementedError("_get_facet_chempots() not implemented for PatchedPhaseDiagram")

    def _get_simplex_intersections(self):
        """Not Implemented - See PhaseDiagram."""
        raise NotImplementedError("_get_simplex_intersections() not implemented for PatchedPhaseDiagram")

    def get_composition_chempots(self):
        """Not Implemented - See PhaseDiagram."""
        raise NotImplementedError("get_composition_chempots() not implemented for PatchedPhaseDiagram")

    def get_all_chempots(self):
        """Not Implemented - See PhaseDiagram."""
        raise NotImplementedError("get_all_chempots() not implemented for PatchedPhaseDiagram")

    def get_transition_chempots(self):
        """Not Implemented - See PhaseDiagram."""
        raise NotImplementedError("get_transition_chempots() not implemented for PatchedPhaseDiagram")

    def get_critical_compositions(self):
        """Not Implemented - See PhaseDiagram."""
        raise NotImplementedError("get_critical_compositions() not implemented for PatchedPhaseDiagram")

    def get_element_profile(self):
        """Not Implemented - See PhaseDiagram."""
        raise NotImplementedError("get_element_profile() not implemented for PatchedPhaseDiagram")

    def get_chempot_range_map(self):
        """Not Implemented - See PhaseDiagram."""
        raise NotImplementedError("get_chempot_range_map() not implemented for PatchedPhaseDiagram")

    def getmu_vertices_stability_phase(self):
        """Not Implemented - See PhaseDiagram."""
        raise NotImplementedError("getmu_vertices_stability_phase() not implemented for PatchedPhaseDiagram")

    def get_chempot_range_stability_phase(self):
        """Not Implemented - See PhaseDiagram."""
        raise NotImplementedError("get_chempot_range_stability_phase() not implemented for PatchedPhaseDiagram")


class ReactionDiagram:
    """
    Analyzes the possible reactions between a pair of compounds, e.g.
    an electrolyte and an electrode.
    """

    def __init__(self, entry1, entry2, all_entries, tol: float = 1e-4, float_fmt="%.4f"):
        """
        Args:
            entry1 (ComputedEntry): Entry for 1st component. Note that
                corrections, if any, must already be pre-applied. This is to
                give flexibility for different kinds of corrections, e.g.
                if a particular entry is fitted to an experimental data (such
                as EC molecule).
            entry2 (ComputedEntry): Entry for 2nd component. Note that
                corrections must already be pre-applied. This is to
                give flexibility for different kinds of corrections, e.g.
                if a particular entry is fitted to an experimental data (such
                as EC molecule).
            all_entries ([ComputedEntry]): All other entries to be
                considered in the analysis. Note that corrections, if any,
                must already be pre-applied.
            tol (float): Tolerance to be used to determine validity of reaction. Defaults to 1e-4.
            float_fmt (str): Formatting string to be applied to all floats. Determines
                number of decimal places in reaction string. Defaults to "%.4f".
        """
        elem_set = set()
        for entry in [entry1, entry2]:
            elem_set |= {el.symbol for el in entry.elements}

        elements = tuple(elem_set)  # Fix elements to ensure order.

        comp_vec1 = np.array([entry1.composition.get_atomic_fraction(el) for el in elements])
        comp_vec2 = np.array([entry2.composition.get_atomic_fraction(el) for el in elements])
        r1 = entry1.composition.reduced_composition
        r2 = entry2.composition.reduced_composition

        logger.debug(f"{len(all_entries)} total entries.")

        pd = PhaseDiagram([*all_entries, entry1, entry2])
        terminal_formulas = [entry1.reduced_formula, entry2.reduced_formula]

        logger.debug(f"{len(pd.stable_entries)} stable entries")
        logger.debug(f"{len(pd.facets)} facets")
        logger.debug(f"{len(pd.qhull_entries)} qhull_entries")

        rxn_entries = []
        done: list[tuple[float, float]] = []

        def fmt(fl):
            return float_fmt % fl

        for facet in pd.facets:
            for face in itertools.combinations(facet, len(facet) - 1):
                face_entries = [pd.qhull_entries[idx] for idx in face]

                if any(entry.reduced_formula in terminal_formulas for entry in face_entries):
                    continue

                try:
                    mat = []
                    for entry in face_entries:
                        mat.append([entry.composition.get_atomic_fraction(el) for el in elements])
                    mat.append(comp_vec2 - comp_vec1)
                    matrix = np.array(mat).T
                    coeffs = np.linalg.solve(matrix, comp_vec2)

                    x = coeffs[-1]

                    if all(c >= -tol for c in coeffs) and (abs(sum(coeffs[:-1]) - 1) < tol) and (tol < x < 1 - tol):
                        c1 = x / r1.num_atoms
                        c2 = (1 - x) / r2.num_atoms
                        factor = 1 / (c1 + c2)

                        c1 *= factor
                        c2 *= factor

                        # Avoid duplicate reactions.
                        if any(np.allclose([c1, c2], cc) for cc in done):
                            continue

                        done.append((c1, c2))

                        rxn_str = f"{fmt(c1)} {r1.reduced_formula} + {fmt(c2)} {r2.reduced_formula} -> "
                        products = []
                        product_entries = []

                        energy = -(x * entry1.energy_per_atom + (1 - x) * entry2.energy_per_atom)

                        for c, entry in zip(coeffs[:-1], face_entries, strict=True):
                            if c > tol:
                                redu_comp = entry.composition.reduced_composition
                                products.append(f"{fmt(c / redu_comp.num_atoms * factor)} {redu_comp.reduced_formula}")
                                product_entries.append((c, entry))
                                energy += c * entry.energy_per_atom

                        rxn_str += " + ".join(products)
                        comp = x * comp_vec1 + (1 - x) * comp_vec2
                        entry = PDEntry(
                            Composition(dict(zip(elements, comp, strict=True))),
                            energy=energy,
                            attribute=rxn_str,
                        )
                        entry.decomposition = product_entries
                        rxn_entries.append(entry)
                except np.linalg.LinAlgError:
                    form_1 = entry1.reduced_formula
                    form_2 = entry2.reduced_formula
                    logger.debug(f"Reactants = {form_1}, {form_2}")
                    logger.debug(f"Products = {', '.join([entry.reduced_formula for entry in face_entries])}")

        rxn_entries = sorted(rxn_entries, key=lambda e: e.name, reverse=True)

        self.entry1 = entry1
        self.entry2 = entry2
        self.rxn_entries = rxn_entries
        self.labels = {}
        for idx, entry in enumerate(rxn_entries, start=1):
            self.labels[str(idx)] = entry.attribute
            entry.name = str(idx)
        self.all_entries = all_entries
        self.pd = pd

    def get_compound_pd(self):
        """Get the CompoundPhaseDiagram object, which can then be used for
        plotting.

        Returns:
            CompoundPhaseDiagram
        """
        # For this plot, since the reactions are reported in formation
        # energies, we need to set the energies of the terminal compositions
        # to 0. So we make create copies with 0 energy.
        entry1 = PDEntry(self.entry1.composition, 0)
        entry2 = PDEntry(self.entry2.composition, 0)

        return CompoundPhaseDiagram(
            [*self.rxn_entries, entry1, entry2],
            [
                Composition(entry1.reduced_formula),
                Composition(entry2.reduced_formula),
            ],
            normalize_terminal_compositions=False,
        )


class PhaseDiagramError(Exception):
    """An exception class for Phase Diagram generation."""


def get_facets(qhull_data: ArrayLike, joggle: bool = False) -> ConvexHull:
    """Get the simplex facets for the Convex hull.

    Args:
        qhull_data (np.ndarray): The data from which to construct the convex
            hull as a Nxd array (N being number of data points and d being the
            dimension)
        joggle (bool): Whether to joggle the input to avoid precision
            errors.

    Returns:
        scipy.spatial.ConvexHull: with list of simplices of the convex hull.
    """
    if joggle:
        return ConvexHull(qhull_data, qhull_options="QJ i").simplices
    return ConvexHull(qhull_data, qhull_options="Qt i").simplices


def _get_slsqp_decomp(
    comp,
    competing_entries,
    tols=(1e-8,),
    maxiter=1000,
):
    """Find the amounts of competing compositions that minimize the energy of a
    given composition.

    The algorithm is based on the work in the following paper:

    1. Bartel, C., Trewartha, A., Wang, Q., Dunn, A., Jain, A., Ceder, G.,
        A critical examination of compound stability predictions from
        machine-learned formation energies, npj Computational Materials 6, 97 (2020)

    Args:
        comp (Composition): A Composition to analyze
        competing_entries ([PDEntry]): List of entries to consider for decomposition
        tols (list): tolerances to try for SLSQP convergence. Issues observed for
            tol > 1e-7 in the fractional composition (default 1e-8)
        maxiter (int): maximum number of SLSQP iterations

    Returns:
            decomposition as a dict of {PDEntry: amount} where amount
            is the amount of the fractional composition.
    """
    # Elemental amount present in given entry
    amounts = comp.get_el_amt_dict()
    chemical_space = tuple(amounts)
    b = np.array([amounts[el] for el in chemical_space])

    # Elemental amounts present in competing entries
    A_transpose = np.zeros((len(chemical_space), len(competing_entries)))
    for ii, comp_entry in enumerate(competing_entries):
        amounts = comp_entry.composition.get_el_amt_dict()
        for jj, el in enumerate(chemical_space):
            A_transpose[jj, ii] = amounts.get(el, 0)

    # NOTE normalize arrays to avoid calls to fractional_composition
    b /= np.sum(b)
    A_transpose /= np.sum(A_transpose, axis=0)

    # Energies of competing entries
    Es = np.array([comp_entry.energy_per_atom for comp_entry in competing_entries])

    molar_constraint = {"type": "eq", "fun": lambda x: np.dot(A_transpose, x) - b, "jac": lambda x: A_transpose}

    options = {"maxiter": maxiter, "disp": False}

    # NOTE max_bound needs to be larger than 1
    max_bound = comp.num_atoms
    bounds = [(0, max_bound)] * len(competing_entries)
    x0 = [1 / len(competing_entries)] * len(competing_entries)

    # NOTE the tolerance needs to be tight to stop the optimization
    # from exiting before convergence is reached. Issues observed for
    # tol > 1e-7 in the fractional composition (default 1e-8).
    for tol in sorted(tols):
        solution = minimize(
            fun=lambda x: np.dot(x, Es),
            x0=x0,
            method="SLSQP",
            jac=lambda x: Es,
            bounds=bounds,
            constraints=[molar_constraint],
            tol=tol,
            options=options,
        )

        if solution.success:
            decomp_amts = solution.x
            return {
                c: amt  # NOTE this is the amount of the fractional composition.
                for c, amt in zip(competing_entries, decomp_amts, strict=True)
                if amt > PhaseDiagram.numerical_tol
            }

    raise ValueError(f"No valid decomp found for {comp}!")


class PDPlotter:
    """
    A plotting class for compositional phase diagrams.

    To use, initialize this class with a PhaseDiagram object containing 1-4 components
    and call get_plot() or show().
    """

    def __init__(
        self,
        phasediagram: PhaseDiagram,
        show_unstable: float = 0.2,
        backend: Literal["plotly", "matplotlib"] = "plotly",
        ternary_style: Literal["2d", "3d"] = "2d",
        **plotkwargs,
    ):
        """
        Args:
            phasediagram (PhaseDiagram): PhaseDiagram object (must be 1-4 components).
            show_unstable (float): Whether unstable (above the hull) phases will be
                plotted. If a number > 0 is entered, all phases with
                e_hull < show_unstable (eV/atom) will be shown.
            backend ("plotly" | "matplotlib"): Python package to use for plotting.
                Defaults to "plotly".
            ternary_style ("2d" | "3d"): Ternary phase diagrams are typically plotted in
                two-dimensions (2d), but can be plotted in three dimensions (3d) to visualize
                the depth of the hull. This argument only applies when backend="plotly".
                Defaults to "2d".
            **plotkwargs (dict): Keyword args passed to matplotlib.pyplot.plot (only
                applies when backend="matplotlib"). Can be used to customize markers
                etc. If not set, the default is:
                    {
                        "markerfacecolor": "#4daf4a",
                        "markersize": 10,
                        "linewidth": 3
                    }.
        """
        dim = len(phasediagram.elements)
        if dim >= 5:
            raise ValueError("Only 1-4 components supported!")

        self._pd = phasediagram
        self.show_unstable = show_unstable
        self.backend = backend
        self.ternary_style = ternary_style.lower()

        self.lines = uniquelines(self._pd.facets) if dim > 1 else [[self._pd.facets[0][0], self._pd.facets[0][0]]]
        self._min_energy = min(self._pd.get_form_energy_per_atom(entry) for entry in self._pd.stable_entries)
        self._dim = dim

        self.plotkwargs = plotkwargs or {
            "markerfacecolor": "#4daf4a",
            "markersize": 10,
            "linewidth": 3,
        }

    def get_plot(
        self,
        label_stable: bool = True,
        label_unstable: bool = True,
        ordering: Sequence[str] | None = None,
        energy_colormap=None,
        process_attributes: bool = False,
        ax: plt.Axes = None,
        label_uncertainties: bool = False,
        fill: bool = True,
        highlight_entries: Collection[PDEntry] | None = None,
    ) -> go.Figure | plt.Axes:
        """
        Args:
            label_stable: Whether to label stable compounds.
            label_unstable: Whether to label unstable compounds.
            ordering: Ordering of vertices, given as a list ['Up',
                'Left','Right'] (matplotlib only).
            energy_colormap: Colormap for coloring energy (matplotlib only).
            process_attributes: Whether to process the attributes (matplotlib only).
            ax: Existing matplotlib Axes object if plotting multiple phase diagrams
                (matplotlib only).
            label_uncertainties: Whether to add error bars to the hull.
                For binaries, this also shades the hull with the uncertainty window.
                (plotly only).
            fill: Whether to shade the hull. For ternary_2d and quaternary plots, this
                colors facets arbitrarily for visual clarity. For ternary_3d plots, this
                shades the hull by formation energy (plotly only).
            highlight_entries: Entries to highlight in the plot (plotly only). This will
                create a new marker trace that is separate from the other entries.

        Returns:
            go.Figure | plt.Axes: Plotly figure or matplotlib axes object depending on backend.
        """
        fig = None
        data = []

        if self.backend == "plotly":
            if self._dim != 1:
                data.append(self._create_plotly_lines())

            stable_marker_plot, unstable_marker_plot, highlight_plot = self._create_plotly_markers(
                highlight_entries,
                label_uncertainties,
            )

            if self._dim == 2 and label_uncertainties:
                data.append(self._create_plotly_uncertainty_shading(stable_marker_plot))

            if self._dim == 3 and self.ternary_style == "3d":
                data.append(self._create_plotly_ternary_support_lines())

            if self._dim != 1 and not (self._dim == 3 and self.ternary_style == "2d"):
                data.append(self._create_plotly_stable_labels(label_stable))

            if fill and self._dim in [3, 4]:
                data.extend(self._create_plotly_fill())

            data.extend([stable_marker_plot, unstable_marker_plot])

            if highlight_plot is not None:
                data.append(highlight_plot)

            fig = go.Figure(data=data)
            fig.layout = self._create_plotly_figure_layout()
            fig.update_layout(coloraxis_colorbar={"yanchor": "top", "y": 0.05, "x": 1})

        elif self.backend == "matplotlib":
            if self._dim <= 3:
                fig = self._get_matplotlib_2d_plot(
                    label_stable,
                    label_unstable,
                    ordering,
                    energy_colormap,
                    ax=ax,
                    process_attributes=process_attributes,
                )
            elif self._dim == 4:
                fig = self._get_matplotlib_3d_plot(label_stable, ax=ax)

        return fig

    def show(self, *args, **kwargs) -> None:
        """
        Draw the phase diagram with the provided arguments and display it. This shows
        the figure but does not return it.

        Args:
            *args: Passed to get_plot.
            **kwargs: Passed to get_plot.
        """
        plot = self.get_plot(*args, **kwargs)
        if self.backend == "matplotlib":
            plot.get_figure().show()
        else:
            plot.show()

    def write_image(self, stream: str | StringIO, image_format: str = "svg", **kwargs) -> None:
        """
        Directly save the plot to a file. This is a wrapper for calling plt.savefig() or
        fig.write_image(), depending on the backend. For more customization, it is
        recommended to call those methods directly.

        Args:
            stream (str | StringIO): Filename or StringIO stream.
            image_format (str): Can be any supported image format for the plotting backend.
                Defaults to 'svg' (vector graphics).
            **kwargs: Optinoal kwargs passed to the get_plot function.
        """
        if self.backend == "matplotlib":
            ax = self.get_plot(**kwargs)
            ax.figure.set_size_inches((12, 10))
            ax.figure.savefig(stream, format=image_format)
        elif self.backend == "plotly":
            fig = self.get_plot(**kwargs)
            fig.write_image(stream, format=image_format)

    def plot_element_profile(self, element, comp, show_label_index=None, xlim=5):
        """
        Draw the element profile plot for a composition varying different
        chemical potential of an element.

        X value is the negative value of the chemical potential reference to
        elemental chemical potential. For example, if choose Element("Li"),
        X= -(μLi-μLi0), which corresponds to the voltage versus metal anode.
        Y values represent for the number of element uptake in this composition
        (unit: per atom). All reactions are printed to help choosing the
        profile steps you want to show label in the plot.

        Args:
            element (Element): An element of which the chemical potential is
                considered. It also must be in the phase diagram.
            comp (Composition): A composition.
            show_label_index (list of integers): The labels for reaction products
                you want to show in the plot. Default to None (not showing any
                annotation for reaction products). For the profile steps you want
                to show the labels, just add it to the show_label_index. The
                profile step counts from zero. For example, you can set
                show_label_index=[0, 2, 5] to label profile step 0,2,5.
            xlim (float): The max x value. x value is from 0 to xlim. Default to
                5 eV.

        Returns:
            Plot of element profile evolution by varying the chemical potential
            of an element.
        """
        ax = pretty_plot(12, 8)
        pd = self._pd
        evolution = pd.get_element_profile(element, comp)
        num_atoms = evolution[0]["reaction"].reactants[0].num_atoms
        element_energy = evolution[0]["chempot"]
        x1, x2, y1 = None, None, None
        for idx, dct in enumerate(evolution):
            v = -(dct["chempot"] - element_energy)
            if idx != 0:
                ax.plot([x2, x2], [y1, dct["evolution"] / num_atoms], "k", linewidth=2.5)
            x1 = v
            y1 = dct["evolution"] / num_atoms

            x2 = -(evolution[idx + 1]["chempot"] - element_energy) if idx != len(evolution) - 1 else 5.0
            if show_label_index is not None and idx in show_label_index:
                products = [
                    re.sub(r"(\d+)", r"$_{\1}$", p.reduced_formula)
                    for p in dct["reaction"].products
                    if p.reduced_formula != element.symbol
                ]
                ax.annotate(
                    ", ".join(products),
                    xy=(v + 0.05, y1 + 0.05),
                    fontsize=24,
                    color="r",
                )
                ax.plot([x1, x2], [y1, y1], "r", linewidth=3)
            else:
                ax.plot([x1, x2], [y1, y1], "k", linewidth=2.5)

        ax.set_xlim((0, xlim))
        ax.set_xlabel("-$\\Delta{\\mu}$ (eV)")
        ax.set_ylabel("Uptake per atom")

        return ax

    def plot_chempot_range_map(self, elements, referenced=True) -> None:
        """
        Plot the chemical potential range _map using matplotlib. Currently works only for
        3-component PDs. This shows the plot but does not return it.

        Note: this functionality is now included in the ChemicalPotentialDiagram
        class (pymatgen.analysis.chempot_diagram).

        Args:
            elements: Sequence of elements to be considered as independent
                variables. e.g. if you want to show the stability ranges of
                all Li-Co-O phases w.r.t. to uLi and uO, you will supply
                [Element("Li"), Element("O")]
            referenced: if True, gives the results with a reference being the
                        energy of the elemental phase. If False, gives absolute values.
        """
        self.get_chempot_range_map_plot(elements, referenced=referenced).show()

    def get_chempot_range_map_plot(self, elements, referenced=True):
        """Get a plot of the chemical potential range _map. Currently works
        only for 3-component PDs.

        Note: this functionality is now included in the ChemicalPotentialDiagram
        class (pymatgen.analysis.chempot_diagram).

        Args:
            elements: Sequence of elements to be considered as independent
                variables. e.g. if you want to show the stability ranges of
                all Li-Co-O phases w.r.t. to uLi and uO, you will supply
                [Element("Li"), Element("O")]
            referenced: if True, gives the results with a reference being the
                energy of the elemental phase. If False, gives absolute values.

        Returns:
            plt.Axes: matplotlib axes object.
        """
        ax = pretty_plot(12, 8)
        chempot_ranges = self._pd.get_chempot_range_map(elements, referenced=referenced)
        missing_lines = {}
        excluded_region = []

        for entry, lines in chempot_ranges.items():
            comp = entry.composition
            center_x = center_y = 0
            coords = []
            contain_zero = any(comp.get_atomic_fraction(el) == 0 for el in elements)
            is_boundary = (not contain_zero) and sum(comp.get_atomic_fraction(el) for el in elements) == 1
            for line in lines:
                x, y = line.coords.transpose()
                plt.plot(x, y, "k-")

                for coord in line.coords:
                    if not in_coord_list(coords, coord):
                        coords.append(coord.tolist())
                        center_x += coord[0]
                        center_y += coord[1]
                if is_boundary:
                    excluded_region.extend(line.coords)

            if coords and contain_zero:
                missing_lines[entry] = coords
            else:
                xy = (center_x / len(coords), center_y / len(coords))
                plt.annotate(latexify(entry.name), xy, fontsize=22)

        xlim = ax.get_xlim()
        ylim = ax.get_ylim()

        # Shade the forbidden chemical potential regions.
        excluded_region.append([xlim[1], ylim[1]])
        excluded_region = sorted(excluded_region, key=lambda c: c[0])
        x, y = np.transpose(excluded_region)
        plt.fill(x, y, "0.80")

        # The hull does not generate the missing horizontal and vertical lines.
        # The following code fixes this.
        el0 = elements[0]
        el1 = elements[1]

        for entry, coords in missing_lines.items():
            center_x = sum(c[0] for c in coords)
            center_y = sum(c[1] for c in coords)
            comp = entry.composition
            is_x = comp.get_atomic_fraction(el0) < 0.01
            is_y = comp.get_atomic_fraction(el1) < 0.01
            n_coords = len(coords)
            if not (is_x and is_y):
                if is_x:
                    coords = sorted(coords, key=lambda c: c[1])
                    for idx in [0, -1]:
                        x = [min(xlim), coords[idx][0]]
                        y = [coords[idx][1], coords[idx][1]]
                        plt.plot(x, y, "k")
                        center_x += min(xlim)
                        center_y += coords[idx][1]
                elif is_y:
                    coords = sorted(coords, key=lambda c: c[0])
                    for idx in [0, -1]:
                        x = [coords[idx][0], coords[idx][0]]
                        y = [coords[idx][1], min(ylim)]
                        plt.plot(x, y, "k")
                        center_x += coords[idx][0]
                        center_y += min(ylim)
                xy = (center_x / (n_coords + 2), center_y / (n_coords + 2))
            else:
                center_x = sum(coord[0] for coord in coords) + xlim[0]
                center_y = sum(coord[1] for coord in coords) + ylim[0]
                xy = (center_x / (n_coords + 1), center_y / (n_coords + 1))

            ax.annotate(
                latexify(entry.name),
                xy,
                horizontalalignment="center",
                verticalalignment="center",
                fontsize=22,
            )

        ax.set_xlabel(f"$\\mu_{{{el0.symbol}}} - \\mu_{{{el0.symbol}}}^0$ (eV)")
        ax.set_ylabel(f"$\\mu_{{{el1.symbol}}} - \\mu_{{{el1.symbol}}}^0$ (eV)")
        plt.tight_layout()
        return ax

    def get_contour_pd_plot(self):
        """
        Plot a contour phase diagram plot, where phase triangles are colored
        according to degree of instability by interpolation. Currently only
        works for 3-component phase diagrams.

        Returns:
            A matplotlib plot object.
        """
        pd = self._pd
        entries = pd.qhull_entries
        data = np.array(pd.qhull_data)

        ax = self._get_matplotlib_2d_plot()
        data[:, 0:2] = triangular_coord(data[:, 0:2]).transpose()
        for idx, entry in enumerate(entries):
            data[idx, 2] = self._pd.get_e_above_hull(entry)

        gridsize = 0.005
        xnew = np.arange(0, 1.0, gridsize)
        ynew = np.arange(0, 1, gridsize)

        f = interpolate.LinearNDInterpolator(data[:, 0:2], data[:, 2])
        znew = np.zeros((len(ynew), len(xnew)))
        for idx, xval in enumerate(xnew):
            for j, yval in enumerate(ynew):
                znew[j, idx] = f(xval, yval)

        contourf = ax.contourf(xnew, ynew, znew, 1000, cmap=cm.autumn_r)

        plt.colorbar(contourf)

        return ax

    @property
    @lru_cache(1)  # noqa: B019
    def pd_plot_data(self):
        """
        Plotting data for phase diagram. Cached for repetitive calls.

        2-comp - Full hull with energies
        3/4-comp - Projection into 2D or 3D Gibbs triangles

        Returns:
            A tuple containing three objects (lines, stable_entries, unstable_entries):
            - lines is a list of list of coordinates for lines in the PD.
            - stable_entries is a dict of {coordinates : entry} for each stable node
                in the phase diagram. (Each coordinate can only have one
                stable phase)
            - unstable_entries is a dict of {entry: coordinates} for all unstable
                nodes in the phase diagram.
        """
        pd = self._pd
        entries = pd.qhull_entries
        data = np.array(pd.qhull_data)
        lines = []
        stable_entries = {}

        for line in self.lines:
            entry1 = entries[line[0]]
            entry2 = entries[line[1]]
            if self._dim < 3:
                x = [data[line[0]][0], data[line[1]][0]]
                y = [
                    pd.get_form_energy_per_atom(entry1),
                    pd.get_form_energy_per_atom(entry2),
                ]
                coord = [x, y]
            elif self._dim == 3:
                coord = triangular_coord(data[line, 0:2])
            else:
                coord = tet_coord(data[line, 0:3])
            lines.append(coord)
            label_coord = list(zip(*coord, strict=True))
            stable_entries[label_coord[0]] = entry1
            stable_entries[label_coord[1]] = entry2

        all_entries = pd.all_entries
        all_data = np.array(pd.all_entries_hulldata)
        unstable_entries = {}
        stable = pd.stable_entries

        for idx, entry in enumerate(all_entries):
            if entry not in stable:
                if self._dim < 3:
                    x = [all_data[idx][0], all_data[idx][0]]
                    y = [pd.get_form_energy_per_atom(entry), pd.get_form_energy_per_atom(entry)]
                    coord = [x, y]
                elif self._dim == 3:
                    coord = triangular_coord([all_data[idx, 0:2], all_data[idx, 0:2]])
                else:
                    coord = tet_coord([all_data[idx, 0:3], all_data[idx, 0:3], all_data[idx, 0:3]])
                label_coord = list(zip(*coord, strict=True))
                unstable_entries[entry] = label_coord[0]

        return lines, stable_entries, unstable_entries

    def _create_plotly_figure_layout(self, label_stable=True):
        """
        Creates layout for plotly phase diagram figure and updates with
        figure annotations.

        Args:
            label_stable (bool): Whether to label stable compounds

        Returns:
            Dictionary with Plotly figure layout settings.
        """
        annotations = None
        layout = {}

        if label_stable:
            annotations = self._create_plotly_element_annotations()

        if self._dim == 1:
            layout = plotly_layouts["default_unary_layout"].copy()
        if self._dim == 2:
            layout = plotly_layouts["default_binary_layout"].copy()
            layout["xaxis"]["title"] = f"Composition (Fraction {self._pd.elements[1]})"
            layout["annotations"] = annotations
        elif self._dim == 3 and self.ternary_style == "2d":
            layout = plotly_layouts["default_ternary_2d_layout"].copy()
            for el, axis in zip(self._pd.elements, ["a", "b", "c"], strict=True):
                el_ref = self._pd.el_refs[el]
                clean_formula = str(el_ref.elements[0])
                if hasattr(el_ref, "original_entry"):  # for grand potential PDs, etc.
                    clean_formula = htmlify(el_ref.original_entry.reduced_formula)

                layout["ternary"][f"{axis}axis"]["title"] = {
                    "text": clean_formula,
                    "font": {"size": 24},
                }
        elif self._dim == 3 and self.ternary_style == "3d":
            layout = plotly_layouts["default_ternary_3d_layout"].copy()
            layout["scene"]["annotations"] = annotations
        elif self._dim == 4:
            layout = plotly_layouts["default_quaternary_layout"].copy()
            layout["scene"]["annotations"] = annotations

        return layout

    def _create_plotly_lines(self):
        """
        Create Plotly scatter plots containing line traces of phase diagram facets.

        Returns:
            Either a go.Scatter (binary), go.Scatterternary (ternary_2d), or
            go.Scatter3d plot (ternary_3d, quaternary)
        """
        line_plot = None
        x, y, z, energies = [], [], [], []

        pd = self._pd

        plot_args = {
            "mode": "lines",
            "hoverinfo": "none",
            "line": {"color": "black", "width": 4.0},
            "showlegend": False,
        }

        if self._dim == 3 and self.ternary_style == "2d":
            plot_args["line"]["width"] = 1.5
            el_a, el_b, el_c = pd.elements
            for line in uniquelines(pd.facets):
                e0 = pd.qhull_entries[line[0]]
                e1 = pd.qhull_entries[line[1]]

                x += [e0.composition[el_a], e1.composition[el_a], None]
                y += [e0.composition[el_b], e1.composition[el_b], None]
                z += [e0.composition[el_c], e1.composition[el_c], None]
        else:
            for line in self.pd_plot_data[0]:
                x += [*line[0], None]
                y += [*line[1], None]

                if self._dim == 3:
                    form_enes = [
                        self._pd.get_form_energy_per_atom(self.pd_plot_data[1][coord])
                        for coord in zip(line[0], line[1], strict=True)
                    ]
                    z += [*form_enes, None]

                elif self._dim == 4:
                    form_enes = [
                        self._pd.get_form_energy_per_atom(self.pd_plot_data[1][coord])
                        for coord in zip(line[0], line[1], line[2], strict=True)
                    ]
                    energies += [*form_enes, None]
                    z += [*line[2], None]

        if self._dim == 2:
            line_plot = go.Scatter(x=x, y=y, **plot_args)
        elif self._dim == 3 and self.ternary_style == "2d":
            line_plot = go.Scatterternary(a=x, b=y, c=z, **plot_args)
        elif self._dim == 3 and self.ternary_style == "3d":
            line_plot = go.Scatter3d(x=y, y=x, z=z, **plot_args)
        elif self._dim == 4:
            plot_args["line"]["width"] = 1.5
            line_plot = go.Scatter3d(x=x, y=y, z=z, **plot_args)

        return line_plot

    def _create_plotly_fill(self):
        """
        Creates shaded mesh traces for coloring the hull.

        For tenrary_3d plots, the color shading is based on formation energy.

        Returns:
            go.Mesh3d plot
        """
        traces = []

        pd = self._pd
        if self._dim == 3 and self.ternary_style == "2d":
            fillcolors = itertools.cycle(plotly_layouts["default_fill_colors"])
            el_a, el_b, el_c = pd.elements

            for _idx, facet in enumerate(pd.facets):
                a = []
                b = []
                c = []

                e0, e1, e2 = sorted((pd.qhull_entries[facet[idx]] for idx in range(3)), key=lambda x: x.reduced_formula)
                a = [e0.composition[el_a], e1.composition[el_a], e2.composition[el_a]]
                b = [e0.composition[el_b], e1.composition[el_b], e2.composition[el_b]]
                c = [e0.composition[el_c], e1.composition[el_c], e2.composition[el_c]]

                e_strs = []
                for entry in (e0, e1, e2):
                    if hasattr(entry, "original_entry"):
                        entry = entry.original_entry
                    e_strs.append(htmlify(entry.reduced_formula))

                name = f"{e_strs[0]}—{e_strs[1]}—{e_strs[2]}"

                traces += [
                    go.Scatterternary(
                        a=a,
                        b=b,
                        c=c,
                        mode="lines",
                        fill="toself",
                        line={"width": 0},
                        fillcolor=next(fillcolors),
                        opacity=0.15,
                        hovertemplate="<extra></extra>",  # removes secondary hover box
                        name=name,
                        showlegend=False,
                    )
                ]
        elif self._dim == 3 and self.ternary_style == "3d":
            facets = np.array(self._pd.facets)
            coords = np.array(
                [
                    triangular_coord(c)
                    for c in zip(self._pd.qhull_data[:-1, 0], self._pd.qhull_data[:-1, 1], strict=True)
                ]
            )
            energies = np.array([self._pd.get_form_energy_per_atom(entry) for entry in self._pd.qhull_entries])

            traces.append(
                go.Mesh3d(
                    x=list(coords[:, 1]),
                    y=list(coords[:, 0]),
                    z=list(energies),
                    i=list(facets[:, 1]),
                    j=list(facets[:, 0]),
                    k=list(facets[:, 2]),
                    opacity=0.7,
                    intensity=list(energies),
                    colorscale=plotly_layouts["stable_colorscale"],
                    colorbar={
                        "title": "Formation energy<br>(eV/atom)",
                        "x": 0.9,
                        "y": 1,
                        "yanchor": "top",
                        "xpad": 0,
                        "ypad": 0,
                        "thickness": 0.02,
                        "thicknessmode": "fraction",
                        "len": 0.5,
                    },
                    hoverinfo="none",
                    lighting={"diffuse": 0.0, "ambient": 1.0},
                    name="Convex Hull (shading)",
                    flatshading=True,
                    showlegend=True,
                )
            )
        elif self._dim == 4:
            all_data = np.array(pd.qhull_data)
            fillcolors = itertools.cycle(plotly_layouts["default_fill_colors"])
            for _idx, facet in enumerate(pd.facets):
                xs, ys, zs = [], [], []
                for v in facet:
                    x, y, z = tet_coord(all_data[v, 0:3])
                    xs.append(x)
                    ys.append(y)
                    zs.append(z)

                if _idx == 1:
                    traces += [
                        go.Mesh3d(
                            x=xs,
                            y=ys,
                            z=zs,
                            opacity=0.05,
                            alphahull=-1,
                            flatshading=True,
                            hoverinfo="skip",
                            color=next(fillcolors),
                            legendgroup="facets",
                            showlegend=True,
                            name="Hull Surfaces (toggle to access points easier)",
                        )
                    ]
                else:
                    traces += [
                        go.Mesh3d(
                            x=xs,
                            y=ys,
                            z=zs,
                            opacity=0.05,
                            alphahull=-1,
                            flatshading=True,
                            hoverinfo="skip",
                            color=next(fillcolors),
                            legendgroup="facets",
                        )
                    ]

        return traces

    def _create_plotly_stable_labels(self, label_stable=True):
        """
        Creates a (hidable) scatter trace containing labels of stable phases.
        Contains some functionality for creating sensible label positions. This method
        does not apply to 2D ternary plots (stable labels are turned off).

        Returns:
            go.Scatter (or go.Scatter3d) plot
        """
        x, y, z, text, textpositions = [], [], [], [], []
        stable_labels_plot = min_energy_x = None
        offset_2d = 0.008  # extra distance to offset label position for clarity
        offset_3d = 0.01

        energy_offset = -0.05 * self._min_energy  # 5% above points

        if self._dim == 2:
            min_energy_x = min(list(self.pd_plot_data[1]), key=lambda c: c[1])[0]

        for coords, entry in self.pd_plot_data[1].items():
            if entry.composition.is_element:  # taken care of by other function
                continue
            x_coord = coords[0]
            y_coord = coords[1]
            textposition = None

            if self._dim == 2:
                textposition = "bottom left"
                if x_coord >= min_energy_x:
                    textposition = "bottom right"
                    x_coord += offset_2d
                else:
                    x_coord -= offset_2d
                y_coord -= offset_2d + 0.005
            elif self._dim == 3 and self.ternary_style == "3d":
                textposition = "middle center"
                if coords[0] > 0.5:  # right half of plot
                    x_coord += offset_3d
                else:
                    x_coord -= offset_3d
                if coords[1] > 0.866 / 2:  # top half of plot (highest point is 0.866)
                    y_coord -= offset_3d
                else:
                    y_coord += offset_3d

                z.append(self._pd.get_form_energy_per_atom(entry) + energy_offset)
            elif self._dim == 4:
                x_coord -= offset_3d
                y_coord -= offset_3d
                textposition = "bottom right"
                z.append(coords[2])

            x.append(x_coord)
            y.append(y_coord)
            textpositions.append(textposition)

            comp = entry.composition
            if hasattr(entry, "original_entry"):
                comp = entry.original_entry.composition

            formula = comp.reduced_formula
            text.append(htmlify(formula))

        visible = True
        if not label_stable or self._dim == 4:
            visible = "legendonly"

        plot_args = {
            "text": text,
            "textposition": textpositions,
            "mode": "text",
            "name": "Labels (stable)",
            "hoverinfo": "skip",
            "opacity": 1.0,
            "visible": visible,
            "showlegend": True,
        }

        if self._dim == 2:
            stable_labels_plot = go.Scatter(x=x, y=y, **plot_args)
        elif self._dim == 3 and self.ternary_style == "3d":
            stable_labels_plot = go.Scatter3d(x=y, y=x, z=z, **plot_args)
        elif self._dim == 4:
            stable_labels_plot = go.Scatter3d(x=x, y=y, z=z, **plot_args)

        return stable_labels_plot

    def _create_plotly_element_annotations(self):
        """
        Creates terminal element annotations for Plotly phase diagrams. This method does
        not apply to ternary_2d plots.

        Functionality is included for phase diagrams with non-elemental endmembers
        (as is true for grand potential phase diagrams).

        Returns:
            List of annotation dicts.
        """
        annotations_list = []
        x, y, z = None, None, None

        if self._dim == 3 and self.ternary_style == "2d":
            return None

        for coords, entry in self.pd_plot_data[1].items():
            if not entry.composition.is_element:
                continue

            x, y = coords[0], coords[1]

            if self._dim == 3:
                z = self._pd.get_form_energy_per_atom(entry)
            elif self._dim == 4:
                z = coords[2]

            if entry.composition.is_element:
                clean_formula = str(entry.elements[0])
                if hasattr(entry, "original_entry"):
                    orig_comp = entry.original_entry.composition
                    clean_formula = htmlify(orig_comp.reduced_formula)

                font_dict = {"color": "#000000", "size": 24.0}
                opacity = 1.0

            else:
                clean_formula = ""
                font_dict = {}
                opacity = 0

            offset = 0.03 if self._dim == 2 else 0.06

            if x < 0.4:
                x -= offset
            elif x > 0.6:
                x += offset
            if y < 0.1:
                y -= offset
            elif y > 0.8:
                y += offset

            if self._dim == 4 and z > 0.8:
                z += offset

            annotation = plotly_layouts["default_annotation_layout"].copy()
            annotation.update(x=x, y=y, font=font_dict, text=clean_formula, opacity=opacity)

            if self._dim in (3, 4):
                for d in ["xref", "yref"]:
                    annotation.pop(d)  # Scatter3d cannot contain xref, yref
                    if self._dim == 3:
                        annotation.update(x=y, y=x)
                        if entry.composition.is_element:
                            z = 0.9 * self._min_energy  # place label 10% above base

                annotation["z"] = z

            annotations_list.append(annotation)

        # extra point ensures equilateral triangular scaling is displayed
        if self._dim == 3:
            annotations_list.append({"x": 1, "y": 1, "z": 0, "opacity": 0, "text": ""})

        return annotations_list

    def _create_plotly_markers(self, highlight_entries=None, label_uncertainties=False):
        """
        Creates stable and unstable marker plots for overlaying on the phase diagram.

        Returns:
            tuple[go.Scatter]: Plotly Scatter objects (unary, binary), go.Scatterternary(ternary_2d),
            or go.Scatter3d (ternary_3d, quaternary) objects in order: (stable markers, unstable markers)
        """

        def get_marker_props(coords, entries):
            """Get marker locations, hovertext, and error bars from pd_plot_data."""
            x, y, z, texts, energies, uncertainties = [], [], [], [], [], []

            is_stable = [entry in self._pd.stable_entries for entry in entries]
            for coord, entry, stable in zip(coords, entries, is_stable, strict=True):
                energy = round(self._pd.get_form_energy_per_atom(entry), 3)

                entry_id = getattr(entry, "entry_id", "no ID")
                comp = entry.composition

                if hasattr(entry, "original_entry"):
                    orig_entry = entry.original_entry
                    comp = orig_entry.composition
                    entry_id = getattr(orig_entry, "entry_id", "no ID")

                formula = comp.reduced_formula
                clean_formula = htmlify(formula)
                label = f"{clean_formula} ({entry_id}) <br>  Formation energy: {energy} eV/atom <br> "
                if not stable:
                    e_above_hull = round(self._pd.get_e_above_hull(entry), 3)
                    if e_above_hull > self.show_unstable:
                        continue
                    label += f" Energy Above Hull: ({e_above_hull:+} eV/atom)"
                    energies.append(e_above_hull)
                else:
                    uncertainty = 0
                    label += " (Stable)"
                    if hasattr(entry, "correction_uncertainty_per_atom") and label_uncertainties:
                        uncertainty = round(entry.correction_uncertainty_per_atom, 4)
                        label += f"<br> (Error: +/- {uncertainty} eV/atom)"
                    uncertainties.append(uncertainty)
                    energies.append(energy)

                if self._dim == 3 and self.ternary_style == "2d":
                    label += "<br>"
                    total_sum_el = sum(
                        entry.composition[el] for el, _axis in zip(self._pd.elements, range(self._dim), strict=True)
                    )
                    for el, axis in zip(self._pd.elements, range(self._dim), strict=True):
                        _cartesian_positions = [x, y, z]
                        _cartesian_positions[axis].append(entry.composition[el])
                        label += f"<br> {el}: {round(entry.composition[el] / total_sum_el, 6)}"
                elif self._dim == 3 and self.ternary_style == "3d":
                    x.append(coord[0])
                    y.append(coord[1])
                    z.append(energy)

                    label += "<br>"
                    total_sum_el = sum(
                        entry.composition[el] for el, _axis in zip(self._pd.elements, range(self._dim), strict=True)
                    )
                    for el, _axis in zip(self._pd.elements, range(self._dim), strict=True):
                        label += f"<br> {el}: {round(entry.composition[el] / total_sum_el, 6)}"
                elif self._dim == 4:
                    x.append(coord[0])
                    y.append(coord[1])
                    z.append(coord[2])

                    label += "<br>"
                    total_sum_el = sum(
                        entry.composition[el] for el, _axis in zip(self._pd.elements, range(self._dim), strict=True)
                    )
                    for el, _axis in zip(self._pd.elements, range(self._dim), strict=True):
                        label += f"<br> {el}: {round(entry.composition[el] / total_sum_el, 6)}"
                else:
                    x.append(coord[0])
                    y.append(coord[1])

                texts.append(label)

            return {"x": x, "y": y, "z": z, "texts": texts, "energies": energies, "uncertainties": uncertainties}

        if highlight_entries is None:
            highlight_entries = []

        stable_coords, stable_entries = [], []
        unstable_coords, unstable_entries = [], []
        highlight_coords, highlight_ents = [], []

        for coord, entry in zip(self.pd_plot_data[1], self.pd_plot_data[1].values(), strict=True):
            if entry in highlight_entries:
                highlight_coords.append(coord)
                highlight_ents.append(entry)
            else:
                stable_coords.append(coord)
                stable_entries.append(entry)

        for coord, entry in zip(self.pd_plot_data[2].values(), self.pd_plot_data[2], strict=True):
            if entry in highlight_entries:
                highlight_coords.append(coord)
                highlight_ents.append(entry)
            else:
                unstable_coords.append(coord)
                unstable_entries.append(entry)

        stable_props = get_marker_props(stable_coords, stable_entries)
        unstable_props = get_marker_props(unstable_coords, unstable_entries)
        highlight_props = get_marker_props(highlight_coords, highlight_entries)

        stable_markers, unstable_markers, highlight_markers = {}, {}, {}

        if self._dim == 1:
            stable_markers = plotly_layouts["default_unary_marker_settings"].copy()
            unstable_markers = plotly_layouts["default_unary_marker_settings"].copy()

            stable_markers.update(
                x=[0] * len(stable_props["y"]),
                y=list(stable_props["x"]),
                name="Stable",
                marker={
                    "color": "darkgreen",
                    "size": 20,
                    "line": {"color": "black", "width": 2},
                    "symbol": "star",
                },
                opacity=0.9,
                hovertext=stable_props["texts"],
                error_y={
                    "array": list(stable_props["uncertainties"]),
                    "type": "data",
                    "color": "gray",
                    "thickness": 2.5,
                    "width": 5,
                },
            )
            plotly_layouts["unstable_colorscale"].copy()
            unstable_markers.update(
                x=[0] * len(unstable_props["y"]),
                y=list(unstable_props["x"]),
                name="Above Hull",
                marker={
                    "color": unstable_props["energies"],
                    "colorscale": plotly_layouts["unstable_colorscale"],
                    "size": 16,
                    "symbol": "diamond-wide",
                    "line": {"color": "black", "width": 2},
                },
                hovertext=unstable_props["texts"],
                opacity=0.9,
            )

            if highlight_entries:
                highlight_markers = plotly_layouts["default_unary_marker_settings"].copy()
                highlight_markers.update(
                    x=[0] * len(highlight_props["y"]),
                    y=list(highlight_props["x"]),
                    name="Highlighted",
                    marker={
                        "color": "mediumvioletred",
                        "size": 22,
                        "line": {"color": "black", "width": 2},
                        "symbol": "square",
                    },
                    opacity=0.9,
                    hovertext=highlight_props["texts"],
                    error_y={
                        "array": list(highlight_props["uncertainties"]),
                        "type": "data",
                        "color": "gray",
                        "thickness": 2.5,
                        "width": 5,
                    },
                )

        if self._dim == 2:
            stable_markers = plotly_layouts["default_binary_marker_settings"].copy()
            unstable_markers = plotly_layouts["default_binary_marker_settings"].copy()

            stable_markers.update(
                x=list(stable_props["x"]),
                y=list(stable_props["y"]),
                name="Stable",
                marker={"color": "darkgreen", "size": 16, "line": {"color": "black", "width": 2}},
                opacity=0.99,
                hovertext=stable_props["texts"],
                error_y={
                    "array": list(stable_props["uncertainties"]),
                    "type": "data",
                    "color": "gray",
                    "thickness": 2.5,
                    "width": 5,
                },
            )
            unstable_markers |= {
                "x": list(unstable_props["x"]),
                "y": list(unstable_props["y"]),
                "name": "Above Hull",
                "marker": {
                    "color": unstable_props["energies"],
                    "colorscale": plotly_layouts["unstable_colorscale"],
                    "size": 7,
                    "symbol": "diamond",
                    "line": {"color": "black", "width": 1},
                    "opacity": 0.8,
                },
                "hovertext": unstable_props["texts"],
            }
            if highlight_entries:
                highlight_markers = plotly_layouts["default_binary_marker_settings"].copy()
                highlight_markers.update(
                    x=list(highlight_props["x"]),
                    y=list(highlight_props["y"]),
                    name="Highlighted",
                    marker={
                        "color": "mediumvioletred",
                        "size": 16,
                        "line": {"color": "black", "width": 2},
                        "symbol": "square",
                    },
                    opacity=0.99,
                    hovertext=highlight_props["texts"],
                    error_y={
                        "array": list(highlight_props["uncertainties"]),
                        "type": "data",
                        "color": "gray",
                        "thickness": 2.5,
                        "width": 5,
                    },
                )

        elif self._dim == 3 and self.ternary_style == "2d":
            stable_markers = plotly_layouts["default_ternary_2d_marker_settings"].copy()
            unstable_markers = plotly_layouts["default_ternary_2d_marker_settings"].copy()

            stable_markers |= {
                "a": list(stable_props["x"]),
                "b": list(stable_props["y"]),
                "c": list(stable_props["z"]),
                "name": "Stable",
                "hovertext": stable_props["texts"],
                "marker": {
                    "color": "green",
                    "line": {"width": 2.0, "color": "black"},
                    "symbol": "circle",
                    "size": 15,
                },
            }
            unstable_markers |= {
                "a": unstable_props["x"],
                "b": unstable_props["y"],
                "c": unstable_props["z"],
                "name": "Above Hull",
                "hovertext": unstable_props["texts"],
                "marker": {
                    "color": unstable_props["energies"],
                    "opacity": 0.8,
                    "colorscale": plotly_layouts["unstable_colorscale"],
                    "line": {"width": 1, "color": "black"},
                    "size": 7,
                    "symbol": "diamond",
                    "colorbar": {
                        "title": "Energy Above Hull<br>(eV/atom)",
                        "x": 0,
                        "y": 1,
                        "yanchor": "top",
                        "xpad": 0,
                        "ypad": 0,
                        "thickness": 0.02,
                        "thicknessmode": "fraction",
                        "len": 0.5,
                    },
                },
            }
            if highlight_entries:
                highlight_markers = plotly_layouts["default_ternary_2d_marker_settings"].copy()
                highlight_markers |= {
                    "a": list(highlight_props["x"]),
                    "b": list(highlight_props["y"]),
                    "c": list(highlight_props["z"]),
                    "name": "Highlighted",
                    "hovertext": highlight_props["texts"],
                    "marker": {
                        "color": "mediumvioletred",
                        "line": {"width": 2.0, "color": "black"},
                        "symbol": "square",
                        "size": 16,
                    },
                }

        elif self._dim == 3 and self.ternary_style == "3d":
            stable_markers = plotly_layouts["default_ternary_3d_marker_settings"].copy()
            unstable_markers = plotly_layouts["default_ternary_3d_marker_settings"].copy()

            stable_markers |= {
                "x": list(stable_props["y"]),
                "y": list(stable_props["x"]),
                "z": list(stable_props["z"]),
                "name": "Stable",
                "marker": {
                    "color": "#1e1e1f",
                    "size": 11,
                    "opacity": 0.99,
                },
                "hovertext": stable_props["texts"],
                "error_z": {
                    "array": list(stable_props["uncertainties"]),
                    "type": "data",
                    "color": "darkgray",
                    "width": 10,
                    "thickness": 5,
                },
            }
            unstable_markers |= {
                "x": unstable_props["y"],
                "y": unstable_props["x"],
                "z": unstable_props["z"],
                "name": "Above Hull",
                "hovertext": unstable_props["texts"],
                "marker": {
                    "color": unstable_props["energies"],
                    "colorscale": plotly_layouts["unstable_colorscale"],
                    "size": 5,
                    "line": {"color": "black", "width": 1},
                    "symbol": "diamond",
                    "opacity": 0.7,
                    "colorbar": {
                        "title": "Energy Above Hull<br>(eV/atom)",
                        "x": 0,
                        "y": 1,
                        "yanchor": "top",
                        "xpad": 0,
                        "ypad": 0,
                        "thickness": 0.02,
                        "thicknessmode": "fraction",
                        "len": 0.5,
                    },
                },
            }
            if highlight_entries:
                highlight_markers = plotly_layouts["default_ternary_3d_marker_settings"].copy()
                highlight_markers |= {
                    "x": list(highlight_props["y"]),
                    "y": list(highlight_props["x"]),
                    "z": list(highlight_props["z"]),
                    "name": "Highlighted",
                    "marker": {
                        "size": 12,
                        "opacity": 0.99,
                        "symbol": "square",
                        "color": "mediumvioletred",
                    },
                    "hovertext": highlight_props["texts"],
                    "error_z": {
                        "array": list(highlight_props["uncertainties"]),
                        "type": "data",
                        "color": "darkgray",
                        "width": 10,
                        "thickness": 5,
                    },
                }

        elif self._dim == 4:
            stable_markers = plotly_layouts["default_quaternary_marker_settings"].copy()
            unstable_markers = plotly_layouts["default_quaternary_marker_settings"].copy()
            stable_markers |= {
                "x": stable_props["x"],
                "y": stable_props["y"],
                "z": stable_props["z"],
                "name": "Stable",
                "marker": {
                    "size": 7,
                    "opacity": 0.99,
                    "color": "darkgreen",
                    "line": {"color": "black", "width": 1},
                },
                "hovertext": stable_props["texts"],
            }
            unstable_markers |= {
                "x": unstable_props["x"],
                "y": unstable_props["y"],
                "z": unstable_props["z"],
                "name": "Above Hull",
                "marker": {
                    "color": unstable_props["energies"],
                    "colorscale": plotly_layouts["unstable_colorscale"],
                    "size": 5,
                    "symbol": "diamond",
                    "line": {"color": "black", "width": 1},
                    "colorbar": {
                        "title": "Energy Above Hull<br>(eV/atom)",
                        "x": 0,
                        "y": 1,
                        "yanchor": "top",
                        "xpad": 0,
                        "ypad": 0,
                        "thickness": 0.02,
                        "thicknessmode": "fraction",
                        "len": 0.5,
                    },
                },
                "hovertext": unstable_props["texts"],
                "visible": "legendonly",
            }
            if highlight_entries:
                highlight_markers = plotly_layouts["default_quaternary_marker_settings"].copy()
                highlight_markers |= {
                    "x": highlight_props["x"],
                    "y": highlight_props["y"],
                    "z": highlight_props["z"],
                    "name": "Highlighted",
                    "marker": {
                        "size": 9,
                        "opacity": 0.99,
                        "symbol": "square",
                        "color": "mediumvioletred",
                        "line": {"color": "black", "width": 1},
                    },
                    "hovertext": highlight_props["texts"],
                }

        highlight_marker_plot = None

        if self._dim in [1, 2]:
            stable_marker_plot, unstable_marker_plot = (
                go.Scatter(**markers) for markers in [stable_markers, unstable_markers]
            )

            if highlight_entries:
                highlight_marker_plot = go.Scatter(**highlight_markers)
        elif self._dim == 3 and self.ternary_style == "2d":
            stable_marker_plot, unstable_marker_plot = (
                go.Scatterternary(**markers) for markers in [stable_markers, unstable_markers]
            )
            if highlight_entries:
                highlight_marker_plot = go.Scatterternary(**highlight_markers)
        else:
            stable_marker_plot, unstable_marker_plot = (
                go.Scatter3d(**markers) for markers in [stable_markers, unstable_markers]
            )
            if highlight_entries:
                highlight_marker_plot = go.Scatter3d(**highlight_markers)

        return stable_marker_plot, unstable_marker_plot, highlight_marker_plot

    def _create_plotly_uncertainty_shading(self, stable_marker_plot):
        """
        Creates shaded uncertainty region for stable entries. Currently only works
        for binary (dim=2) phase diagrams.

        Args:
            stable_marker_plot: go.Scatter object with stable markers and their
            error bars.

        Returns:
            Plotly go.Scatter object with uncertainty window shading.
        """
        uncertainty_plot = None

        x = stable_marker_plot.x
        y = stable_marker_plot.y

        transformed = False
        if hasattr(self._pd, "original_entries") or hasattr(self._pd, "chempots"):
            transformed = True

        if self._dim == 2:
            error = stable_marker_plot.error_y["array"]

            points = np.append(x, [y, error]).reshape(3, -1).T
            points = points[points[:, 0].argsort()]  # sort by composition

            # these steps trace out the boundary pts of the uncertainty window
            outline = points[:, :2].copy()
            outline[:, 1] += points[:, 2]

            last = -1
            if transformed:
                last = None  # allows for uncertainty in terminal compounds

            flipped_points = np.flip(points[:last, :].copy(), axis=0)
            flipped_points[:, 1] -= flipped_points[:, 2]
            outline = np.vstack((outline, flipped_points[:, :2]))

            uncertainty_plot = go.Scatter(
                x=outline[:, 0],
                y=outline[:, 1],
                name="Uncertainty (window)",
                fill="toself",
                mode="lines",
                line={"width": 0},
                fillcolor="lightblue",
                hoverinfo="skip",
                opacity=0.4,
            )

        return uncertainty_plot

    def _create_plotly_ternary_support_lines(self):
        """
        Creates support lines which aid in seeing the ternary hull in three
        dimensions.

        Returns:
            go.Scatter3d plot of support lines for ternary phase diagram.
        """
        stable_entry_coords = dict(map(reversed, self.pd_plot_data[1].items()))

        elem_coords = [stable_entry_coords[entry] for entry in self._pd.el_refs.values()]

        # add top and bottom triangle guidelines
        x, y, z = [], [], []
        for line in itertools.combinations(elem_coords, 2):
            x.extend([line[0][0], line[1][0], None] * 2)
            y.extend([line[0][1], line[1][1], None] * 2)
            z.extend([0, 0, None, self._min_energy, self._min_energy, None])

        # add vertical guidelines
        for elem in elem_coords:
            x.extend([elem[0], elem[0], None])
            y.extend([elem[1], elem[1], None])
            z.extend([0, self._min_energy, None])

        return go.Scatter3d(
            x=list(y),
            y=list(x),
            z=list(z),
            mode="lines",
            hoverinfo="none",
            line={"color": "rgba (0, 0, 0, 0.4)", "dash": "solid", "width": 1.0},
            showlegend=False,
        )

    @no_type_check
    def _get_matplotlib_2d_plot(
        self,
        label_stable=True,
        label_unstable=True,
        ordering=None,
        energy_colormap=None,
        vmin_mev=-60.0,
        vmax_mev=60.0,
        show_colorbar=True,
        process_attributes=False,
        ax: plt.Axes = None,
    ):
        """Show the plot using matplotlib.

        Imports are done within the function as matplotlib is no longer the default.
        """
        ax = ax or pretty_plot(8, 6)

        if ordering is None:
            lines, labels, unstable = self.pd_plot_data
        else:
            _lines, _labels, _unstable = self.pd_plot_data
            lines, labels, unstable = order_phase_diagram(_lines, _labels, _unstable, ordering)
        if energy_colormap is None:
            if process_attributes:
                for x, y in lines:
                    plt.plot(x, y, "k-", linewidth=3, markeredgecolor="k")
                # One should think about a clever way to have "complex"
                # attributes with complex processing options but with a clear
                # logic. At this moment, I just use the attributes to know
                # whether an entry is a new compound or an existing (from the
                #  ICSD or from the MP) one.
                for x, y in labels:
                    if labels[x, y].attribute is None or labels[x, y].attribute == "existing":
                        plt.plot(x, y, "ko", **self.plotkwargs)
                    else:
                        plt.plot(x, y, "k*", **self.plotkwargs)
            else:
                for x, y in lines:
                    plt.plot(x, y, "ko-", **self.plotkwargs)
        else:
            for x, y in lines:
                plt.plot(x, y, "k-", markeredgecolor="k")
            vmin = vmin_mev / 1000.0
            vmax = vmax_mev / 1000.0
            if energy_colormap == "default":
                mid = -vmin / (vmax - vmin)
                cmap = LinearSegmentedColormap.from_list(
                    "custom_colormap",
                    [(0.0, "#005500"), (mid, "#55FF55"), (mid, "#FFAAAA"), (1.0, "#FF0000")],
                )
            else:
                cmap = energy_colormap
            norm = Normalize(vmin=vmin, vmax=vmax)
            _map = ScalarMappable(norm=norm, cmap=cmap)
            _energies = [self._pd.get_equilibrium_reaction_energy(entry) for coord, entry in labels.items()]
            energies = [en if en < 0 else -0.000_000_01 for en in _energies]
            vals_stable = _map.to_rgba(energies)
            ii = 0
            if process_attributes:
                for x, y in labels:
                    if labels[x, y].attribute is None or labels[x, y].attribute == "existing":
                        plt.plot(x, y, "o", markerfacecolor=vals_stable[ii], markersize=12)
                    else:
                        plt.plot(x, y, "*", markerfacecolor=vals_stable[ii], markersize=18)
                    ii += 1
            else:
                for x, y in labels:
                    plt.plot(x, y, "o", markerfacecolor=vals_stable[ii], markersize=15)
                    ii += 1

        font = FontProperties()
        font.set_weight("bold")
        font.set_size(24)

        # Sets a nice layout depending on the type of PD. Also defines a
        # "center" for the PD, which then allows the annotations to be spread
        # out in a nice manner.
        if len(self._pd.elements) == 3:
            plt.axis("equal")
            plt.xlim((-0.1, 1.2))
            plt.ylim((-0.1, 1.0))
            plt.axis("off")
            center = (0.5, math.sqrt(3) / 6)
        else:
            miny = min(c[1] for c in labels)
            ybuffer = max(abs(miny) * 0.1, 0.1)
            plt.xlim((-0.1, 1.1))
            plt.ylim((miny - ybuffer, ybuffer))
            center = (0.5, miny / 2)
            plt.xlabel("Fraction", fontsize=28, fontweight="bold")
            plt.ylabel("Formation energy (eV/atom)", fontsize=28, fontweight="bold")

        for coords in sorted(labels, key=lambda x: -x[1]):
            entry = labels[coords]
            label = entry.name

            # The follow defines an offset for the annotation text emanating
            # from the center of the PD. Results in fairly nice layouts for the
            # most part.
            vec = np.array(coords) - center
            vec = vec / np.linalg.norm(vec) * 10 if np.linalg.norm(vec) != 0 else vec
            valign = "bottom" if vec[1] > 0 else "top"
            if vec[0] < -0.01:
                halign = "right"
            elif vec[0] > 0.01:
                halign = "left"
            else:
                halign = "center"
            if label_stable:
                if process_attributes and entry.attribute == "new":
                    plt.annotate(
                        latexify(label),
                        coords,
                        xytext=vec,
                        textcoords="offset points",
                        horizontalalignment=halign,
                        verticalalignment=valign,
                        fontproperties=font,
                        color="g",
                    )
                else:
                    plt.annotate(
                        latexify(label),
                        coords,
                        xytext=vec,
                        textcoords="offset points",
                        horizontalalignment=halign,
                        verticalalignment=valign,
                        fontproperties=font,
                    )

        if self.show_unstable:
            font = FontProperties()
            font.set_size(16)
            energies_unstable = [self._pd.get_e_above_hull(entry) for entry, coord in unstable.items()]
            if energy_colormap is not None:
                energies.extend(energies_unstable)
                vals_unstable = _map.to_rgba(energies_unstable)
            ii = 0
            for entry, coords in unstable.items():
                ehull = self._pd.get_e_above_hull(entry)
                if ehull < self.show_unstable:
                    vec = np.array(coords) - center
                    vec = vec / np.linalg.norm(vec) * 10 if np.linalg.norm(vec) != 0 else vec
                    label = entry.name
                    if energy_colormap is None:
                        plt.plot(
                            coords[0],
                            coords[1],
                            "ks",
                            linewidth=3,
                            markeredgecolor="k",
                            markerfacecolor="r",
                            markersize=8,
                        )
                    else:
                        plt.plot(
                            coords[0],
                            coords[1],
                            "s",
                            linewidth=3,
                            markeredgecolor="k",
                            markerfacecolor=vals_unstable[ii],
                            markersize=8,
                        )
                    if label_unstable:
                        plt.annotate(
                            latexify(label),
                            coords,
                            xytext=vec,
                            textcoords="offset points",
                            horizontalalignment=halign,
                            color="b",
                            verticalalignment=valign,
                            fontproperties=font,
                        )
                    ii += 1
        if energy_colormap is not None and show_colorbar:
            _map.set_array(energies)
            cbar = plt.colorbar(_map)
            cbar.set_label(
                "Energy [meV/at] above hull (positive values)\nInverse energy [meV/at] above hull (negative values)",
                rotation=-90,
                ha="center",
                va="bottom",
            )
        fig = plt.gcf()
        fig.set_size_inches((8, 6))
        plt.subplots_adjust(left=0.09, right=0.98, top=0.98, bottom=0.07)
        return ax

    @no_type_check
    def _get_matplotlib_3d_plot(self, label_stable=True, ax: plt.Axes = None):
        """Show the plot using matplotlib.

        Args:
            label_stable (bool): Whether to label stable compounds.
            ax (plt.Axes): An existing axes object (optional). If not provided, a new one will be created.

        Returns:
            plt.Axes: The axes object with the plot.
        """
        ax = ax or plt.figure().add_subplot(111, projection="3d")

        font = FontProperties(weight="bold", size=13)
        lines, labels, _ = self.pd_plot_data
        count = 1
        newlabels = []
        for x, y, z in lines:
            ax.plot(
                x,
                y,
                z,
                "bo-",
                linewidth=3,
                markeredgecolor="b",
                markerfacecolor="r",
                markersize=10,
            )
        for coords in sorted(labels):
            entry = labels[coords]
            label = entry.name
            if label_stable:
                if len(entry.elements) == 1:
                    ax.text(coords[0], coords[1], coords[2], label, fontproperties=font)
                else:
                    ax.text(coords[0], coords[1], coords[2], str(count), fontsize=12)
                    newlabels.append(f"{count} : {latexify(label)}")
                    count += 1
        plt.figtext(0.01, 0.01, "\n".join(newlabels), fontproperties=font)
        ax.axis("off")
        ax.set(xlim=(-0.1, 0.72), ylim=(0, 0.66), zlim=(0, 0.56))
        return ax


def uniquelines(q):
    """
    Given all the facets, convert it into a set of unique lines. Specifically
    used for converting convex hull facets into line pairs of coordinates.

    Args:
        q: A 2-dim sequence, where each row represents a facet. e.g.
            [[1,2,3],[3,6,7],...]

    Returns:
        setoflines:
            A set of tuple of lines. e.g. ((1,2), (1,3), (2,3), ....)
    """
    return {tuple(sorted(line)) for facets in q for line in itertools.combinations(facets, 2)}


def triangular_coord(coord):
    """
    Convert a 2D coordinate into a triangle-based coordinate system for a
    prettier phase diagram.

    Args:
        coord: coordinate used in the convex hull computation.

    Returns:
        coordinates in a triangular-based coordinate system.
    """
    unit_vec = np.array([[1, 0], [0.5, math.sqrt(3) / 2]])

    result = np.dot(np.array(coord), unit_vec)
    return result.transpose()


def tet_coord(coord):
    """
    Convert a 3D coordinate into a tetrahedron based coordinate system for a
    prettier phase diagram.

    Args:
        coord: coordinate used in the convex hull computation.

    Returns:
        coordinates in a tetrahedron-based coordinate system.
    """
    unitvec = np.array(
        [
            [1, 0, 0],
            [0.5, math.sqrt(3) / 2, 0],
            [0.5, 1 / 3 * math.sqrt(3) / 2, math.sqrt(6) / 3],
        ]
    )
    result = np.dot(np.array(coord), unitvec)
    return result.transpose()


def order_phase_diagram(lines, stable_entries, unstable_entries, ordering):
    """
    Orders the entries (their coordinates) in a phase diagram plot according
    to the user specified ordering.
    Ordering should be given as ['Up', 'Left', 'Right'], where Up,
    Left and Right are the names of the entries in the upper, left and right
    corners of the triangle respectively.

    Args:
        lines: list of list of coordinates for lines in the PD.
        stable_entries: {coordinate : entry} for each stable node in the
            phase diagram. (Each coordinate can only have one stable phase)
        unstable_entries: {entry: coordinates} for all unstable nodes in the
            phase diagram.
        ordering: Ordering of the phase diagram, given as a list ['Up',
            'Left','Right']

    Returns:
        tuple[list, dict, dict]:
            - new_lines is a list of list of coordinates for lines in the PD.
            - new_stable_entries is a {coordinate: entry} for each stable node
            in the phase diagram. (Each coordinate can only have one
            stable phase)
            - new_unstable_entries is a {entry: coordinates} for all unstable
            nodes in the phase diagram.
    """
    yup = -1000.0
    xleft = 1000.0
    xright = -1000.0

    nameup = ""
    nameleft = ""
    nameright = ""
    for coord in stable_entries:
        if coord[0] > xright:
            xright = coord[0]
            nameright = stable_entries[coord].name
        if coord[0] < xleft:
            xleft = coord[0]
            nameleft = stable_entries[coord].name
        if coord[1] > yup:
            yup = coord[1]
            nameup = stable_entries[coord].name

    if (nameup not in ordering) or (nameright not in ordering) or (nameleft not in ordering):
        raise ValueError(
            "Error in ordering_phase_diagram :\n"
            f"{nameup!r}, {nameleft!r} and {nameright!r} should be in ordering : {ordering}"
        )

    cc = np.array([0.5, np.sqrt(3.0) / 6.0], float)

    if nameup == ordering[0]:
        if nameleft == ordering[1]:
            # The coordinates were already in the user ordering
            return lines, stable_entries, unstable_entries

        new_lines = [[np.array(1 - x), y] for x, y in lines]
        new_stable_entries = {(1 - c[0], c[1]): entry for c, entry in stable_entries.items()}
        new_unstable_entries = {entry: (1 - c[0], c[1]) for entry, c in unstable_entries.items()}
        return new_lines, new_stable_entries, new_unstable_entries
    if nameup == ordering[1]:
        if nameleft == ordering[2]:
            c120 = np.cos(2 * np.pi / 3.0)
            s120 = np.sin(2 * np.pi / 3.0)
            new_lines = []
            for x, y in lines:
                newx = np.zeros_like(x)
                newy = np.zeros_like(y)
                for ii, xx in enumerate(x):
                    newx[ii] = c120 * (xx - cc[0]) - s120 * (y[ii] - cc[1]) + cc[0]
                    newy[ii] = s120 * (xx - cc[0]) + c120 * (y[ii] - cc[1]) + cc[1]
                new_lines.append([newx, newy])
            new_stable_entries = {
                (
                    c120 * (c[0] - cc[0]) - s120 * (c[1] - cc[1]) + cc[0],
                    s120 * (c[0] - cc[0]) + c120 * (c[1] - cc[1]) + cc[1],
                ): entry
                for c, entry in stable_entries.items()
            }
            new_unstable_entries = {
                entry: (
                    c120 * (c[0] - cc[0]) - s120 * (c[1] - cc[1]) + cc[0],
                    s120 * (c[0] - cc[0]) + c120 * (c[1] - cc[1]) + cc[1],
                )
                for entry, c in unstable_entries.items()
            }
            return new_lines, new_stable_entries, new_unstable_entries
        c120 = np.cos(2 * np.pi / 3.0)
        s120 = np.sin(2 * np.pi / 3.0)
        new_lines = []
        for x, y in lines:
            newx = np.zeros_like(x)
            newy = np.zeros_like(y)
            for ii, xx in enumerate(x):
                newx[ii] = -c120 * (xx - 1.0) - s120 * y[ii] + 1.0
                newy[ii] = -s120 * (xx - 1.0) + c120 * y[ii]
            new_lines.append([newx, newy])
        new_stable_entries = {
            (
                -c120 * (c[0] - 1.0) - s120 * c[1] + 1.0,
                -s120 * (c[0] - 1.0) + c120 * c[1],
            ): entry
            for c, entry in stable_entries.items()
        }
        new_unstable_entries = {
            entry: (
                -c120 * (c[0] - 1.0) - s120 * c[1] + 1.0,
                -s120 * (c[0] - 1.0) + c120 * c[1],
            )
            for entry, c in unstable_entries.items()
        }
        return new_lines, new_stable_entries, new_unstable_entries
    if nameup == ordering[2]:
        if nameleft == ordering[0]:
            c240 = np.cos(4 * np.pi / 3.0)
            s240 = np.sin(4 * np.pi / 3.0)
            new_lines = []
            for x, y in lines:
                newx = np.zeros_like(x)
                newy = np.zeros_like(y)
                for ii, xx in enumerate(x):
                    newx[ii] = c240 * (xx - cc[0]) - s240 * (y[ii] - cc[1]) + cc[0]
                    newy[ii] = s240 * (xx - cc[0]) + c240 * (y[ii] - cc[1]) + cc[1]
                new_lines.append([newx, newy])
            new_stable_entries = {
                (
                    c240 * (c[0] - cc[0]) - s240 * (c[1] - cc[1]) + cc[0],
                    s240 * (c[0] - cc[0]) + c240 * (c[1] - cc[1]) + cc[1],
                ): entry
                for c, entry in stable_entries.items()
            }
            new_unstable_entries = {
                entry: (
                    c240 * (c[0] - cc[0]) - s240 * (c[1] - cc[1]) + cc[0],
                    s240 * (c[0] - cc[0]) + c240 * (c[1] - cc[1]) + cc[1],
                )
                for entry, c in unstable_entries.items()
            }
            return new_lines, new_stable_entries, new_unstable_entries
        c240 = np.cos(4 * np.pi / 3.0)
        s240 = np.sin(4 * np.pi / 3.0)
        new_lines = []
        for x, y in lines:
            newx = np.zeros_like(x)
            newy = np.zeros_like(y)
            for ii, xx in enumerate(x):
                newx[ii] = -c240 * xx - s240 * y[ii]
                newy[ii] = -s240 * xx + c240 * y[ii]
            new_lines.append([newx, newy])
        new_stable_entries = {
            (-c240 * c[0] - s240 * c[1], -s240 * c[0] + c240 * c[1]): entry for c, entry in stable_entries.items()
        }
        new_unstable_entries = {
            entry: (-c240 * c[0] - s240 * c[1], -s240 * c[0] + c240 * c[1]) for entry, c in unstable_entries.items()
        }
        return new_lines, new_stable_entries, new_unstable_entries
    raise ValueError("Invalid ordering.")
