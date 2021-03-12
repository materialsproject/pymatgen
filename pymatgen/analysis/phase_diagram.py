# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module defines tools to generate and analyze phase diagrams.
"""

import collections
import itertools
import json
import logging
import math
import os
import re
from functools import lru_cache

import numpy as np
import plotly.graph_objs as go
from monty.json import MontyDecoder, MSONable
from scipy.optimize import minimize
from scipy.spatial import ConvexHull

from pymatgen.analysis.reaction_calculator import Reaction, ReactionError
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import DummySpecies, Element, get_el_sp
from pymatgen.entries import Entry
from pymatgen.util.coord import Simplex, in_coord_list
from pymatgen.util.plotting import pretty_plot
from pymatgen.util.string import latexify

logger = logging.getLogger(__name__)

with open(os.path.join(os.path.dirname(__file__), "..", "util", "plotly_pd_layouts.json")) as f:
    plotly_layouts = json.load(f)


class PDEntry(Entry):
    """
    An object encompassing all relevant data for phase diagrams.

    Attributes:
        composition (Composition): The composition associated with the PDEntry.
        energy (float): The energy associated with the entry.
        name (str):  A name for the entry. This is the string shown in the phase diagrams.
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
        name: str = None,
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
        self.name = name if name else self.composition.reduced_formula
        self.attribute = attribute

    @property
    def energy(self) -> float:
        """
        :return: the energy of the entry.
        """
        return self._energy

    def as_dict(self):
        """
        :return: MSONable dict.
        """
        return_dict = super().as_dict()
        return_dict.update({"name": self.name, "attribute": self.attribute})
        return return_dict

    @classmethod
    def from_dict(cls, d):
        """
        :param d: Dict representation
        :return: PDEntry
        """
        return cls(
            Composition(d["composition"]),
            d["energy"],
            d["name"] if "name" in d else None,
            d["attribute"] if "attribute" in d else None,
        )


class GrandPotPDEntry(PDEntry):
    """
    A grand potential pd entry object encompassing all relevant data for phase
    diagrams.  Chemical potentials are given as a element-chemical potential
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
            name if name else entry.name,
            entry.attribute if hasattr(entry, "attribute") else None,
        )
        # NOTE if we init GrandPotPDEntry from ComputedEntry _energy is the
        # corrected energy of the ComputedEntry hence the need to keep
        # the original entry to not lose data.
        self.original_entry = entry
        self.original_comp = self._composition
        self.chempots = chempots

    @property
    def composition(self) -> Composition:
        """The composition after removing free species

        Returns:
            Composition
        """
        return Composition({el: self._composition[el] for el in self._composition.elements if el not in self.chempots})

    @property
    def chemical_energy(self):
        """The chemical energy term mu*N in the grand potential

        Returns:
            The chemical energy term mu*N in the grand potential
        """
        return sum([self._composition[el] * pot for el, pot in self.chempots.items()])

    @property
    def energy(self):
        """
        Returns:
            The grand potential energy
        """
        return self._energy - self.chemical_energy

    def __repr__(self):
        chempot_str = " ".join(["mu_%s = %.4f" % (el, mu) for el, mu in self.chempots.items()])
        return "GrandPotPDEntry with original composition " + "{}, energy = {:.4f}, {}".format(
            self.original_entry.composition, self.original_entry.energy, chempot_str
        )

    def as_dict(self):
        """
        :return: MSONAble dict
        """
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "entry": self.original_entry.as_dict(),
            "chempots": {el.symbol: u for el, u in self.chempots.items()},
            "name": self.name,
        }

    @classmethod
    def from_dict(cls, d):
        """
        :param d: Dict representation
        :return: PDStructureEntry
        """
        chempots = {Element(symbol): u for symbol, u in d["chempots"].items()}
        entry = MontyDecoder().process_decoded(d["entry"])
        return cls(entry, chempots, d["name"])


class TransformedPDEntry(PDEntry):
    """
    This class repesents a TransformedPDEntry, which allows for a PDEntry to be
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
            sp_mapping ({Composition: DummySpecies}): dictionary
                mapping Terminal Compositions to Dummy Species

        """
        super().__init__(
            entry.composition,
            entry.energy,
            name if name else entry.name,
            entry.attribute if hasattr(entry, "attribute") else None,
        )
        self.original_entry = entry
        self.sp_mapping = sp_mapping

        self.rxn = Reaction(list(self.sp_mapping.keys()), [self._composition])
        self.rxn.normalize_to(self.original_entry.composition)

        # NOTE We only allow reactions that have positive amounts of reactants.
        if not all(self.rxn.get_coeff(comp) <= TransformedPDEntry.amount_tol for comp in self.sp_mapping.keys()):
            raise TransformedPDEntryError("Only reactions with positive amounts of reactants allowed")

    @property
    def composition(self) -> Composition:
        """The composition in the dummy species space

        Returns:
            Composition
        """
        # NOTE this is not infallable as the original entry is mutable and an
        # end user could choose to normalize or change the original entry.
        # However, the risk of this seems low.
        factor = self._composition.num_atoms / self.original_entry.composition.num_atoms

        trans_comp = {self.sp_mapping[comp]: -self.rxn.get_coeff(comp) for comp in self.sp_mapping}

        trans_comp = {k: v * factor for k, v in trans_comp.items() if v > TransformedPDEntry.amount_tol}

        return Composition(trans_comp)

    def __repr__(self):
        output = [
            "TransformedPDEntry {}".format(self.composition),
            " with original composition {}".format(self.original_entry.composition),
            ", E = {:.4f}".format(self.original_entry.energy),
        ]
        return "".join(output)

    def as_dict(self):
        """
        :return: MSONable dict
        """
        d = {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "sp_mapping": self.sp_mapping,
        }
        d.update(self.original_entry.as_dict())
        return d

    @classmethod
    def from_dict(cls, d):
        """
        :param d: Dict representation
        :return: TransformedPDEntry
        """
        sp_mapping = d["sp_mapping"]
        del d["sp_mapping"]
        entry = MontyDecoder().process_decoded(d)
        return cls(entry, sp_mapping)


class TransformedPDEntryError(Exception):
    """
    An exception class for TransformedPDEntry.
    """

    pass


class BasePhaseDiagram(MSONable):
    """
    BasePhaseDiagram is not intended to be used directly, and PhaseDiagram should be preferred.

    When constructing a PhaseDiagram, a lot of heavy processing is performed to calculate the
    phase diagram information such as facets, simplexes, etc. The BasePhaseDiagram offers a way to
    store this information so that a phase diagram can be re-constructed without doing this heavy
    processing. It is primarily intended for database applications.
    """

    # Tolerance for determining if formation energy is positive.
    formation_energy_tol = 1e-11
    numerical_tol = 1e-8

    def __init__(
        self,
        facets,
        simplexes,
        all_entries,
        qhull_data,
        dim,
        el_refs,
        elements,
        qhull_entries,
    ):
        """
        This class uses casting to bypass the init, so this constructor should only be
        called by as_dict and from_dict functions. Prefer the PhaseDiagram class for
        typical use cases.
        """
        self.facets = facets
        self.simplexes = simplexes
        self.all_entries = all_entries
        self.qhull_data = qhull_data
        self.dim = dim
        self.el_refs = el_refs
        self.elements = elements
        self.qhull_entries = qhull_entries
        self._stable_entries = set(self.qhull_entries[i] for i in set(itertools.chain(*self.facets)))

    @classmethod
    def from_entries(cls, entries, elements=None):
        """
        Construct the PhaseDiagram object and recast it as a BasePhaseDiagram

        Args:
            entries ([PDEntry]): A list of PDEntry-like objects having an
                energy, energy_per_atom and composition.
            elements ([Element]): Optional list of elements in the phase
                diagram. If set to None, the elements are determined from
                the the entries themselves and are sorted alphabetically.
                If specified, element ordering (e.g. for pd coordinates)
                is preserved.
        """
        return cls(**cls._kwargs_from_entries(entries, elements))

    @classmethod
    def _kwargs_from_entries(cls, entries, elements):
        if elements is None:
            elements = sorted({els for e in entries for els in e.composition.elements})

        elements = list(elements)
        dim = len(elements)

        entries = sorted(entries, key=lambda e: e.composition.reduced_composition)

        el_refs = {}
        min_entries = []
        all_entries = []
        for c, g in itertools.groupby(entries, key=lambda e: e.composition.reduced_composition):
            g = list(g)
            min_entry = min(g, key=lambda e: e.energy_per_atom)
            if c.is_element:
                el_refs[c.elements[0]] = min_entry
            min_entries.append(min_entry)
            all_entries.extend(g)

        if len(el_refs) != dim:
            missing = set(elements).difference(el_refs.keys())
            raise ValueError(f"There are no entries for the terminal elements: {missing}")

        data = np.array(
            [[e.composition.get_atomic_fraction(el) for el in elements] + [e.energy_per_atom] for e in min_entries]
        )

        # Use only entries with negative formation energy
        vec = [el_refs[el].energy_per_atom for el in elements] + [-1]
        form_e = -np.dot(data, vec)
        inds = np.where(form_e < -cls.formation_energy_tol)[0].tolist()

        # Add the elemental references
        inds.extend([min_entries.index(el) for el in el_refs.values()])

        qhull_entries = [min_entries[i] for i in inds]
        qhull_data = data[inds][:, 1:]

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
                m = qhull_data[facet]
                m[:, -1] = 1
                if abs(np.linalg.det(m)) > 1e-14:
                    final_facets.append(facet)
            facets = final_facets

        simplexes = [Simplex(qhull_data[f, :-1]) for f in facets]

        return dict(
            facets=facets,
            simplexes=simplexes,
            all_entries=all_entries,
            qhull_data=qhull_data,
            dim=dim,
            el_refs=el_refs,
            elements=elements,
            qhull_entries=qhull_entries,
        )

    def pd_coords(self, comp):
        """
        The phase diagram is generated in a reduced dimensional space
        (n_elements - 1). This function returns the coordinates in that space.
        These coordinates are compatible with the stored simplex objects.

        Args:
            comp (Composition): A composition

        Returns:
            The coordinates for a given composition in the PhaseDiagram's basis

        """
        if set(comp.elements).difference(self.elements):
            raise ValueError("{} has elements not in the phase diagram {}" "".format(comp, self.elements))
        return np.array([comp.get_atomic_fraction(el) for el in self.elements[1:]])

    @property
    def all_entries_hulldata(self):
        """
        :return: The actual ndarray used to construct the convex hull.
        """
        data = []
        for entry in self.all_entries:
            comp = entry.composition
            row = [comp.get_atomic_fraction(el) for el in self.elements]
            row.append(entry.energy_per_atom)
            data.append(row)
        return np.array(data)[:, 1:]

    @property
    def unstable_entries(self):
        """
        Returns a list of Entries that are unstable in the phase diagram.
        Includes positive formation energy entries.
        """
        return [e for e in self.all_entries if e not in self.stable_entries]

    @property
    def stable_entries(self):
        """
        Returns the set of stable entries in the phase diagram.
        """
        return self._stable_entries

    def get_form_energy(self, entry):
        """
        Returns the formation energy for an entry (NOT normalized) from the
        elemental references.

        Args:
            entry (PDEntry): A PDEntry-like object.

        Returns:
            Formation energy from the elemental references.
        """
        c = entry.composition
        return entry.energy - sum([c[el] * self.el_refs[el].energy_per_atom for el in c.elements])

    def get_form_energy_per_atom(self, entry):
        """
        Returns the formation energy per atom for an entry from the
        elemental references.

        Args:
            entry (PDEntry): An PDEntry-like object

        Returns:
            Formation energy **per atom** from the elemental references.
        """
        return self.get_form_energy(entry) / entry.composition.num_atoms

    def __repr__(self):
        symbols = [el.symbol for el in self.elements]
        output = [
            "{} phase diagram".format("-".join(symbols)),
            "{} stable phases: ".format(len(self.stable_entries)),
            ", ".join([entry.name for entry in self.stable_entries]),
        ]
        return "\n".join(output)

    @lru_cache(1)
    def _get_facet_and_simplex(self, comp):
        """
        Get any facet that a composition falls into. Cached so successive
        calls at same composition are fast.

        Args:
            comp (Composition): A composition

        """
        c = self.pd_coords(comp)
        for f, s in zip(self.facets, self.simplexes):
            if s.in_simplex(c, PhaseDiagram.numerical_tol / 10):
                return f, s
        raise RuntimeError("No facet found for comp = {}".format(comp))

    def _get_facet_chempots(self, facet):
        """
        Calculates the chemical potentials for each element within a facet.

        Args:
            facet: Facet of the phase diagram.

        Returns:
            {element: chempot} for all elements in the phase diagram.
        """
        complist = [self.qhull_entries[i].composition for i in facet]
        energylist = [self.qhull_entries[i].energy_per_atom for i in facet]
        m = [[c.get_atomic_fraction(e) for e in self.elements] for c in complist]
        chempots = np.linalg.solve(m, energylist)
        return dict(zip(self.elements, chempots))

    def get_decomposition(self, comp):
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
            self.qhull_entries[f]: amt for f, amt in zip(facet, decomp_amts) if abs(amt) > PhaseDiagram.numerical_tol
        }

    def get_hull_energy(self, comp):
        """
        Args:
            comp (Composition): Input composition

        Returns:
            Energy of lowest energy equilibrium at desired composition. Not
                normalized by atoms, i.e. E(Li4O2) = 2 * E(Li2O)
        """
        decomp = self.get_decomposition(comp)
        return comp.num_atoms * sum([e.energy_per_atom * n for e, n in decomp.items()])

    def get_decomp_and_e_above_hull(self, entry, allow_negative=False):
        """
        Provides the decomposition and energy above convex hull for an entry.
        Due to caching, can be much faster if entries with the same composition
        are processed together.

        Args:
            entry (PDEntry): A PDEntry like object
            allow_negative (bool): Whether to allow negative e_above_hulls. Used to
                calculate equilibrium reaction energies. Defaults to False.

        Returns:
            (decomp, energy_above_hull). The decomposition is provided
                as a dict of {PDEntry: amount} where amount is the amount of the
                fractional composition. Stable entries should have energy above
                convex hull of 0. The energy is given per atom.
        """
        # Avoid computation for stable_entries.
        # NOTE scaled duplicates of stable_entries will not be caught.
        if entry in list(self.stable_entries):
            return {entry: 1}, 0

        decomp = self.get_decomposition(entry.composition)
        e_above_hull = entry.energy_per_atom - sum([e.energy_per_atom * n for e, n in decomp.items()])

        if allow_negative or e_above_hull >= -PhaseDiagram.numerical_tol:
            return decomp, e_above_hull

        raise ValueError("No valid decomp found for {}! (e {})".format(entry, e_above_hull))

    def get_e_above_hull(self, entry, **kwargs):
        """
        Provides the energy above convex hull for an entry

        Args:
            entry (PDEntry): A PDEntry like object

        Returns:
            Energy above convex hull of entry. Stable entries should have
            energy above hull of 0. The energy is given per atom.
        """
        return self.get_decomp_and_e_above_hull(entry, **kwargs)[1]

    def get_equilibrium_reaction_energy(self, entry):
        """
        Provides the reaction energy of a stable entry from the neighboring
        equilibrium stable entries (also known as the inverse distance to
        hull).

        Args:
            entry (PDEntry): A PDEntry like object

        Returns:
            Equilibrium reaction energy of entry. Stable entries should have
            equilibrium reaction energy <= 0. The energy is given per atom.
        """
        # NOTE scaled duplicates of stable_entries will not be caught.
        if entry not in self.stable_entries:
            raise ValueError(
                "{} is unstable, the equilibrium reaction energy is" "available only for stable entries.".format(entry)
            )

        if entry.is_element:
            return 0

        entries = [e for e in self.stable_entries if e != entry]
        modpd = PhaseDiagram(entries, self.elements)
        return modpd.get_decomp_and_e_above_hull(entry, allow_negative=True)[1]

    def get_decomp_and_phase_separation_energy(
        self, entry, space_limit=200, stable_only=False, tol=1e-10, maxiter=1000
    ):
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
        optimisation rather than recalculation of the convex hull making it
        algorithmically cheaper. However, if `tol` is too loose there is potential
        for this algorithm to converge to a different solution.

        Args:
            entry (PDEntry): A PDEntry like object.
            space_limit (int): The maximum number of competing entries to consider
                before calculating a second convex hull to reducing the complexity
                of the optimization.
            stable_only (bool): Only use stable materials as competing entries.
            tol (float): The tolerence for convergence of the SLSQP optimization
                when finding the equilibrium reaction.
            maxiter (int): The maximum number of iterations of the SLSQP optimizer
                when finding the equilibrium reaction.

        Returns:
            (decomp, energy). The decompostion  is given as a dict of {PDEntry, amount}
            for all entries in the decomp reaction where amount is the amount of the
            fractional composition. The phase separation energy is given per atom.
        """

        # For unstable or novel materials use simplex approach
        if entry.composition.fractional_composition not in [
            e.composition.fractional_composition for e in self.stable_entries
        ]:
            return self.get_decomp_and_e_above_hull(entry, allow_negative=True)

        # Handle elemental materials
        if entry.is_element:
            return self.get_decomp_and_e_above_hull(entry, allow_negative=True)

        # Select space to compare against
        if stable_only:
            compare_entries = self.stable_entries
        else:
            compare_entries = self.qhull_entries

        # take entries with negative formation enthalpies as competing entries
        competing_entries = [
            c
            for c in compare_entries
            if (c.composition.fractional_composition != entry.composition.fractional_composition)
            if set(c.composition.elements).issubset(entry.composition.elements)
        ]

        # NOTE SLSQP optimizer doesn't scale well for > 300 competing entries. As a
        # result in phase diagrams where we have too many competing entries we can
        # reduce the number by looking at the first and second convex hulls. This
        # requires computing the convex hull of a second (hopefully smallish) space
        # and so is not done by default
        if len(competing_entries) > space_limit and not stable_only:
            inner_hull = PhaseDiagram(
                list(
                    set.intersection(
                        set(competing_entries),  # same chemical space
                        set(self.qhull_entries),  # negative E_f
                        set(self.unstable_entries),  # not already on hull
                    )
                )
                + list(self.el_refs.values())
            )  # terminal points

            competing_entries = list(self.stable_entries.union(inner_hull.stable_entries))
            competing_entries = [c for c in competing_entries if c != entry]

        solution = _get_slsqp_decomp(entry, competing_entries, tol, maxiter)

        if solution.success:
            decomp_amts = solution.x
            decomp = {c: amt for c, amt in zip(competing_entries, decomp_amts) if amt > PhaseDiagram.numerical_tol}

            # find the minimum alternative formation energy for the decomposition
            decomp_enthalpy = np.sum([c.energy_per_atom * amt for c, amt in decomp.items()])

            decomp_enthalpy = entry.energy_per_atom - decomp_enthalpy

            return decomp, decomp_enthalpy

        raise ValueError("No valid decomp found for {}!".format(entry))

    def get_phase_separation_energy(self, entry, **kwargs):
        """
        Provides the energy to the convex hull for the given entry. For stable entries
        already in the phase diagram the algorithm provides the phase separation energy
        which is refered to as the decomposition enthalpy in:

        1. Bartel, C., Trewartha, A., Wang, Q., Dunn, A., Jain, A., Ceder, G.,
            A critical examination of compound stability predictions from
            machine-learned formation energies, npj Computational Materials 6, 97 (2020)

        Args:
            entry (PDEntry): A PDEntry like object
            **kwargs: Keyword args passed to `get_decomp_and_decomp_energy`
                space_limit (int): The maximum number of competing entries to consider.
                stable_only (bool): Only use stable materials as competing entries
                tol (float): The tolerence for convergence of the SLSQP optimization
                    when finding the equilibrium reaction.
                maxiter (int): The maximum number of iterations of the SLSQP optimizer
                    when finding the equilibrium reaction.

        Returns:
            phase separation energy per atom of entry. Stable entries should have
            energies <= 0, Stable elemental entries should have energies = 0 and
            unstable entries should have energies > 0. Entries that have the same
            composition as a stable energy may have postive or negative phase
            separation energies depending on their own energy.
        """
        return self.get_decomp_and_phase_separation_energy(entry, **kwargs)[1]

    def get_composition_chempots(self, comp):
        """
        Get the chemical potentials for all elements at a given composition.

        :param comp: Composition
        :return: Dict of chemical potentials.
        """
        facet = self._get_facet_and_simplex(comp)[0]
        return self._get_facet_chempots(facet)

    def get_all_chempots(self, comp):
        """
        Get chemical potentials at a given compositon.

        :param comp: Composition
        :return: Chemical potentials.
        """
        # NOTE the top part takes from format of _get_facet_and_simplex,
        # but wants to return all facets rather than the first one that
        # meets this criteria
        c = self.pd_coords(comp)
        all_facets = []
        for f, s in zip(self.facets, self.simplexes):
            if s.in_simplex(c, PhaseDiagram.numerical_tol / 10):
                all_facets.append(f)

        if not len(all_facets):
            raise RuntimeError("No facets found for comp = {}".format(comp))

        chempots = {}
        for facet in all_facets:
            facet_name = "-".join([self.qhull_entries[j].name for j in facet])
            chempots[facet_name] = self._get_facet_chempots(facet)
        return chempots

    def get_transition_chempots(self, element):
        """
        Get the critical chemical potentials for an element in the Phase
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
            if len(clean_pots) == 0:
                clean_pots.append(c)
            else:
                if abs(c - clean_pots[-1]) > PhaseDiagram.numerical_tol:
                    clean_pots.append(c)
        clean_pots.reverse()
        return tuple(clean_pots)

    def get_critical_compositions(self, comp1, comp2):
        """
        Get the critical compositions along the tieline between two
        compositions. I.e. where the decomposition products change.
        The endpoints are also returned.
        Args:
            comp1, comp2 (Composition): compositions that define the tieline
        Returns:
            [(Composition)]: list of critical compositions. All are of
                the form x * comp1 + (1-x) * comp2
        """

        n1 = comp1.num_atoms
        n2 = comp2.num_atoms
        pd_els = self.elements

        # the reduced dimensionality Simplexes don't use the
        # first element in the PD
        c1 = self.pd_coords(comp1)
        c2 = self.pd_coords(comp2)

        # none of the projections work if c1 == c2, so just return *copies*
        # of the inputs
        if np.all(c1 == c2):
            return [comp1.copy(), comp2.copy()]

        intersections = [c1, c2]
        for sc in self.simplexes:
            intersections.extend(sc.line_intersection(c1, c2))
        intersections = np.array(intersections)

        # find position along line
        l = c2 - c1
        l /= np.sum(l ** 2) ** 0.5
        proj = np.dot(intersections - c1, l)

        # only take compositions between endpoints
        proj = proj[
            np.logical_and(proj > -self.numerical_tol, proj < proj[1] + self.numerical_tol)  # proj[1] is |c2-c1|
        ]
        proj.sort()

        # only unique compositions
        valid = np.ones(len(proj), dtype=np.bool)
        valid[1:] = proj[1:] > proj[:-1] + self.numerical_tol
        proj = proj[valid]

        ints = c1 + l * proj[:, None]
        # reconstruct full-dimensional composition array
        cs = np.concatenate([np.array([1 - np.sum(ints, axis=-1)]).T, ints], axis=-1)
        # mixing fraction when compositions are normalized
        x = proj / np.dot(c2 - c1, l)
        # mixing fraction when compositions are not normalized
        x_unnormalized = x * n1 / (n2 + x * (n1 - n2))
        num_atoms = n1 + (n2 - n1) * x_unnormalized
        cs *= num_atoms[:, None]
        return [Composition((c, v) for c, v in zip(pd_els, m)) for m in cs]

    def get_element_profile(self, element, comp, comp_tol=1e-5):
        """
        Provides the element evolution data for a composition.
        For example, can be used to analyze Li conversion voltages by varying
        uLi and looking at the phases formed. Also can be used to analyze O2
        evolution by varying uO2.

        Args:
            element: An element. Must be in the phase diagram.
            comp: A Composition
            comp_tol: The tolerance to use when calculating decompositions.
                Phases with amounts less than this tolerance are excluded.
                Defaults to 1e-5.

        Returns:
            Evolution data as a list of dictionaries of the following format:
            [ {'chempot': -10.487582010000001, 'evolution': -2.0,
            'reaction': Reaction Object], ...]
        """
        element = get_el_sp(element)

        if element not in self.elements:
            raise ValueError("get_transition_chempots can only be called with" " elements in the phase diagram.")
        gccomp = Composition({el: amt for el, amt in comp.items() if el != element})
        elref = self.el_refs[element]
        elcomp = Composition(element.symbol)
        evolution = []

        for cc in self.get_critical_compositions(elcomp, gccomp)[1:]:
            decomp_entries = self.get_decomposition(cc).keys()
            decomp = [k.composition for k in decomp_entries]
            rxn = Reaction([comp], decomp + [elcomp])
            rxn.normalize_to(comp)
            c = self.get_composition_chempots(cc + elcomp * 1e-5)[element]
            amt = -rxn.coeffs[rxn.all_comp.index(elcomp)]
            evolution.append(
                {
                    "chempot": c,
                    "evolution": amt,
                    "element_reference": elref,
                    "reaction": rxn,
                    "entries": decomp_entries,
                }
            )
        return evolution

    def get_chempot_range_map(self, elements, referenced=True, joggle=True):
        """
        Returns a chemical potential range map for each stable entry.

        Args:
            elements: Sequence of elements to be considered as independent
                variables. E.g., if you want to show the stability ranges
                of all Li-Co-O phases wrt to uLi and uO, you will supply
                [Element("Li"), Element("O")]
            referenced: If True, gives the results with a reference being the
                energy of the elemental phase. If False, gives absolute values.
            joggle (boolean): Whether to joggle the input to avoid precision
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
            el_energies = {el: 0.0 for el in elements}
        chempot_ranges = collections.defaultdict(list)
        vertices = [list(range(len(self.elements)))]
        if len(all_chempots) > len(self.elements):
            vertices = get_facets(all_chempots, joggle=joggle)
        for ufacet in vertices:
            for combi in itertools.combinations(ufacet, 2):
                data1 = self.facets[combi[0]]
                data2 = self.facets[combi[1]]
                common_ent_ind = set(data1).intersection(set(data2))
                if len(common_ent_ind) == len(elements):
                    common_entries = [self.qhull_entries[i] for i in common_ent_ind]
                    data = np.array([[all_chempots[i][j] - el_energies[self.elements[j]] for j in inds] for i in combi])
                    sim = Simplex(data)
                    for entry in common_entries:
                        chempot_ranges[entry].append(sim)

        return chempot_ranges

    def getmu_vertices_stability_phase(self, target_comp, dep_elt, tol_en=1e-2):
        """
        returns a set of chemical potentials corresponding to the vertices of
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
             [{Element:mu}]: An array of conditions on simplex vertices for
             which each element has a chemical potential set to a given
             value. "absolute" values (i.e., not referenced to element energies)
        """
        muref = np.array([self.el_refs[e].energy_per_atom for e in self.elements if e != dep_elt])
        chempot_ranges = self.get_chempot_range_map([e for e in self.elements if e != dep_elt])

        for e in self.elements:
            if e not in target_comp.elements:
                target_comp = target_comp + Composition({e: 0.0})
        coeff = [-target_comp[e] for e in self.elements if e != dep_elt]
        for e in chempot_ranges.keys():
            if e.composition.reduced_composition == target_comp.reduced_composition:
                multiplicator = e.composition[dep_elt] / target_comp[dep_elt]
                ef = e.energy / multiplicator
                all_coords = []
                for s in chempot_ranges[e]:
                    for v in s._coords:
                        elts = [e for e in self.elements if e != dep_elt]
                        res = {}
                        for i, el in enumerate(elts):
                            res[el] = v[i] + muref[i]
                        res[dep_elt] = (np.dot(v + muref, coeff) + ef) / target_comp[dep_elt]
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

    def get_chempot_range_stability_phase(self, target_comp, open_elt):
        """
        returns a set of chemical potentials corresponding to the max and min
        chemical potential of the open element for a given composition. It is
        quite common to have for instance a ternary oxide (e.g., ABO3) for
        which you want to know what are the A and B chemical potential leading
        to the highest and lowest oxygen chemical potential (reducing and
        oxidizing conditions). This is useful for defect computations.

        Args:
            target_comp: A Composition object
            open_elt: Element that you want to constrain to be max or min

        Returns:
             {Element:(mu_min,mu_max)}: Chemical potentials are given in
             "absolute" values (i.e., not referenced to 0)
        """
        muref = np.array([self.el_refs[e].energy_per_atom for e in self.elements if e != open_elt])
        chempot_ranges = self.get_chempot_range_map([e for e in self.elements if e != open_elt])
        for e in self.elements:
            if e not in target_comp.elements:
                target_comp = target_comp + Composition({e: 0.0})
        coeff = [-target_comp[e] for e in self.elements if e != open_elt]
        max_open = -float("inf")
        min_open = float("inf")
        max_mus = None
        min_mus = None
        for e in chempot_ranges.keys():
            if e.composition.reduced_composition == target_comp.reduced_composition:
                multiplicator = e.composition[open_elt] / target_comp[open_elt]
                ef = e.energy / multiplicator
                all_coords = []
                for s in chempot_ranges[e]:
                    for v in s._coords:
                        all_coords.append(v)
                        test_open = (np.dot(v + muref, coeff) + ef) / target_comp[open_elt]
                        if test_open > max_open:
                            max_open = test_open
                            max_mus = v
                        if test_open < min_open:
                            min_open = test_open
                            min_mus = v
        elts = [e for e in self.elements if e != open_elt]
        res = {}
        for i, el in enumerate(elts):
            res[el] = (min_mus[i] + muref[i], max_mus[i] + muref[i])
        res[open_elt] = (min_open, max_open)
        return res


class PhaseDiagram(BasePhaseDiagram):
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

    .. attribute: elements:

        Elements in the phase diagram.

    ..attribute: all_entries

        All entries provided for Phase Diagram construction. Note that this
        does not mean that all these entries are actually used in the phase
        diagram. For example, this includes the positive formation energy
        entries that are filtered out before Phase Diagram construction.

    .. attribute: qhull_data

        Data used in the convex hull operation. This is essentially a matrix of
        composition data and energy per atom values created from qhull_entries.

    .. attribute: qhull_entries:

        Actual entries used in convex hull. Excludes all positive formation
        energy entries.

    .. attribute: dim

        The dimensionality of the phase diagram.

    .. attribute: facets

        Facets of the phase diagram in the form of  [[1,2,3],[4,5,6]...].
        For a ternary, it is the indices (references to qhull_entries and
        qhull_data) for the vertices of the phase triangles. Similarly
        extended to higher D simplices for higher dimensions.

    .. attribute: el_refs:

        List of elemental references for the phase diagrams. These are
        entries corresponding to the lowest energy element entries for simple
        compositional phase diagrams.

    .. attribute: simplices:

        The simplices of the phase diagram as a list of np.ndarray, i.e.,
        the list of stable compositional coordinates in the phase diagram.
    """

    def __init__(self, entries, elements=None):
        """
        Standard constructor for phase diagram.

        Args:
            entries ([PDEntry]): A list of PDEntry-like objects having an
                energy, energy_per_atom and composition.
            elements ([Element]): Optional list of elements in the phase
                diagram. If set to None, the elements are determined from
                the the entries themselves and are sorted alphabetically.
                If specified, element ordering (e.g. for pd coordinates)
                is preserved.
        """
        super().__init__(**BasePhaseDiagram._kwargs_from_entries(entries, elements))

    def as_dict(self):
        """
        :return: MSONAble dict
        """
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "all_entries": [e.as_dict() for e in self.all_entries],
            "elements": [e.as_dict() for e in self.elements],
        }

    @classmethod
    def from_dict(cls, d):
        """
        :param d: Dict representation
        :return: PhaseDiagram
        """
        entries = [MontyDecoder().process_decoded(dd) for dd in d["all_entries"]]
        elements = [Element.from_dict(dd) for dd in d["elements"]]
        return cls(entries, elements)


class GrandPotentialPhaseDiagram(PhaseDiagram):
    """
    A class representing a Grand potential phase diagram. Grand potential phase
    diagrams are essentially phase diagrams that are open to one or more
    components. To construct such phase diagrams, the relevant free energy is
    the grand potential, which can be written as the Legendre transform of the
    Gibbs free energy as follows

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

    def __init__(self, entries, chempots, elements=None):
        """
        Standard constructor for grand potential phase diagram.

        Args:
            entries ([PDEntry]): A list of PDEntry-like objects having an
                energy, energy_per_atom and composition.
            chempots ({Element: float}): Specify the chemical potentials
                of the open elements.
            elements ([Element]): Optional list of elements in the phase
                diagram. If set to None, the elements are determined from
                the the entries themselves.
        """
        if elements is None:
            elements = {els for e in entries for els in e.composition.elements}

        self.chempots = {get_el_sp(el): u for el, u in chempots.items()}
        elements = set(elements).difference(self.chempots.keys())

        all_entries = [
            GrandPotPDEntry(e, self.chempots) for e in entries if len(elements.intersection(e.composition.elements)) > 0
        ]

        super().__init__(all_entries, elements)

    def __repr__(self):
        chemsys = "-".join([el.symbol for el in self.elements])
        chempots = ", ".join(["u{}={}".format(el, v) for el, v in self.chempots.items()])

        output = [
            "{} grand potential phase diagram with {}".format(chemsys, chempots),
            "{} stable phases: ".format(len(self.stable_entries)),
            ", ".join([entry.name for entry in self.stable_entries]),
        ]
        return "\n".join(output)

    def as_dict(self):
        """
        :return: MSONable dict
        """
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "all_entries": [e.as_dict() for e in self.all_entries],
            "chempots": self.chempots,
            "elements": [e.as_dict() for e in self.elements],
        }

    @classmethod
    def from_dict(cls, d):
        """
        :param d: Dict representation
        :return: GrandPotentialPhaseDiagram
        """
        entries = MontyDecoder().process_decoded(d["all_entries"])
        elements = MontyDecoder().process_decoded(d["elements"])
        return cls(entries, d["chempots"], elements)


class CompoundPhaseDiagram(PhaseDiagram):
    """
    Generates phase diagrams from compounds as terminations instead of
    elements.
    """

    # Tolerance for determining if amount of a composition is positive.
    amount_tol = 1e-5

    def __init__(self, entries, terminal_compositions, normalize_terminal_compositions=True):
        """
        Initializes a CompoundPhaseDiagram.

        Args:
            entries ([PDEntry]): Sequence of input entries. For example,
               if you want a Li2O-P2O5 phase diagram, you might have all
               Li-P-O entries as an input.
            terminal_compositions ([Composition]): Terminal compositions of
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
        (pentries, species_mapping) = self.transform_entries(entries, terminal_compositions)
        self.species_mapping = species_mapping
        super().__init__(pentries, elements=species_mapping.values())

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
        sp_mapping = collections.OrderedDict()
        for i, comp in enumerate(terminal_compositions):
            sp_mapping[comp] = DummySpecies("X" + chr(102 + i))

        for entry in entries:
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
        """
        :return: MSONable dict
        """
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "original_entries": [e.as_dict() for e in self.original_entries],
            "terminal_compositions": [c.as_dict() for c in self.terminal_compositions],
            "normalize_terminal_compositions": self.normalize_terminals,
        }

    @classmethod
    def from_dict(cls, d):
        """
        :param d: Dict Representation
        :return: CompoundPhaseDiagram
        """
        dec = MontyDecoder()
        entries = dec.process_decoded(d["original_entries"])
        terminal_compositions = dec.process_decoded(d["terminal_compositions"])
        return cls(entries, terminal_compositions, d["normalize_terminal_compositions"])


class ReactionDiagram:
    """
    Analyzes the possible reactions between a pair of compounds, e.g.,
    an electrolyte and an electrode.
    """

    def __init__(self, entry1, entry2, all_entries, tol=1e-4, float_fmt="%.4f"):
        """
        Args:
            entry1 (ComputedEntry): Entry for 1st component. Note that
                corrections, if any, must already be pre-applied. This is to
                give flexibility for different kinds of corrections, e.g.,
                if a particular entry is fitted to an experimental data (such
                as EC molecule).
            entry2 (ComputedEntry): Entry for 2nd component. Note that
                corrections must already be pre-applied. This is to
                give flexibility for different kinds of corrections, e.g.,
                if a particular entry is fitted to an experimental data (such
                as EC molecule).
            all_entries ([ComputedEntry]): All other entries to be
                considered in the analysis. Note that corrections, if any,
                must already be pre-applied.
            tol (float): Tolerance to be used to determine validity of reaction.
            float_fmt (str): Formatting string to be applied to all floats.
                Determines number of decimal places in reaction string.
        """
        elements = set()
        for e in [entry1, entry2]:
            elements.update([el.symbol for el in e.composition.elements])

        elements = tuple(elements)  # Fix elements to ensure order.

        comp_vec1 = np.array([entry1.composition.get_atomic_fraction(el) for el in elements])
        comp_vec2 = np.array([entry2.composition.get_atomic_fraction(el) for el in elements])
        r1 = entry1.composition.reduced_composition
        r2 = entry2.composition.reduced_composition

        logger.debug("%d total entries." % len(all_entries))

        pd = PhaseDiagram(all_entries + [entry1, entry2])
        terminal_formulas = [
            entry1.composition.reduced_formula,
            entry2.composition.reduced_formula,
        ]

        logger.debug("%d stable entries" % len(pd.stable_entries))
        logger.debug("%d facets" % len(pd.facets))
        logger.debug("%d qhull_entries" % len(pd.qhull_entries))

        rxn_entries = []
        done = []

        def fmt(fl):
            return float_fmt % fl

        for facet in pd.facets:
            for face in itertools.combinations(facet, len(facet) - 1):
                face_entries = [pd.qhull_entries[i] for i in face]

                if any(e.composition.reduced_formula in terminal_formulas for e in face_entries):
                    continue

                try:

                    m = []
                    for e in face_entries:
                        m.append([e.composition.get_atomic_fraction(el) for el in elements])
                    m.append(comp_vec2 - comp_vec1)
                    m = np.array(m).T
                    coeffs = np.linalg.solve(m, comp_vec2)

                    x = coeffs[-1]
                    # pylint: disable=R1716
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

                        rxn_str = "%s %s + %s %s -> " % (
                            fmt(c1),
                            r1.reduced_formula,
                            fmt(c2),
                            r2.reduced_formula,
                        )
                        products = []
                        product_entries = []

                        energy = -(x * entry1.energy_per_atom + (1 - x) * entry2.energy_per_atom)

                        for c, e in zip(coeffs[:-1], face_entries):
                            if c > tol:
                                r = e.composition.reduced_composition
                                products.append("%s %s" % (fmt(c / r.num_atoms * factor), r.reduced_formula))
                                product_entries.append((c, e))
                                energy += c * e.energy_per_atom

                        rxn_str += " + ".join(products)
                        comp = x * comp_vec1 + (1 - x) * comp_vec2
                        entry = PDEntry(
                            Composition(dict(zip(elements, comp))),
                            energy=energy,
                            attribute=rxn_str,
                        )
                        entry.decomposition = product_entries
                        rxn_entries.append(entry)
                except np.linalg.LinAlgError:
                    logger.debug(
                        "Reactants = %s"
                        % (
                            ", ".join(
                                [
                                    entry1.composition.reduced_formula,
                                    entry2.composition.reduced_formula,
                                ]
                            )
                        )
                    )
                    logger.debug("Products = %s" % (", ".join([e.composition.reduced_formula for e in face_entries])))

        rxn_entries = sorted(rxn_entries, key=lambda e: e.name, reverse=True)

        self.entry1 = entry1
        self.entry2 = entry2
        self.rxn_entries = rxn_entries
        self.labels = collections.OrderedDict()
        for i, e in enumerate(rxn_entries):
            self.labels[str(i + 1)] = e.attribute
            e.name = str(i + 1)
        self.all_entries = all_entries
        self.pd = pd

    def get_compound_pd(self):
        """
        Get the CompoundPhaseDiagram object, which can then be used for
        plotting.

        Returns:
            (CompoundPhaseDiagram)
        """
        # For this plot, since the reactions are reported in formation
        # energies, we need to set the energies of the terminal compositions
        # to 0. So we make create copies with 0 energy.
        entry1 = PDEntry(self.entry1.composition, 0)
        entry2 = PDEntry(self.entry2.composition, 0)

        cpd = CompoundPhaseDiagram(
            self.rxn_entries + [entry1, entry2],
            [
                Composition(entry1.composition.reduced_formula),
                Composition(entry2.composition.reduced_formula),
            ],
            normalize_terminal_compositions=False,
        )
        return cpd


class PhaseDiagramError(Exception):
    """
    An exception class for Phase Diagram generation.
    """

    pass


def get_facets(qhull_data, joggle=False):
    """
    Get the simplex facets for the Convex hull.

    Args:
        qhull_data (np.ndarray): The data from which to construct the convex
            hull as a Nxd array (N being number of data points and d being the
            dimension)
        joggle (boolean): Whether to joggle the input to avoid precision
            errors.

    Returns:
        List of simplices of the Convex Hull.
    """
    if joggle:
        return ConvexHull(qhull_data, qhull_options="QJ i").simplices
    return ConvexHull(qhull_data, qhull_options="Qt i").simplices


def _get_slsqp_decomp(comp, competing_entries, tol=1e-10, maxiter=1000):
    """
    Finds the amounts of competing compositions that minimize the energy of a
    given composition

    The algorithm is based on the work in the following paper:

    1. Bartel, C., Trewartha, A., Wang, Q., Dunn, A., Jain, A., Ceder, G.,
        A critical examination of compound stability predictions from
        machine-learned formation energies, npj Computational Materials 6, 97 (2020)

    Args:
        entry (PDEntry): A PDEntry like entry to analyze
        competing_entries ([PDEntry]): List of entries to consider for decomposition

    Returns:
        scipy.optimize.minimize result. If sucessful this gives the linear combination of
            competing entrys that minimizes the competing formation energy
    """
    if not isinstance(comp, Composition):
        comp = comp.composition

    # Elemental amount present in given entry
    amts = comp.fractional_composition.get_el_amt_dict()
    chemical_space = tuple(amts.keys())
    b = np.array([amts[el] for el in chemical_space])

    # Elemental amounts present in competing entries
    A_transpose = np.zeros((len(chemical_space), len(competing_entries)))
    for j, comp_entry in enumerate(competing_entries):
        amts = comp_entry.composition.fractional_composition.get_el_amt_dict()
        for i, el in enumerate(chemical_space):
            A_transpose[i, j] = amts[el]

    # Energies of competing entries
    Es = np.array([comp_entry.energy_per_atom for comp_entry in competing_entries])

    molar_constraint = {
        "type": "eq",
        "fun": lambda x: np.dot(A_transpose, x) - b,
        "jac": lambda x: A_transpose,
    }

    options = {"maxiter": maxiter, "disp": False}

    max_bound = comp.num_atoms
    bounds = [(0, max_bound)] * len(competing_entries)
    x0 = [1 / len(competing_entries)] * len(competing_entries)

    # NOTE the tolerence needs to be tight to stop the optimization
    # from exiting before convergence is reached. Issues observed for
    # tol > 1e-7 in the fractional composition (default 1e-10).
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

    return solution


class PDPlotter:
    """
    A plotter class for compositional phase diagrams.
    """

    def __init__(
        self,
        phasediagram: PhaseDiagram,
        show_unstable: float = 0.2,
        backend: str = "plotly",
        **plotkwargs,
    ):
        """
        Args:
            phasediagram (PhaseDiagram): PhaseDiagram object.
            show_unstable (float): Whether unstable (above the hull) phases will be
                plotted. If a number > 0 is entered, all phases with
                e_hull < show_unstable (eV/atom) will be shown.
            backend (str): Python package used for plotting ("matplotlib" or
                "plotly"). Defaults to "plotly".
            **plotkwargs (dict): Keyword args passed to matplotlib.pyplot.plot. Can
                be used to customize markers etc. If not set, the default is
                {
                    "markerfacecolor": (0.2157, 0.4941, 0.7216),
                    "markersize": 10,
                    "linewidth": 3
                }
        """
        # note: palettable imports matplotlib
        from palettable.colorbrewer.qualitative import Set1_3

        self._pd = phasediagram
        self._dim = len(self._pd.elements)
        if self._dim > 4:
            raise ValueError("Only 1-4 components supported!")
        self.lines = uniquelines(self._pd.facets) if self._dim > 1 else [[self._pd.facets[0][0], self._pd.facets[0][0]]]
        self.show_unstable = show_unstable
        self.backend = backend
        self._min_energy = min([self._pd.get_form_energy_per_atom(e) for e in self._pd.stable_entries])
        colors = Set1_3.mpl_colors
        self.plotkwargs = plotkwargs or {
            "markerfacecolor": colors[2],
            "markersize": 10,
            "linewidth": 3,
        }

    @property  # type: ignore
    @lru_cache(1)
    def pd_plot_data(self):
        """
        Plotting data for phase diagram. Cached for repetitive calls.
        2-comp - Full hull with energies
        3/4-comp - Projection into 2D or 3D Gibbs triangle.

        Returns:
            (lines, stable_entries, unstable_entries):
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
            labelcoord = list(zip(*coord))
            stable_entries[labelcoord[0]] = entry1
            stable_entries[labelcoord[1]] = entry2

        all_entries = pd.all_entries
        all_data = np.array(pd.all_entries_hulldata)
        unstable_entries = dict()
        stable = pd.stable_entries
        for i, entry in enumerate(all_entries):
            if entry not in stable:
                if self._dim < 3:
                    x = [all_data[i][0], all_data[i][0]]
                    y = [
                        pd.get_form_energy_per_atom(entry),
                        pd.get_form_energy_per_atom(entry),
                    ]
                    coord = [x, y]
                elif self._dim == 3:
                    coord = triangular_coord([all_data[i, 0:2], all_data[i, 0:2]])
                else:
                    coord = tet_coord([all_data[i, 0:3], all_data[i, 0:3], all_data[i, 0:3]])
                labelcoord = list(zip(*coord))
                unstable_entries[entry] = labelcoord[0]

        return lines, stable_entries, unstable_entries

    def get_plot(
        self,
        label_stable=True,
        label_unstable=True,
        ordering=None,
        energy_colormap=None,
        process_attributes=False,
        plt=None,
        label_uncertainties=False,
    ):
        """
        :param label_stable: Whether to label stable compounds.
        :param label_unstable: Whether to label unstable compounds.
        :param ordering: Ordering of vertices (matplotlib backend only).
        :param energy_colormap: Colormap for coloring energy (matplotlib backend only).
        :param process_attributes: Whether to process the attributes (matplotlib
            backend only).
        :param plt: Existing plt object if plotting multiple phase diagrams (
            matplotlib backend only).
        :param label_uncertainties: Whether to add error bars to the hull (plotly
            backend only). For binaries, this also shades the hull with the
            uncertainty window.
        :return: go.Figure (plotly) or matplotlib.pyplot (matplotlib)
        """
        fig = None

        if self.backend == "plotly":
            data = [self._create_plotly_lines()]

            if self._dim == 3:
                data.append(self._create_plotly_ternary_support_lines())
                data.append(self._create_plotly_ternary_hull())

            stable_labels_plot = self._create_plotly_stable_labels(label_stable)
            stable_marker_plot, unstable_marker_plot = self._create_plotly_markers(label_uncertainties)

            if self._dim == 2 and label_uncertainties:
                data.append(self._create_plotly_uncertainty_shading(stable_marker_plot))

            data.append(stable_labels_plot)
            data.append(unstable_marker_plot)
            data.append(stable_marker_plot)

            fig = go.Figure(data=data)
            fig.layout = self._create_plotly_figure_layout()

        elif self.backend == "matplotlib":
            if self._dim <= 3:
                fig = self._get_2d_plot(
                    label_stable,
                    label_unstable,
                    ordering,
                    energy_colormap,
                    plt=plt,
                    process_attributes=process_attributes,
                )
            elif self._dim == 4:
                fig = self._get_3d_plot(label_stable)

        return fig

    def plot_element_profile(self, element, comp, show_label_index=None, xlim=5):
        """
        Draw the element profile plot for a composition varying different
        chemical potential of an element.
        X value is the negative value of the chemical potential reference to
        elemental chemical potential. For example, if choose Element("Li"),
        X= -(µLi-µLi0), which corresponds to the voltage versus metal anode.
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
        plt = pretty_plot(12, 8)
        pd = self._pd
        evolution = pd.get_element_profile(element, comp)
        num_atoms = evolution[0]["reaction"].reactants[0].num_atoms
        element_energy = evolution[0]["chempot"]
        x1, x2, y1 = None, None, None
        for i, d in enumerate(evolution):
            v = -(d["chempot"] - element_energy)
            if i != 0:
                plt.plot([x2, x2], [y1, d["evolution"] / num_atoms], "k", linewidth=2.5)
            x1 = v
            y1 = d["evolution"] / num_atoms

            if i != len(evolution) - 1:
                x2 = -(evolution[i + 1]["chempot"] - element_energy)
            else:
                x2 = 5.0
            if show_label_index is not None and i in show_label_index:
                products = [
                    re.sub(r"(\d+)", r"$_{\1}$", p.reduced_formula)
                    for p in d["reaction"].products
                    if p.reduced_formula != element.symbol
                ]
                plt.annotate(
                    ", ".join(products),
                    xy=(v + 0.05, y1 + 0.05),
                    fontsize=24,
                    color="r",
                )
                plt.plot([x1, x2], [y1, y1], "r", linewidth=3)
            else:
                plt.plot([x1, x2], [y1, y1], "k", linewidth=2.5)

        plt.xlim((0, xlim))
        plt.xlabel("-$\\Delta{\\mu}$ (eV)")
        plt.ylabel("Uptake per atom")

        return plt

    def show(self, *args, **kwargs):
        r"""
        Draw the phase diagram using Plotly (or Matplotlib) and show it.

        Args:
            *args: Passed to get_plot.
            **kwargs: Passed to get_plot.
        """
        self.get_plot(*args, **kwargs).show()

    def _get_2d_plot(
        self,
        label_stable=True,
        label_unstable=True,
        ordering=None,
        energy_colormap=None,
        vmin_mev=-60.0,
        vmax_mev=60.0,
        show_colorbar=True,
        process_attributes=False,
        plt=None,
    ):
        """
        Shows the plot using pylab. Contains import statements since matplotlib is a
        fairly extensive library to load.
        """
        if plt is None:
            plt = pretty_plot(8, 6)
        from matplotlib.font_manager import FontProperties

        if ordering is None:
            (lines, labels, unstable) = self.pd_plot_data
        else:
            (_lines, _labels, _unstable) = self.pd_plot_data
            (lines, labels, unstable) = order_phase_diagram(_lines, _labels, _unstable, ordering)
        if energy_colormap is None:
            if process_attributes:
                for x, y in lines:
                    plt.plot(x, y, "k-", linewidth=3, markeredgecolor="k")
                # One should think about a clever way to have "complex"
                # attributes with complex processing options but with a clear
                # logic. At this moment, I just use the attributes to know
                # whether an entry is a new compound or an existing (from the
                #  ICSD or from the MP) one.
                for x, y in labels.keys():
                    if labels[(x, y)].attribute is None or labels[(x, y)].attribute == "existing":
                        plt.plot(x, y, "ko", **self.plotkwargs)
                    else:
                        plt.plot(x, y, "k*", **self.plotkwargs)
            else:
                for x, y in lines:
                    plt.plot(x, y, "ko-", **self.plotkwargs)
        else:
            from matplotlib.cm import ScalarMappable
            from matplotlib.colors import LinearSegmentedColormap, Normalize

            for x, y in lines:
                plt.plot(x, y, "k-", markeredgecolor="k")
            vmin = vmin_mev / 1000.0
            vmax = vmax_mev / 1000.0
            if energy_colormap == "default":
                mid = -vmin / (vmax - vmin)
                cmap = LinearSegmentedColormap.from_list(
                    "my_colormap",
                    [
                        (0.0, "#005500"),
                        (mid, "#55FF55"),
                        (mid, "#FFAAAA"),
                        (1.0, "#FF0000"),
                    ],
                )
            else:
                cmap = energy_colormap
            norm = Normalize(vmin=vmin, vmax=vmax)
            _map = ScalarMappable(norm=norm, cmap=cmap)
            _energies = [self._pd.get_equilibrium_reaction_energy(entry) for coord, entry in labels.items()]
            energies = [en if en < 0.0 else -0.00000001 for en in _energies]
            vals_stable = _map.to_rgba(energies)
            ii = 0
            if process_attributes:
                for x, y in labels.keys():
                    if labels[(x, y)].attribute is None or labels[(x, y)].attribute == "existing":
                        plt.plot(x, y, "o", markerfacecolor=vals_stable[ii], markersize=12)
                    else:
                        plt.plot(x, y, "*", markerfacecolor=vals_stable[ii], markersize=18)
                    ii += 1
            else:
                for x, y in labels.keys():
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
            all_coords = labels.keys()
            miny = min([c[1] for c in all_coords])
            ybuffer = max(abs(miny) * 0.1, 0.1)
            plt.xlim((-0.1, 1.1))
            plt.ylim((miny - ybuffer, ybuffer))
            center = (0.5, miny / 2)
            plt.xlabel("Fraction", fontsize=28, fontweight="bold")
            plt.ylabel("Formation energy (eV/atom)", fontsize=28, fontweight="bold")

        for coords in sorted(labels.keys(), key=lambda x: -x[1]):
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
                "Energy [meV/at] above hull (in red)\nInverse energy [" "meV/at] above hull (in green)",
                rotation=-90,
                ha="left",
                va="center",
            )
        f = plt.gcf()
        f.set_size_inches((8, 6))
        plt.subplots_adjust(left=0.09, right=0.98, top=0.98, bottom=0.07)
        return plt

    def _get_3d_plot(self, label_stable=True):
        """
        Shows the plot using pylab.  Usually I won"t do imports in methods,
        but since plotting is a fairly expensive library to load and not all
        machines have matplotlib installed, I have done it this way.
        """
        import matplotlib.pyplot as plt
        import mpl_toolkits.mplot3d.axes3d as p3
        from matplotlib.font_manager import FontProperties

        fig = plt.figure()
        ax = p3.Axes3D(fig)
        font = FontProperties(weight="bold", size=13)
        (lines, labels, unstable) = self.pd_plot_data
        count = 1
        newlabels = list()
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
        for coords in sorted(labels.keys()):
            entry = labels[coords]
            label = entry.name
            if label_stable:
                if len(entry.composition.elements) == 1:
                    ax.text(coords[0], coords[1], coords[2], label, fontproperties=font)
                else:
                    ax.text(coords[0], coords[1], coords[2], str(count), fontsize=12)
                    newlabels.append("{} : {}".format(count, latexify(label)))
                    count += 1
        plt.figtext(0.01, 0.01, "\n".join(newlabels), fontproperties=font)
        ax.axis("off")
        ax.set_xlim(-0.1, 0.72)
        ax.set_ylim(0, 0.66)
        ax.set_zlim(0, 0.56)  # pylint: disable=E1101
        return plt

    def write_image(self, stream, image_format="svg", **kwargs):
        r"""
        Writes the phase diagram to an image in a stream.

        Args:
            stream:
                stream to write to. Can be a file stream or a StringIO stream.
            image_format
                format for image. Can be any of matplotlib supported formats.
                Defaults to svg for best results for vector graphics.
            **kwargs: Pass through to get_plot functino.
        """
        plt = self.get_plot(**kwargs)

        f = plt.gcf()
        f.set_size_inches((12, 10))

        plt.savefig(stream, format=image_format)

    def plot_chempot_range_map(self, elements, referenced=True):
        """
        Plot the chemical potential range _map. Currently works only for
        3-component PDs.

        Args:
            elements: Sequence of elements to be considered as independent
                variables. E.g., if you want to show the stability ranges of
                all Li-Co-O phases wrt to uLi and uO, you will supply
                [Element("Li"), Element("O")]
            referenced: if True, gives the results with a reference being the
                        energy of the elemental phase. If False, gives absolute values.
        """
        self.get_chempot_range_map_plot(elements, referenced=referenced).show()

    def get_chempot_range_map_plot(self, elements, referenced=True):
        """
        Returns a plot of the chemical potential range _map. Currently works
        only for 3-component PDs.

        Args:
            elements: Sequence of elements to be considered as independent
                variables. E.g., if you want to show the stability ranges of
                all Li-Co-O phases wrt to uLi and uO, you will supply
                [Element("Li"), Element("O")]
            referenced: if True, gives the results with a reference being the
                        energy of the elemental phase. If False, gives absolute values.

        Returns:
            A matplotlib plot object.
        """

        plt = pretty_plot(12, 8)
        chempot_ranges = self._pd.get_chempot_range_map(elements, referenced=referenced)
        missing_lines = {}
        excluded_region = []
        for entry, lines in chempot_ranges.items():
            comp = entry.composition
            center_x = 0
            center_y = 0
            coords = []
            contain_zero = any(comp.get_atomic_fraction(el) == 0 for el in elements)
            is_boundary = (not contain_zero) and sum([comp.get_atomic_fraction(el) for el in elements]) == 1
            for line in lines:
                (x, y) = line.coords.transpose()
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

        ax = plt.gca()
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()

        # Shade the forbidden chemical potential regions.
        excluded_region.append([xlim[1], ylim[1]])
        excluded_region = sorted(excluded_region, key=lambda c: c[0])
        (x, y) = np.transpose(excluded_region)
        plt.fill(x, y, "0.80")

        # The hull does not generate the missing horizontal and vertical lines.
        # The following code fixes this.
        el0 = elements[0]
        el1 = elements[1]
        for entry, coords in missing_lines.items():
            center_x = sum([c[0] for c in coords])
            center_y = sum([c[1] for c in coords])
            comp = entry.composition
            is_x = comp.get_atomic_fraction(el0) < 0.01
            is_y = comp.get_atomic_fraction(el1) < 0.01
            n = len(coords)
            if not (is_x and is_y):
                if is_x:
                    coords = sorted(coords, key=lambda c: c[1])
                    for i in [0, -1]:
                        x = [min(xlim), coords[i][0]]
                        y = [coords[i][1], coords[i][1]]
                        plt.plot(x, y, "k")
                        center_x += min(xlim)
                        center_y += coords[i][1]
                elif is_y:
                    coords = sorted(coords, key=lambda c: c[0])
                    for i in [0, -1]:
                        x = [coords[i][0], coords[i][0]]
                        y = [coords[i][1], min(ylim)]
                        plt.plot(x, y, "k")
                        center_x += coords[i][0]
                        center_y += min(ylim)
                xy = (center_x / (n + 2), center_y / (n + 2))
            else:
                center_x = sum(coord[0] for coord in coords) + xlim[0]
                center_y = sum(coord[1] for coord in coords) + ylim[0]
                xy = (center_x / (n + 1), center_y / (n + 1))

            plt.annotate(
                latexify(entry.name),
                xy,
                horizontalalignment="center",
                verticalalignment="center",
                fontsize=22,
            )

        plt.xlabel("$\\mu_{{{0}}} - \\mu_{{{0}}}^0$ (eV)".format(el0.symbol))
        plt.ylabel("$\\mu_{{{0}}} - \\mu_{{{0}}}^0$ (eV)".format(el1.symbol))
        plt.tight_layout()
        return plt

    def get_contour_pd_plot(self):
        """
        Plot a contour phase diagram plot, where phase triangles are colored
        according to degree of instability by interpolation. Currently only
        works for 3-component phase diagrams.

        Returns:
            A matplotlib plot object.
        """
        from matplotlib import cm
        from scipy import interpolate

        pd = self._pd
        entries = pd.qhull_entries
        data = np.array(pd.qhull_data)

        plt = self._get_2d_plot()
        data[:, 0:2] = triangular_coord(data[:, 0:2]).transpose()
        for i, e in enumerate(entries):
            data[i, 2] = self._pd.get_e_above_hull(e)

        gridsize = 0.005
        xnew = np.arange(0, 1.0, gridsize)
        ynew = np.arange(0, 1, gridsize)

        f = interpolate.LinearNDInterpolator(data[:, 0:2], data[:, 2])
        znew = np.zeros((len(ynew), len(xnew)))
        for (i, xval) in enumerate(xnew):
            for (j, yval) in enumerate(ynew):
                znew[j, i] = f(xval, yval)

        # pylint: disable=E1101
        plt.contourf(xnew, ynew, znew, 1000, cmap=cm.autumn_r)

        plt.colorbar()
        return plt

    def _create_plotly_lines(self):
        """
        Creates Plotly scatter (line) plots for all phase diagram facets.

        :return: go.Scatter (or go.Scatter3d) plot
        """
        line_plot = None
        x, y, z, energies = [], [], [], []

        for line in self.pd_plot_data[0]:
            x.extend(list(line[0]) + [None])
            y.extend(list(line[1]) + [None])

            if self._dim == 3:
                z.extend(
                    [self._pd.get_form_energy_per_atom(self.pd_plot_data[1][coord]) for coord in zip(line[0], line[1])]
                    + [None]
                )

            elif self._dim == 4:
                energies.extend(
                    [
                        self._pd.get_form_energy_per_atom(self.pd_plot_data[1][coord])
                        for coord in zip(line[0], line[1], line[2])
                    ]
                    + [None]
                )
                z.extend(list(line[2]) + [None])

        plot_args = dict(
            mode="lines",
            hoverinfo="none",
            line={"color": "rgba(0,0,0,1.0)", "width": 7.0},
            showlegend=False,
        )

        if self._dim == 2:
            line_plot = go.Scatter(x=x, y=y, **plot_args)
        elif self._dim == 3:
            line_plot = go.Scatter3d(x=y, y=x, z=z, **plot_args)
        elif self._dim == 4:
            line_plot = go.Scatter3d(x=x, y=y, z=z, **plot_args)

        return line_plot

    def _create_plotly_stable_labels(self, label_stable=True):
        """
        Creates a (hidable) scatter trace containing labels of stable phases.
        Contains some functionality for creating sensible label positions.

        :return: go.Scatter (or go.Scatter3d) plot
        """
        x, y, z, text, textpositions = [], [], [], [], []
        stable_labels_plot = None
        min_energy_x = None
        offset_2d = 0.005  # extra distance to offset label position for clarity
        offset_3d = 0.01

        energy_offset = -0.1 * self._min_energy

        if self._dim == 2:
            min_energy_x = min(list(self.pd_plot_data[1].keys()), key=lambda c: c[1])[0]

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
                y_coord -= offset_2d
            elif self._dim == 3:
                textposition = "middle center"
                if coords[0] > 0.5:
                    x_coord += offset_3d
                else:
                    x_coord -= offset_3d
                if coords[1] > 0.866 / 2:
                    y_coord -= offset_3d
                else:
                    y_coord += offset_3d

                z.append(self._pd.get_form_energy_per_atom(entry) + energy_offset)

            elif self._dim == 4:
                x_coord = x_coord - offset_3d
                y_coord = y_coord - offset_3d
                textposition = "bottom right"
                z.append(coords[2])

            x.append(x_coord)
            y.append(y_coord)
            textpositions.append(textposition)

            comp = entry.composition
            if hasattr(entry, "original_entry"):
                comp = entry.original_entry.composition

            formula = list(comp.reduced_formula)
            text.append(self._htmlize_formula(formula))

        visible = True
        if not label_stable or self._dim == 4:
            visible = "legendonly"

        plot_args = dict(
            text=text,
            textposition=textpositions,
            mode="text",
            name="Labels (stable)",
            hoverinfo="skip",
            opacity=1.0,
            visible=visible,
            showlegend=True,
        )

        if self._dim == 2:
            stable_labels_plot = go.Scatter(x=x, y=y, **plot_args)
        elif self._dim == 3:
            stable_labels_plot = go.Scatter3d(x=y, y=x, z=z, **plot_args)
        elif self._dim == 4:
            stable_labels_plot = go.Scatter3d(x=x, y=y, z=z, **plot_args)

        return stable_labels_plot

    def _create_plotly_element_annotations(self):
        """
        Creates terminal element annotations for Plotly phase diagrams.

        :return: list of annotation dicts.
        """
        annotations_list = []
        x, y, z = None, None, None

        for coords, entry in self.pd_plot_data[1].items():
            if not entry.composition.is_element:
                continue

            x, y = coords[0], coords[1]

            if self._dim == 3:
                z = self._pd.get_form_energy_per_atom(entry)
            elif self._dim == 4:
                z = coords[2]

            if entry.composition.is_element:
                clean_formula = str(entry.composition.elements[0])
                if hasattr(entry, "original_entry"):
                    orig_comp = entry.original_entry.composition
                    clean_formula = self._htmlize_formula(orig_comp.reduced_formula)

                font_dict = {"color": "#000000", "size": 24.0}
                opacity = 1.0

            annotation = plotly_layouts["default_annotation_layout"].copy()
            annotation.update(
                {
                    "x": x,
                    "y": y,
                    "font": font_dict,
                    "text": clean_formula,
                    "opacity": opacity,
                }
            )

            if self._dim == 3 or self._dim == 4:
                for d in ["xref", "yref"]:
                    annotation.pop(d)  # Scatter3d cannot contain xref, yref
                    if self._dim == 3:
                        annotation.update({"x": y, "y": x})
                        if entry.composition.is_element:
                            z = 0.9 * self._min_energy  # place label 10% above base

                annotation.update({"z": z})

            annotations_list.append(annotation)

        # extra point ensures equilateral triangular scaling is displayed
        if self._dim == 3:
            annotations_list.append(dict(x=1, y=1, z=0, opacity=0, text=""))

        return annotations_list

    def _create_plotly_figure_layout(self, label_stable=True):
        """
        Creates layout for plotly phase diagram figure and updates with
        figure annotations.

        :return: Dictionary with Plotly figure layout settings.
        """
        annotations_list = None
        layout = dict()

        if label_stable:
            annotations_list = self._create_plotly_element_annotations()

        if self._dim == 2:
            layout = plotly_layouts["default_binary_layout"].copy()
            layout["annotations"] = annotations_list
        elif self._dim == 3:
            layout = plotly_layouts["default_ternary_layout"].copy()
            layout["scene"].update({"annotations": annotations_list})
        elif self._dim == 4:
            layout = plotly_layouts["default_quaternary_layout"].copy()
            layout["scene"].update({"annotations": annotations_list})

        return layout

    def _create_plotly_markers(self, label_uncertainties=False):
        """
        Creates stable and unstable marker plots for overlaying on the phase diagram.

        :return: Tuple of Plotly go.Scatter (or go.Scatter3d) objects in order: (
            stable markers, unstable markers)
        """

        def get_marker_props(coords, entries, stable=True):
            """Method for getting marker locations, hovertext, and error bars
            from pd_plot_data"""
            x, y, z, texts, energies, uncertainties = [], [], [], [], [], []

            for coord, entry in zip(coords, entries):
                energy = round(self._pd.get_form_energy_per_atom(entry), 3)

                entry_id = getattr(entry, "entry_id", "no ID")
                comp = entry.composition

                if hasattr(entry, "original_entry"):
                    comp = entry.original_entry.composition

                formula = comp.reduced_formula
                clean_formula = self._htmlize_formula(formula)
                label = f"{clean_formula} ({entry_id}) <br> " f"{energy} eV/atom"

                if not stable:
                    e_above_hull = round(self._pd.get_e_above_hull(entry), 3)
                    if e_above_hull > self.show_unstable:
                        continue
                    label += f" (+{e_above_hull} eV/atom)"
                    energies.append(e_above_hull)
                else:
                    uncertainty = 0
                    if hasattr(entry, "correction_uncertainty_per_atom") and label_uncertainties:
                        uncertainty = round(entry.correction_uncertainty_per_atom, 4)
                        label += f"<br> (Error: +/- {uncertainty} eV/atom)"

                    uncertainties.append(uncertainty)
                    energies.append(energy)

                texts.append(label)

                x.append(coord[0])
                y.append(coord[1])

                if self._dim == 3:
                    z.append(energy)
                elif self._dim == 4:
                    z.append(coord[2])

            return {
                "x": x,
                "y": y,
                "z": z,
                "texts": texts,
                "energies": energies,
                "uncertainties": uncertainties,
            }

        stable_coords, stable_entries = (
            self.pd_plot_data[1].keys(),
            self.pd_plot_data[1].values(),
        )
        unstable_entries, unstable_coords = (
            self.pd_plot_data[2].keys(),
            self.pd_plot_data[2].values(),
        )

        stable_props = get_marker_props(stable_coords, stable_entries)

        unstable_props = get_marker_props(unstable_coords, unstable_entries, stable=False)

        stable_markers, unstable_markers = dict(), dict()

        if self._dim == 2:
            stable_markers = plotly_layouts["default_binary_marker_settings"].copy()
            stable_markers.update(
                dict(
                    x=list(stable_props["x"]),
                    y=list(stable_props["y"]),
                    name="Stable",
                    marker=dict(color="darkgreen", size=11, line=dict(color="black", width=2)),
                    opacity=0.9,
                    hovertext=stable_props["texts"],
                    error_y=dict(
                        array=list(stable_props["uncertainties"]),
                        type="data",
                        color="gray",
                        thickness=2.5,
                        width=5,
                    ),
                )
            )

            unstable_markers = plotly_layouts["default_binary_marker_settings"].copy()
            unstable_markers.update(
                dict(
                    x=list(unstable_props["x"]),
                    y=list(unstable_props["y"]),
                    name="Above Hull",
                    marker=dict(
                        color=unstable_props["energies"],
                        colorscale=plotly_layouts["unstable_colorscale"],
                        size=6,
                        symbol="diamond",
                    ),
                    hovertext=unstable_props["texts"],
                )
            )

        elif self._dim == 3:
            stable_markers = plotly_layouts["default_ternary_marker_settings"].copy()
            stable_markers.update(
                dict(
                    x=list(stable_props["y"]),
                    y=list(stable_props["x"]),
                    z=list(stable_props["z"]),
                    name="Stable",
                    marker=dict(
                        color="black",
                        size=12,
                        opacity=0.8,
                        line=dict(color="black", width=3),
                    ),
                    hovertext=stable_props["texts"],
                    error_z=dict(
                        array=list(stable_props["uncertainties"]),
                        type="data",
                        color="darkgray",
                        width=10,
                        thickness=5,
                    ),
                )
            )

            unstable_markers = plotly_layouts["default_ternary_marker_settings"].copy()
            unstable_markers.update(
                dict(
                    x=unstable_props["y"],
                    y=unstable_props["x"],
                    z=unstable_props["z"],
                    name="Above Hull",
                    marker=dict(
                        color=unstable_props["energies"],
                        colorscale=plotly_layouts["unstable_colorscale"],
                        size=6,
                        symbol="diamond",
                        colorbar=dict(title="Energy Above Hull<br>(eV/atom)", x=0.05, len=0.75),
                    ),
                    hovertext=unstable_props["texts"],
                )
            )

        elif self._dim == 4:
            stable_markers = plotly_layouts["default_quaternary_marker_settings"].copy()
            stable_markers.update(
                dict(
                    x=stable_props["x"],
                    y=stable_props["y"],
                    z=stable_props["z"],
                    name="Stable",
                    marker=dict(
                        color=stable_props["energies"],
                        colorscale=plotly_layouts["stable_markers_colorscale"],
                        size=8,
                        opacity=0.9,
                    ),
                    hovertext=stable_props["texts"],
                )
            )

            unstable_markers = plotly_layouts["default_quaternary_marker_settings"].copy()
            unstable_markers.update(
                dict(
                    x=unstable_props["x"],
                    y=unstable_props["y"],
                    z=unstable_props["z"],
                    name="Above Hull",
                    marker=dict(
                        color=unstable_props["energies"],
                        colorscale=plotly_layouts["unstable_colorscale"],
                        size=5,
                        symbol="diamond",
                        colorbar=dict(title="Energy Above Hull<br>(eV/atom)", x=0.05, len=0.75),
                    ),
                    hovertext=unstable_props["texts"],
                    visible="legendonly",
                )
            )

        stable_marker_plot = go.Scatter(**stable_markers) if self._dim == 2 else go.Scatter3d(**stable_markers)
        unstable_marker_plot = go.Scatter(**unstable_markers) if self._dim == 2 else go.Scatter3d(**unstable_markers)

        return stable_marker_plot, unstable_marker_plot

    def _create_plotly_uncertainty_shading(self, stable_marker_plot):
        """
        Creates shaded uncertainty region for stable entries. Currently only works
        for binary (dim=2) phase diagrams.

        :param stable_marker_plot: go.Scatter object with stable markers and their
            error bars.
        :return: Plotly go.Scatter object with uncertainty window shading.
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
            points = points[points[:, 0].argsort()]  # sort by composition  # pylint: disable=E1136

            # these steps trace out the boundary pts of the uncertainty window
            outline = points[:, :2].copy()
            outline[:, 1] = outline[:, 1] + points[:, 2]

            last = -1
            if transformed:
                last = None  # allows for uncertainty in terminal compounds

            flipped_points = np.flip(points[:last, :].copy(), axis=0)
            flipped_points[:, 1] = flipped_points[:, 1] - flipped_points[:, 2]
            outline = np.vstack((outline, flipped_points[:, :2]))

            uncertainty_plot = go.Scatter(
                x=outline[:, 0],
                y=outline[:, 1],
                name="Uncertainty (window)",
                fill="toself",
                mode="lines",
                line=dict(width=0),
                fillcolor="lightblue",
                hoverinfo="skip",
                opacity=0.4,
            )

        return uncertainty_plot

    def _create_plotly_ternary_support_lines(self):
        """
        Creates support lines which aid in seeing the ternary hull in three
        dimensions.

        :return: go.Scatter3d plot of support lines for ternary phase diagram.
        """
        stable_entry_coords = dict(map(reversed, self.pd_plot_data[1].items()))

        elem_coords = [stable_entry_coords[e] for e in self._pd.el_refs.values()]

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
            line=dict(color="rgba (0, 0, 0, 0.4)", dash="solid", width=1.0),
            showlegend=False,
        )

    def _create_plotly_ternary_hull(self):
        """
        Creates shaded mesh plot for coloring the ternary hull by formation energy.

        :return: go.Mesh3d plot
        """
        facets = np.array(self._pd.facets)
        coords = np.array([triangular_coord(c) for c in zip(self._pd.qhull_data[:-1, 0], self._pd.qhull_data[:-1, 1])])
        energies = np.array([self._pd.get_form_energy_per_atom(e) for e in self._pd.qhull_entries])

        return go.Mesh3d(
            x=list(coords[:, 1]),
            y=list(coords[:, 0]),
            z=list(energies),
            i=list(facets[:, 1]),
            j=list(facets[:, 0]),
            k=list(facets[:, 2]),
            opacity=0.8,
            intensity=list(energies),
            colorscale=plotly_layouts["stable_colorscale"],
            colorbar=dict(title="Formation energy<br>(eV/atom)", x=0.9, len=0.75),
            hoverinfo="none",
            lighting=dict(diffuse=0.0, ambient=1.0),
            name="Convex Hull (shading)",
            flatshading=True,
            showlegend=True,
        )

    @staticmethod
    def _htmlize_formula(formula: str):
        """
        Adds HTML tags for displaying chemical formula in Plotly figure annotations.

        :param formula: chemical formula
        :return: clean chemical formula with necessary HTML tags
        """
        s = []
        for char in formula:
            if char.isdigit():
                s.append(f"<sub>{char}</sub>")
            else:
                s.append(char)

        return "".join(s)


def uniquelines(q):
    """
    Given all the facets, convert it into a set of unique lines.  Specifically
    used for converting convex hull facets into line pairs of coordinates.

    Args:
        q: A 2-dim sequence, where each row represents a facet. E.g.,
            [[1,2,3],[3,6,7],...]

    Returns:
        setoflines:
            A set of tuple of lines.  E.g., ((1,2), (1,3), (2,3), ....)
    """
    setoflines = set()
    for facets in q:
        for line in itertools.combinations(facets, 2):
            setoflines.add(tuple(sorted(line)))
    return setoflines


def triangular_coord(coord):
    """
    Convert a 2D coordinate into a triangle-based coordinate system for a
    prettier phase diagram.

    Args:
        coord: coordinate used in the convex hull computation.

    Returns:
        coordinates in a triangular-based coordinate system.
    """
    unitvec = np.array([[1, 0], [0.5, math.sqrt(3) / 2]])

    result = np.dot(np.array(coord), unitvec)
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
            [0.5, 1.0 / 3.0 * math.sqrt(3) / 2, math.sqrt(6) / 3],
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
        (newlines, newstable_entries, newunstable_entries):
        - newlines is a list of list of coordinates for lines in the PD.
        - newstable_entries is a {coordinate : entry} for each stable node
        in the phase diagram. (Each coordinate can only have one
        stable phase)
        - newunstable_entries is a {entry: coordinates} for all unstable
        nodes in the phase diagram.
    """
    yup = -1000.0
    xleft = 1000.0
    xright = -1000.0

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
            'Error in ordering_phase_diagram : \n"{up}", "{left}" and "{'
            'right}"'
            " should be in ordering : {ord}".format(up=nameup, left=nameleft, right=nameright, ord=ordering)
        )

    cc = np.array([0.5, np.sqrt(3.0) / 6.0], np.float_)

    if nameup == ordering[0]:
        if nameleft == ordering[1]:
            # The coordinates were already in the user ordering
            return lines, stable_entries, unstable_entries

        newlines = [[np.array(1.0 - x), y] for x, y in lines]
        newstable_entries = {(1.0 - c[0], c[1]): entry for c, entry in stable_entries.items()}
        newunstable_entries = {entry: (1.0 - c[0], c[1]) for entry, c in unstable_entries.items()}
        return newlines, newstable_entries, newunstable_entries
    if nameup == ordering[1]:
        if nameleft == ordering[2]:
            c120 = np.cos(2.0 * np.pi / 3.0)
            s120 = np.sin(2.0 * np.pi / 3.0)
            newlines = []
            for x, y in lines:
                newx = np.zeros_like(x)
                newy = np.zeros_like(y)
                for ii, xx in enumerate(x):
                    newx[ii] = c120 * (xx - cc[0]) - s120 * (y[ii] - cc[1]) + cc[0]
                    newy[ii] = s120 * (xx - cc[0]) + c120 * (y[ii] - cc[1]) + cc[1]
                newlines.append([newx, newy])
            newstable_entries = {
                (
                    c120 * (c[0] - cc[0]) - s120 * (c[1] - cc[1]) + cc[0],
                    s120 * (c[0] - cc[0]) + c120 * (c[1] - cc[1]) + cc[1],
                ): entry
                for c, entry in stable_entries.items()
            }
            newunstable_entries = {
                entry: (
                    c120 * (c[0] - cc[0]) - s120 * (c[1] - cc[1]) + cc[0],
                    s120 * (c[0] - cc[0]) + c120 * (c[1] - cc[1]) + cc[1],
                )
                for entry, c in unstable_entries.items()
            }
            return newlines, newstable_entries, newunstable_entries
        c120 = np.cos(2.0 * np.pi / 3.0)
        s120 = np.sin(2.0 * np.pi / 3.0)
        newlines = []
        for x, y in lines:
            newx = np.zeros_like(x)
            newy = np.zeros_like(y)
            for ii, xx in enumerate(x):
                newx[ii] = -c120 * (xx - 1.0) - s120 * y[ii] + 1.0
                newy[ii] = -s120 * (xx - 1.0) + c120 * y[ii]
            newlines.append([newx, newy])
        newstable_entries = {
            (
                -c120 * (c[0] - 1.0) - s120 * c[1] + 1.0,
                -s120 * (c[0] - 1.0) + c120 * c[1],
            ): entry
            for c, entry in stable_entries.items()
        }
        newunstable_entries = {
            entry: (
                -c120 * (c[0] - 1.0) - s120 * c[1] + 1.0,
                -s120 * (c[0] - 1.0) + c120 * c[1],
            )
            for entry, c in unstable_entries.items()
        }
        return newlines, newstable_entries, newunstable_entries
    if nameup == ordering[2]:
        if nameleft == ordering[0]:
            c240 = np.cos(4.0 * np.pi / 3.0)
            s240 = np.sin(4.0 * np.pi / 3.0)
            newlines = []
            for x, y in lines:
                newx = np.zeros_like(x)
                newy = np.zeros_like(y)
                for ii, xx in enumerate(x):
                    newx[ii] = c240 * (xx - cc[0]) - s240 * (y[ii] - cc[1]) + cc[0]
                    newy[ii] = s240 * (xx - cc[0]) + c240 * (y[ii] - cc[1]) + cc[1]
                newlines.append([newx, newy])
            newstable_entries = {
                (
                    c240 * (c[0] - cc[0]) - s240 * (c[1] - cc[1]) + cc[0],
                    s240 * (c[0] - cc[0]) + c240 * (c[1] - cc[1]) + cc[1],
                ): entry
                for c, entry in stable_entries.items()
            }
            newunstable_entries = {
                entry: (
                    c240 * (c[0] - cc[0]) - s240 * (c[1] - cc[1]) + cc[0],
                    s240 * (c[0] - cc[0]) + c240 * (c[1] - cc[1]) + cc[1],
                )
                for entry, c in unstable_entries.items()
            }
            return newlines, newstable_entries, newunstable_entries
        c240 = np.cos(4.0 * np.pi / 3.0)
        s240 = np.sin(4.0 * np.pi / 3.0)
        newlines = []
        for x, y in lines:
            newx = np.zeros_like(x)
            newy = np.zeros_like(y)
            for ii, xx in enumerate(x):
                newx[ii] = -c240 * xx - s240 * y[ii]
                newy[ii] = -s240 * xx + c240 * y[ii]
            newlines.append([newx, newy])
        newstable_entries = {
            (-c240 * c[0] - s240 * c[1], -s240 * c[0] + c240 * c[1]): entry for c, entry in stable_entries.items()
        }
        newunstable_entries = {
            entry: (-c240 * c[0] - s240 * c[1], -s240 * c[0] + c240 * c[1]) for entry, c in unstable_entries.items()
        }
        return newlines, newstable_entries, newunstable_entries
    raise ValueError("Invalid ordering.")
