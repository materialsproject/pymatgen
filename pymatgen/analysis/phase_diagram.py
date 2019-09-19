# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module defines tools to generate and analyze phase diagrams.
"""

import re
import collections
import itertools
import math
import logging

from monty.json import MSONable, MontyDecoder
from functools import lru_cache

import numpy as np
from scipy.spatial import ConvexHull

from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element, DummySpecie, get_el_sp
from pymatgen.util.coord import Simplex, in_coord_list
from pymatgen.util.string import latexify
from pymatgen.util.plotting import pretty_plot
from pymatgen.analysis.reaction_calculator import Reaction, \
    ReactionError

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "May 16, 2011"

logger = logging.getLogger(__name__)


class PDEntry(MSONable):
    """
    An object encompassing all relevant data for phase diagrams.

    .. attribute:: composition

        The composition associated with the PDEntry.

    .. attribute:: energy

        The energy associated with the entry.

    .. attribute:: name

        A name for the entry. This is the string shown in the phase diagrams.
        By default, this is the reduced formula for the composition, but can be
        set to some other string for display purposes.

    .. attribute:: attribute

        A arbitrary attribute.
    """

    def __init__(self, composition: Composition, energy: float,
                 name: str = None, attribute: object = None):
        """
        Args:
            composition (Composition): Composition
            energy (float): Energy for composition.
            name (str): Optional parameter to name the entry. Defaults to the
                reduced chemical formula.
            attribute: Optional attribute of the entry. This can be used to
                specify that the entry is a newly found compound, or to specify a
                particular label for the entry, or else ... Used for further
                analysis and plotting purposes. An attribute can be anything
                but must be MSONable.
        """
        self.energy = energy
        self.composition = Composition(composition)
        self.name = name if name else self.composition.reduced_formula
        self.attribute = attribute

    @property
    def energy_per_atom(self):
        """
        Returns the final energy per atom.
        """
        return self.energy / self.composition.num_atoms

    @property
    def is_element(self):
        """
        True if the entry is an element.
        """
        return self.composition.is_element

    def __repr__(self):
        return "PDEntry : {} with energy = {:.4f}".format(self.composition,
                                                          self.energy)

    def __str__(self):
        return self.__repr__()

    def as_dict(self):
        """
        :return: MSONable dict.
        """
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "composition": self.composition.as_dict(),
                "energy": self.energy,
                "name": self.name,
                "attribute": self.attribute}

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.as_dict() == other.as_dict()
        else:
            return False

    def __hash__(self):
        return id(self)

    @classmethod
    def from_dict(cls, d):
        """
        :param d: Dict representation
        :return: PDEntry
        """
        return cls(Composition(d["composition"]), d["energy"],
                   d["name"] if "name" in d else None,
                   d["attribute"] if "attribute" in d else None)


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
        comp = entry.composition
        self.original_entry = entry
        self.original_comp = comp
        grandpot = entry.energy - sum([comp[el] * pot
                                       for el, pot in chempots.items()])
        self.chempots = chempots
        new_comp_map = {el: comp[el] for el in comp.elements
                        if el not in chempots}
        super().__init__(new_comp_map, grandpot, entry.name)
        self.name = name if name else entry.name

    @property
    def is_element(self):
        """
        True if the entry is an element.
        """
        return self.original_comp.is_element

    def __repr__(self):
        chempot_str = " ".join(["mu_%s = %.4f" % (el, mu)
                                for el, mu in self.chempots.items()])
        return "GrandPotPDEntry with original composition " + \
               "{}, energy = {:.4f}, {}".format(self.original_entry.composition,
                                                self.original_entry.energy,
                                                chempot_str)

    def __str__(self):
        return self.__repr__()

    def as_dict(self):
        """
        :return: MSONAble dict
        """
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "entry": self.original_entry.as_dict(),
                "chempots": {el.symbol: u for el, u in self.chempots.items()},
                "name": self.name}

    @classmethod
    def from_dict(cls, d):
        """
        :param d: Dict representation
        :return: PDStructureEntry
        """
        chempots = {Element(symbol): u for symbol, u in d["chempots"].items()}
        entry = MontyDecoder().process_decoded(d["entry"])
        return cls(entry, chempots, d["name"])

    def __getattr__(self, a):
        """
        Delegate attribute to original entry if available.
        """
        if hasattr(self.original_entry, a):
            return getattr(self.original_entry, a)
        raise AttributeError(a)


class TransformedPDEntry(PDEntry):
    """
    This class repesents a TransformedPDEntry, which allows for a PDEntry to be
    transformed to a different composition coordinate space. It is used in the
    construction of phase diagrams that do not have elements as the terminal
    compositions.
    """

    def __init__(self, comp, original_entry):
        """
        Args:
            comp (Composition): Transformed composition as a Composition.
            original_entry (PDEntry): Original entry that this entry arose from.
        """
        super().__init__(comp, original_entry.energy)
        self.original_entry = original_entry
        self.name = original_entry.name

    def __getattr__(self, a):
        """
        Delegate attribute to original entry if available.
        """
        if hasattr(self.original_entry, a):
            return getattr(self.original_entry, a)
        raise AttributeError(a)

    def __repr__(self):
        output = ["TransformedPDEntry {}".format(self.composition),
                  " with original composition {}".format(self.original_entry.composition),
                  ", E = {:.4f}".format(self.original_entry.energy)]
        return "".join(output)

    def __str__(self):
        return self.__repr__()

    def as_dict(self):
        """
        :return: MSONable dict
        """
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "entry": self.original_entry.as_dict(),
                "composition": self.composition}

    @classmethod
    def from_dict(cls, d):
        """
        :param d: Dict representation
        :return: TransformedPDEntry
        """
        entry = MontyDecoder().process_decoded(d["entry"])
        return cls(d["composition"], entry)


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

    # Tolerance for determining if formation energy is positive.
    formation_energy_tol = 1e-11
    numerical_tol = 1e-8

    def __init__(self, entries, elements=None):
        """
        Standard constructor for phase diagram.

        Args:
            entries ([PDEntry]): A list of PDEntry-like objects having an
                energy, energy_per_atom and composition.
            elements ([Element]): Optional list of elements in the phase
                diagram. If set to None, the elements are determined from
                the the entries themselves.
        """
        if elements is None:
            elements = set()
            for entry in entries:
                elements.update(entry.composition.elements)
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
            raise PhaseDiagramError(
                "There are no entries associated with a terminal element!.")

        data = np.array([
            [e.composition.get_atomic_fraction(el) for el in elements] + [
                e.energy_per_atom]
            for e in min_entries
        ])

        # Use only entries with negative formation energy
        vec = [el_refs[el].energy_per_atom for el in elements] + [-1]
        form_e = -np.dot(data, vec)
        inds = np.where(form_e < -self.formation_energy_tol)[0].tolist()

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
            self.facets = [qhull_data.argmin(axis=0)]
        else:
            facets = get_facets(qhull_data)
            finalfacets = []
            for facet in facets:
                # Skip facets that include the extra point
                if max(facet) == len(qhull_data) - 1:
                    continue
                m = qhull_data[facet]
                m[:, -1] = 1
                if abs(np.linalg.det(m)) > 1e-14:
                    finalfacets.append(facet)
            self.facets = finalfacets

        self.simplexes = [Simplex(qhull_data[f, :-1]) for f in self.facets]
        self.all_entries = all_entries
        self.qhull_data = qhull_data
        self.dim = dim
        self.el_refs = el_refs
        self.elements = elements
        self.qhull_entries = qhull_entries
        self._stable_entries = set(self.qhull_entries[i] for i in
                                   set(itertools.chain(*self.facets)))

    def pd_coords(self, comp):
        """
        The phase diagram is generated in a reduced dimensional space
        (n_elements - 1). This function returns the coordinates in that space.
        These coordinates are compatible with the stored simplex objects.
        """
        if set(comp.elements).difference(self.elements):
            raise ValueError('{} has elements not in the phase diagram {}'
                             ''.format(comp, self.elements))
        return np.array(
            [comp.get_atomic_fraction(el) for el in self.elements[1:]])

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
        Entries that are unstable in the phase diagram. Includes positive
        formation energy entries.
        """
        return [e for e in self.all_entries if e not in self.stable_entries]

    @property
    def stable_entries(self):
        """
        Returns the stable entries in the phase diagram.
        """
        return self._stable_entries

    def get_form_energy(self, entry):
        """
        Returns the formation energy for an entry (NOT normalized) from the
        elemental references.

        Args:
            entry: A PDEntry-like object.

        Returns:
            Formation energy from the elemental references.
        """
        c = entry.composition
        return entry.energy - sum([c[el] * self.el_refs[el].energy_per_atom
                                   for el in c.elements])

    def get_form_energy_per_atom(self, entry):
        """
        Returns the formation energy per atom for an entry from the
        elemental references.

        Args:
            entry: An PDEntry-like object

        Returns:
            Formation energy **per atom** from the elemental references.
        """
        return self.get_form_energy(entry) / entry.composition.num_atoms

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        symbols = [el.symbol for el in self.elements]
        output = ["{} phase diagram".format("-".join(symbols)),
                  "{} stable phases: ".format(len(self.stable_entries)),
                  ", ".join([entry.name
                             for entry in self.stable_entries])]
        return "\n".join(output)

    def as_dict(self):
        """
        :return: MSONAble dict
        """
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "all_entries": [e.as_dict() for e in self.all_entries],
                "elements": [e.as_dict() for e in self.elements]}

    @classmethod
    def from_dict(cls, d):
        """
        :param d: Dict representation
        :return: PhaseDiagram
        """
        entries = [MontyDecoder().process_decoded(dd) for dd in d["all_entries"]]
        elements = [Element.from_dict(dd) for dd in d["elements"]]
        return cls(entries, elements)

    @lru_cache(1)
    def _get_facet_and_simplex(self, comp):
        """
        Get any facet that a composition falls into. Cached so successive
        calls at same composition are fast.
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
            { element: chempot } for all elements in the phase diagram.
        """
        complist = [self.qhull_entries[i].composition for i in facet]
        energylist = [self.qhull_entries[i].energy_per_atom for i in facet]
        m = [[c.get_atomic_fraction(e) for e in self.elements] for c in
             complist]
        chempots = np.linalg.solve(m, energylist)
        return dict(zip(self.elements, chempots))

    def get_decomposition(self, comp):
        """
        Provides the decomposition at a particular composition.

        Args:
            comp: A composition

        Returns:
            Decomposition as a dict of {Entry: amount}
        """
        facet, simplex = self._get_facet_and_simplex(comp)
        decomp_amts = simplex.bary_coords(self.pd_coords(comp))
        return {self.qhull_entries[f]: amt
                for f, amt in zip(facet, decomp_amts)
                if abs(amt) > PhaseDiagram.numerical_tol}

    def get_hull_energy(self, comp):
        """
        Args:
            comp (Composition): Input composition

        Returns:
            Energy of lowest energy equilibrium at desired composition. Not
            normalized by atoms, i.e. E(Li4O2) = 2 * E(Li2O)
        """
        e = 0
        for k, v in self.get_decomposition(comp).items():
            e += k.energy_per_atom * v
        return e * comp.num_atoms

    def get_decomp_and_e_above_hull(self, entry, allow_negative=False):
        """
        Provides the decomposition and energy above convex hull for an entry.
        Due to caching, can be much faster if entries with the same composition
        are processed together.

        Args:
            entry: A PDEntry like object
            allow_negative: Whether to allow negative e_above_hulls. Used to
                calculate equilibrium reaction energies. Defaults to False.

        Returns:
            (decomp, energy above convex hull)  Stable entries should have
            energy above hull of 0. The decomposition is provided as a dict of
            {Entry: amount}.
        """
        if entry in self.stable_entries:
            return {entry: 1}, 0

        comp = entry.composition
        facet, simplex = self._get_facet_and_simplex(comp)
        decomp_amts = simplex.bary_coords(self.pd_coords(comp))
        decomp = {self.qhull_entries[f]: amt
                  for f, amt in zip(facet, decomp_amts)
                  if abs(amt) > PhaseDiagram.numerical_tol}
        energies = [self.qhull_entries[i].energy_per_atom for i in facet]
        ehull = entry.energy_per_atom - np.dot(decomp_amts, energies)
        if allow_negative or ehull >= -PhaseDiagram.numerical_tol:
            return decomp, ehull
        raise ValueError("No valid decomp found!")

    def get_e_above_hull(self, entry):
        """
        Provides the energy above convex hull for an entry

        Args:
            entry: A PDEntry like object

        Returns:
            Energy above convex hull of entry. Stable entries should have
            energy above hull of 0.
        """
        return self.get_decomp_and_e_above_hull(entry)[1]

    def get_equilibrium_reaction_energy(self, entry):
        """
        Provides the reaction energy of a stable entry from the neighboring
        equilibrium stable entries (also known as the inverse distance to
        hull).

        Args:
            entry: A PDEntry like object

        Returns:
            Equilibrium reaction energy of entry. Stable entries should have
            equilibrium reaction energy <= 0.
        """
        if entry not in self.stable_entries:
            raise ValueError("Equilibrium reaction energy is available only "
                             "for stable entries.")
        if entry.is_element:
            return 0
        entries = [e for e in self.stable_entries if e != entry]
        modpd = PhaseDiagram(entries, self.elements)
        return modpd.get_decomp_and_e_above_hull(entry,
                                                 allow_negative=True)[1]

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
        # note the top part takes from format of _get_facet_and_simplex,
        #   but wants to return all facets rather than the first one that meets this criteria
        c = self.pd_coords(comp)
        allfacets = []
        for f, s in zip(self.facets, self.simplexes):
            if s.in_simplex(c, PhaseDiagram.numerical_tol / 10):
                allfacets.append(f)

        if not len(allfacets):
            raise RuntimeError("No facets found for comp = {}".format(comp))
        else:
            chempots = {}
            for facet in allfacets:
                facet_elt_list = [self.qhull_entries[j].name for j in facet]
                facet_name = '-'.join(facet_elt_list)
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
            raise ValueError("get_transition_chempots can only be called with "
                             "elements in the phase diagram.")

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
        l = (c2 - c1)
        l /= np.sum(l ** 2) ** 0.5
        proj = np.dot(intersections - c1, l)

        # only take compositions between endpoints
        proj = proj[np.logical_and(proj > -self.numerical_tol,
                                   proj < proj[1] + self.numerical_tol)]
        proj.sort()

        # only unique compositions
        valid = np.ones(len(proj), dtype=np.bool)
        valid[1:] = proj[1:] > proj[:-1] + self.numerical_tol
        proj = proj[valid]

        ints = c1 + l * proj[:, None]
        # reconstruct full-dimensional composition array
        cs = np.concatenate([np.array([1 - np.sum(ints, axis=-1)]).T,
                             ints], axis=-1)
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
            raise ValueError("get_transition_chempots can only be called with"
                             " elements in the phase diagram.")
        gccomp = Composition({el: amt for el, amt in comp.items()
                              if el != element})
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
            evolution.append({'chempot': c,
                              'evolution': amt,
                              'element_reference': elref,
                              'reaction': rxn, 'entries': decomp_entries})
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
        pd = self
        facets = pd.facets
        for facet in facets:
            chempots = self._get_facet_chempots(facet)
            all_chempots.append([chempots[el] for el in pd.elements])
        inds = [pd.elements.index(el) for el in elements]
        el_energies = {el: 0.0 for el in elements}
        if referenced:
            el_energies = {el: pd.el_refs[el].energy_per_atom
                           for el in elements}
        chempot_ranges = collections.defaultdict(list)
        vertices = [list(range(len(self.elements)))]
        if len(all_chempots) > len(self.elements):
            vertices = get_facets(all_chempots, joggle=joggle)
        for ufacet in vertices:
            for combi in itertools.combinations(ufacet, 2):
                data1 = facets[combi[0]]
                data2 = facets[combi[1]]
                common_ent_ind = set(data1).intersection(set(data2))
                if len(common_ent_ind) == len(elements):
                    common_entries = [pd.qhull_entries[i]
                                      for i in common_ent_ind]
                    data = np.array([[all_chempots[i][j]
                                      - el_energies[pd.elements[j]]
                                      for j in inds] for i in combi])
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
                from the energy of
            the stable phase at the target composition
            tol_en: a tolerance on the energy to set

        Returns:
             [{Element:mu}]: An array of conditions on simplex vertices for
             which each element has a chemical potential set to a given
             value. "absolute" values (i.e., not referenced to element energies)
        """
        muref = np.array([self.el_refs[e].energy_per_atom
                          for e in self.elements if e != dep_elt])
        chempot_ranges = self.get_chempot_range_map(
            [e for e in self.elements if e != dep_elt])

        for e in self.elements:
            if e not in target_comp.elements:
                target_comp = target_comp + Composition({e: 0.0})
        coeff = [-target_comp[e] for e in self.elements if e != dep_elt]
        for e in chempot_ranges.keys():
            if e.composition.reduced_composition == \
                    target_comp.reduced_composition:
                multiplicator = e.composition[dep_elt] / target_comp[dep_elt]
                ef = e.energy / multiplicator
                all_coords = []
                for s in chempot_ranges[e]:
                    for v in s._coords:
                        elts = [e for e in self.elements if e != dep_elt]
                        res = {}
                        for i in range(len(elts)):
                            res[elts[i]] = v[i] + muref[i]
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
        muref = np.array([self.el_refs[e].energy_per_atom
                          for e in self.elements if e != open_elt])
        chempot_ranges = self.get_chempot_range_map(
            [e for e in self.elements if e != open_elt])
        for e in self.elements:
            if e not in target_comp.elements:
                target_comp = target_comp + Composition({e: 0.0})
        coeff = [-target_comp[e] for e in self.elements if e != open_elt]
        max_open = -float('inf')
        min_open = float('inf')
        max_mus = None
        min_mus = None
        for e in chempot_ranges.keys():
            if e.composition.reduced_composition == \
                    target_comp.reduced_composition:
                multiplicator = e.composition[open_elt] / target_comp[open_elt]
                ef = e.energy / multiplicator
                all_coords = []
                for s in chempot_ranges[e]:
                    for v in s._coords:
                        all_coords.append(v)
                        if (np.dot(v + muref, coeff) + ef) / target_comp[open_elt] > max_open:
                            max_open = (np.dot(v + muref, coeff) + ef) / target_comp[open_elt]
                            max_mus = v
                        if (np.dot(v + muref, coeff) + ef) / target_comp[open_elt] < min_open:
                            min_open = (np.dot(v + muref, coeff) + ef) / target_comp[open_elt]
                            min_mus = v
        elts = [e for e in self.elements if e != open_elt]
        res = {}
        for i in range(len(elts)):
            res[elts[i]] = (min_mus[i] + muref[i], max_mus[i] + muref[i])
        res[open_elt] = (min_open, max_open)
        return res


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
            chempots {Element: float}: Specify the chemical potentials
                of the open elements.
            elements ([Element]): Optional list of elements in the phase
                diagram. If set to None, the elements are determined from
                the the entries themselves.
        """
        if elements is None:
            elements = set()
            for entry in entries:
                elements.update(entry.composition.elements)
        self.chempots = {get_el_sp(el): u for el, u in chempots.items()}
        elements = set(elements).difference(self.chempots.keys())
        all_entries = []
        for e in entries:
            if len(set(e.composition.elements).intersection(set(elements))) > 0:
                all_entries.append(GrandPotPDEntry(e, self.chempots))
        super().__init__(all_entries, elements)

    def __str__(self):
        output = []
        chemsys = "-".join([el.symbol for el in self.elements])
        output.append("{} grand potential phase diagram with ".format(chemsys))
        output[-1] += ", ".join(["u{}={}".format(el, v)
                                 for el, v in self.chempots.items()])
        output.append("{} stable phases: ".format(len(self.stable_entries)))
        output.append(", ".join([entry.name
                                 for entry in self.stable_entries]))
        return "\n".join(output)

    def as_dict(self):
        """
        :return: MSONable dict
        """
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "all_entries": [e.as_dict() for e in self.all_entries],
                "chempots": self.chempots,
                "elements": [e.as_dict() for e in self.elements]}

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

    def __init__(self, entries, terminal_compositions,
                 normalize_terminal_compositions=True):
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
        (pentries, species_mapping) = \
            self.transform_entries(entries, terminal_compositions)
        self.species_mapping = species_mapping
        super().__init__(
            pentries, elements=species_mapping.values())

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
            fractional_comp = [c.fractional_composition
                               for c in terminal_compositions]
        else:
            fractional_comp = terminal_compositions

        # Map terminal compositions to unique dummy species.
        sp_mapping = collections.OrderedDict()
        for i, comp in enumerate(fractional_comp):
            sp_mapping[comp] = DummySpecie("X" + chr(102 + i))

        for entry in entries:
            try:
                rxn = Reaction(fractional_comp, [entry.composition])
                rxn.normalize_to(entry.composition)
                # We only allow reactions that have positive amounts of
                # reactants.
                if all([rxn.get_coeff(comp) <= CompoundPhaseDiagram.amount_tol
                        for comp in fractional_comp]):
                    newcomp = {sp_mapping[comp]: -rxn.get_coeff(comp)
                               for comp in fractional_comp}
                    newcomp = {k: v for k, v in newcomp.items()
                               if v > CompoundPhaseDiagram.amount_tol}
                    transformed_entry = \
                        TransformedPDEntry(Composition(newcomp), entry)
                    new_entries.append(transformed_entry)
            except ReactionError:
                # If the reaction can't be balanced, the entry does not fall
                # into the phase space. We ignore them.
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
            "terminal_compositions": [c.as_dict()
                                      for c in self.terminal_compositions],
            "normalize_terminal_compositions":
                self.normalize_terminals}

    @classmethod
    def from_dict(cls, d):
        """
        :param d: Dict Representation
        :return: CompoundPhaseDiagram
        """
        dec = MontyDecoder()
        entries = dec.process_decoded(d["original_entries"])
        terminal_compositions = dec.process_decoded(d["terminal_compositions"])
        return cls(entries, terminal_compositions,
                   d["normalize_terminal_compositions"])


class ReactionDiagram:
    """
    Analyzes the possible reactions between a pair of compounds, e.g.,
    an electrolyte and an electrode.
    """

    def __init__(self, entry1, entry2, all_entries, tol=1e-4,
                 float_fmt="%.4f"):
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

        comp_vec1 = np.array([entry1.composition.get_atomic_fraction(el)
                              for el in elements])
        comp_vec2 = np.array([entry2.composition.get_atomic_fraction(el)
                              for el in elements])
        r1 = entry1.composition.reduced_composition
        r2 = entry2.composition.reduced_composition

        logger.debug("%d total entries." % len(all_entries))

        pd = PhaseDiagram(all_entries + [entry1, entry2])
        terminal_formulas = [entry1.composition.reduced_formula,
                             entry2.composition.reduced_formula]

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

                if any([e.composition.reduced_formula in terminal_formulas
                        for e in face_entries]):
                    continue

                try:

                    m = []
                    for e in face_entries:
                        m.append([e.composition.get_atomic_fraction(el)
                                  for el in elements])
                    m.append(comp_vec2 - comp_vec1)
                    m = np.array(m).T
                    coeffs = np.linalg.solve(m, comp_vec2)

                    x = coeffs[-1]

                    if all([c >= -tol for c in coeffs]) and \
                            (abs(sum(coeffs[:-1]) - 1) < tol) and \
                            (tol < x < 1 - tol):

                        c1 = x / r1.num_atoms
                        c2 = (1 - x) / r2.num_atoms
                        factor = 1 / (c1 + c2)

                        c1 *= factor
                        c2 *= factor

                        # Avoid duplicate reactions.
                        if any([np.allclose([c1, c2], cc) for cc in done]):
                            continue

                        done.append((c1, c2))

                        rxn_str = "%s %s + %s %s -> " % (
                            fmt(c1), r1.reduced_formula,
                            fmt(c2), r2.reduced_formula)
                        products = []

                        energy = - (x * entry1.energy_per_atom +
                                    (1 - x) * entry2.energy_per_atom)
                        for c, e in zip(coeffs[:-1], face_entries):
                            if c > tol:
                                r = e.composition.reduced_composition
                                products.append("%s %s" % (
                                    fmt(c / r.num_atoms * factor),
                                    r.reduced_formula))
                                energy += c * e.energy_per_atom

                        rxn_str += " + ".join(products)
                        comp = x * comp_vec1 + (1 - x) * comp_vec2
                        entry = PDEntry(
                            Composition(dict(zip(elements, comp))),
                            energy=energy, attribute=rxn_str)
                        rxn_entries.append(entry)
                except np.linalg.LinAlgError:
                    logger.debug("Reactants = %s" % (", ".join([
                        entry1.composition.reduced_formula,
                        entry2.composition.reduced_formula])))
                    logger.debug("Products = %s" % (
                        ", ".join([e.composition.reduced_formula
                                   for e in face_entries])))

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
            [Composition(entry1.composition.reduced_formula),
             Composition(entry2.composition.reduced_formula)],
            normalize_terminal_compositions=False)
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
    else:
        return ConvexHull(qhull_data, qhull_options="Qt i").simplices


class PDPlotter:
    """
    A plotter class for phase diagrams.
    """

    def __init__(self, phasediagram, show_unstable=0, **plotkwargs):
        r"""

        Args:
            phasediagram: PhaseDiagram object.
            show_unstable (float): Whether unstable phases will be plotted as
                well as red crosses. If a number > 0 is entered, all phases with
                ehull < show_unstable will be shown.
            **plotkwargs: Keyword args passed to matplotlib.pyplot.plot. Can
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
        self.lines = uniquelines(self._pd.facets) if self._dim > 1 else \
            [[self._pd.facets[0][0], self._pd.facets[0][0]]]
        self.show_unstable = show_unstable
        colors = Set1_3.mpl_colors
        self.plotkwargs = plotkwargs or {
            "markerfacecolor": colors[2],
            "markersize": 10,
            "linewidth": 3
        }

    @property
    def pd_plot_data(self):
        """
        Plot data for phase diagram.
        2-comp - Full hull with energies
        3/4-comp - Projection into 2D or 3D Gibbs triangle.

        Returns:
            (lines, stable_entries, unstable_entries):
            - lines is a list of list of coordinates for lines in the PD.
            - stable_entries is a {coordinate : entry} for each stable node
            in the phase diagram. (Each coordinate can only have one
            stable phase)
            - unstable_entries is a {entry: coordinates} for all unstable
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
                y = [pd.get_form_energy_per_atom(entry1),
                     pd.get_form_energy_per_atom(entry2)]
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
        for i in range(0, len(all_entries)):
            entry = all_entries[i]
            if entry not in stable:
                if self._dim < 3:
                    x = [all_data[i][0], all_data[i][0]]
                    y = [pd.get_form_energy_per_atom(entry),
                         pd.get_form_energy_per_atom(entry)]
                    coord = [x, y]
                elif self._dim == 3:
                    coord = triangular_coord([all_data[i, 0:2],
                                              all_data[i, 0:2]])
                else:
                    coord = tet_coord([all_data[i, 0:3], all_data[i, 0:3],
                                       all_data[i, 0:3]])
                labelcoord = list(zip(*coord))
                unstable_entries[entry] = labelcoord[0]

        return lines, stable_entries, unstable_entries

    def get_plot(self, label_stable=True, label_unstable=True, ordering=None,
                 energy_colormap=None, process_attributes=False, plt=None):
        """
        :param label_stable: Whether to label stable compounds.
        :param label_unstable: Whether to label unstable compounds.
        :param ordering: Ordering of vertices.
        :param energy_colormap: Colormap for coloring energy.
        :param process_attributes: Whether to process the attributes.
        :param plt: Existing plt object if plotting multiple phase diagrams.
        :return: matplotlib.pyplot.
        """
        if self._dim < 4:
            plt = self._get_2d_plot(label_stable, label_unstable, ordering,
                                    energy_colormap, plt=plt,
                                    process_attributes=process_attributes)
        elif self._dim == 4:
            plt = self._get_3d_plot(label_stable)

        return plt

    def plot_element_profile(self, element, comp, show_label_index=None,
                             xlim=5):
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
        element_energy = evolution[0]['chempot']
        x1, x2, y1 = None, None, None
        for i, d in enumerate(evolution):
            v = -(d["chempot"] - element_energy)
            print("index= %s, -\u0394\u03BC=%.4f(eV)," % (i, v), d["reaction"])
            if i != 0:
                plt.plot([x2, x2], [y1, d["evolution"] / num_atoms],
                         'k', linewidth=2.5)
            x1 = v
            y1 = d["evolution"] / num_atoms

            if i != len(evolution) - 1:
                x2 = - (evolution[i + 1]["chempot"] - element_energy)
            else:
                x2 = 5.0
            if show_label_index is not None and i in show_label_index:
                products = [re.sub(r"(\d+)", r"$_{\1}$", p.reduced_formula)
                            for p in d["reaction"].products
                            if p.reduced_formula != element.symbol]
                plt.annotate(", ".join(products), xy=(v + 0.05, y1 + 0.05),
                             fontsize=24, color='r')
                plt.plot([x1, x2], [y1, y1], 'r', linewidth=3)
            else:
                plt.plot([x1, x2], [y1, y1], 'k', linewidth=2.5)

        plt.xlim((0, xlim))
        plt.xlabel("-$\\Delta{\\mu}$ (eV)")
        plt.ylabel("Uptake per atom")

        return plt

    def show(self, *args, **kwargs):
        r"""
        Draws the phase diagram using Matplotlib and show it.

        Args:
            *args: Passed to get_plot.
            **kwargs: Passed to get_plot.
        """
        self.get_plot(*args, **kwargs).show()

    def _get_2d_plot(self, label_stable=True, label_unstable=True,
                     ordering=None, energy_colormap=None, vmin_mev=-60.0,
                     vmax_mev=60.0, show_colorbar=True,
                     process_attributes=False, plt=None):
        """
        Shows the plot using pylab.  Usually I won't do imports in methods,
        but since plotting is a fairly expensive library to load and not all
        machines have matplotlib installed, I have done it this way.
        """
        if plt is None:
            plt = pretty_plot(8, 6)
        from matplotlib.font_manager import FontProperties
        if ordering is None:
            (lines, labels, unstable) = self.pd_plot_data
        else:
            (_lines, _labels, _unstable) = self.pd_plot_data
            (lines, labels, unstable) = order_phase_diagram(
                _lines, _labels, _unstable, ordering)
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
                    if labels[(x, y)].attribute is None or \
                            labels[(x, y)].attribute == "existing":
                        plt.plot(x, y, "ko", **self.plotkwargs)
                    else:
                        plt.plot(x, y, "k*", **self.plotkwargs)
            else:
                for x, y in lines:
                    plt.plot(x, y, "ko-", **self.plotkwargs)
        else:
            from matplotlib.colors import Normalize, LinearSegmentedColormap
            from matplotlib.cm import ScalarMappable
            for x, y in lines:
                plt.plot(x, y, "k-", markeredgecolor="k")
            vmin = vmin_mev / 1000.0
            vmax = vmax_mev / 1000.0
            if energy_colormap == 'default':
                mid = - vmin / (vmax - vmin)
                cmap = LinearSegmentedColormap.from_list(
                    'my_colormap', [(0.0, '#005500'), (mid, '#55FF55'),
                                    (mid, '#FFAAAA'), (1.0, '#FF0000')])
            else:
                cmap = energy_colormap
            norm = Normalize(vmin=vmin, vmax=vmax)
            _map = ScalarMappable(norm=norm, cmap=cmap)
            _energies = [self._pd.get_equilibrium_reaction_energy(entry)
                         for coord, entry in labels.items()]
            energies = [en if en < 0.0 else -0.00000001 for en in _energies]
            vals_stable = _map.to_rgba(energies)
            ii = 0
            if process_attributes:
                for x, y in labels.keys():
                    if labels[(x, y)].attribute is None or \
                            labels[(x, y)].attribute == "existing":
                        plt.plot(x, y, "o", markerfacecolor=vals_stable[ii],
                                 markersize=12)
                    else:
                        plt.plot(x, y, "*", markerfacecolor=vals_stable[ii],
                                 markersize=18)
                    ii += 1
            else:
                for x, y in labels.keys():
                    plt.plot(x, y, "o", markerfacecolor=vals_stable[ii],
                             markersize=15)
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
            plt.xlabel("Fraction", fontsize=28, fontweight='bold')
            plt.ylabel("Formation energy (eV/fu)", fontsize=28,
                       fontweight='bold')

        for coords in sorted(labels.keys(), key=lambda x: -x[1]):
            entry = labels[coords]
            label = entry.name

            # The follow defines an offset for the annotation text emanating
            # from the center of the PD. Results in fairly nice layouts for the
            # most part.
            vec = (np.array(coords) - center)
            vec = vec / np.linalg.norm(vec) * 10 if np.linalg.norm(vec) != 0 \
                else vec
            valign = "bottom" if vec[1] > 0 else "top"
            if vec[0] < -0.01:
                halign = "right"
            elif vec[0] > 0.01:
                halign = "left"
            else:
                halign = "center"
            if label_stable:
                if process_attributes and entry.attribute == 'new':
                    plt.annotate(latexify(label), coords, xytext=vec,
                                 textcoords="offset points",
                                 horizontalalignment=halign,
                                 verticalalignment=valign,
                                 fontproperties=font,
                                 color='g')
                else:
                    plt.annotate(latexify(label), coords, xytext=vec,
                                 textcoords="offset points",
                                 horizontalalignment=halign,
                                 verticalalignment=valign,
                                 fontproperties=font)

        if self.show_unstable:
            font = FontProperties()
            font.set_size(16)
            energies_unstable = [self._pd.get_e_above_hull(entry)
                                 for entry, coord in unstable.items()]
            if energy_colormap is not None:
                energies.extend(energies_unstable)
                vals_unstable = _map.to_rgba(energies_unstable)
            ii = 0
            for entry, coords in unstable.items():
                ehull = self._pd.get_e_above_hull(entry)
                if ehull < self.show_unstable:
                    vec = (np.array(coords) - center)
                    vec = vec / np.linalg.norm(vec) * 10 \
                        if np.linalg.norm(vec) != 0 else vec
                    label = entry.name
                    if energy_colormap is None:
                        plt.plot(coords[0], coords[1], "ks", linewidth=3,
                                 markeredgecolor="k", markerfacecolor="r",
                                 markersize=8)
                    else:
                        plt.plot(coords[0], coords[1], "s", linewidth=3,
                                 markeredgecolor="k",
                                 markerfacecolor=vals_unstable[ii],
                                 markersize=8)
                    if label_unstable:
                        plt.annotate(latexify(label), coords, xytext=vec,
                                     textcoords="offset points",
                                     horizontalalignment=halign, color="b",
                                     verticalalignment=valign,
                                     fontproperties=font)
                    ii += 1
        if energy_colormap is not None and show_colorbar:
            _map.set_array(energies)
            cbar = plt.colorbar(_map)
            cbar.set_label(
                'Energy [meV/at] above hull (in red)\nInverse energy ['
                'meV/at] above hull (in green)',
                rotation=-90, ha='left', va='center')
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
        font = FontProperties()
        font.set_weight("bold")
        font.set_size(20)
        (lines, labels, unstable) = self.pd_plot_data
        count = 1
        newlabels = list()
        for x, y, z in lines:
            ax.plot(x, y, z, "bo-", linewidth=3, markeredgecolor="b",
                    markerfacecolor="r", markersize=10)
        for coords in sorted(labels.keys()):
            entry = labels[coords]
            label = entry.name
            if label_stable:
                if len(entry.composition.elements) == 1:
                    ax.text(coords[0], coords[1], coords[2], label)
                else:
                    ax.text(coords[0], coords[1], coords[2], str(count))
                    newlabels.append("{} : {}".format(count, latexify(label)))
                    count += 1
        plt.figtext(0.01, 0.01, "\n".join(newlabels))
        ax.axis("off")
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
        chempot_ranges = self._pd.get_chempot_range_map(
            elements, referenced=referenced)
        missing_lines = {}
        excluded_region = []
        for entry, lines in chempot_ranges.items():
            comp = entry.composition
            center_x = 0
            center_y = 0
            coords = []
            contain_zero = any([comp.get_atomic_fraction(el) == 0
                                for el in elements])
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

            plt.annotate(latexify(entry.name), xy,
                         horizontalalignment="center",
                         verticalalignment="center", fontsize=22)

        plt.xlabel("$\\mu_{{{0}}} - \\mu_{{{0}}}^0$ (eV)"
                   .format(el0.symbol))
        plt.ylabel("$\\mu_{{{0}}} - \\mu_{{{0}}}^0$ (eV)"
                   .format(el1.symbol))
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
        from scipy import interpolate
        from matplotlib import cm

        pd = self._pd
        entries = pd.qhull_entries
        data = np.array(pd.qhull_data)

        plt = self._get_2d_plot()
        data[:, 0:2] = triangular_coord(data[:, 0:2]).transpose()
        for i, e in enumerate(entries):
            data[i, 2] = self._pd.get_e_above_hull(e)

        gridsize = 0.005
        xnew = np.arange(0, 1., gridsize)
        ynew = np.arange(0, 1, gridsize)

        f = interpolate.LinearNDInterpolator(data[:, 0:2], data[:, 2])
        znew = np.zeros((len(ynew), len(xnew)))
        for (i, xval) in enumerate(xnew):
            for (j, yval) in enumerate(ynew):
                znew[j, i] = f(xval, yval)

        plt.contourf(xnew, ynew, znew, 1000, cmap=cm.autumn_r)

        plt.colorbar()
        return plt


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
        coordinate: coordinate used in the convex hull computation.

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
        coordinate: coordinate used in the convex hull computation.

    Returns:
        coordinates in a tetrahedron-based coordinate system.
    """
    unitvec = np.array([[1, 0, 0], [0.5, math.sqrt(3) / 2, 0],
                        [0.5, 1.0 / 3.0 * math.sqrt(3) / 2, math.sqrt(6) / 3]])
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
            ' should be in ordering : {ord}'.format(up=nameup, left=nameleft,
                                                    right=nameright,
                                                    ord=ordering))

    cc = np.array([0.5, np.sqrt(3.0) / 6.0], np.float)

    if nameup == ordering[0]:
        if nameleft == ordering[1]:
            # The coordinates were already in the user ordering
            return lines, stable_entries, unstable_entries
        else:
            newlines = [[np.array(1.0 - x), y] for x, y in lines]
            newstable_entries = {(1.0 - c[0], c[1]): entry
                                 for c, entry in stable_entries.items()}
            newunstable_entries = {entry: (1.0 - c[0], c[1])
                                   for entry, c in
                                   unstable_entries.items()}
            return newlines, newstable_entries, newunstable_entries
    elif nameup == ordering[1]:
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
                (c120 * (c[0] - cc[0]) - s120 * (c[1] - cc[1]) + cc[0],
                 s120 * (c[0] - cc[0]) + c120 * (c[1] - cc[1]) + cc[1]): entry
                for c, entry in stable_entries.items()}
            newunstable_entries = {
                entry: (c120 * (c[0] - cc[0]) - s120 * (c[1] - cc[1]) + cc[0],
                        s120 * (c[0] - cc[0]) + c120 * (c[1] - cc[1]) + cc[1])
                for entry, c in unstable_entries.items()}
            return newlines, newstable_entries, newunstable_entries
        else:
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
            newstable_entries = {(-c120 * (c[0] - 1.0) - s120 * c[1] + 1.0,
                                  -s120 * (c[0] - 1.0) + c120 * c[1]): entry
                                 for c, entry in stable_entries.items()}
            newunstable_entries = {
                entry: (-c120 * (c[0] - 1.0) - s120 * c[1] + 1.0,
                        -s120 * (c[0] - 1.0) + c120 * c[1])
                for entry, c in unstable_entries.items()}
            return newlines, newstable_entries, newunstable_entries
    elif nameup == ordering[2]:
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
                (c240 * (c[0] - cc[0]) - s240 * (c[1] - cc[1]) + cc[0],
                 s240 * (c[0] - cc[0]) + c240 * (c[1] - cc[1]) + cc[1]): entry
                for c, entry in stable_entries.items()}
            newunstable_entries = {
                entry: (c240 * (c[0] - cc[0]) - s240 * (c[1] - cc[1]) + cc[0],
                        s240 * (c[0] - cc[0]) + c240 * (c[1] - cc[1]) + cc[1])
                for entry, c in unstable_entries.items()}
            return newlines, newstable_entries, newunstable_entries
        else:
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
            newstable_entries = {(-c240 * c[0] - s240 * c[1],
                                  -s240 * c[0] + c240 * c[1]): entry
                                 for c, entry in stable_entries.items()}
            newunstable_entries = {entry: (-c240 * c[0] - s240 * c[1],
                                           -s240 * c[0] + c240 * c[1])
                                   for entry, c in unstable_entries.items()}
            return newlines, newstable_entries, newunstable_entries
