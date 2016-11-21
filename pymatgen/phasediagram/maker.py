# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
This module provides classes to create phase diagrams.
"""

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "2.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Nov 25, 2012"

import collections
import numpy as np
import itertools
from monty.json import MSONable, MontyDecoder

from scipy.spatial import ConvexHull

from pymatgen.core.periodic_table import get_el_sp
from pymatgen.core.composition import Composition
from pymatgen.phasediagram.entries import GrandPotPDEntry, TransformedPDEntry
from pymatgen.entries.computed_entries import ComputedEntry

from pymatgen.core.periodic_table import DummySpecie, Element
from pymatgen.analysis.reaction_calculator import Reaction, ReactionError


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

        get_reduced_comp = lambda e: e.composition.reduced_composition

        entries = sorted(entries, key=get_reduced_comp)

        el_refs = {}
        min_entries = []
        all_entries = []
        for c, g in itertools.groupby(entries, key=get_reduced_comp):
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
            [e.composition.get_atomic_fraction(el) for el in elements] + [e.energy_per_atom]
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
                if max(facet) == len(qhull_data)-1:
                    continue
                m = qhull_data[facet]
                m[:, -1] = 1
                if abs(np.linalg.det(m)) > 1e-14:
                    finalfacets.append(facet)
            self.facets = finalfacets

        self.simplices = [qhull_data[f, :-1] for f in self.facets]
        self.all_entries = all_entries
        self.qhull_data = qhull_data
        self.dim = dim
        self.el_refs = el_refs
        self.elements = elements
        self.qhull_entries = qhull_entries

    @property
    def all_entries_hulldata(self):
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
        return set((self.qhull_entries[i]
                    for i in itertools.chain(*self.facets)))

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
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "all_entries": [e.as_dict() for e in self.all_entries],
                "elements": [e.as_dict() for e in self.elements]}

    @classmethod
    def from_dict(cls, d):
        entries = [ComputedEntry.from_dict(dd) for dd in d["all_entries"]]
        elements = [Element.from_dict(dd) for dd in d["elements"]]
        return cls(entries, elements)


class GrandPotentialPhaseDiagram(PhaseDiagram):
    """
    A class representing a Grand potential phase diagram. Grand potential phase
    diagrams are essentially phase diagrams that are open to one or more
    components. To construct such phase diagrams, the relevant free energy is
    the grand potential, which can be written as the Legendre transform of the
    Gibbs free energy as follows

    Grand potential = G - u\ :sub:`X` N\ :sub:`X`\

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
        super(GrandPotentialPhaseDiagram, self).__init__(all_entries, elements)

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
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "all_entries": [e.as_dict() for e in self.all_entries],
                "chempots": self.chempots,
                "elements": [e.as_dict() for e in self.elements]}

    @classmethod
    def from_dict(cls, d):
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
        super(CompoundPhaseDiagram, self).__init__(
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
        dec = MontyDecoder()
        entries = dec.process_decoded(d["original_entries"])
        terminal_compositions = dec.process_decoded(d["terminal_compositions"])
        return cls(entries, terminal_compositions,
                   d["normalize_terminal_compositions"])


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
