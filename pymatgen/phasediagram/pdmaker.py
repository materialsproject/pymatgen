#!/usr/bin/env python

"""
This module provides classes to create phase diagrams.
"""

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "2.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__status__ = "Production"
__date__ = "Jul 11, 2012"

import numpy as np
import collections
import logging

from pymatgen.core.structure import Composition
from pymatgen.phasediagram.entries import GrandPotPDEntry, TransformedPDEntry
from pymatgen.util.coord_utils import get_convex_hull
from pymatgen.core.periodic_table import DummySpecie
from pymatgen.analysis.reaction_calculator import Reaction, ReactionError


logger = logging.getLogger(__name__)


class PhaseDiagram (object):
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
    """

    """
    Tolerance for determining if formation energy is positive.
    """
    formation_energy_tol = 1e-11

    def __init__(self, entries, elements=None, use_external_qhull=False):
        """
        Standard constructor for phase diagram.

        Args:
            entries:
                A list of PDEntry-like objects having an energy,
                energy_per_atom and composition.
            elements:
                Optional list of elements in the phase diagram. If set to None,
                the elements are determined from the the entries themselves.
            use_external_qhull:
                If set to True, the code will use an external command line call
                to Qhull to calculate the convex hull data. This requires the
                user to have qhull installed and the executables "qconvex"
                available in his path. By default, the code uses the
                scipy.spatail.Delaunay.

                The benefit of using external qhull is that a) it is much
                faster, especially for higher-dimensional hulls with many
                entries, and b) it is more robustly tested. The scipy Delaunay
                class is relatively new, and we have found some issues with
                higher-dimensional hulls, with non-sensical results given in
                some instances. Nonetheless, scipy Delaunay seems to work well
                enough for phase diagrams with < 4 components.
        """
        if elements == None:
            elements = set()
            map(elements.update, [entry.composition.elements
                                  for entry in entries])
        self._all_entries = entries
        self._elements = tuple(elements)
        self._qhull_data = None
        self._facets = None
        self._qhull_entries = None
        self._stable_entries = None
        self._use_external_qhull = use_external_qhull
        self._make_phasediagram()

    @property
    def all_entries(self):
        """
        All entries provided for Phase Diagram construction. Note that this
        does not mean that all these entries are actually used in the phase
        diagram. For example, this includes the positive formation energy
        entries that are filtered out before Phase Diagram construction.
        """
        return self._all_entries

    @property
    def dim(self):
        """
        The dimensionality of the phase diagram
        """
        return len(self._elements)

    @property
    def elements(self):
        """
        Elements in the phase diagram.
        """
        return self._elements

    @property
    def facets(self):
        """
        Facets of the phase diagram in the form of  [[1,2,3],[4,5,6]...]
        """
        return self._facets

    @property
    def qhull_data(self):
        """
        Data used in the convex hull operation. This is essentially a matrix of
        composition data and energy per atom values created from qhull_entries.
        """
        return self._qhull_data

    @property
    def qhull_entries(self):
        """
        Actual entries used in convex hull. Excludes all positive formation
        energy entries.
        """
        return self._qhull_entries

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

    @property
    def all_entries_hulldata(self):
        """
        Same as qhull_data, but for all entries rather than just positive
        formation energy ones.
        """
        return self._process_entries_qhulldata(self._all_entries)

    @property
    def el_refs(self):
        """
        List of elemental references for the phase diagrams. These are
        entries corresponding to the lowest energy element entries for simple
        compositional phase diagrams.
        """
        return self._el_refs

    def get_form_energy(self, entry):
        """
        Returns the formation energy for an entry (NOT normalized) from the
        elemental references.

        Args:
            entry:
                A PDEntry-like object.

        Returns:
            Formation energy from the elemental references.
        """
        comp = entry.composition
        energy = entry.energy - sum([comp[el] *
                                     self._el_refs[el].energy_per_atom
                                     for el in comp.elements])
        return energy

    def get_form_energy_per_atom(self, entry):
        """
        Returns the formation energy per atom for an entry from the
        elemental references.

        Args:
            entry:
                An PDEntry-like object

        Returns:
            Formation energy **per atom** from the elemental references.
        """
        comp = entry.composition
        return self.get_form_energy(entry) / comp.num_atoms

    def _process_entries_qhulldata(self, entries_to_process):
        """
        From a sequence of entries, generate the necessary for the convex hull.
        Using the Li-Fe-O phase diagram as an example, this is of the form:
        [[ Fe_fraction_entry_1, O_fraction_entry_1, Energy_per_atom_entry_1],
         [ Fe_fraction_entry_2, O_fraction_entry_2, Energy_per_atom_entry_2],
         ...]]

        Note that there are only two independent variables, since the third
        elemental fraction is fixed by the constraint that all compositions sum
        to 1. The choice of the elements is arbitrary.
        """
        data = []
        for entry in entries_to_process:
            comp = entry.composition
            energy_per_atom = entry.energy_per_atom
            row = []
            for i in xrange(1, len(self._elements)):
                row.append(comp.get_atomic_fraction(self._elements[i]))
            row.append(energy_per_atom)
            data.append(row)
        return data

    def _create_convhull_data(self):
        """
        Make data suitable for convex hull procedure from the list of entries.
        The procedure is as follows:

        1. First find the elemental references, i.e., the lowest energy entry
           for the vertices of the phase diagram. Using the Li-Fe-O phase
           diagram as an example, this means the lowest energy Li, Fe, and O
           phases.
        2. Calculate the formation energies from these elemental references for
           all entries. Exclude all positive formation energy ones from the
           data for convex hull.
        3. Generate the convex hull data.
        """
        logger.debug("Creating convex hull data...")
        #Determine the elemental references based on lowest energy for each.
        self._el_refs = dict()
        for entry in self._all_entries:
            if entry.composition.is_element:
                for el in entry.composition.elements:
                    if entry.composition[el] > Composition.amount_tolerance:
                        break
                e_per_atom = entry.energy_per_atom
                if el not in self._el_refs:
                    self._el_refs[el] = entry
                elif self._el_refs[el].energy_per_atom > e_per_atom:
                    self._el_refs[el] = entry
        # Remove positive formation energy entries
        entries_to_process = list()
        for entry in self._all_entries:
            if self.get_form_energy(entry) <= -self.formation_energy_tol:
                entries_to_process.append(entry)
            else:
                logger.debug("Removing positive formation energy entry " +
                             "{}".format(entry))
        entries_to_process.extend([entry for entry in self._el_refs.values()])

        self._qhull_entries = entries_to_process
        return self._process_entries_qhulldata(entries_to_process)

    def _make_phasediagram(self):
        """
        Make the phase diagram.
        """
        stable_entries = set()
        dim = len(self._elements)
        self._qhull_data = self._create_convhull_data()
        if len(self._qhull_data) == dim:
            self._facets = [range(len(self._elements))]
        else:
            self._facets = get_convex_hull(self._qhull_data,
                                           self._use_external_qhull)
            logger.debug("Final facets are\n{}".format(self._facets))

            logger.debug("Removing vertical facets...")
            finalfacets = list()
            for facet in self._facets:
                facetmatrix = np.zeros((len(facet), len(facet)))
                count = 0
                is_element_facet = True
                for vertex in facet:
                    facetmatrix[count] = np.array(self._qhull_data[vertex])
                    facetmatrix[count, dim - 1] = 1
                    count += 1
                    if len(self._qhull_entries[vertex].composition) > 1:
                        is_element_facet = False
                if abs(np.linalg.det(facetmatrix)) > 1e-8 and \
                   (not is_element_facet):
                    finalfacets.append(facet)
                else:
                    logger.debug("Removing vertical facet : {}".format(facet))
            self._facets = finalfacets

        for facet in self._facets:
            for vertex in facet:
                stable_entries.add(self._qhull_entries[vertex])
        self._stable_entries = stable_entries

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        symbols = [el.symbol for el in self._elements]
        output = []
        output.append("{} phase diagram".format("-".join(symbols)))
        output.append("{} stable phases: ".format(len(self._stable_entries)))
        output.append(", ".join([entry.name
                                 for entry in self._stable_entries]))
        return "\n".join(output)


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

    def __init__(self, entries, chempots, elements=None,
                 use_external_qhull=False):
        """
        Standard constructor for grand potential phase diagram.

        Args:
            entries:
                A list of PDEntry-like objects having an energy,
                energy_per_atom and composition.
            chempots:
                A dict of {element: float} to specify the chemical potentials
                of the open elements.
            elements:
                Optional list of elements in the phase diagram. If set to None,
                the elements are determined from the entries themselves.
            use_external_qhull:
                If set to True, the code will use an external command line call
                to Qhull to calculate the convex hull data instead of
                scipy.spatial.Delaunay. See the doc for the PhaseDiagram class
                for an explanation of the pros and cons.
        """
        if elements == None:
            elements = set()
            map(elements.update, [entry.composition.elements
                                  for entry in entries])
        allentries = list()
        for entry in entries:
            if not (entry.is_element and
                    (entry.composition.elements[0] in chempots)):
                allentries.append(GrandPotPDEntry(entry, chempots))
        self.chempots = chempots
        filteredels = list()
        for el in elements:
            if el not in chempots:
                filteredels.append(el)
        elements = sorted(filteredels)
        super(GrandPotentialPhaseDiagram, self).__init__(allentries, elements,
                                                         use_external_qhull)

    def __str__(self):
        output = []
        chemsys = "-".join([el.symbol for el in self._elements])
        output.append("{} grand potential phase diagram with ".format(chemsys))
        output[-1] += ", ".join(["u{}={}".format(el, v)
                                 for el, v in self.chempots.items()])
        output.append("{} stable phases: ".format(len(self._stable_entries)))
        output.append(", ".join([entry.name
                                 for entry in self._stable_entries]))
        return "\n".join(output)


class CompoundPhaseDiagram(PhaseDiagram):
    """
    Generates phase diagrams from compounds as terminations instead of
    elements.
    """

    """
    Tolerance for determining if amount of a composition is positive.
    """
    amount_tol = 1e-5

    def __init__(self, entries, terminal_compositions,
                 use_external_qhull=False):
        """
        Args:
            entries:
                Sequence of input entries. For example, if you want a Li2O-P2O5
                phase diagram, you might have all Li-P-O entries as an input.
            terminal_compositions:
                Terminal compositions of phase space. In the Li2O-P2O5 example,
                these will be the Li2O and P2O5 compositions.
            use_external_qhull:
                Similar to PhaseDiagram, set to True if you wish to use
                external qhull command.
        """
        self.original_entries = entries
        self.terminal_compositions = terminal_compositions
        (pentries, species_mapping) = \
            self.transform_entries(entries, terminal_compositions)
        PhaseDiagram.__init__(self, pentries,
                              elements=species_mapping.values(),
                              use_external_qhull=use_external_qhull)

    def transform_entries(self, entries, terminal_compositions):
        """
        Method to transform all entries to the composition coordinate in the
        terminal compositions. If the entry does not fall within the space
        defined by the terminal compositions, they are excluded. For example,
        Li3PO4 is mapped into a Li2O:1.5, P2O5:0.5 composition. The terminal
        compositions are represented by DummySpecies.

        Args:
            entries:
                Sequence of all input entries
            terminal_compositions:
                Terminal compositions of phase space.

        Returns:
            Sequence of TransformedPDEntries falling within the phase space.
        """
        new_entries = []

        fractional_comp = [c.get_fractional_composition()
                           for c in terminal_compositions]

        #Map terminal compositions to unique dummy species.
        sp_mapping = collections.OrderedDict()
        for i, comp in enumerate(fractional_comp):
            sp_mapping[comp] = DummySpecie("X" + chr(102 + i))

        for entry in entries:
            try:
                rxn = Reaction(fractional_comp, [entry.composition])
                rxn.normalize_to(entry.composition)
                #We only allow reactions that have positive amounts of
                #reactants.
                if all([rxn.get_coeff(comp) <= CompoundPhaseDiagram.amount_tol
                        for comp in fractional_comp]):
                    newcomp = {sp_mapping[comp]:-rxn.get_coeff(comp)
                               for comp in fractional_comp}
                    newcomp = {k: v for k, v in newcomp.items()
                               if v > CompoundPhaseDiagram.amount_tol}
                    transformed_entry = \
                        TransformedPDEntry(Composition(newcomp), entry)
                    new_entries.append(transformed_entry)
            except ReactionError:
                #If the reaction can't be balanced, the entry does not fall
                #into the phase space. We ignore them.
                pass
        return (new_entries, sp_mapping)
