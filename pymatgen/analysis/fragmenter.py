# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import logging

from monty.json import MSONable
from pymatgen.analysis.graphs import MoleculeGraph, MolGraphSplitError
from pymatgen.analysis.local_env import OpenBabelNN
from pymatgen.io.babel import BabelMolAdaptor


__author__ = "Samuel Blau"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"
__status__ = "Alpha"
__date__ = "9/7/18"


logger = logging.getLogger(__name__)


class Fragmenter(MSONable):

    def __init__(self, molecule, edges=None, depth=1, open_rings=True, opt_steps=10000):
        """
        Standard constructor for molecule fragmentation

        Args:
            molecule (Molecule): The molecule to fragment
            edges (list): List of index pairs that define graph edges, aka molecule bonds. If not
                          set, edges will be determined with OpenBabel.
            depth (int): The number of levels of iterative fragmentation to perform, where each
                         level will include fragments obtained by breaking one bond of a fragment
                         one level up. Defaults to 1. However, if set to 0, instead all possible
                         fragments are generated using an alternative, non-iterative scheme.
            open_rings (bool): Whether or not to open any rings encountered during fragmentation.
                               Defaults to False. If true, any bond that fails to yield disconnected
                               graphs when broken is instead removed and the entire structure is
                               optimized with OpenBabel in order to obtain a good initial guess for
                               an opened geometry that can then be put back into QChem to be
                               optimized without the ring just reforming.
            opt_steps (int): Number of optimization steps when opening rings. Defaults to 1000.

        """

        self.open_rings = open_rings
        self.opt_steps = opt_steps

        if edges is None:
            self.mol_graph = MoleculeGraph.with_local_env_strategy(molecule, OpenBabelNN(),
                                                                   reorder=False,
                                                                   extend_structure=False)
        else:
            edges = {(e[0], e[1]): None for e in edges}
            self.mol_graph = MoleculeGraph.with_edges(molecule, edges)

        self.unique_fragments = []
        self.unique_fragments_from_ring_openings = []

        if depth == 0: # Non-iterative, find all possible fragments:

            # Find all unique fragments besides those involving ring opening
            self.unique_fragments = self.mol_graph.build_unique_fragments()

            # Then, if self.open_rings is True, open all rings present in self.unique_fragments
            # in order to capture all unique fragments that require ring opening.
            if self.open_rings:
                self._open_all_rings()

        else: # Iterative fragment generation:
            self.fragments_by_level = {}

            # Loop through the number of levels,
            for level in range(depth):
                # If on the first level, perform one level of fragmentation on the principle molecule graph:
                if level == 0:
                    self.fragments_by_level["0"] = self._fragment_one_level([self.mol_graph])
                else:
                    if len(self.fragments_by_level[str(level-1)]) == 0:
                        # Nothing left to fragment, so exit the loop:
                        break
                    else: # If not on the first level, and there are fragments present in the previous level, then
                          # perform one level of fragmentation on all fragments present in the previous level:
                        self.fragments_by_level[str(level)] = self._fragment_one_level(self.fragments_by_level[str(level-1)])

    def _fragment_one_level(self, mol_graphs):
        """
        Perform one step of iterative fragmentation on a list of molecule graphs. Loop through the graphs,
        then loop through each graph's edges and attempt to remove that edge in order to obtain two
        disconnected subgraphs, aka two new fragments. If successful, check to see if the new fragments
        are already present in self.unique_fragments, and append them if not. If unsucessful, we know
        that edge belongs to a ring. If we are opening rings, do so with that bond, and then again
        check if the resulting fragment is present in self.unique_fragments and add it if it is not.
        """
        unique_fragments_on_this_level = []
        for mol_graph in mol_graphs:
            for edge in mol_graph.graph.edges:
                bond = [(edge[0],edge[1])]
                try:
                    fragments = mol_graph.split_molecule_subgraphs(bond, allow_reverse=True)
                    for fragment in fragments:
                        found = False
                        for unique_fragment in self.unique_fragments:
                            if unique_fragment.isomorphic_to(fragment):
                                found = True
                                break
                        if not found:
                            self.unique_fragments.append(fragment)
                            unique_fragments_on_this_level.append(fragment)
                except MolGraphSplitError:
                    if self.open_rings:
                        fragment = open_ring(mol_graph, bond, self.opt_steps)
                        found = False
                        for unique_fragment in self.unique_fragments:
                            if unique_fragment.isomorphic_to(fragment):
                                found = True
                                break
                        if not found:
                            self.unique_fragments.append(fragment)
                            self.unique_fragments_from_ring_openings.append(fragment)
                            unique_fragments_on_this_level.append(fragment)
        return unique_fragments_on_this_level

    def _open_all_rings(self):
        """
        Having already generated all unique fragments that did not require ring opening,
        now we want to also obtain fragments that do require opening. We achieve this by
        looping through all unique fragments and opening each bond present in any ring
        we find. We also temporarily add the principle molecule graph to self.unique_fragments
        so that its rings are opened as well.
        """
        self.unique_fragments.insert(0, self.mol_graph)
        for fragment in self.unique_fragments:
            ring_edges = fragment.find_rings()
            if ring_edges != []:
                for bond in ring_edges[0]:
                    new_fragment = open_ring(fragment, [bond], self.opt_steps)
                    found = False
                    for unique_fragment in self.unique_fragments:
                        if unique_fragment.isomorphic_to(new_fragment):
                            found = True
                            break
                    if not found:
                        self.unique_fragments_from_ring_openings.append(new_fragment)
                        self.unique_fragments.append(new_fragment)
        # Finally, remove the principle molecule graph:
        self.unique_fragments.pop(0)

def open_ring(mol_graph, bond, opt_steps):
    """
    Function to actually open a ring using OpenBabel's local opt. Given a molecule
    graph and a bond, convert the molecule graph into an OpenBabel molecule, remove
    the given bond, perform the local opt with the number of steps determined by
    self.steps, and then convert the resulting structure back into a molecule graph
    to be returned.
    """
    obmol = BabelMolAdaptor.from_molecule_graph(mol_graph)
    obmol.remove_bond(bond[0][0]+1, bond[0][1]+1)
    obmol.localopt(steps=opt_steps)
    return MoleculeGraph.with_local_env_strategy(obmol.pymatgen_mol, OpenBabelNN(), reorder=False, extend_structure=False)
