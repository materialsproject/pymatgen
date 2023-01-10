# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Perform fragmentation of molecules.
"""

from __future__ import annotations

import copy
import logging

from monty.json import MSONable

from pymatgen.analysis.graphs import MoleculeGraph, MolGraphSplitError
from pymatgen.analysis.local_env import OpenBabelNN, metal_edge_extender
from pymatgen.io.babel import BabelMolAdaptor

__author__ = "Samuel Blau"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "2.0"
__maintainer__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"
__status__ = "Beta"
__date__ = "8/21/19"

logger = logging.getLogger(__name__)


class Fragmenter(MSONable):
    """
    Molecule fragmenter class.
    """

    def __init__(
        self,
        molecule,
        edges=None,
        depth=1,
        open_rings=False,
        use_metal_edge_extender=False,
        opt_steps=10000,
        prev_unique_frag_dict=None,
        assume_previous_thoroughness=True,
    ):
        """
        Standard constructor for molecule fragmentation

        Args:
            molecule (Molecule): The molecule to fragment.
            edges (list): List of index pairs that define graph edges, aka molecule bonds. If not set,
                edges will be determined with OpenBabel. Defaults to None.
            depth (int): The number of levels of iterative fragmentation to perform, where each level
                will include fragments obtained by breaking one bond of a fragment one level up.
                Defaults to 1. However, if set to 0, instead all possible fragments are generated
                using an alternative, non-iterative scheme.
            open_rings (bool): Whether or not to open any rings encountered during fragmentation.
                Defaults to False. If true, any bond that fails to yield disconnected graphs when
                broken is instead removed and the entire structure is optimized with OpenBabel in
                order to obtain a good initial guess for an opened geometry that can then be put
                back into QChem to be optimized without the ring just reforming.
            use_metal_edge_extender (bool): Whether or not to attempt to add additional edges from
                O, N, F, or Cl to any Li or Mg atoms present that OpenBabel may have missed. Defaults
                to False. Most important for ionic bonding. Note that additional metal edges may yield
                new "rings" (e.g. -C-O-Li-O- in LiEC) that will not play nicely with ring opening.
            opt_steps (int): Number of optimization steps when opening rings. Defaults to 10000.
            prev_unique_frag_dict (dict): A dictionary of previously identified unique fragments.
                Defaults to None. Typically only used when trying to find the set of unique fragments
                that come from multiple molecules.
            assume_previous_thoroughness (bool): Whether or not to assume that a molecule / fragment
                provided in prev_unique_frag_dict has all of its unique subfragments also provided in
                prev_unique_frag_dict. Defaults to True. This is an essential optimization when trying
                to find the set of unique fragments that come from multiple molecules if all of those
                molecules are being fully iteratively fragmented. However, if you're passing a
                prev_unique_frag_dict which includes a molecule and its fragments that were generated
                at insufficient depth to find all possible subfragments to a fragmentation calculation
                of a different molecule that you aim to find all possible subfragments of and which has
                common subfragments with the previous molecule, this optimization will cause you to
                miss some unique subfragments.
        """
        self.assume_previous_thoroughness = assume_previous_thoroughness
        self.open_rings = open_rings
        self.opt_steps = opt_steps

        if edges is None:
            self.mol_graph = MoleculeGraph.with_local_env_strategy(molecule, OpenBabelNN())
        else:
            edges = {(e[0], e[1]): None for e in edges}
            self.mol_graph = MoleculeGraph.with_edges(molecule, edges)

        if ("Li" in molecule.composition or "Mg" in molecule.composition) and use_metal_edge_extender:
            self.mol_graph = metal_edge_extender(self.mol_graph)

        self.prev_unique_frag_dict = prev_unique_frag_dict or {}
        self.new_unique_frag_dict = {}  # new fragments from the given molecule not contained in prev_unique_frag_dict
        self.all_unique_frag_dict = {}  # all fragments from just the given molecule
        self.unique_frag_dict = {}  # all fragments from both the given molecule and prev_unique_frag_dict

        if depth == 0:  # Non-iterative, find all possible fragments:

            # Find all unique fragments besides those involving ring opening
            self.all_unique_frag_dict = self.mol_graph.build_unique_fragments()

            # Then, if self.open_rings is True, open all rings present in self.unique_fragments
            # in order to capture all unique fragments that require ring opening.
            if self.open_rings:
                self._open_all_rings()

        else:  # Iterative fragment generation:
            self.fragments_by_level = {}

            # Loop through the number of levels,
            for level in range(depth):
                # If on the first level, perform one level of fragmentation on the principle molecule graph:
                if level == 0:
                    self.fragments_by_level["0"] = self._fragment_one_level(
                        {
                            str(self.mol_graph.molecule.composition.alphabetical_formula)
                            + " E"
                            + str(len(self.mol_graph.graph.edges())): [self.mol_graph]
                        }
                    )
                else:
                    num_frags_prev_level = 0
                    for key in self.fragments_by_level[str(level - 1)]:
                        num_frags_prev_level += len(self.fragments_by_level[str(level - 1)][key])
                    if num_frags_prev_level == 0:
                        # Nothing left to fragment, so exit the loop:
                        break
                    # If not on the first level, and there are fragments present in the previous level, then
                    # perform one level of fragmentation on all fragments present in the previous level:
                    self.fragments_by_level[str(level)] = self._fragment_one_level(
                        self.fragments_by_level[str(level - 1)]
                    )

        if self.prev_unique_frag_dict == {}:
            self.new_unique_frag_dict = copy.deepcopy(self.all_unique_frag_dict)
        else:
            for frag_key in self.all_unique_frag_dict:
                if frag_key not in self.prev_unique_frag_dict:
                    self.new_unique_frag_dict[frag_key] = copy.deepcopy(self.all_unique_frag_dict[frag_key])
                else:
                    for fragment in self.all_unique_frag_dict[frag_key]:
                        found = False
                        for prev_frag in self.prev_unique_frag_dict[frag_key]:
                            if fragment.isomorphic_to(prev_frag):
                                found = True
                        if not found:
                            if frag_key not in self.new_unique_frag_dict:
                                self.new_unique_frag_dict[frag_key] = [fragment]
                            else:
                                self.new_unique_frag_dict[frag_key].append(fragment)

        self.new_unique_fragments = 0
        for frag_key in self.new_unique_frag_dict:
            self.new_unique_fragments += len(self.new_unique_frag_dict[frag_key])

        if self.prev_unique_frag_dict == {}:
            self.unique_frag_dict = self.new_unique_frag_dict
            self.total_unique_fragments = self.new_unique_fragments
        else:
            self.unique_frag_dict = copy.deepcopy(self.prev_unique_frag_dict)
            for frag_key in self.new_unique_frag_dict:
                if frag_key in self.unique_frag_dict:
                    for new_frag in self.new_unique_frag_dict[frag_key]:
                        self.unique_frag_dict[frag_key].append(new_frag)
                else:
                    self.unique_frag_dict[frag_key] = copy.deepcopy(self.new_unique_frag_dict[frag_key])

            self.total_unique_fragments = 0
            for frag_key in self.unique_frag_dict:
                self.total_unique_fragments += len(self.unique_frag_dict[frag_key])

    def _fragment_one_level(self, old_frag_dict):
        """
        Perform one step of iterative fragmentation on a list of molecule graphs. Loop through the graphs,
        then loop through each graph's edges and attempt to remove that edge in order to obtain two
        disconnected subgraphs, aka two new fragments. If successful, check to see if the new fragments
        are already present in self.unique_fragments, and append them if not. If unsuccessful, we know
        that edge belongs to a ring. If we are opening rings, do so with that bond, and then again
        check if the resulting fragment is present in self.unique_fragments and add it if it is not.
        """
        new_frag_dict = {}
        for old_frag_key in old_frag_dict:
            for old_frag in old_frag_dict[old_frag_key]:
                for edge in old_frag.graph.edges:
                    bond = [(edge[0], edge[1])]
                    fragments = []
                    try:
                        fragments = old_frag.split_molecule_subgraphs(bond, allow_reverse=True)
                    except MolGraphSplitError:
                        if self.open_rings:
                            fragments = [open_ring(old_frag, bond, self.opt_steps)]
                    for fragment in fragments:
                        new_frag_key = (
                            str(fragment.molecule.composition.alphabetical_formula)
                            + " E"
                            + str(len(fragment.graph.edges()))
                        )
                        proceed = True
                        if self.assume_previous_thoroughness and self.prev_unique_frag_dict != {}:
                            if new_frag_key in self.prev_unique_frag_dict:
                                for unique_fragment in self.prev_unique_frag_dict[new_frag_key]:
                                    if unique_fragment.isomorphic_to(fragment):
                                        proceed = False
                                        break
                        if proceed:
                            if new_frag_key not in self.all_unique_frag_dict:
                                self.all_unique_frag_dict[new_frag_key] = [fragment]
                                new_frag_dict[new_frag_key] = [fragment]
                            else:
                                found = False
                                for unique_fragment in self.all_unique_frag_dict[new_frag_key]:
                                    if unique_fragment.isomorphic_to(fragment):
                                        found = True
                                        break
                                if not found:
                                    self.all_unique_frag_dict[new_frag_key].append(fragment)
                                    if new_frag_key in new_frag_dict:
                                        new_frag_dict[new_frag_key].append(fragment)
                                    else:
                                        new_frag_dict[new_frag_key] = [fragment]
        return new_frag_dict

    def _open_all_rings(self):
        """
        Having already generated all unique fragments that did not require ring opening,
        now we want to also obtain fragments that do require opening. We achieve this by
        looping through all unique fragments and opening each bond present in any ring
        we find. We also temporarily add the principle molecule graph to self.unique_fragments
        so that its rings are opened as well.
        """
        mol_key = (
            str(self.mol_graph.molecule.composition.alphabetical_formula)
            + " E"
            + str(len(self.mol_graph.graph.edges()))
        )
        self.all_unique_frag_dict[mol_key] = [self.mol_graph]
        new_frag_keys = {"0": []}
        new_frag_key_dict = {}
        for key in self.all_unique_frag_dict:
            for fragment in self.all_unique_frag_dict[key]:
                ring_edges = fragment.find_rings()
                if ring_edges != []:
                    for bond in ring_edges[0]:
                        new_fragment = open_ring(fragment, [bond], self.opt_steps)
                        frag_key = (
                            str(new_fragment.molecule.composition.alphabetical_formula)
                            + " E"
                            + str(len(new_fragment.graph.edges()))
                        )
                        if frag_key not in self.all_unique_frag_dict:
                            if frag_key not in new_frag_keys["0"]:
                                new_frag_keys["0"].append(copy.deepcopy(frag_key))
                                new_frag_key_dict[frag_key] = copy.deepcopy([new_fragment])
                            else:
                                found = False
                                for unique_fragment in new_frag_key_dict[frag_key]:
                                    if unique_fragment.isomorphic_to(new_fragment):
                                        found = True
                                        break
                                if not found:
                                    new_frag_key_dict[frag_key].append(copy.deepcopy(new_fragment))
                        else:
                            found = False
                            for unique_fragment in self.all_unique_frag_dict[frag_key]:
                                if unique_fragment.isomorphic_to(new_fragment):
                                    found = True
                                    break
                            if not found:
                                self.all_unique_frag_dict[frag_key].append(copy.deepcopy(new_fragment))
        for key, value in new_frag_key_dict.items():
            self.all_unique_frag_dict[key] = copy.deepcopy(value)
        idx = 0
        while len(new_frag_keys[str(idx)]) != 0:
            new_frag_key_dict = {}
            idx += 1
            new_frag_keys[str(idx)] = []
            for key in new_frag_keys[str(idx - 1)]:
                for fragment in self.all_unique_frag_dict[key]:
                    ring_edges = fragment.find_rings()
                    if ring_edges != []:
                        for bond in ring_edges[0]:
                            new_fragment = open_ring(fragment, [bond], self.opt_steps)
                            frag_key = (
                                str(new_fragment.molecule.composition.alphabetical_formula)
                                + " E"
                                + str(len(new_fragment.graph.edges()))
                            )
                            if frag_key not in self.all_unique_frag_dict:
                                if frag_key not in new_frag_keys[str(idx)]:
                                    new_frag_keys[str(idx)].append(copy.deepcopy(frag_key))
                                    new_frag_key_dict[frag_key] = copy.deepcopy([new_fragment])
                                else:
                                    found = False
                                    for unique_fragment in new_frag_key_dict[frag_key]:
                                        if unique_fragment.isomorphic_to(new_fragment):
                                            found = True
                                            break
                                    if not found:
                                        new_frag_key_dict[frag_key].append(copy.deepcopy(new_fragment))
                            else:
                                found = False
                                for unique_fragment in self.all_unique_frag_dict[frag_key]:
                                    if unique_fragment.isomorphic_to(new_fragment):
                                        found = True
                                        break
                                if not found:
                                    self.all_unique_frag_dict[frag_key].append(copy.deepcopy(new_fragment))
            for key, value in new_frag_key_dict.items():
                self.all_unique_frag_dict[key] = copy.deepcopy(value)
        self.all_unique_frag_dict.pop(mol_key)


def open_ring(mol_graph, bond, opt_steps):
    """
    Function to actually open a ring using OpenBabel's local opt. Given a molecule
    graph and a bond, convert the molecule graph into an OpenBabel molecule, remove
    the given bond, perform the local opt with the number of steps determined by
    self.steps, and then convert the resulting structure back into a molecule graph
    to be returned.
    """
    obmol = BabelMolAdaptor.from_molecule_graph(mol_graph)
    obmol.remove_bond(bond[0][0] + 1, bond[0][1] + 1)
    obmol.localopt(steps=opt_steps, forcefield="uff")
    return MoleculeGraph.with_local_env_strategy(obmol.pymatgen_mol, OpenBabelNN())
