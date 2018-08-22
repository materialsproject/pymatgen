# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import logging

from monty.json import MSONable

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import build_MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN
import networkx as nx


__author__ = "Samuel Blau"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"
__status__ = "Alpha"
__date__ = "7/26/18"


logger = logging.getLogger(__name__)


class BondDissociationEnergies(MSONable):

    def __init__(self, molecule_entry, fragment_entries):
        """
        Standard constructor for bond dissociation energies.

        Args:
            

        Get list of all bonds in target molecule
        Build molecule graph of target molecule
        Loop through bonds aka edges and attempt to split into disconnected subgraphs
            On success:
                Search through fragments: find entries with initial OR final molecules which are isomorphic
                Should be three entires per fragment: original charge, OC+1, OC-1
                If there are too few entries:
                    If one of the fragments is an H atom, you should only find two. H+ should be set (to what value???)
                    If fragment matches initial, but not final structure:
                        Just accept the rearrangement as part of the energetic cost of dissociation?
                        At the very least, flag this structure change in the eventual entry
                    If fragment matches final, but not initial structure:
                        Ignore unless no other fragment matches either?
                    If fragment is unstable:
                        Find its subfragments? Is that even sensical?
                        Just pass?
                    If fragment truly cannot be found, aka has not been calculated or calculate failed to complete:
                        Tell user to calculate it and pass
                If there are too many entries:
                    Get rid of duplicates
                        Going to be annoying given time stamps, rounding errors, etc
                    Remove those whose initial rems don't conform to current defaults - Happens before passing
                    If still more than three, for now, isolate differences, report to user, and grab most recent
                    Just take more negative?
                Grab energies of each fragment at all three charges
                If target molecule is neutral:
                    Calculate E_frag1(0)  + E_frag2(0)
                    Calculate E_frag1(+1) + E_frag2(-1)
                    Calculate E_frag1(-1) + E_frag2(+1)
                If target molecule is +1:
                    Calculate E_frag1(0)  + E_frag2(+1)
                    Calculate E_frag1(+1) + E_frag2(0)
                    Calculate E_frag1(+2) + E_frag2(-1)?
                    Calculate E_frag1(-1) + E_frag2(+2)?
                If target molecule is -1:
                    Calculate E_frag1(0)  + E_frag2(-1)
                    Calculate E_frag1(-1) + E_frag2(0)
                    Calculate E_frag1(-2) + E_frag2(+1)?
                    Calculate E_frag1(+1) + E_frag2(-2)?
                Don't limit - just make final charge add to original charge
                Save E_molecule - calculated values associated with bond indices & charges in some data structure
            On failure:
                Track bonds that do not break to yield subgraphs as ring bonds
        Make all possible pairs of ring bonds
        Loop through ring bond pairs and attempt to split into disconnected subgraphs
            Same success/failure procedure as above. Should be wrapped in a function.
        Resolve individual rings
        Attempt to solve for single ring bond dissociation energies, careful with charge

        """

        self.molecule_entry = molecule_entry
        self.filter_fragment_entries(fragment_entries)
        print(str(len(self.filtered_entries)) + " filtered entries")
        self.bond_dissociation_energies = []
        self.done_frag_pairs = []
        self.ring_bonds = []
        self.bad_pairs = []
        self.mol_graph = build_MoleculeGraph(Molecule.from_dict(molecule_entry["final_molecule"]),
                                             strategy=OpenBabelNN,
                                             reorder=False,
                                             extend_structure=False)
        for bond in self.mol_graph.graph.edges:
            bonds = [(bond[0],bond[1])]
            self.fragment_and_process(bonds)
        self.bond_pairs = []
        for ii,bond in enumerate(self.ring_bonds):
            for jj in range(ii+1,len(self.ring_bonds)):
                bond_pair = [bond, self.ring_bonds[jj]]
                self.bond_pairs += [bond_pair]
        if len(self.bond_pairs) > 0:
            print("Ring bonds detected!")
        for bond_pair in self.bond_pairs:
            self.fragment_and_process(bond_pair)
        self.solve_ring_bonds()

    def fragment_and_process(self, bonds):
        try:
            frags = self.mol_graph.split_molecule_subgraphs(bonds,allow_reverse=True)
            frag_success = True
        except RuntimeError:
            if len(bonds) == 1:
                self.ring_bonds += bonds
            elif len(bonds) == 2:
                self.bad_pairs += bonds
            else:
                print('No reason to try and break more than two bonds at once! Exiting...')
                raise ValueError
            frag_success = False
        if frag_success:
            frags_done = False
            for frag_pair in self.done_frag_pairs:
                if frag_pair[0].isomorphic_to(frags[0]):
                    if frag_pair[1].isomorphic_to(frags[1]):
                        frags_done = True
                        break
                elif frag_pair[1].isomorphic_to(frags[0]):
                    if frag_pair[0].isomorphic_to(frags[1]):
                        frags_done = True
                        break
            if not frags_done:
                self.done_frag_pairs += [frags]
                frag1_entries = self.search_fragment_entries(frags[0])
                frag2_entries = self.search_fragment_entries(frags[1])
                for frag1 in frag1_entries[0]:
                    for frag2 in frag2_entries[0]:
                        if frag1["initial_molecule"]["charge"] + frag2["initial_molecule"]["charge"] == self.molecule_entry["final_molecule"]["charge"]:
                            coords = nx.get_node_attributes(self.mol_graph.graph, "coords")
                            specie = nx.get_node_attributes(self.mol_graph.graph, "specie")
                            frag1_charge = frag1["initial_molecule"]["charge"]
                            frag2_charge = frag2["initial_molecule"]["charge"]
                            if frag1["final_energy"] == "unstable" or frag2["final_energy"] == "unstable":
                                new_entry = ["unstable", bonds, specie[bonds[0][0]], specie[bonds[0][1]], coords[bonds[0][0]], coords[bonds[0][1]], frag1["smiles"], frag1_charge, frag1["final_energy"], frag2["smiles"], frag2_charge, frag2["final_energy"]]
                            else:
                                new_entry = [self.molecule_entry["final_energy"] - (frag1["final_energy"] + frag2["final_energy"]), bonds, specie[bonds[0][0]], specie[bonds[0][1]], coords[bonds[0][0]], coords[bonds[0][1]], frag1["smiles"], frag1_charge, frag1["final_energy"], frag2["smiles"], frag2_charge, frag2["final_energy"]]
                            self.bond_dissociation_energies += [new_entry]

    def search_fragment_entries(self, frag):
        entries = []
        initial_entries = []
        final_entries = []
        for entry in self.filtered_entries:
            if frag.isomorphic_to(entry["initial_molgraph"]) and frag.isomorphic_to(entry["final_molgraph"]):
                entries += [entry]
            elif frag.isomorphic_to(entry["initial_molgraph"]):
                entries += [entry]
                initial_entries += [entry]
            elif frag.isomorphic_to(entry["final_molgraph"]):
                final_entries += [entry]
        return [entries, initial_entries, final_entries]

    def filter_fragment_entries(self,fragment_entries):
        self.filtered_entries = []
        for entry in fragment_entries:
            print(len(self.filtered_entries))
            entry["initial_molgraph"] = build_MoleculeGraph(Molecule.from_dict(entry["initial_molecule"]),
                                          strategy=OpenBabelNN,
                                          reorder=False,
                                          extend_structure=False)
            entry["final_molgraph"] = build_MoleculeGraph(Molecule.from_dict(entry["final_molecule"]),
                                        strategy=OpenBabelNN,
                                        reorder=False,
                                        extend_structure=False)
            found_similar_entry = False
            for ii,filtered_entry in enumerate(self.filtered_entries):
                if filtered_entry["formula_pretty"] == entry["formula_pretty"]:
                    if filtered_entry["initial_molgraph"].isomorphic_to(entry["initial_molgraph"]) and filtered_entry["final_molgraph"].isomorphic_to(entry["final_molgraph"]) and filtered_entry["initial_molecule"]["charge"] == entry["initial_molecule"]["charge"]:
                        found_similar_entry = True
                        if entry["final_energy"] < filtered_entry["final_energy"]:
                            self.filtered_entries[ii] = entry
                        break
            if not found_similar_entry:
                self.filtered_entries += [entry]

    def solve_ring_bonds(self):
        pass
