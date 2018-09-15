# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import logging

from monty.json import MSONable

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import build_MoleculeGraph, MolGraphSplitError
from pymatgen.analysis.local_env import OpenBabelNN
from pymatgen.analysis.fragmenter import open_ring
from pymatgen.io.babel import BabelMolAdaptor
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

    def __init__(self, molecule_entry, fragment_entries, allow_additional_charge_separation, multibreak):
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

        """

        self.molecule_entry = molecule_entry
        self.filter_fragment_entries(fragment_entries)
        print(str(len(self.filtered_entries)) + " filtered entries")
        self.bond_dissociation_energies = []
        self.done_frag_pairs = []
        self.ring_bonds = []

        if not allow_additional_charge_separation:
            if molecule_entry["final_molecule"]["charge"] == 0:
                self.expected_charges = [-1, 0, 1]
            elif molecule_entry["final_molecule"]["charge"] < 0:
                self.expected_charges = [molecule_entry["final_molecule"]["charge"], molecule_entry["final_molecule"]["charge"]+1]
            else:
                self.expected_charges = [molecule_entry["final_molecule"]["charge"]-1, molecule_entry["final_molecule"]["charge"]]
        else:
            if molecule_entry["final_molecule"]["charge"] == 0:
                self.expected_charges = [-2, -1, 0, 1, 2]
            elif molecule_entry["final_molecule"]["charge"] < 0:
                self.expected_charges = [molecule_entry["final_molecule"]["charge"]-1, molecule_entry["final_molecule"]["charge"], molecule_entry["final_molecule"]["charge"]+1, molecule_entry["final_molecule"]["charge"]+2]
            else:
                self.expected_charges = [molecule_entry["final_molecule"]["charge"]-2, molecule_entry["final_molecule"]["charge"]-1, molecule_entry["final_molecule"]["charge"], molecule_entry["final_molecule"]["charge"]+1]

        self.mol_graph = build_MoleculeGraph(Molecule.from_dict(molecule_entry["final_molecule"]),
                                             strategy=OpenBabelNN,
                                             reorder=False,
                                             extend_structure=False)
        for bond in self.mol_graph.graph.edges:
            bonds = [(bond[0],bond[1])]
            self.fragment_and_process(bonds)
        if multibreak:
            print("Breaking pairs of ring bonds. WARNING: Structure changes much more likely, meaning dissociation values are less reliable!")
            self.bond_pairs = []
            for ii,bond in enumerate(self.ring_bonds):
                for jj in range(ii+1,len(self.ring_bonds)):
                    bond_pair = [bond, self.ring_bonds[jj]]
                    self.bond_pairs += [bond_pair]
            for bond_pair in self.bond_pairs:
                self.fragment_and_process(bond_pair)

    def fragment_and_process(self, bonds):
        try:
            frags = self.mol_graph.split_molecule_subgraphs(bonds,allow_reverse=True)
            frag_success = True
        except MolGraphSplitError:
            if len(bonds) == 1:
                self.ring_bonds += bonds
                opened_frag = open_ring(self.mol_graph, bonds, 1000)
                # print(opened_frag)
                opened_entries = self.search_fragment_entries(opened_frag)
                if len(opened_entries) == 0:
                    print("Missing ring opening fragment resulting from the breakage of bond " + str(bonds[0][0]) + " " + str(bonds[0][1]))
                else:
                    print(len(opened_entries))
            elif len(bonds) == 2:
                if not multibreak:
                    raise RuntimeError("Should only be trying to break two bonds if multibreak is true! Exiting...")
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
                num_entries_for_this_frag_pair = 0
                frag1_entries = self.search_fragment_entries(frags[0])
                frag2_entries = self.search_fragment_entries(frags[1])
                frag1_charges_found = []
                frag2_charges_found = []
                for frag1 in frag1_entries[0] + frag1_entries[1]:
                    if frag1["initial_molecule"]["charge"] not in frag1_charges_found:
                        frag1_charges_found += [frag1["initial_molecule"]["charge"]]
                for frag2 in frag2_entries[0] + frag2_entries[1]:
                    if frag2["initial_molecule"]["charge"] not in frag2_charges_found:
                        frag2_charges_found += [frag2["initial_molecule"]["charge"]]
                if len(frag1_charges_found) < len(self.expected_charges):
                    bb = BabelMolAdaptor(frags[0].molecule)
                    pbmol = bb.pybel_mol
                    smiles = pbmol.write(str("smi")).split()[0]
                    for charge in self.expected_charges:
                        if charge not in frag1_charges_found:
                            print("Missing charge " + str(charge) + " for fragment " + smiles)
                if len(frag2_charges_found) < len(self.expected_charges):
                    bb = BabelMolAdaptor(frags[1].molecule)
                    pbmol = bb.pybel_mol
                    smiles = pbmol.write(str("smi")).split()[0]
                    for charge in self.expected_charges:
                        if charge not in frag2_charges_found:
                            print("Missing charge " + str(charge) + " for fragment " + smiles)
                for frag1 in frag1_entries[0]:
                    for frag2 in frag2_entries[0]:
                        if frag1["initial_molecule"]["charge"] + frag2["initial_molecule"]["charge"] == self.molecule_entry["final_molecule"]["charge"]:
                            self.bond_dissociation_energies += [self.build_new_entry(frag1, frag2, bonds)]
                            num_entries_for_this_frag_pair += 1
                if num_entries_for_this_frag_pair < len(self.expected_charges):
                    for frag1 in frag1_entries[0]:
                        for frag2 in frag2_entries[1]:
                            if frag1["initial_molecule"]["charge"] + frag2["initial_molecule"]["charge"] == self.molecule_entry["final_molecule"]["charge"]:
                                self.bond_dissociation_energies += [self.build_new_entry(frag1, frag2, bonds)]
                                num_entries_for_this_frag_pair += 1
                    for frag1 in frag1_entries[1]:
                        for frag2 in frag2_entries[0]:
                            if frag1["initial_molecule"]["charge"] + frag2["initial_molecule"]["charge"] == self.molecule_entry["final_molecule"]["charge"]:
                                self.bond_dissociation_energies += [self.build_new_entry(frag1, frag2, bonds)]
                                num_entries_for_this_frag_pair += 1

    def search_fragment_entries(self, frag):
        entries = []
        initial_entries = []
        final_entries = []
        for entry in self.filtered_entries:
            if frag.isomorphic_to(entry["initial_molgraph"]) and frag.isomorphic_to(entry["final_molgraph"]):
                entries += [entry]
            elif frag.isomorphic_to(entry["initial_molgraph"]):
                initial_entries += [entry]
            elif frag.isomorphic_to(entry["final_molgraph"]):
                final_entries += [entry]
        return [entries, initial_entries, final_entries]

    def filter_fragment_entries(self,fragment_entries):
        self.filtered_entries = []
        for entry in fragment_entries:
            entry["initial_molgraph"] = build_MoleculeGraph(Molecule.from_dict(entry["initial_molecule"]),
                                          strategy=OpenBabelNN,
                                          reorder=False,
                                          extend_structure=False)
            entry["final_molgraph"] = build_MoleculeGraph(Molecule.from_dict(entry["final_molecule"]),
                                        strategy=OpenBabelNN,
                                        reorder=False,
                                        extend_structure=False)
            if entry["initial_molgraph"].isomorphic_to(entry["final_molgraph"]):
                entry["structure_change"] = "no_change"
            else:
                initial_graph = entry["initial_molgraph"].graph
                final_graph = entry["final_molgraph"].graph
                if nx.is_connected(initial_graph.to_undirected()) and not nx.is_connected(final_graph.to_undirected()):
                    entry["structure_change"] = "unconnected_fragments"
                elif final_graph.number_of_edges() < initial_graph.number_of_edges():
                    entry["structure_change"] = "fewer_bonds"
                elif final_graph.number_of_edges() > initial_graph.number_of_edges():
                    entry["structure_change"] = "more_bonds"
                else:
                    entry["structure_change"] = "bond_change"
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

    def build_new_entry(self, frag1, frag2, bonds):
        specie = nx.get_node_attributes(self.mol_graph.graph, "specie")
        new_entry = [self.molecule_entry["final_energy"] - (frag1["final_energy"] + frag2["final_energy"]), bonds, specie[bonds[0][0]], specie[bonds[0][1]], frag1["smiles"], frag1["structure_change"], frag1["initial_molecule"]["charge"], frag1["initial_molecule"]["spin_multiplicity"], frag1["final_energy"], frag2["smiles"], frag2["structure_change"], frag2["initial_molecule"]["charge"], frag2["initial_molecule"]["spin_multiplicity"], frag2["final_energy"]]
        return new_entry
