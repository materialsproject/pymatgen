# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Module for BondDissociationEnergies.
"""

from __future__ import annotations

import logging

import networkx as nx
from monty.json import MSONable

from pymatgen.analysis.fragmenter import open_ring
from pymatgen.analysis.graphs import MoleculeGraph, MolGraphSplitError
from pymatgen.analysis.local_env import OpenBabelNN
from pymatgen.core.structure import Molecule
from pymatgen.io.babel import BabelMolAdaptor

__author__ = "Samuel Blau"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"
__status__ = "Alpha"
__date__ = "7/26/18"

logger = logging.getLogger(__name__)


class BondDissociationEnergies(MSONable):
    """
    Standard constructor for bond dissociation energies. All bonds in the principle molecule are
    looped through and their dissociation energies are calculated given the energies of the resulting
    fragments, or, in the case of a ring bond, from the energy of the molecule obtained from breaking
    the bond and opening the ring. This class should only be called after the energies of the optimized
    principle molecule and all relevant optimized fragments have been determined, either from quantum
    chemistry or elsewhere. It was written to provide the analysis after running an Atomate fragmentation
    workflow.
    """

    def __init__(
        self,
        molecule_entry,
        fragment_entries,
        allow_additional_charge_separation=False,
        multibreak=False,
    ):
        """
        Note that the entries passed by the user must have the following keys: formula_pretty, initial_molecule,
        final_molecule. If a PCM is present, all entries should also have a pcm_dielectric key.

        Args:
            molecule_entry (dict): Entry for the principle molecule. Should have the keys mentioned above.
            fragment_entries (list of dicts): List of fragment entries. Each should have the keys mentioned above.
            allow_additional_charge_separation (bool): If True, consider larger than normal charge separation
                                                       among fragments. Defaults to False. See the definition
                                                       of self.expected_charges below for more specific information.
            multibreak (bool): If True, additionally attempt to break pairs of bonds. Defaults to False.
        """

        self.molecule_entry = molecule_entry
        self.filter_fragment_entries(fragment_entries)
        print(str(len(self.filtered_entries)) + " filtered entries")
        self.bond_dissociation_energies = []
        self.done_frag_pairs = []
        self.done_RO_frags = []
        self.ring_bonds = []

        required_keys = ["formula_pretty", "initial_molecule", "final_molecule"]
        if "pcm_dielectric" in self.molecule_entry:
            required_keys.append("pcm_dielectric")
        for key in required_keys:
            if key not in self.molecule_entry:
                raise RuntimeError(key + " must be present in molecule entry! Exiting...")
            for entry in self.filtered_entries:
                if key not in entry:
                    raise RuntimeError(key + " must be present in all fragment entries! Exiting...")

        # Define expected charges
        if not allow_additional_charge_separation:
            if molecule_entry["final_molecule"]["charge"] == 0:
                self.expected_charges = [-1, 0, 1]
            elif molecule_entry["final_molecule"]["charge"] < 0:
                self.expected_charges = [
                    molecule_entry["final_molecule"]["charge"],
                    molecule_entry["final_molecule"]["charge"] + 1,
                ]
            else:
                self.expected_charges = [
                    molecule_entry["final_molecule"]["charge"] - 1,
                    molecule_entry["final_molecule"]["charge"],
                ]
        else:
            if molecule_entry["final_molecule"]["charge"] == 0:
                self.expected_charges = [-2, -1, 0, 1, 2]
            elif molecule_entry["final_molecule"]["charge"] < 0:
                self.expected_charges = [
                    molecule_entry["final_molecule"]["charge"] - 1,
                    molecule_entry["final_molecule"]["charge"],
                    molecule_entry["final_molecule"]["charge"] + 1,
                    molecule_entry["final_molecule"]["charge"] + 2,
                ]
            else:
                self.expected_charges = [
                    molecule_entry["final_molecule"]["charge"] - 2,
                    molecule_entry["final_molecule"]["charge"] - 1,
                    molecule_entry["final_molecule"]["charge"],
                    molecule_entry["final_molecule"]["charge"] + 1,
                ]

        # Build principle molecule graph
        self.mol_graph = MoleculeGraph.with_local_env_strategy(
            Molecule.from_dict(molecule_entry["final_molecule"]), OpenBabelNN()
        )
        # Loop through bonds, aka graph edges, and fragment and process:
        for bond in self.mol_graph.graph.edges:
            bonds = [(bond[0], bond[1])]
            self.fragment_and_process(bonds)
        # If mulitbreak, loop through pairs of ring bonds.
        if multibreak:
            print(
                "Breaking pairs of ring bonds. WARNING: Structure changes much more likely, meaning dissociation values"
                " are less reliable! This is a bad idea!"
            )
            self.bond_pairs = []
            for ii, bond in enumerate(self.ring_bonds):
                for jj in range(ii + 1, len(self.ring_bonds)):
                    bond_pair = [bond, self.ring_bonds[jj]]
                    self.bond_pairs += [bond_pair]
            for bond_pair in self.bond_pairs:
                self.fragment_and_process(bond_pair)

    def fragment_and_process(self, bonds):
        """
        Fragment and process bonds.

        :param bonds: Bonds to process.
        :return:
        """
        # Try to split the principle:
        try:
            frags = self.mol_graph.split_molecule_subgraphs(bonds, allow_reverse=True)
            frag_success = True
        except MolGraphSplitError:
            # If split is unsuccessful, then we have encountered a ring bond
            if len(bonds) == 1:
                self.ring_bonds += bonds
                # So we open the ring and make sure we haven't already encountered an identically opened fragment:
                RO_frag = open_ring(self.mol_graph, bonds, 1000)
                frag_done = False
                for done_RO_frag in self.done_RO_frags:
                    if RO_frag.isomorphic_to(done_RO_frag):
                        frag_done = True
                if not frag_done:
                    # If this is a new fragment, save the record and then search for relevant fragment entries:
                    self.done_RO_frags.append(RO_frag)
                    opened_entries = self.search_fragment_entries(RO_frag)
                    good_entries = []
                    # We will start by looking at entries with no structure change
                    for frag in opened_entries[0]:  # 0 -> no structural change
                        # Since a ring opening still yields a single molecule, it should have the same charge as the
                        # principle:
                        if frag["initial_molecule"]["charge"] == self.molecule_entry["final_molecule"]["charge"]:
                            good_entries.append(frag)
                    # If we didn't find any good entries, let's also look at those that exhibit structural changes:
                    if len(good_entries) == 0:
                        for frag in opened_entries[1]:  # 1 -> YES structural change
                            if frag["initial_molecule"]["charge"] == self.molecule_entry["final_molecule"]["charge"]:
                                good_entries.append(frag)
                    # If we still have no good entries, something must have gone wrong with the calculations:
                    if len(good_entries) == 0:
                        bb = BabelMolAdaptor.from_molecule_graph(RO_frag)
                        pbmol = bb.pybel_mol
                        smiles = pbmol.write("smi").split()[0]
                        specie = nx.get_node_attributes(self.mol_graph.graph, "specie")
                        print(
                            "Missing ring opening fragment resulting from the breakage of "
                            + specie[bonds[0][0]]
                            + " "
                            + specie[bonds[0][1]]
                            + " bond "
                            + str(bonds[0][0])
                            + " "
                            + str(bonds[0][1])
                            + " which would yield a molecule with this SMILES string: "
                            + smiles
                        )
                    elif len(good_entries) == 1:
                        # If we have only one good entry, format it and add it to the list that will eventually return:
                        self.bond_dissociation_energies += [self.build_new_entry(good_entries, bonds)]
                    else:
                        # We shouldn't ever encounter more than one good entry.
                        raise RuntimeError("There should only be one valid ring opening fragment! Exiting...")
            elif len(bonds) == 2:
                raise RuntimeError("Should only be trying to break two bonds if multibreak is true! Exiting...")
            else:
                print("No reason to try and break more than two bonds at once! Exiting...")
                raise ValueError
            frag_success = False
        if frag_success:
            # If the principle did successfully split, then we aren't dealing with a ring bond.
            # As above, we begin by making sure we haven't already encountered an identical pair of fragments:
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
                # If we haven't, we save this pair and search for the relevant fragment entries:
                self.done_frag_pairs += [frags]
                num_entries_for_this_frag_pair = 0
                frag1_entries = self.search_fragment_entries(frags[0])
                frag2_entries = self.search_fragment_entries(frags[1])
                frag1_charges_found = []
                frag2_charges_found = []
                # We then check for our expected charges of each fragment:
                for frag1 in frag1_entries[0] + frag1_entries[1]:
                    if frag1["initial_molecule"]["charge"] not in frag1_charges_found:
                        frag1_charges_found += [frag1["initial_molecule"]["charge"]]
                for frag2 in frag2_entries[0] + frag2_entries[1]:
                    if frag2["initial_molecule"]["charge"] not in frag2_charges_found:
                        frag2_charges_found += [frag2["initial_molecule"]["charge"]]
                # If we're missing some of either, tell the user:
                if len(frag1_charges_found) < len(self.expected_charges):
                    bb = BabelMolAdaptor(frags[0].molecule)
                    pbmol = bb.pybel_mol
                    smiles = pbmol.write("smi").split()[0]
                    for charge in self.expected_charges:
                        if charge not in frag1_charges_found:
                            print("Missing charge " + str(charge) + " for fragment " + smiles)
                if len(frag2_charges_found) < len(self.expected_charges):
                    bb = BabelMolAdaptor(frags[1].molecule)
                    pbmol = bb.pybel_mol
                    smiles = pbmol.write("smi").split()[0]
                    for charge in self.expected_charges:
                        if charge not in frag2_charges_found:
                            print("Missing charge " + str(charge) + " for fragment " + smiles)
                # Now we attempt to pair fragments with the right total charge, starting with only fragments with no
                # structural change:
                for frag1 in frag1_entries[0]:  # 0 -> no structural change
                    for frag2 in frag2_entries[0]:  # 0 -> no structural change
                        if (
                            frag1["initial_molecule"]["charge"] + frag2["initial_molecule"]["charge"]
                            == self.molecule_entry["final_molecule"]["charge"]
                        ):
                            self.bond_dissociation_energies += [self.build_new_entry([frag1, frag2], bonds)]
                            num_entries_for_this_frag_pair += 1
                # If we haven't found the number of fragment pairs that we expect, we expand our search to include
                # fragments that do exhibit structural change:
                if num_entries_for_this_frag_pair < len(self.expected_charges):
                    for frag1 in frag1_entries[0]:  # 0 -> no structural change
                        for frag2 in frag2_entries[1]:  # 1 -> YES structural change
                            if (
                                frag1["initial_molecule"]["charge"] + frag2["initial_molecule"]["charge"]
                                == self.molecule_entry["final_molecule"]["charge"]
                            ):
                                self.bond_dissociation_energies += [self.build_new_entry([frag1, frag2], bonds)]
                                num_entries_for_this_frag_pair += 1
                    for frag1 in frag1_entries[1]:  # 1 -> YES structural change
                        for frag2 in frag2_entries[0]:  # 0 -> no structural change
                            if (
                                frag1["initial_molecule"]["charge"] + frag2["initial_molecule"]["charge"]
                                == self.molecule_entry["final_molecule"]["charge"]
                            ):
                                self.bond_dissociation_energies += [self.build_new_entry([frag1, frag2], bonds)]
                                num_entries_for_this_frag_pair += 1

    def search_fragment_entries(self, frag):
        """
        Search all fragment entries for those isomorphic to the given fragment.
        We distinguish between entries where both initial and final molgraphs are isomorphic to the
        given fragment (entries) vs those where only the initial molgraph is isomorphic to the given
        fragment (initial_entries) vs those where only the final molgraph is isomorphic (final_entries)

        Args:
            frag: Fragment
        """
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

    def filter_fragment_entries(self, fragment_entries):
        """
        Filter the fragment entries.

        :param fragment_entries:
        :return:
        """
        self.filtered_entries = []
        for entry in fragment_entries:
            # Check and make sure that PCM dielectric is consistent with principle:
            if "pcm_dielectric" in self.molecule_entry:
                if "pcm_dielectric" not in entry:
                    raise RuntimeError(
                        "Principle molecule has a PCM dielectric of "
                        + str(self.molecule_entry["pcm_dielectric"])
                        + " but a fragment entry has no PCM dielectric! Please only pass fragment entries"
                        " with PCM details consistent with the principle entry. Exiting..."
                    )
                if entry["pcm_dielectric"] != self.molecule_entry["pcm_dielectric"]:
                    raise RuntimeError(
                        "Principle molecule has a PCM dielectric of "
                        + str(self.molecule_entry["pcm_dielectric"])
                        + " but a fragment entry has a different PCM dielectric! Please only pass"
                        " fragment entries with PCM details consistent with the principle entry."
                        " Exiting..."
                    )
            # Build initial and final molgraphs:
            entry["initial_molgraph"] = MoleculeGraph.with_local_env_strategy(
                Molecule.from_dict(entry["initial_molecule"]), OpenBabelNN()
            )
            entry["final_molgraph"] = MoleculeGraph.with_local_env_strategy(
                Molecule.from_dict(entry["final_molecule"]), OpenBabelNN()
            )
            # Classify any potential structural change that occurred during optimization:
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
            # Check for uniqueness
            for ii, filtered_entry in enumerate(self.filtered_entries):
                if filtered_entry["formula_pretty"] == entry["formula_pretty"]:
                    if (
                        filtered_entry["initial_molgraph"].isomorphic_to(entry["initial_molgraph"])
                        and filtered_entry["final_molgraph"].isomorphic_to(entry["final_molgraph"])
                        and filtered_entry["initial_molecule"]["charge"] == entry["initial_molecule"]["charge"]
                    ):
                        found_similar_entry = True
                        # If two entries are found that pass the above similarity check, take the one with the lower
                        # energy:
                        if entry["final_energy"] < filtered_entry["final_energy"]:
                            self.filtered_entries[ii] = entry
                        # Note that this will essentially choose between singlet and triplet entries assuming both have
                        # the same structural details
                        break
            if not found_similar_entry:
                self.filtered_entries += [entry]

    def build_new_entry(self, frags, bonds):
        """
        Simple function to format a bond dissociation entry that will eventually be returned to the user.

        :param frags:
        :param bonds:
        :return:
        """
        specie = nx.get_node_attributes(self.mol_graph.graph, "specie")
        if len(frags) == 2:
            new_entry = [
                self.molecule_entry["final_energy"] - (frags[0]["final_energy"] + frags[1]["final_energy"]),
                bonds,
                specie[bonds[0][0]],
                specie[bonds[0][1]],
                frags[0]["smiles"],
                frags[0]["structure_change"],
                frags[0]["initial_molecule"]["charge"],
                frags[0]["initial_molecule"]["spin_multiplicity"],
                frags[0]["final_energy"],
                frags[1]["smiles"],
                frags[1]["structure_change"],
                frags[1]["initial_molecule"]["charge"],
                frags[1]["initial_molecule"]["spin_multiplicity"],
                frags[1]["final_energy"],
            ]
        else:
            new_entry = [
                self.molecule_entry["final_energy"] - frags[0]["final_energy"],
                bonds,
                specie[bonds[0][0]],
                specie[bonds[0][1]],
                frags[0]["smiles"],
                frags[0]["structure_change"],
                frags[0]["initial_molecule"]["charge"],
                frags[0]["initial_molecule"]["spin_multiplicity"],
                frags[0]["final_energy"],
            ]
        return new_entry
