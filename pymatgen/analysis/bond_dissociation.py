# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import re
import csv
import collections
import itertools
from io import open
import math
from six.moves import zip
import logging

from monty.json import MSONable, MontyDecoder

import numpy as np

from pymatgen.core.structure import Molecule
from pymatgen.io.qchem.outputs import edges_from_babel, build_MoleculeGraph, is_isomorphic
import networkx as nx

# have_babel = True
# try:
#     from pymatgen.io.babel import BabelMolAdaptor
#     import openbabel as ob
# except ImportError:
#     print("Cannot find OpenBabel! Thus, bonds must be provided by the user.")
#     have_babel = False


"""
This module defines tools to calculate bond dissocation energies.
"""

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
        self.fragment_entries = fragment_entries
        self.bond_dissociation_energies = []
        # self.all_bonds = edges_from_babel(target_molecule)
        self.ring_bonds = []
        self.bad_pairs = []
        self.mol_graph = build_MoleculeGraph(Molecule.from_dict(molecule_entry["output"]["optimized_molecule"]))#, self.all_bonds)
        for bond in self.mol_graph.graph.edges: #switch to all_bonds?
            bonds = [(bond[0],bond[1])]
            self.fragment_and_process(bonds)
        self.bond_pairs = []
        for ii,bond in enumerate(self.ring_bonds):
            for jj in range(ii+1,len(self.ring_bonds)):
                bond_pair = [bond, self.ring_bonds[jj]]
                self.bond_pairs += [bond_pair]
        for bond_pair in self.bond_pairs:
            self.fragment_and_process(bond_pair)
        self.solve_ring_bonds()
        return self.bond_dissociation_energies

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
                frag1_entries = self.search_fragment_entries(frags[0])[0]
                frag2_entries = self.search_fragment_entries(frags[1])[0]
                for frag1 in frag1_entries:
                    for frag2 in frag2_entries:
                        if frag1["output"]["optimized_molecule"].charge + frag2["output"]["optimized_molecule"].charge == self.molecule_entry["output"]["optimized_molecule"].charge:
                            new_entry = [bonds, self.molecule_entry["output"]["final_energy"] - (frag1["output"]["final_energy"] + frag2["output"]["final_energy"]), frag1["output"]["final_energy"], frag1["output"]["optimized_molecule"].charge, frag2["output"]["final_energy"], frag2["output"]["optimized_molecule"].charge]
                            self.bond_dissociation_energies += new_entry

        def search_fragment_entries(self, frag):
            entries = []
            initial_entries = []
            final_entries = []
            for entry in self.fragment_entries:
                initial_graph = build_MoleculeGraph(Molecule.from_dict(entry["input"]["initial_molecule"])).graph
                final_graph = build_MoleculeGraph(Molecule.from_dict(entry["output"]["optimized_molecule"])).graph
                if is_isomorphic(frag.graph, initial_graph) and is_isomorphic(frag.graph, final_graph):
                    entries += entry
                elif is_isomorphic(frag.graph, initial_graph):
                    initial_entries += entry
                elif is_isomorphic(frag.graph, final_graph):
                    final_entries += entry
            return [entries, initial_entries, final_entries]

        def solve_ring_bonds(self):
            pass
