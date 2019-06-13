# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from pymatgen.analysis.graphs import StructureGraph
from pymatgen.analysis.local_env import MinimumDistanceNN
from pymatgen.analysis.magnetism import MagneticStructureEnumerator, CollinearMagneticStructureAnalyzer
from pymatgen.analysis.structure_matcher import StructureMatcher, ElementComparator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from pymatgen import Structure

import pandas as pd
import numpy as np

"""
This module implements a simple algorithm for extracting nearest neighbor
exchange parameters by mapping low energy magnetic orderings to a Heisenberg
model.
"""

__author__ = "ncfrey"
__version__ = "0.1"
__maintainer__ = "Nathan C. Frey"
__email__ = "ncfrey@lbl.gov"
__status__ = "Pre-Beta"
__date__ = "6/12/19"

class HeisenbergMapper:
   
    def __init__(self, ordered_structures, energies):        
        """Compute exchange parameters from low energy magnetic orderings.
        
        Args:
        	ordered_structures (list): Structure objects with magmoms.
        	energies (list): Energies of each relaxed magnetic structure.

        """
        
        self.ordered_structures = ordered_structures
        self.energies = energies
    
    def fit_graphs(self, strategy=MinimumDistanceNN()):
        """
        Generate graph representations of magnetic structures with nearest 
        neighbor bonds. Right now this only works for MinimumDistanceNN.

        Args:
        	strategy (object): Class from pymatgen.analysis.local_env

        Returns:
        	sgraphs (list): StructureGraph objects.

        
        Todo:
            * When supersizing everything, make site labels based on one
            convential structure and enforce that all other structs follow
            that labeling scheme.
            * Implement other strategies that capture NNN bonds, etc.
        """

        # Find a maximally sized supercell to reference site labels to
        imax = np.argmax([len(s) for s in self.ordered_structures])
        s1 = self.ordered_structs[imax]

        # Match all structures to reference so site labels are consistent
        sm = StructureMatcher(primitive_cell=False, attempt_supercell=True, 
        	comparator=ElementComparator())
        matched_structures = self.ordered_structures[:]

        for i, s in enumerate(self.ordered_structures, 0):
        	if len(s) < len(s1):
        		s2 = sm.get_s2_like_s1(s1, s)
        		matched_structures[i] = s2
        		s_index = self.ordered_structures.index(s)
        		scale = len(s2) / len(s)
        		self.energies[s_index] *= scale  # scale the energy by supercell

        self.ordered_structures = matched_structures

        # Generate structure graphs
        sgraphs = [StructureGraph.with_local_env_strategy(s, strategy=strategy)
        for s in self.ordered_structures]

        return sgraphs


    def _get_unique_sites(self, sgraphs):
        """
        Get dict that maps site indices to unique identifiers. If there is only
        one unique site, get the exchange <J> by averaging.
        
        Args:
        	sgraphs (list): StructureGraph objects from fit_graphs() method.

        Returns:
        	unique_site_ids (dict): Maps each site to its unique identifier.
        """

        # Find the unique sites in a structure graph
        s = sgraphs[0].structure
        symm_struct = SpacegroupAnalyzer(s).get_symmetrized_structure()
        unique_sites =[]
        unique_site_ids = {}
        i = 0

        for site in s:
        	if site not in (elem for sublist in unique_sites for elem in sublist):
        		equiv_sites = symm_struct.find_equivalent_sites(site)
        		unique_sites.append(equiv_sites)
        		equiv_sites_ids = [s.index(site) for site in equiv_sites]
        		unique_site_ids[tuple(equiv_site_ids)] = i
        		i += 1

        return unique_site_ids


    def _get_nn_dict(self, sgraphs, unique_site_ids):
        """Get dict of unique nearest neighbor interactions.

        Args:
        	sgraphs (list): StructureGraph objects.
        	unique_site_ids (dict): Site to unique id map.
        
        Returns:
        	nn_interacations (dict): {i: j} pairs of NN interactions between 
        	unique sites.
        """

        nn_ids = {}
        sgraph = sgraphs[0]

        # Loop over unique sites and get index of NN
        for k in unique_site_ids:
        	i = k[0]
        	i_key = unique_site_ids[k]
        	nn = sgraph.get_connected_sites(i)[0]
        	j = nn[2]  # index of NN
        	for key in unique_site_ids.keys():
        		if j in key:
        			j_key = unique_site_ids[key]
        	nn_interactions[i_key] = j_key

        return nn_interactions

 
    def _get_exchange_df(self, sgraphs, energies, unique_site_ids, 
    	nn_interactions):
        """
        Loop over all sites in a graph and count the number and types of
        nearest neighbor interactions, computing +-|S_i . S_j| to construct
        a Heisenberg Hamiltonian for each graph.
        
        Args:
        	sgraphs (list): StructureGraph objects.
        	energies (list): Energies of StructureGraph.structure
        	unique_site_ids (dict): Maps each site to its unique identifier.
        	nn_interactions (dict): Pairs of NN interactions between 
        	unique sites.

        Returns:
        	ex_mat (DataFrame): Invertible Heisenberg Hamiltonian for each graph.
        """

        # Total energy and nonmagnetic energy contribution
        columns = ['E', 'E0']   

        # Get labels of unique NN interactions
        for k0, v0 in nn_interactions.items():
        	c = str(k0) + '-' + str(v0)
        	c_rev = str(v0) + '-' + str(k0)
        	if c not in columns and c_rev not in columns:
        		columns.append(c)

        num_nn_j = len(columns) - 1  # ignore total energy
        ex_mat_empty = pd.DataFrame(columns=columns)

        ex_mat = ex_mat_empty.copy()
        sgraphs_copy = sgraphs[:]  # deep copy of sgraphs
        sgraph_index = 0

        # Loop over all sites in each graph and compute |S_i . S_j|
        # for n+1 unique graphs to compute n exchange params
        while sgraph_index < len(num_nn_j):
        	sgraph = sgraphs_copy.pop(0)
        	ex_row = pd.DataFrame(np.zeros((1, num_nn_j)), index=[sgraph_index],
        	 columns=columns)

        	for i in range(len(sgraph)):
        		try:
        			S_i = sgraph.as_dict()['structure']['sites'][i]['species'][0]['properties']['spin']
        		except:  # no spin/magmom property
        			S_i = 0
        		for k in unique_site_ids.keys():
        			if i in k:
        				i_index = unique_site_ids[k]
        				break

        		# Get all connections for ith site and compute |S_i . S_j|
        		connections = sgraph.get_connected_sites(i)
        		for j in range(len(connections)):
        			try:
        				S_j = connections[j][0].specie.spin
        			except:  # no spin/magmom
        				S_j = 0
        			for k in unique_site_ids.keys():
        				if j in k:
        					j_index = unique_site_ids[k]
        					break

        			J_ij = str(i_index) + '-' + str(j_index)
        			J_ji = str(j_index) + '-' + str(i_index)
        			if J_ij in ex_mat.columns:
        				ex_row.at[sgraph_index, J_ij] -= S_i * S_j
        			elif J_ji in ex_mat.columns:
        				ex_row.at[sgraph_index, J_ji] -= S_i * S_j

        		# Ignore the row if it is a duplicate to avoid singular matrix
        		if ex_mat.append(ex_row).equals(ex_mat.append(ex_row).drop_duplicates(keep='first')):
        			e_index = self.ordered_structures.index(sgraph.structure)
        			ex_row.at[sgraph_index, 'E'] = self.energies[e_index]
        			sgraph_index += 1
        			ex_mat.append(ex_row)

        j_columns = columns.remove
        ex_mat[[]]






    
    def get_exchange(self, exchange_df):
        
        """
        Take Heisenberg Hamiltonian and corresponding energy for each row and
        solve for the exchange parameters.
        
        :return (dict): {J_ij: value (meV/atom)}.
        """
