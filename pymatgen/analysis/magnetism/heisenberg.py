# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from pymatgen.analysis.magnetism import CollinearMagneticStructureAnalyzer
from pymatgen.analysis.graphs import StructureGraph
from pymatgen.analysis.local_env import MinimumDistanceNN
from pymatgen.analysis.structure_matcher import ElementComparator, StructureMatcher
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from monty.serialization import dumpfn

import numpy as np
import pandas as pd
import copy

import warnings

"""
This module implements a simple algorithm for extracting nearest neighbor
exchange parameters by mapping low energy magnetic orderings to a Heisenberg
model.
"""

__author__ = "ncfrey"
__version__ = "0.1"
__maintainer__ = "Nathan C. Frey"
__email__ = "ncfrey@lbl.gov"
__status__ = "Development"
__date__ = "June 2019"


class HeisenbergMapper:

    def __init__(self, ordered_structures, energies,
                 strategy=MinimumDistanceNN()):
        """Compute exchange parameters from low energy magnetic orderings.

        Exchange parameters are computed by mapping to a classical Heisenberg
        model. Strategy is the scheme for generating neighbors. Currently only
        MinimumDistanceNN is implemented.
        n+1 unique orderings are required to compute n exchange
        parameters.

        Args:
            ordered_structures (list): Structure objects with magmoms.
            energies (list): Energies of each relaxed magnetic structure.
            strategy (object): Class from pymatgen.analysis.local_env

        Parameters:
            sgraphs (list): StructureGraph objects.
            unique_site_ids (dict): Maps each site to its unique identifier
            nn_interacations (dict): {i: j} pairs of NN interactions
                between unique sites.
            ex_mat (DataFrame): Invertible Heisenberg Hamiltonian for each
                graph.
            ex_params (dict): Exchange parameter values (meV/atom) 

        """

        # Get only magnetic ions & give all structures site_properties['magmom']
        # threshold set to 0.01 uB so that magnetic ions with small moments
        # are preserved
        ordered_structures = [CollinearMagneticStructureAnalyzer(s, 
            threshold=0.01).get_structure_with_only_magnetic_atoms()
                              for s in ordered_structures]

        # Normalize supercells and get graph representations
        self.sgraphs = None
        self.ordered_structures = ordered_structures
        self.energies = energies
        self.strategy = strategy
        self._get_graphs()

        # Sort by energy if not already sorted
        self.ordered_structures = [s for _, s in
                              sorted(zip(self.energies, self.ordered_structures), reverse=False)]

        self.energies = sorted(self.energies, reverse=False)

        # These attributes are set by internal methods
        self.unique_site_ids = None
        self.nn_interactions = None
        self.ex_mat = None
        self.ex_params = None

        # Set attributes
        self._get_unique_sites()  # xx - change to static method
        self._get_nn_dict()
        self._get_exchange_df()


    def _get_graphs(self):
        """
        Generate graph representations of magnetic structures with nearest
        neighbor bonds. Right now this only works for MinimumDistanceNN.

        Returns:
            None (sets self.sgraphs instance variable)

        Todo:
            * Don't supersize / structure match; instead, find symmetrically
            equiv sites and assign unique labels
            * Structure matcher is throwing out structures with incommensurate
            supercells...need to fix.
            * Implement other strategies that capture NNN bonds, etc.
        """

        strategy = self.strategy

        # Check if supercells are present
        if any(len(s) != len(self.ordered_structures[0]) for s in 
            self.ordered_structures):

            # Find a maximally sized supercell to reference site labels to
            imax = np.argmax([s.volume for s in self.ordered_structures])
            s1 = self.ordered_structures[imax]

            # Match all structures to reference so site labels are consistent
            sm = StructureMatcher(primitive_cell=False, attempt_supercell=True,
                                  comparator=ElementComparator())
            matched_structures = []
            matched_energies = []

            for i, s in enumerate(self.ordered_structures, 0):
                s2 = sm.get_s2_like_s1(s1, s)
                if s2 is not None:
                    matched_structures.append(s2)
                    scale = len(s2) / len(s)
                    matched_energies.append(self.energies[i] * scale)

            self.ordered_structures = matched_structures
            self.energies = matched_energies

        # Generate structure graphs
        sgraphs = [StructureGraph.with_local_env_strategy(s, strategy=strategy)
                   for s in self.ordered_structures]

        self.sgraphs = sgraphs

    def _get_unique_sites(self):
        """
        Get dict that maps site indices to unique identifiers. If there is only
        one unique site, get the exchange <J> by averaging.

        Returns:
            None: (sets self.unique_site_ids instance variable)

        """

        sgraphs = self.sgraphs
        # Find the unique sites in a structure graph
        s = sgraphs[0].structure
        symm_struct = SpacegroupAnalyzer(s).get_symmetrized_structure()
        unique_sites = []
        unique_site_ids = {}
        i = 0

        for site in s:
            if site not in (elem for sublist in unique_sites for elem in sublist):
                equiv_sites = symm_struct.find_equivalent_sites(site)
                unique_sites.append(equiv_sites)
                equiv_site_ids = [s.index(site) for site in equiv_sites]
                unique_site_ids[tuple(equiv_site_ids)] = i
                i += 1

        self.unique_site_ids = unique_site_ids

    def _get_nn_dict(self):
        """Get dict of unique nearest neighbor interactions.

        Returns:
            None: (sets self.nn_interactions instance variable)

        Todo:
            * Add NNN interactions, etc.

        """

        sgraphs = self.sgraphs
        unique_site_ids = self.unique_site_ids

        nn_interactions = {}
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

        self.nn_interactions = nn_interactions

    def _get_exchange_df(self):
        """
        Loop over all sites in a graph and count the number and types of
        nearest neighbor interactions, computing +-|S_i . S_j| to construct
        a Heisenberg Hamiltonian for each graph.

        Returns:
            None: (sets self.ex_mat instance variable)

        """

        sgraphs = self.sgraphs
        unique_site_ids = self.unique_site_ids
        nn_interactions = self.nn_interactions

        # Total energy and nonmagnetic energy contribution
        columns = ['E', 'E0']

        # Get labels of unique NN interactions
        for k0, v0 in nn_interactions.items():
            c = str(k0) + '-' + str(v0)
            c_rev = str(v0) + '-' + str(k0)
            if c not in columns and c_rev not in columns:
                columns.append(c)

        num_nn_j = len(columns) - 1  # ignore total energy
        j_columns = [name for name in columns if name not in ['E', 'E0']]
        ex_mat_empty = pd.DataFrame(columns=columns)

        ex_mat = ex_mat_empty.copy()

        # Only 1 unique NN interaction
        if len(j_columns) < 2:
            self.ex_mat = ex_mat  # Just use the J_ij column labels
        else:
            sgraphs_copy = copy.deepcopy(sgraphs)
            sgraph_index = 0

            # Loop over all sites in each graph and compute |S_i . S_j|
            # for n+1 unique graphs to compute n exchange params
            for graph in sgraphs:
                sgraph = sgraphs_copy.pop(0)
                ex_row = pd.DataFrame(np.zeros((1, num_nn_j + 1)),
                                      index=[sgraph_index], columns=columns)

                for i, node in enumerate(sgraph.graph.nodes):
                    # try:
                    s_i = sgraph.structure.site_properties['magmom'][i]
                    # except:  # no spin/magmom property
                    #     s_i = 0
                    for k in unique_site_ids.keys():
                        if i in k:
                            i_index = unique_site_ids[k]
                            break

                    # Get all connections for ith site and compute |S_i . S_j|
                    connections = sgraph.get_connected_sites(i)
                    for j, connection in enumerate(connections):
                        j_site = connection[2]
                        # try:
                        s_j = sgraph.structure.site_properties['magmom'][j]
                        # except:  # no spin/magmom
                        #     s_j = 0
                        for k in unique_site_ids.keys():
                            if j in k:
                                j_index = unique_site_ids[k]
                                break

                        j_ij = str(i_index) + '-' + str(j_index)
                        j_ji = str(j_index) + '-' + str(i_index)
                        if j_ij in ex_mat.columns:
                            ex_row.at[sgraph_index, j_ij] -= s_i * s_j
                        elif j_ji in ex_mat.columns:
                            ex_row.at[sgraph_index, j_ji] -= s_i * s_j

                # Ignore the row if it is a duplicate to avoid singular matrix
                if ex_mat.append(ex_row)[j_columns].equals(ex_mat.append(ex_row)[j_columns].drop_duplicates(keep='first')):          
                    e_index = self.ordered_structures.index(sgraph.structure)
                    ex_row.at[sgraph_index, 'E'] = self.energies[e_index]
                    sgraph_index += 1
                    ex_mat = ex_mat.append(ex_row)
                    if sgraph_index == num_nn_j:  # check for zero columns
                        zeros = [b for b in (ex_mat[j_columns] == 0).all(axis=0)]
                        if True in zeros:
                            sgraph_index -= 1  # keep looking
        
            ex_mat[j_columns] = ex_mat[j_columns].div(2.)  # 1/2 factor in Heisenberg Hamiltonian
            ex_mat[['E0']] = 1  # Nonmagnetic contribution

            # Check for singularities and delete columns with all zeros
            zeros = [b for b in (ex_mat == 0).all(axis=0)]
            if True in zeros:
                c = ex_mat.columns[zeros.index(True)]
                len_zeros = len(c)
                ex_mat = ex_mat.drop(columns=[c], axis=1)
                ex_mat = ex_mat.drop(ex_mat.tail(len_zeros).index)

            self.ex_mat = ex_mat

    def get_exchange(self):
        """
        Take Heisenberg Hamiltonian and corresponding energy for each row and
        solve for the exchange parameters.

        Returns:
            ex_params (dict): Exchange parameter values (meV/atom).

        """

        ex_mat = self.ex_mat
        # Solve the matrix equation for J_ij values
        E = ex_mat[['E']]
        j_names = [j for j in ex_mat.columns if j not in ['E']]

        # Only 1 NN interaction
        if len(j_names) < 3:

            # Find lowest energy FM and AFM configs
            fm_struct, afm_struct, fm_e, afm_e = self.get_low_energy_orderings()

            # Estimate exchange by J ~ E_AFM - E_FM
            j_avg = self.estimate_exchange(fm_struct, afm_struct, fm_e, afm_e)
            j_name = [j for j in self.ex_mat.columns if j not in ['E', 'E0']][0]
            ex_params = {j_name: j_avg}
            self.ex_params = ex_params

            return ex_params

        # Solve eigenvalue problem for more than 1 NN interaction
        H = ex_mat.loc[:, ex_mat.columns != 'E'].values

        H_inv = np.linalg.inv(H)
        j_ij = np.dot(H_inv, E)
        j_ij[1:] *= 1000  # meV
        j_ij[1:] /= len(self.ordered_structures[0])  # / atom
        j_ij = j_ij.tolist()
        ex_params = {j_name: j for j_name, j in zip(j_names, j_ij)}

        self.ex_params = ex_params

        return ex_params

    def get_low_energy_orderings(self):
        """
        Find lowest energy FM and AFM orderings to compute E_AFM - E_FM.

        Returns:
            fm_struct (Structure): fm structure with 'magmom' site property
            afm_struct (Structure): afm structure with 'magmom' site property
            fm_e (float): fm energy
            afm_e (float): afm energy

        """

        mag_min = np.inf
        mag_max = 0.01
        fm_e_min = 0
        afm_e_min = 0
        afm_threshold = 0.1  # total magnetization < threshold -> AFM, not FiM
            
        for s, e in zip(self.ordered_structures, self.energies):
            magmoms = s.site_properties['magmom']
            if abs(sum(magmoms)) > mag_max and e < fm_e_min:
                fm_struct = s  # FM config
                fm_e = e
                mag_max = abs(sum(magmoms))
                fm_e_min = e
            if abs(sum(magmoms)) < mag_min or abs(sum(magmoms)) < afm_threshold and e < afm_e_min:
                afm_struct = s
                afm_e = e
                mag_min = abs(sum(magmoms))
                afm_e_min = e

        # Convert to magnetic structures with 'magmom' site property
        fm_struct = CollinearMagneticStructureAnalyzer(fm_struct, 
            threshold=0.01).get_structure_with_only_magnetic_atoms()
        afm_struct = CollinearMagneticStructureAnalyzer(afm_struct, 
            threshold=0.01).get_structure_with_only_magnetic_atoms()

        return fm_struct, afm_struct, fm_e, afm_e

    def estimate_exchange(self, fm_struct=None, afm_struct=None, fm_e=None, afm_e=None):
        """
        Estimate <J> for a structure based on low energy FM and AFM orderings.

        Args:
            fm_struct (Structure): fm structure with 'magmom' site property
            afm_struct (Structure): afm structure with 'magmom' site property
            fm_e (float): fm energy
            afm_e (float): afm energy
        
        Returns:
            j_avg (float): Average exchange parameter (meV/atom)

        """

        # Get low energy orderings if not supplied
        if any(arg is None for arg in [fm_struct, afm_struct, fm_e, afm_e]):
            fm_struct, afm_struct, fm_e, afm_e = self.get_low_energy_orderings()

        n = len(fm_struct)
        magmoms = fm_struct.site_properties['magmom']
        afm_magmoms = afm_struct.site_properties['magmom']
        m_avg = np.mean([np.sqrt(m**2) for m in magmoms])
        afm_m_avg = np.mean([np.sqrt(m**2) for m in magmoms])
        delta_e = afm_e - fm_e  # J > 0 -> FM
        j_avg = delta_e / (n*m_avg**2)
        j_avg *= 1000  # meV

        return j_avg

    def get_mft_temperature(self, j_avg):
        """
        Crude mean field estimate of critical temperature based on <J> for
        one sublattice, or solving the coupled equations for a multisublattice
        material.

        Args:
            j_avg (float): j_avg (float): Average exchange parameter (meV/atom)

        Returns:
            mft_t (float): Critical temperature (K)

        """

        num_sublattices = len(self.unique_site_ids)
        k_boltzmann = 0.0861733  # meV/K

        # Only 1 magnetic sublattice
        if num_sublattices == 1:
            
            mft_t = 2 * abs(j_avg) / 3 / k_boltzmann

            return mft_t

        else:  # multiple magnetic sublattices
            omega = np.zeros((num_sublattices, num_sublattices))
            for k in self.ex_params:
                # split into i, j unique site identifiers
                sites = [int(num) for num in k.split('-')]
                i, j = sites[0], sites[1]
                omega[i, j] = self.ex_params[k]
                omega[j, i] = self.ex_params[k]

            omega = omega * 2 / 3 / k_boltzmann
            eigenvals, eigenvecs = np.linalg.eig(omega)
            mft_t = max(eigenvals)

            return mft_t

    def get_interaction_graph(self):
        """
        Get a StructureGraph with edges and weights that correspond to exchange
        interactions and J_ij values, respectively.

        Returns:
            igraph (StructureGraph): Exchange interaction graph.
        """

        structure = self.ordered_structures[0]
        sgraph = self.sgraphs[0]
        unique_site_ids = self.unique_site_ids

        igraph = StructureGraph.with_empty_graph(structure, edge_weight_name=
            'exchange_constant', edge_weight_units='meV/atom')

        # J_ij exchange interaction matrix
        for i, node in enumerate(sgraph.graph.nodes):
            connections = sgraph.get_connected_sites(i)
            for c in connections:
                jimage = c[1]  # relative integer coordinates of atom j
                dx = jimage[0]
                dy = jimage[1]
                dz = jimage[2]
                j = c[2]  # index of neighbor

                # uniqe site identifiers
                for key in unique_site_ids.keys():
                    if i in key:
                        i_site = unique_site_ids[key]
                    if j in key:
                        j_site = unique_site_ids[key]

                j_ij = str(i_site) + '-' + str(j_site)
                j_ji = str(j_site) + '-' + str(i_site)

                if j_ij in self.ex_params:
                    j_value = self.ex_params[j_ij]
                elif j_ji in self.ex_params:
                    j_value = self.ex_params[j_ji]
                else:
                    j_value = 0

                igraph.add_edge(i, j, to_jimage=jimage, weight=j_value,
                    warn_duplicates=False)

        dumpfn(igraph, 'interaction_graph.json')

        return igraph


