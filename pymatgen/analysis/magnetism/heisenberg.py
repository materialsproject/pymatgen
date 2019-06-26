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
import sys

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

    def __init__(self, ordered_structures, energies, cutoff=0.0):
        """Compute exchange parameters from low energy magnetic orderings.

        Exchange parameters are computed by mapping to a classical Heisenberg
        model. Strategy is the scheme for generating neighbors. Currently only
        MinimumDistanceNN is implemented.
        n+1 unique orderings are required to compute n exchange
        parameters.

        First run a MagneticOrderingsWF to obtain low energy collinear magnetic
        orderings and find the magnetic ground state. Then enumerate magnetic
        states with the ground state as the input structure, find the subset
        of supercells that map to the ground state, and do static calculations
        for these orderings.  

        Args:
            ordered_structures (list): Structure objects with magmoms.
            energies (list): Energies of each relaxed magnetic structure.
            cutoff (float): Cutoff in Angstrom for nearest neighbor search.
                Defaults to 0 (only NN, no NNN, etc.)

        Parameters:
            strategy (object): Class from pymatgen.analysis.local_env for
                constructing graphs.
            sgraphs (list): StructureGraph objects.
            unique_site_ids (dict): Maps each site to its unique numerical
                identifier.
            wyckoff_ids (dict): Maps unique numerical identifier to wyckoff
                position.
            nn_interacations (dict): {i: j} pairs of NN interactions
                between unique sites.
            ex_mat (DataFrame): Invertible Heisenberg Hamiltonian for each
                graph.
            ex_params (dict): Exchange parameter values (meV/atom) 

        """

        # Save original copies of inputs
        self.ordered_structures_ = ordered_structures
        self.energies_ = energies

        # Get energies / supercell if needed
        # energies = [e * len(s) for (e, s) in zip(energies, ordered_structures)]

        # Get only magnetic ions & give all structures site_properties['magmom']
        # zero threshold so that magnetic ions with small moments
        # are preserved
        ordered_structures = [CollinearMagneticStructureAnalyzer(s, make_primitive=False, 
            threshold=0.0).get_structure_with_only_magnetic_atoms(make_primitive=False)
                              for s in ordered_structures]

        # Check for duplicate / degenerate states (sometimes different initial
        # configs relax to the same state)
        remove_list = []
        for i, e in enumerate(energies):
            e = round(e, 6)  # 10^-6 eV tol on equal energies
            if i not in remove_list:
                for i_check, e_check in enumerate(energies):
                    e_check = round(e_check, 6)
                    if i != i_check and i_check not in remove_list and e == e_check:
                        remove_list.append(i_check)

        if len(remove_list):
            ordered_structures = [s for i, s in enumerate(ordered_structures) 
            if i not in remove_list]
            energies = [e for i, e in enumerate(energies) if i not in remove_list]

        # Sort by energy if not already sorted
        ordered_structures = [s for _, s in
                              sorted(zip(energies, ordered_structures), reverse=False)]
        energies = sorted(energies, reverse=False)

        self.ordered_structures = ordered_structures
        self.energies = energies
        self.cutoff = cutoff

        # These attributes are set by internal methods
        self.sgraphs = None
        self.wyckoff_ids = None
        self.unique_site_ids = None
        self.nn_interactions = None
        self.ex_mat = None
        self.ex_params = None

        # Get number of unique sites
        self._get_unique_sites()
        
        # Get graph representations
        self._get_graphs(self.cutoff)

        # Check how many commensurate graphs we found
        if len(self.sgraphs) < 2:
            print('We need at least 2 unique orderings.')
            sys.exit(1)
        else:  # Set attributes
            self._get_nn_dict()
            self._get_exchange_df()


    def _get_graphs(self, cutoff):
        """
        Generate graph representations of magnetic structures with nearest
        neighbor bonds. Right now this only works for MinimumDistanceNN.

        Args:
            cutoff (float): Cutoff in Angstrom for nearest neighbor search.

        Returns:
            None (sets self.sgraphs instance variable)

        Todo:
            * Implement other strategies that capture NNN bonds, etc.
        """

        # Strategy for finding neighbors
        if cutoff:
            strategy = MinimumDistanceNN(cutoff=cutoff, get_all_sites=True)
        else:
            strategy = MinimumDistanceNN()  # only NN

        # Generate structure graphs
        sgraphs = [StructureGraph.with_local_env_strategy(s, strategy=strategy)
                   for s in self.ordered_structures]

        self.sgraphs = sgraphs

    def _get_unique_sites(self):
        """
        Get dict that maps site indices to unique identifiers.

        Returns:
            None: (sets self.unique_site_ids and self.wyckoff_ids 
            instance variables)

        """

        # Get a nonmagnetic representation of the supercell geometry
        s0 = CollinearMagneticStructureAnalyzer(self.ordered_structures[0], 
            make_primitive=False,
            threshold=0.0).get_nonmagnetic_structure(make_primitive=False)

        # Get unique sites and wyckoff positions
        if 'wyckoff' in s0.site_properties:
            s0.remove_site_property('wyckoff')

        symm_s0 = SpacegroupAnalyzer(s0).get_symmetrized_structure()
        wyckoff = ['n/a'] * len(symm_s0)
        equivalent_indices = symm_s0.equivalent_indices
        wyckoff_symbols = symm_s0.wyckoff_symbols

        # Construct dictionaries that map sites to numerical and wyckoff
        # identifiers
        unique_site_ids = {}
        wyckoff_ids = {}
        
        i = 0
        for indices, symbol in zip(equivalent_indices, wyckoff_symbols):
            unique_site_ids[tuple(indices)] = i
            wyckoff_ids[i] = symbol
            i += 1
            for index in indices:
                wyckoff[index] = symbol

        self.unique_site_ids = unique_site_ids
        self.wyckoff_ids = wyckoff_ids

    def _get_nn_dict(self):
        """Get dict of unique nearest neighbor interactions.

        Returns:
            None: (sets self.nn_interactions instance variable)

        """

        sgraph = self.sgraphs[0]
        unique_site_ids = self.unique_site_ids

        nn_dict = {}
        nnn_dict = {}
        nnnn_dict = {}

        # Loop over unique sites and get neighbor distances
        for k in unique_site_ids:
            i = k[0]
            i_key = unique_site_ids[k]
            connected_sites = sgraph.get_connected_sites(i)
            dists = [round(cs[-1], 2) for cs in connected_sites]  # i<->j distances
            dists = sorted(list(set(dists)))  # NN, NNN, NNNN, etc.

            # Loop over sites and determine unique NN, NNN, etc. interactions 
            for cs in connected_sites:
                dist = round(cs[-1], 2)  # i_j distance
                j = cs[2]  # j index
                for key in unique_site_ids.keys():
                    if j in key:
                        j_key = unique_site_ids[key]
                if dist == dists[0]:
                    nn_dict[i_key] = j_key
                elif dist == dists[1]:
                    nnn_dict[i_key] = j_key
                elif dist == dists[2]:
                    nnnn_dict[i_key] = j_key

        nn_interactions = {'nn': nn_dict, 'nnn': nnn_dict, 'nnnn': nnnn_dict}

        self.nn_interactions = nn_interactions

    def _get_exchange_df(self):
        """
        Loop over all sites in a graph and count the number and types of
        nearest neighbor interactions, computing +-|S_i . S_j| to construct
        a Heisenberg Hamiltonian for each graph.

        Returns:
            None: (sets self.ex_mat instance variable)

        TODO:
            * Account for beyond NN if possible

        """

        sgraphs = self.sgraphs
        unique_site_ids = self.unique_site_ids
        nn_interactions = self.nn_interactions

        # Total energy and nonmagnetic energy contribution
        columns = ['E', 'E0']

        # Get labels of unique NN interactions
        for k0, v0 in nn_interactions.items():
            for i, j in v0.items():  # i and j indices
                c = str(i) + '-' + str(j) + '-' + str(k0)
                c_rev = str(j) + '-' + str(i) + '-' + str(k0)
                if c not in columns and c_rev not in columns:
                    columns.append(c)

        num_sgraphs = len(sgraphs)

        # Keep n interactions (not counting 'E') for n+1 structure graphs
        columns = columns[:num_sgraphs]

        num_nn_j = len(columns) - 1  # ignore total energy
        j_columns = [name for name in columns if name not in ['E', 'E0']]
        ex_mat_empty = pd.DataFrame(columns=columns)
        ex_mat = ex_mat_empty.copy()

        if len(j_columns) < 2:
            self.ex_mat = ex_mat  # Only <J> can be calculated here
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
                    s_i = sgraph.structure.site_properties['magmom'][i]

                    for k in unique_site_ids.keys():
                        if i in k:
                            i_index = unique_site_ids[k]

                    # Get all connections for ith site and compute |S_i . S_j|
                    connections = sgraph.get_connected_sites(i)
                    dists = [round(cs[-1], 2) for cs in connections]  # i<->j distances
                    dists = sorted(list(set(dists)))  # NN, NNN, NNNN, etc.

                    for j, connection in enumerate(connections):
                        j_site = connection[2]
                        dist = round(connection[-1], 2)  # i_j distance

                        s_j = sgraph.structure.site_properties['magmom'][j_site]

                        for k in unique_site_ids.keys():
                            if j in k:
                                j_index = unique_site_ids[k]

                        # Determine order of connection
                        if dist == dists[0]:
                            order = '-nn'
                        elif dist == dists[1]:
                            order = '-nnn'
                        elif dist == dists[2]:
                            order = '-nnn'
                        j_ij = str(i_index) + '-' + str(j_index) + order
                        j_ji = str(j_index) + '-' + str(i_index) + order

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
                    # if sgraph_index == num_nn_j:  # check for zero columns
                    #     zeros = [b for b in (ex_mat[j_columns] == 0).all(axis=0)]
                    #     if True in zeros:
                    #         sgraph_index -= 1  # keep looking
            
            
            ex_mat[j_columns] = ex_mat[j_columns].div(2.)  # 1/2 factor in Heisenberg Hamiltonian
            ex_mat[['E0']] = 1  # Nonmagnetic contribution

            # Check for singularities and delete columns with all zeros
            zeros = [b for b in (ex_mat == 0).all(axis=0)]
            if True in zeros:
                c = ex_mat.columns[zeros.index(True)]
                len_zeros = len(c)
                ex_mat = ex_mat.drop(columns=[c], axis=1)
                ex_mat = ex_mat.drop(ex_mat.tail(len_zeros).index)

            # Check if ex_mat is square
            ex_mat = ex_mat[:ex_mat.shape[1]-1]                 

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

            # Estimate exchange by J ~ E_AFM - E_FM
            j_avg = self.estimate_exchange()
            ex_params = {'<J>': j_avg}
            self.ex_params = ex_params

            return ex_params

        # Solve eigenvalue problem for more than 1 NN interaction
        H = ex_mat.loc[:, ex_mat.columns != 'E'].values
        np.savetxt('H.json', H)
        H_inv = np.linalg.inv(H)
        j_ij = np.dot(H_inv, E)
        j_ij[1:] *= 1000  # meV
        #j_ij[1:] /= len(self.ordered_structures[0])  # / atom
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

        TODO:
            * Should we just take the energy diff between two lowest energy
            orderings?

        """

        mag_min = np.inf
        mag_max = 0.01
        fm_e_min = 0
        afm_e_min = 0
        afm_threshold = 0.1  # total magnetization < threshold -> AFM, not FiM
            
        for s, e in zip(self.ordered_structures, self.energies):
            magmoms = s.site_properties['magmom']
            if abs(sum(magmoms)) > mag_max and e < fm_e_min:  # FM ground state
                fm_struct = s
                fm_e = e
                mag_max = abs(sum(magmoms))
                fm_e_min = e

            # AFM ground state
            if abs(sum(magmoms)) < mag_min or abs(sum(magmoms)) < afm_threshold and e < afm_e_min:
                afm_struct = s
                afm_e = e
                mag_min = abs(sum(magmoms))
                afm_e_min = e

        # Get two lowest energy configs and determine FM vs AFM
        # s0, s1 = self.ordered_structures[0], self.ordered_structures[1]
        # e0, e1 = self.energies[0], self.energies[1]
        # tot_mag_s0 = abs(sum(s0.site_properties['magmom']))
        # tot_mag_s1 = abs(sum(s1.site_properties['magmom']))

        # if tot_mag_s0 > tot_mag_s1:
        #     fm_struct = s0
        #     fm_e = e0
        #     afm_struct = s1
        #     afm_e = e1
        # else:
        #     fm_struct = s1
        #     fm_e = e1
        #     afm_struct = s0
        #     afm_e = e0

        # Convert to magnetic structures with 'magmom' site property
        fm_struct = CollinearMagneticStructureAnalyzer(fm_struct,
            make_primitive=False, 
            threshold=0.0).get_structure_with_only_magnetic_atoms(make_primitive=False)
        afm_struct = CollinearMagneticStructureAnalyzer(afm_struct,
            make_primitive=False, 
            threshold=0.0).get_structure_with_only_magnetic_atoms(make_primitive=False)

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

        n = max(len(fm_struct), len(afm_struct))  # num of magnetic atoms
        # fm_e *= n  # eV/ion -> eV
        # afm_e *= n

        magmoms = fm_struct.site_properties['magmom']
        afm_magmoms = afm_struct.site_properties['magmom']

        m_avg = np.mean([np.sqrt(m**2) for m in magmoms])
        afm_m_avg = np.mean([np.sqrt(m**2) for m in afm_magmoms])

        delta_e = afm_e - fm_e  # J > 0 -> FM
        #j_avg = delta_e / (n*m_avg**2)
        j_avg = delta_e / (n * m_avg**2)  # eV / magnetic ion
        j_avg *= 1000 # meV / ion

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

    def get_interaction_graph(self, filename=None):
        """
        Get a StructureGraph with edges and weights that correspond to exchange
        interactions and J_ij values, respectively.
        
        Args:
            filename (str): if not None, save interaction graph to filename.
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

        # Save to a json file if desired
        if filename:
            if filename.endswith('.json'):
                dumpfn(igraph, filename)
            else:
                filename += '.json'
                dumpfn(igraph, filename)

        return igraph


