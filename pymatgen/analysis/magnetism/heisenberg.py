# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from pymatgen.analysis.magnetism import CollinearMagneticStructureAnalyzer
from pymatgen.analysis.graphs import StructureGraph
from pymatgen.analysis.local_env import MinimumDistanceNN
from pymatgen.analysis.structure_matcher import ElementComparator, StructureMatcher
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

import numpy as np
import pandas as pd


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
                 strategy=MinimumDistanceNN(), cutoff=10):
        """Compute exchange parameters from low energy magnetic orderings.

        Exchange parameters are computed by mapping to a classical Heisenberg
        model. Strategy is the scheme for generating neighbors. Currently only
        MinimumDistanceNN is implemented.
        n+1 unique orderings are required to compute n exchange
        parameters. The default cutoff is to consider 10 structures.

        Args:
            ordered_structures (list): Structure objects with magmoms.
            energies (list): Energies of each relaxed magnetic structure.
            strategy (object): Class from pymatgen.analysis.local_env
            cutoff (int): Number of structures to consider.

        Parameters:
            sgraphs (list): StructureGraph objects.
            unique_site_ids (dict): Maps each site to its unique identifier
            nn_interacations (dict): {i: j} pairs of NN interactions
                between unique sites.
            ex_mat (DataFrame): Invertible Heisenberg Hamiltonian for each
                graph.
            ex_params (dict): Exchange parameter values (meV/atom) 

        """

        # Sort by energy if not already sorted
        ordered_structures = [s for _, s in
                              sorted(zip(energies, ordered_structures), reverse=True)]

        energies = sorted(energies, reverse=True)

        # Cutoff for structure graphs
        ordered_structures = ordered_structures[:cutoff]
        energies = energies[:cutoff]

        # Get only magnetic ions
        ordered_structures = [CollinearMagneticStructureAnalyzer(s).get_structure_with_only_magnetic_atoms()
                              for s in ordered_structures]

        self.ordered_structures = ordered_structures
        self.energies = energies

        # These attributes are set by internal methods
        self.sgraphs = None
        self.unique_site_ids = None
        self.nn_interactions = None
        self.ex_mat = None
        self.ex_params = None

        # Set attributes
        self._get_graphs()
        self._get_unique_sites()
        self._get_nn_dict()
        self._get_exchange_df()

    def _get_graphs(self):
        """
        Generate graph representations of magnetic structures with nearest
        neighbor bonds. Right now this only works for MinimumDistanceNN.

        Returns:
            None (sets self.sgraphs instance variable)

        Todo:
            * When supersizing everything, make site labels based on one
            convential structure and enforce that all other structs follow
            that labeling scheme.
            * Implement other strategies that capture NNN bonds, etc.
        """

        strategy = self.strategy

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
                # scale the energy by supercell
                self.energies[s_index] *= scale

        self.ordered_structures = matched_structures

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
        s = sgraphs[0].structures
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
        ex_mat_empty = pd.DataFrame(columns=columns)

        ex_mat = ex_mat_empty.copy()
        sgraphs_copy = sgraphs[:]  # deep copy of sgraphs
        sgraph_index = 0

        # Loop over all sites in each graph and compute |S_i . S_j|
        # for n+1 unique graphs to compute n exchange params
        while sgraph_index < len(num_nn_j):
            sgraph = sgraphs_copy.pop(0)
            ex_row = pd.DataFrame(np.zeros((1, num_nn_j)),
                                  index=[sgraph_index], columns=columns)

            for i in range(len(sgraph)):
                try:
                    s_i = sgraph.as_dict()['structure']['sites'][i][
                        'species'][0]['properties']['spin']
                except:  # no spin/magmom property
                    s_i = 0
                for k in unique_site_ids.keys():
                    if i in k:
                        i_index = unique_site_ids[k]
                        break

                # Get all connections for ith site and compute |S_i . S_j|
                connections = sgraph.get_connected_sites(i)
                for j in range(len(connections)):
                    try:
                        s_j = connections[j][0].specie.spin
                    except:  # no spin/magmom
                        s_j = 0
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
                if ex_mat.append(ex_row).equals(ex_mat.append(ex_row).drop_duplicates(keep='first')):
                    e_index = self.ordered_structures.index(sgraph.structure)
                    ex_row.at[sgraph_index, 'E'] = self.energies[e_index]
                    sgraph_index += 1
                    ex_mat.append(ex_row)

        j_columns = [name for name in columns if name not in ['E', 'E0']]
        ex_mat[[j_columns]].div(2)  # 1/2 factor in Heisenberg Hamiltonian

        ex_mat[['E0']] = 1  # Nonmagnetic contribution

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
        H = ex_mat.loc[:, ex_mat.columns != 'E'].values
        H_inv = np.linalg.inv(H)
        j_ij = np.dot(H_inv, E)
        j_ij *= 1000  # meV
        j_ij /= len(self.ordered_structures[0])  # / atom
        ex_params = {j_name: j for j_name, j in zip(j_names, j_ij)}

        self.ex_params = ex_params

        return ex_params

    def estimate_exchange(self, fm_struct, afm_struct, fm_e, afm_e):
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

        n = len(fm_struct)
        magmoms = fm_struct.site_properties['magmom']
        m_avg = np.mean([np.sqrt(m**2) for m in magmoms])
        delta_e = afm_e - fm_e  # J > 0 -> FM
        j_avg = delta_e / (n*m_avg**2)
        j_avg *= 1000  # meV

        return j_avg

    def get_mft_temperature(self, j_avg):
        """Crude mean field estimate of critical temperature based on <J>.

        Args:
            j_avg (float): j_avg (float): Average exchange parameter (meV/atom)

        Returns:
            mft_t (float): Critical temperature (K)

        """

        k_boltzmann = 0.0861733  # meV/K
        mft_t = 2 * j_avg / 3 / k_boltzmann

        return mft_t  


