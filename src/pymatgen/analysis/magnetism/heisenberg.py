"""
This module implements a simple algorithm for extracting nearest neighbor
exchange parameters by mapping low energy magnetic orderings to a Heisenberg
model.
"""

from __future__ import annotations

import copy
import logging
from ast import literal_eval
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from monty.json import MSONable, jsanitize
from monty.serialization import dumpfn

from pymatgen.analysis.graphs import StructureGraph
from pymatgen.analysis.local_env import MinimumDistanceNN
from pymatgen.analysis.magnetism import CollinearMagneticStructureAnalyzer, Ordering
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

if TYPE_CHECKING:
    from typing_extensions import Self

__author__ = "ncfrey"
__version__ = "0.1"
__maintainer__ = "Nathan C. Frey"
__email__ = "ncfrey@lbl.gov"
__status__ = "Development"
__date__ = "June 2019"

logger = logging.getLogger(__name__)


class HeisenbergMapper:
    """Compute exchange parameters from low energy magnetic orderings.

    Attributes:
        strategy (object): Class from pymatgen.analysis.local_env for constructing graphs.
        sgraphs (list): StructureGraph objects.
        unique_site_ids (dict): Maps each site to its unique numerical identifier.
        wyckoff_ids (dict): Maps unique numerical identifier to wyckoff position.
        nn_interactions (dict): {i: j} pairs of NN interactions between unique sites.
        dists (dict): NN, NNN, and NNNN interaction distances
        ex_mat (DataFrame): Invertible Heisenberg Hamiltonian for each graph.
        ex_params (dict): Exchange parameter values (meV/atom)
    """

    def __init__(self, ordered_structures, energies, cutoff=0, tol: float = 0.02):
        """Exchange parameters are computed by mapping to a classical Heisenberg
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
            energies (list): Total energies of each relaxed magnetic structure.
            cutoff (float): Cutoff in Angstrom for nearest neighbor search.
                Defaults to 0 (only NN, no NNN, etc.)
            tol (float): Tolerance (in Angstrom) on nearest neighbor distances
                being equal.
        """
        # Save original copies of inputs
        self.ordered_structures_ = ordered_structures
        self.energies_ = energies

        # Sanitize inputs and optionally order them by energy / magnetic moments
        hs = HeisenbergScreener(ordered_structures, energies, screen=False)
        ordered_structures = hs.screened_structures
        energies = hs.screened_energies

        self.ordered_structures = ordered_structures
        self.energies = energies
        self.cutoff = cutoff
        self.tol = tol

        # Get graph representations
        self.sgraphs = self._get_graphs(cutoff, ordered_structures)

        # Get unique site ids and wyckoff symbols
        self.unique_site_ids, self.wyckoff_ids = self._get_unique_sites(ordered_structures[0])

        # These attributes are set by internal methods
        self.nn_interactions = self.dists = self.ex_mat = self.ex_params = None

        # Check how many commensurate graphs we found
        if len(self.sgraphs) < 2:
            raise SystemExit("We need at least 2 unique orderings.")

        # Set attributes
        self._get_nn_dict()
        self._get_exchange_df()

    @staticmethod
    def _get_graphs(cutoff, ordered_structures):
        """Generate graph representations of magnetic structures with nearest
        neighbor bonds. Right now this only works for MinimumDistanceNN.

        Args:
            cutoff (float): Cutoff in Angstrom for nearest neighbor search.
            ordered_structures (list): Structure objects.

        Returns:
            sgraphs (list): StructureGraph objects.
        """
        # Strategy for finding neighbors
        strategy = MinimumDistanceNN(cutoff=cutoff, get_all_sites=True) if cutoff else MinimumDistanceNN()  # only NN

        # Generate structure graphs
        return [StructureGraph.from_local_env_strategy(s, strategy=strategy) for s in ordered_structures]

    @staticmethod
    def _get_unique_sites(structure):
        """Get dict that maps site indices to unique identifiers.

        Args:
            structure (Structure): ground state Structure object.

        Returns:
            tuple[dict, dict]: unique_site_ids maps tuples of equivalent site indices to a
                unique int identifier.
                wyckoff_ids maps tuples of equivalent site indices to their wyckoff symbols
        """
        # Get a nonmagnetic representation of the supercell geometry
        s0 = CollinearMagneticStructureAnalyzer(
            structure, make_primitive=False, threshold=0.0
        ).get_nonmagnetic_structure(make_primitive=False)

        # Get unique sites and wyckoff positions
        if "wyckoff" in s0.site_properties:
            s0.remove_site_property("wyckoff")

        symm_s0 = SpacegroupAnalyzer(s0).get_symmetrized_structure()
        wyckoff = ["n/a"] * len(symm_s0)
        equivalent_indices = symm_s0.equivalent_indices
        wyckoff_symbols = symm_s0.wyckoff_symbols

        # Construct dictionaries that map sites to numerical and wyckoff
        # identifiers
        unique_site_ids = {}
        wyckoff_ids = {}

        for idx, (indices, symbol) in enumerate(zip(equivalent_indices, wyckoff_symbols, strict=True)):
            unique_site_ids[tuple(indices)] = idx
            wyckoff_ids[idx] = symbol
            for index in indices:
                wyckoff[index] = symbol

        return unique_site_ids, wyckoff_ids

    def _get_nn_dict(self):
        """Set self.nn_interactions and self.dists instance variables describing unique
        nearest neighbor interactions.
        """
        tol = self.tol  # tolerance on NN distances
        sgraph = self.sgraphs[0]
        unique_site_ids = self.unique_site_ids

        nn_dict = {}
        nnn_dict = {}
        nnnn_dict = {}

        all_dists = []

        # Loop over unique sites and get neighbor distances up to NNNN
        for k in unique_site_ids:
            i = k[0]
            i_key = unique_site_ids[k]
            connected_sites = sgraph.get_connected_sites(i)
            dists = [round(cs[-1], 2) for cs in connected_sites]  # i<->j distances
            dists = sorted(set(dists))  # NN, NNN, NNNN, etc.

            dists = dists[:3]  # keep up to NNNN
            all_dists += dists

        # Keep only up to NNNN and call dists equal if they are within tol
        all_dists = sorted(set(all_dists))
        rm_list = []
        for idx, d in enumerate(all_dists[:-1], start=1):
            if abs(d - all_dists[idx]) < tol:
                rm_list.append(idx)

        all_dists = [d for idx, d in enumerate(all_dists) if idx not in rm_list]

        if len(all_dists) < 3:  # pad with zeros
            all_dists += [0] * (3 - len(all_dists))

        all_dists = all_dists[:3]
        labels = ("nn", "nnn", "nnnn")
        dists = dict(zip(labels, all_dists, strict=True))

        # Get dictionary keys for interactions
        for k in unique_site_ids:
            i = k[0]
            i_key = unique_site_ids[k]
            connected_sites = sgraph.get_connected_sites(i)

            # Loop over sites and determine unique NN, NNN, etc. interactions
            for cs in connected_sites:
                dist = round(cs[-1], 2)  # i_j distance

                j = cs[2]  # j index
                j_key = None
                for key, value in unique_site_ids.items():
                    if j in key:
                        j_key = value
                if abs(dist - dists["nn"]) <= tol:
                    nn_dict[i_key] = j_key
                elif abs(dist - dists["nnn"]) <= tol:
                    nnn_dict[i_key] = j_key
                elif abs(dist - dists["nnnn"]) <= tol:
                    nnnn_dict[i_key] = j_key

        nn_interactions = {"nn": nn_dict, "nnn": nnn_dict, "nnnn": nnnn_dict}

        self.dists = dists
        self.nn_interactions = nn_interactions

    def _get_exchange_df(self):
        """
        Loop over all sites in a graph and count the number and types of
        nearest neighbor interactions, computing +-|S_i . S_j| to construct
        a Heisenberg Hamiltonian for each graph. Sets self.ex_mat instance variable.

        TODO Deal with large variance in |S| across configs
        """
        sgraphs = self.sgraphs
        tol = self.tol
        unique_site_ids = self.unique_site_ids
        nn_interactions = self.nn_interactions
        dists = self.dists

        # Get |site magmoms| from FM ordering so that S_i and S_j are consistent?
        # Large S variations is throwing a loop
        # fm_struct = self.get_low_energy_orderings()[0]

        # Total energy and nonmagnetic energy contribution
        columns = ["E", "E0"]

        # Get labels of unique NN interactions
        for k0, v0 in nn_interactions.items():
            for idx, j in v0.items():  # i and j indices
                c = f"{idx}-{j}-{k0}"
                c_rev = f"{j}-{idx}-{k0}"
                if c not in columns and c_rev not in columns:
                    columns.append(c)

        n_sgraphs = len(sgraphs)

        # Keep n interactions (not counting 'E') for n+1 structure graphs
        columns = columns[: n_sgraphs + 1]

        n_nn_j = len(columns) - 1  # ignore total energy
        j_columns = [name for name in columns if name not in ["E", "E0"]]
        ex_mat_empty = pd.DataFrame(columns=columns)
        ex_mat = ex_mat_empty.copy()

        if len(j_columns) < 2:
            self.ex_mat = ex_mat  # Only <J> can be calculated here
        else:
            sgraphs_copy = copy.deepcopy(sgraphs)
            sgraph_index = 0

            # Loop over all sites in each graph and compute |S_i . S_j|
            # for n+1 unique graphs to compute n exchange params
            order = ""
            for _graph in sgraphs:
                sgraph = sgraphs_copy.pop(0)
                ex_row = pd.DataFrame(np.zeros((1, n_nn_j + 1)), index=[sgraph_index], columns=columns)

                for idx, _node in enumerate(sgraph.graph.nodes):
                    # s_i_sign = np.sign(sgraph.structure.site_properties['magmom'][i])
                    s_i = sgraph.structure.site_properties["magmom"][idx]

                    i_index = None
                    for k, v in unique_site_ids.items():
                        if idx in k:
                            i_index = v

                    # Get all connections for ith site and compute |S_i . S_j|
                    connections = sgraph.get_connected_sites(idx)
                    # dists = [round(cs[-1], 2) for cs in connections]  # i<->j distances
                    # dists = sorted(list(set(dists)))  # NN, NNN, NNNN, etc.

                    for connection in connections:
                        j_site = connection[2]
                        dist = round(connection[-1], 2)  # i_j distance

                        # s_j_sign = np.sign(sgraph.structure.site_properties['magmom'][j_site])
                        s_j = sgraph.structure.site_properties["magmom"][j_site]

                        j_index = None
                        for k, v in unique_site_ids.items():
                            if j_site in k:
                                j_index = v

                        # Determine order of connection
                        if abs(dist - dists["nn"]) <= tol:
                            order = "-nn"
                        elif abs(dist - dists["nnn"]) <= tol:
                            order = "-nnn"
                        elif abs(dist - dists["nnnn"]) <= tol:
                            order = "-nnnn"

                        j_ij = f"{i_index}-{j_index}{order}"
                        j_ji = f"{j_index}-{i_index}{order}"

                        if j_ij in ex_mat.columns:
                            ex_row.loc[sgraph_index, j_ij] -= s_i * s_j
                        elif j_ji in ex_mat.columns:
                            ex_row.loc[sgraph_index, j_ji] -= s_i * s_j

                # Ignore the row if it is a duplicate to avoid singular matrix
                # Create a temporary DataFrame with the new row
                ex_mat = ex_mat.dropna(how="all", axis=1)
                ex_row = ex_row.dropna(how="all", axis=1)
                temp_df = pd.concat([ex_mat, ex_row], ignore_index=True)
                if temp_df[j_columns].equals(temp_df[j_columns].drop_duplicates(keep="first")):
                    e_index = self.ordered_structures.index(sgraph.structure)
                    ex_row.loc[sgraph_index, "E"] = self.energies[e_index]
                    sgraph_index += 1
                    ex_mat = pd.concat([ex_mat, ex_row], ignore_index=True)
                    # if sgraph_index == num_nn_j:  # check for zero columns
                    #     zeros = [b for b in (ex_mat[j_columns] == 0).all(axis=0)]
                    #     if True in zeros:
                    #         sgraph_index -= 1  # keep looking

            ex_mat[j_columns] = ex_mat[j_columns].div(2)  # 1/2 factor in Heisenberg Hamiltonian
            ex_mat[["E0"]] = 1  # Nonmagnetic contribution

            # Check for singularities and delete columns with all zeros
            zeros = list((ex_mat == 0).all(axis=0))
            if True in zeros:
                c = ex_mat.columns[zeros.index(True)]
                ex_mat = ex_mat.drop(columns=[c], axis=1)
                # ex_mat = ex_mat.drop(ex_mat.tail(len_zeros).index)

            # Force ex_mat to be square
            ex_mat = ex_mat[: ex_mat.shape[1] - 1]

            self.ex_mat = ex_mat

    def get_exchange(self):
        """
        Take Heisenberg Hamiltonian and corresponding energy for each row and
        solve for the exchange parameters.

        Returns:
            dict[str, float]: Exchange parameters (meV/atom).
        """
        ex_mat = self.ex_mat
        # Solve the matrix equation for J_ij values
        E = ex_mat[["E"]]
        j_names = [j for j in ex_mat.columns if j != "E"]

        # Only 1 NN interaction
        if len(j_names) < 3:
            # Estimate exchange by J ~ E_AFM - E_FM
            j_avg = self.estimate_exchange()
            ex_params = {"<J>": j_avg}
            self.ex_params = ex_params

            return ex_params

        # Solve eigenvalue problem for more than 1 NN interaction
        H = np.array(ex_mat.loc[:, ex_mat.columns != "E"].values).astype(float)
        H_inv = np.linalg.inv(H)
        j_ij = np.dot(H_inv, E)

        # Convert J_ij to meV
        j_ij[1:] *= 1000  # J_ij in meV
        j_ij = j_ij.tolist()
        ex_params = {j_name: j[0] for j_name, j in zip(j_names, j_ij, strict=True)}

        self.ex_params = ex_params

        return ex_params

    def get_low_energy_orderings(self):
        """Find lowest energy FM and AFM orderings to compute E_AFM - E_FM.

        Returns:
            fm_struct (Structure): fm structure with 'magmom' site property
            afm_struct (Structure): afm structure with 'magmom' site property
            fm_e (float): fm energy
            afm_e (float): afm energy
        """
        fm_struct, afm_struct = None, None
        mag_min = np.inf
        mag_max = 0.001
        fm_e = afm_e = fm_e_min = afm_e_min = 0

        # epas = [e / len(s) for (e, s) in zip(self.energies, self.ordered_structures)]

        for s, e in zip(self.ordered_structures, self.energies, strict=True):
            ordering = CollinearMagneticStructureAnalyzer(s, threshold=0, make_primitive=False).ordering
            magmoms = s.site_properties["magmom"]

            # Try to find matching orderings first
            if ordering == Ordering.FM and e < fm_e_min:
                fm_struct = s
                mag_max = abs(sum(magmoms))
                fm_e = e
                fm_e_min = e

            if ordering == Ordering.AFM and e < afm_e_min:
                afm_struct = s
                afm_e = e
                mag_min = abs(sum(magmoms))
                afm_e_min = e

        # Brute force search for closest thing to FM and AFM
        if not fm_struct or not afm_struct:
            for s, e in zip(self.ordered_structures, self.energies, strict=True):
                magmoms = s.site_properties["magmom"]

                if abs(sum(magmoms)) > mag_max:  # FM ground state
                    fm_struct = s
                    fm_e = e
                    mag_max = abs(sum(magmoms))

                # AFM ground state
                if abs(sum(magmoms)) < mag_min:
                    afm_struct = s
                    afm_e = e
                    mag_min = abs(sum(magmoms))
                    afm_e_min = e
                elif abs(sum(magmoms)) == 0 and mag_min == 0 and e < afm_e_min:
                    afm_struct = s
                    afm_e = e
                    afm_e_min = e

        # Convert to magnetic structures with 'magmom' site property
        fm_struct = CollinearMagneticStructureAnalyzer(
            fm_struct, make_primitive=False, threshold=0.0
        ).get_structure_with_only_magnetic_atoms(make_primitive=False)

        afm_struct = CollinearMagneticStructureAnalyzer(
            afm_struct, make_primitive=False, threshold=0.0
        ).get_structure_with_only_magnetic_atoms(make_primitive=False)

        return fm_struct, afm_struct, fm_e, afm_e

    def estimate_exchange(self, fm_struct=None, afm_struct=None, fm_e=None, afm_e=None):
        """Estimate <J> for a structure based on low energy FM and AFM orderings.

        Args:
            fm_struct (Structure): fm structure with 'magmom' site property
            afm_struct (Structure): afm structure with 'magmom' site property
            fm_e (float): fm energy/atom
            afm_e (float): afm energy/atom

        Returns:
            float: Average J exchange parameter (meV/atom)
        """
        # Get low energy orderings if not supplied
        if any(arg is None for arg in [fm_struct, afm_struct, fm_e, afm_e]):
            fm_struct, afm_struct, fm_e, afm_e = self.get_low_energy_orderings()

        magmoms = fm_struct.site_properties["magmom"]

        # Normalize energies by number of magnetic ions
        # fm_e = fm_e / len(magmoms)
        # afm_e = afm_e / len(afm_magmoms)

        m_avg = np.mean([np.sqrt(m**2) for m in magmoms])

        # If m_avg for FM config is < 1 we won't get sensible results.
        if m_avg < 1:
            logger.warning(
                "Local magnetic moments are small (< 1 muB / atom). The exchange parameters may "
                "be wrong, but <J> and the mean field critical temperature estimate may be OK."
            )

        delta_e = afm_e - fm_e  # J > 0 -> FM
        j_avg = delta_e / (m_avg**2)  # eV / magnetic ion
        j_avg *= 1000  # meV / ion

        return j_avg

    def get_mft_temperature(self, j_avg):
        """
        Crude mean field estimate of critical temperature based on <J> for
        one sublattice, or solving the coupled equations for a multi-sublattice
        material.

        Args:
            j_avg (float): j_avg (float): Average exchange parameter (meV/atom)

        Returns:
            float: Critical temperature mft_t (K)
        """
        n_sub_lattices = len(self.unique_site_ids)
        k_boltzmann = 0.0861733  # meV/K

        # Only 1 magnetic sublattice
        if n_sub_lattices == 1:
            mft_t = 2 * abs(j_avg) / 3 / k_boltzmann

        else:  # multiple magnetic sublattices
            omega = np.zeros((n_sub_lattices, n_sub_lattices))
            ex_params = self.ex_params
            ex_params = {k: v for (k, v) in ex_params.items() if k != "E0"}  # ignore E0
            for k in ex_params:
                # split into i, j unique site identifiers
                sites = k.split("-")
                sites = [int(num) for num in sites[:2]]  # cut 'nn' identifier
                i, j = sites[0], sites[1]
                omega[i, j] += ex_params[k]
                omega[j, i] += ex_params[k]

            omega = omega * 2 / 3 / k_boltzmann
            eigen_vals, _eigen_vecs = np.linalg.eig(omega)
            mft_t = max(eigen_vals)

        if mft_t > 1500:  # Not sensible!
            logger.warning(
                "This mean field estimate is too high! Probably the true low energy orderings were not given as inputs."
            )

        return mft_t

    def get_interaction_graph(self, filename=None):
        """Get a StructureGraph with edges and weights that correspond to exchange
        interactions and J_ij values, respectively.

        Args:
            filename (str): if not None, save interaction graph to filename.

        Returns:
            StructureGraph: Exchange interaction graph.
        """
        structure = self.ordered_structures[0]
        sgraph = self.sgraphs[0]

        igraph = StructureGraph.from_empty_graph(
            structure, edge_weight_name="exchange_constant", edge_weight_units="meV"
        )

        if "<J>" in self.ex_params:  # Only <J> is available
            warning_msg = """
                Only <J> is available. The interaction graph will not tell
                you much.
                """
            logger.warning(warning_msg)

        # J_ij exchange interaction matrix
        for idx in range(len(sgraph.graph.nodes)):
            connections = sgraph.get_connected_sites(idx)
            for c in connections:
                jimage = c[1]  # relative integer coordinates of atom j
                j = c[2]  # index of neighbor
                dist = c[-1]  # i <-> j distance

                j_exc = self._get_j_exc(idx, j, dist)

                igraph.add_edge(idx, j, to_jimage=jimage, weight=j_exc, warn_duplicates=False)

        # Save to a JSON file if desired
        if filename:
            if not filename.endswith(".json"):
                filename += ".json"

            dumpfn(igraph, filename)

        return igraph

    def _get_j_exc(self, i, j, dist):
        """
        Convenience method for looking up exchange parameter between two sites.

        Args:
            i (int): index of ith site
            j (int): index of jth site
            dist (float): distance (Angstrom) between sites
                (10E-2 precision)

        Returns:
            float: Exchange parameter J_exc in meV
        """
        # Get unique site identifiers
        i_index = j_index = 0
        for k, v in self.unique_site_ids.items():
            if i in k:
                i_index = v
            if j in k:
                j_index = v

        # Determine order of interaction
        order = ""
        if abs(dist - self.dists["nn"]) <= self.tol:
            order = "-nn"
        elif abs(dist - self.dists["nnn"]) <= self.tol:
            order = "-nnn"
        elif abs(dist - self.dists["nnnn"]) <= self.tol:
            order = "-nnnn"

        j_ij = f"{i_index}-{j_index}{order}"
        j_ji = f"{j_index}-{i_index}{order}"

        if j_ij in self.ex_params:
            j_exc = self.ex_params[j_ij]
        elif j_ji in self.ex_params:
            j_exc = self.ex_params[j_ji]
        else:
            j_exc = 0

        # Check if only averaged NN <J> values are available
        if "<J>" in self.ex_params and order == "-nn":
            j_exc = self.ex_params["<J>"]

        return j_exc

    def get_heisenberg_model(self):
        """Save results of mapping to a HeisenbergModel object.

        Returns:
            HeisenbergModel: MSONable object.
        """
        # Original formula unit with nonmagnetic ions
        hm_formula = str(self.ordered_structures_[0].reduced_formula)

        hm_structures = self.ordered_structures
        hm_energies = self.energies
        hm_cutoff = self.cutoff
        hm_tol = self.tol
        hm_sgraphs = self.sgraphs
        hm_usi = self.unique_site_ids
        hm_wids = self.wyckoff_ids
        hm_nni = self.nn_interactions
        hm_d = self.dists

        # Exchange matrix DataFrame in JSON format
        hm_em = self.ex_mat.to_json()
        hm_ep = self.get_exchange()
        hm_javg = self.estimate_exchange()
        hm_igraph = self.get_interaction_graph()

        return HeisenbergModel(
            hm_formula,
            hm_structures,
            hm_energies,
            hm_cutoff,
            hm_tol,
            hm_sgraphs,
            hm_usi,
            hm_wids,
            hm_nni,
            hm_d,
            hm_em,
            hm_ep,
            hm_javg,
            hm_igraph,
        )


class HeisenbergScreener:
    """Clean and screen magnetic orderings."""

    def __init__(self, structures, energies, screen=False):
        """Pre-processes magnetic orderings and energies for HeisenbergMapper.
        It prioritizes low-energy orderings with large and localized magnetic moments.

        Args:
            structures (list): Structure objects with magnetic moments.
            energies (list): Energies/atom of magnetic orderings.
            screen (bool): Try to screen out high energy and low-spin configurations.

        Attributes:
            screened_structures (list): Sorted structures.
            screened_energies (list): Sorted energies.
        """
        # Cleanup
        structures, energies = self._do_cleanup(structures, energies)

        n_structures = len(structures)

        # If there are more than 2 structures, we want to perform a
        # screening to prioritize well-behaved orderings
        if screen and n_structures > 2:
            structures, energies = self._do_screen(structures, energies)

        self.screened_structures = structures
        self.screened_energies = energies

    @staticmethod
    def _do_cleanup(structures, energies):
        """Sanitize input structures and energies.

        Takes magnetic structures and performs the following operations
        - Erases nonmagnetic ions and gives all ions ['magmom'] site prop
        - Converts total energies -> energy / magnetic ion
        - Checks for duplicate/degenerate orderings
        - Sorts by energy

        Args:
            structures (list): Structure objects with magmoms.
            energies (list): Corresponding energies.

        Returns:
            ordered_structures (list): Sanitized structures.
            ordered_energies (list): Sorted energies.
        """
        # Get only magnetic ions & give all structures site_properties['magmom']
        # zero threshold so that magnetic ions with small moments
        # are preserved
        ordered_structures = [
            CollinearMagneticStructureAnalyzer(
                s, make_primitive=False, threshold=0.0
            ).get_structure_with_only_magnetic_atoms(make_primitive=False)
            for s in structures
        ]

        # Convert to energies / magnetic ion
        energies = [e / len(s) for (e, s) in zip(energies, ordered_structures, strict=True)]

        # Check for duplicate / degenerate states (sometimes different initial
        # configs relax to the same state)
        remove_list = []
        e_tol = 6  # 10^-6 eV/atom tol on energies

        for idx, energy in enumerate(energies):
            energy = round(energy, e_tol)
            if idx not in remove_list:
                for i_check, e_check in enumerate(energies):
                    e_check = round(e_check, e_tol)
                    if idx != i_check and i_check not in remove_list and energy == e_check:
                        remove_list.append(i_check)

        # Also discard structures with small |magmoms| < 0.1 uB
        # xx - get rid of these or just bury them in the list?
        # for idx, struct in enumerate(ordered_structures):
        #     magmoms = struct.site_properties["magmom"]
        #     if idx not in remove_list and any(abs(m) < 0.1 for m in magmoms):
        #         remove_list.append(idx)

        # Remove duplicates
        if remove_list:
            ordered_structures = [struct for idx, struct in enumerate(ordered_structures) if idx not in remove_list]
            energies = [energy for idx, energy in enumerate(energies) if idx not in remove_list]

        # Sort by energy if not already sorted
        ordered_structures = [s for _, s in sorted(zip(energies, ordered_structures, strict=True), reverse=False)]
        ordered_energies = sorted(energies, reverse=False)

        return ordered_structures, ordered_energies

    @staticmethod
    def _do_screen(structures, energies):
        """Screen and sort magnetic orderings based on some criteria.

        Prioritize low energy orderings and large, localized magmoms. do_clean should be run first to sanitize inputs.

        Args:
            structures (list): At least three structure objects.
            energies (list): Energies.

        Returns:
            screened_structures (list): Sorted structures.
            screened_energies (list): Sorted energies.
        """
        magmoms = [struct.site_properties["magmom"] for struct in structures]
        n_below_1ub = [sum(abs(m) < 1 for m in ms) for ms in magmoms]

        df_mag = pd.DataFrame(
            {
                "structure": structures,
                "energy": energies,
                "magmoms": magmoms,
                "n_below_1ub": n_below_1ub,
            }
        )

        # keep the ground and first excited state fixed to capture the
        # low-energy spectrum
        index = list(df_mag.index)[2:]
        df_high_energy = df_mag.iloc[2:]

        # Prioritize structures with fewer magmoms < 1 uB
        df_high_energy = df_high_energy.sort_values(by="n_below_1ub")

        index = [0, 1, *df_high_energy.index]

        # sort
        df_mag = df_mag.reindex(index)
        screened_structures = list(df_mag["structure"].values)
        screened_energies = list(df_mag["energy"].values)

        return screened_structures, screened_energies


class HeisenbergModel(MSONable):
    """
    Store a Heisenberg model fit to low-energy magnetic orderings.
    Intended to be generated by HeisenbergMapper.get_heisenberg_model().
    """

    def __init__(
        self,
        formula=None,
        structures=None,
        energies=None,
        cutoff=None,
        tol=None,
        sgraphs=None,
        unique_site_ids=None,
        wyckoff_ids=None,
        nn_interactions=None,
        dists=None,
        ex_mat=None,
        ex_params=None,
        javg=None,
        igraph=None,
    ):
        """
        Args:
            formula (str): Reduced formula of compound.
            structures (list): Structure objects with magmoms.
            energies (list): Energies of each relaxed magnetic structure.
            cutoff (float): Cutoff in Angstrom for nearest neighbor search.
            tol (float): Tolerance (in Angstrom) on nearest neighbor distances being equal.
            sgraphs (list): StructureGraph objects.
            unique_site_ids (dict): Maps each site to its unique numerical
                identifier.
            wyckoff_ids (dict): Maps unique numerical identifier to wyckoff
                position.
            nn_interactions (dict): {i: j} pairs of NN interactions
                between unique sites.
            dists (dict): NN, NNN, and NNNN interaction distances
            ex_mat (DataFrame): Invertible Heisenberg Hamiltonian for each
                graph.
            ex_params (dict): Exchange parameter values (meV/atom).
            javg (float): <J> exchange param (meV/atom).
            igraph (StructureGraph): Exchange interaction graph.
        """
        self.formula = formula
        self.structures = structures
        self.energies = energies
        self.cutoff = cutoff
        self.tol = tol
        self.sgraphs = sgraphs
        self.unique_site_ids = unique_site_ids
        self.wyckoff_ids = wyckoff_ids
        self.nn_interactions = nn_interactions
        self.dists = dists
        self.ex_mat = ex_mat
        self.ex_params = ex_params
        self.javg = javg
        self.igraph = igraph

    def as_dict(self):
        """Because some dicts have tuple keys, some sanitization is required for JSON compatibility."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "@version": __version__,
            "formula": self.formula,
            "structures": [struct.as_dict() for struct in self.structures],
            "energies": self.energies,
            "cutoff": self.cutoff,
            "tol": self.tol,
            "sgraphs": [sgraph.as_dict() for sgraph in self.sgraphs],
            "dists": self.dists,
            "ex_params": self.ex_params,
            "javg": self.javg,
            "igraph": self.igraph.as_dict(),
            # Sanitize tuple & int keys
            "ex_mat": jsanitize(self.ex_mat),
            "nn_interactions": jsanitize(self.nn_interactions),
            "unique_site_ids": jsanitize(self.unique_site_ids),
            "wyckoff_ids": jsanitize(self.wyckoff_ids),
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Create a HeisenbergModel from a dict."""
        # Reconstitute the site ids
        unique_site_ids = {}
        wyckoff_ids = {}
        nn_interactions = {}

        for k, v in dct["nn_interactions"].items():
            nn_dict = {}
            for k1, v1 in v.items():
                key = literal_eval(k1)
                nn_dict[key] = v1
            nn_interactions[k] = nn_dict

        for k, v in dct["unique_site_ids"].items():
            key = literal_eval(k)
            if isinstance(key, int):
                unique_site_ids[key,] = v
            elif isinstance(key, tuple):
                unique_site_ids[key] = v

        for k, v in dct["wyckoff_ids"].items():
            wyckoff_ids[literal_eval(k)] = v

        # Reconstitute the structure and graph objects
        structures = []
        sgraphs = []
        for v in dct["structures"]:
            structures.append(Structure.from_dict(v))
        for v in dct["sgraphs"]:
            sgraphs.append(StructureGraph.from_dict(v))

        # Interaction graph
        igraph = StructureGraph.from_dict(dct["igraph"])

        # Reconstitute the exchange matrix DataFrame
        try:
            ex_mat = literal_eval(dct["ex_mat"])
            ex_mat = pd.DataFrame.from_dict(ex_mat)
        except SyntaxError:  # if ex_mat is empty
            ex_mat = pd.DataFrame(columns=["E", "E0"])

        return HeisenbergModel(
            formula=dct["formula"],
            structures=structures,
            energies=dct["energies"],
            cutoff=dct["cutoff"],
            tol=dct["tol"],
            sgraphs=sgraphs,
            unique_site_ids=unique_site_ids,
            wyckoff_ids=wyckoff_ids,
            nn_interactions=nn_interactions,
            dists=dct["dists"],
            ex_mat=ex_mat,
            ex_params=dct["ex_params"],
            javg=dct["javg"],
            igraph=igraph,
        )

    def _get_j_exc(self, i, j, dist):
        """
        Convenience method for looking up exchange parameter between two sites.

        Args:
            i (int): index of ith site
            j (int): index of jth site
            dist (float): distance (Angstrom) between sites +- tol

        Returns:
            float: Exchange parameter J_exc in meV
        """
        # Get unique site identifiers
        i_index = j_index = 0
        for k in self.unique_site_ids:
            if i in k:
                i_index = self.unique_site_ids[k]
            if j in k:
                j_index = self.unique_site_ids[k]

        # Determine order of interaction
        order = ""
        if abs(dist - self.dists["nn"]) <= self.tol:
            order = "-nn"
        elif abs(dist - self.dists["nnn"]) <= self.tol:
            order = "-nnn"
        elif abs(dist - self.dists["nnnn"]) <= self.tol:
            order = "-nnnn"

        j_ij = f"{i_index}-{j_index}{order}"
        j_ji = f"{j_index}-{i_index}{order}"

        if j_ij in self.ex_params:
            j_exc = self.ex_params[j_ij]
        elif j_ji in self.ex_params:
            j_exc = self.ex_params[j_ji]
        else:
            j_exc = 0

        # Check if only averaged NN <J> values are available
        if "<J>" in self.ex_params and order == "-nn":
            j_exc = self.ex_params["<J>"]

        return j_exc
