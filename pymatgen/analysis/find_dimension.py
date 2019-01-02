# coding: utf-8


"""
This module provides functions for finding the dimensions of connected subunits in a crystal structure.
This method finds the dimensionality of the material even when the material is not layered along low-index planes, or does not have flat layers/molecular wires.
See details at : Cheon, G.; Duerloo, K.-A. N.; Sendek, A. D.; Porter, C.; Chen, Y.; Reed, E. J. Data Mining for New Two- and One-Dimensional Weakly Bonded Solids and Lattice-Commensurate Heterostructures. Nano Lett. 2017.

"""

__author__ = "Gowoon Cheon"
__version__ = "1.0"
__email__ = "gcheon@stanford.edu"

import numpy as np
import itertools
import copy

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.local_env import JmolNN
from pymatgen.core.periodic_table import Specie
from monty.dev import deprecated


@deprecated(message=("find_connected_atoms has been moved to "
                     "pymatgen.analysis.dimensionality.find_connected_atoms. "
                     "This method will be removed in pymatgen v2019.1.1."))
def find_connected_atoms(struct, tolerance=0.45, ldict=JmolNN().el_radius):
    """
    Finds the list of bonded atoms.

    Args:
        struct (Structure): Input structure
        tolerance: length in angstroms used in finding bonded atoms. Two atoms are considered bonded if (radius of atom 1) + (radius of atom 2) + (tolerance) < (distance between atoms 1 and 2). Default value = 0.45, the value used by JMol and Cheon et al.
        ldict: dictionary of bond lengths used in finding bonded atoms. Values from JMol are used as default
        standardize: works with conventional standard structures if True. It is recommended to keep this as True.

    Returns:
        connected_list: A numpy array of shape (number of bonded pairs, 2); each row of is of the form [atomi, atomj].
        atomi and atomj are the indices of the atoms in the input structure.
        If any image of atomj is bonded to atomi with periodic boundary conditions, [atomi, atomj] is included in the list.
        If atomi is bonded to multiple images of atomj, it is only counted once.
    """
    n_atoms = len(struct.species)
    fc = np.array(struct.frac_coords)
    species = list(map(str, struct.species))
    #in case of charged species
    for i,item in enumerate(species):
        if not item in ldict.keys():
            species[i]=str(Specie.from_string(item).element)
    latmat = struct.lattice.matrix
    connected_list = []

    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            max_bond_length = ldict[species[i]] + ldict[species[j]] + tolerance
            add_ij = False
            for move_cell in itertools.product([0, 1, -1], [0, 1, -1], [0, 1, -1]):
                if not add_ij:
                    frac_diff = fc[j] + move_cell - fc[i]
                    distance_ij = np.dot(latmat.T, frac_diff)
                    if np.linalg.norm(distance_ij) < max_bond_length:
                        add_ij = True
            if add_ij:
                connected_list.append([i, j])
    return np.array(connected_list)


@deprecated(message=("find_clusters has been moved to"
                     "pymatgen.analysis.dimensionality.find_clusters. "
                     "This method will be removed in pymatgen v2019.1.1."))
def find_clusters(struct, connected_list):
    """
    Finds bonded clusters of atoms in the structure with periodic boundary conditions.
    If there are atoms that are not bonded to anything, returns [0,1,0].(For faster computation time in FindDimension())

    Args:
        struct (Structure): Input structure
        connected_list: Must be made from the same structure with FindConnected() function.
            An array of shape (number of bonded pairs, 2); each row of is of the form [atomi, atomj].

    Returns:
        max_cluster: the size of the largest cluster in the crystal structure
        min_cluster: the size of the smallest cluster in the crystal structure
        clusters: list of bonded clusters found here, clusters are formatted as sets of indices of atoms
    """
    n_atoms = len(struct.species)
    if len(np.unique(connected_list)) != n_atoms:
        return [0, 1, 0]
    if n_atoms == 0:
        return [0, 0, 0]
    cluster_sizes = []
    clusters = []
    for atom in range(n_atoms):
        connected_inds = np.where(connected_list == atom)[0]
        atom_cluster = np.unique(connected_list[connected_inds])
        atom_cluster = set(atom_cluster)
        if len(clusters) == 0:
            new_clusters = [atom_cluster]
            new_cluster_sizes = [len(atom_cluster)]
        else:
            clusters_w_atom = [atom_cluster]
            clusters_noatom = []
            clusters_noatom_sizes = []
            for cluster in clusters:
                if len(cluster.intersection(atom_cluster)) > 0:
                    clusters_w_atom.append(cluster)
                else:
                    clusters_noatom.append(cluster)
                    clusters_noatom_sizes.append(len(cluster))
            if len(clusters_w_atom) > 1:
                clusters_w_atom = [set.union(*clusters_w_atom)]
            new_clusters = clusters_noatom + clusters_w_atom
            new_cluster_sizes = clusters_noatom_sizes + [len(clusters_w_atom[0])]
        clusters = list(new_clusters)
        cluster_sizes = list(new_cluster_sizes)
        if n_atoms in cluster_sizes:
            break
    max_cluster = max(cluster_sizes)
    min_cluster = min(cluster_sizes)
    return [max_cluster, min_cluster, clusters]


@deprecated(message=("find_dimension has been moved to"
                     "pymatgen.analysis.dimensionality.get_dimensionality_cheon. "
                     "This method will be removed in pymatgen v2019.1.1."))
def find_dimension(structure_raw, tolerance=0.45, ldict=JmolNN().el_radius, standardize=True):
    """
    Algorithm for finding the dimensions of connected subunits in a crystal structure.
    This method finds the dimensionality of the material even when the material is not layered along low-index planes, or does not have flat layers/molecular wires.
    See details at : Cheon, G.; Duerloo, K.-A. N.; Sendek, A. D.; Porter, C.; Chen, Y.; Reed, E. J. Data Mining for New Two- and One-Dimensional Weakly Bonded Solids and Lattice-Commensurate Heterostructures. Nano Lett. 2017.

    Args:
        structure (Structure): Input structure
        tolerance: length in angstroms used in finding bonded atoms. Two atoms are considered bonded if (radius of atom 1) + (radius of atom 2) + (tolerance) < (distance between atoms 1 and 2). Default value = 0.45, the value used by JMol and Cheon et al.
        ldict: dictionary of bond lengths used in finding bonded atoms. Values from JMol are used as default
        standardize: works with conventional standard structures if True. It is recommended to keep this as True.

    Returns:
        dim: dimension of the largest cluster as a string. If there are ions or molecules it returns 'intercalated ion/molecule'
    """
    if standardize:
        structure = SpacegroupAnalyzer(structure_raw).get_conventional_standard_structure()
    structure_save = copy.copy(structure_raw)
    connected_list1 = find_connected_atoms(structure, tolerance=tolerance, ldict=ldict)
    max1, min1, clusters1 = find_clusters(structure, connected_list1)
    structure.make_supercell([[2, 0, 0], [0, 2, 0], [0, 0, 2]])
    connected_list2 = find_connected_atoms(structure, tolerance=tolerance, ldict=ldict)
    max2, min2, clusters2 = find_clusters(structure, connected_list2)
    if min2 == 1:
        dim = 'intercalated ion'
    elif min2 == min1:
        if max2 == max1:
            dim = '0D'
        else:
            dim = 'intercalated molecule'
    else:
        dim = np.log2(float(max2) / max1)
        if dim == int(dim):
            dim = str(int(dim)) + 'D'
        else:
            structure=copy.copy(structure_save)
            structure.make_supercell([[3, 0, 0], [0, 3, 0], [0, 0, 3]])
            connected_list3 = find_connected_atoms(structure, tolerance=tolerance, ldict=ldict)
            max3, min3, clusters3 = find_clusters(structure, connected_list3)
            if min3 == min2:
                if max3 == max2:
                    dim = '0D'
                else:
                    dim = 'intercalated molecule'
            else:
                dim = np.log2(float(max3) / max1) / np.log2(3)
                if dim == int(dim):
                    dim = str(int(dim)) + 'D'
                else:
                    return
    return dim
