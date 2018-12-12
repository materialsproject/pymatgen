"""
This module provides functions to get the dimensionality of a structure.

A number of different algorithms are implemented. These are based on the
following publications:

get_dimensionality_larsen:
  - P. Larsen, M. Pandey, M. Strange, K. W. Jacobsen, 2018, arXiv:1808.02114

get_dimensionality_cheon:
  - Cheon, G.; Duerloo, K.-A. N.; Sendek, A. D.; Porter, C.; Chen, Y.; Reed,
    E. J. Data Mining for New Two- and One-Dimensional Weakly Bonded Solids
    and Lattice-Commensurate Heterostructures. Nano Lett. 2017.

get_dimensionality_gorai:
  - Gorai, P., Toberer, E. & Stevanovic, V. Computational Identification of
    Promising Thermoelectric Materials Among Known Quasi-2D Binary Compounds.
    J. Mater. Chem. A 2, 4136 (2016).
"""


import itertools
import copy

import numpy as np

from collections import defaultdict

from networkx.readwrite import json_graph

from pymatgen.analysis.graphs import MoleculeGraph, StructureGraph
from pymatgen.core.lattice import get_integer_index
from pymatgen.core.structure import Structure, Molecule
from pymatgen.core.periodic_table import Specie
from pymatgen.core.surface import SlabGenerator
from pymatgen.analysis.local_env import JmolNN
from pymatgen.analysis.structure_analyzer import get_max_bond_lengths
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

__author__ = "Alex Ganose, Gowoon Cheon, Prashun Gorai"


def get_dimensionality_larsen(bonded_structure):
    """
    Gets the dimensionality of a bonded structure.

    The dimensionality of the structure is the highest dimensionality of all
    structure components. This method is very robust and can handle
    many tricky structures, regardless of structure type or improper connections
    due to periodic boundary conditions.

    Requires a StructureGraph object as input. This can be generated using one
    of the NearNeighbor classes. For example, using the CrystalNN class::

        bonded_structure = CrystalNN().get_bonded_structure(structure)

    Based on the modified breadth-first-search algorithm described in:

    P. Larsem, M. Pandey, M. Strange, K. W. Jacobsen, 2018, arXiv:1808.02114

    Args:
        bonded_structure (StructureGraph): A structure with bonds, represented
            as a pymatgen structure graph. For example, generated using the
            CrystalNN.get_bonded_structure() method.

    Returns:
        (int): The dimensionality of the structure.
    """
    return max([c['dimensionality'] for c in
                get_structure_components(bonded_structure)])


def get_structure_components(bonded_structure, inc_orientation=False,
                             inc_site_ids=False, inc_molecule_graph=False):
    """
    Gets information on the components in a bonded structure.

    Correctly determines the dimensionality of all structures, regardless of
    structure type or improper connections due to periodic boundary conditions.

    Requires a StructureGraph object as input. This can be generated using one
    of the NearNeighbor classes. For example, using the CrystalNN class::

        bonded_structure = CrystalNN().get_bonded_structure(structure)

    Based on the modified breadth-first-search algorithm described in:

    P. Larsem, M. Pandey, M. Strange, K. W. Jacobsen, 2018, arXiv:1808.02114

    Args:
        bonded_structure (StructureGraph): A structure with bonds, represented
            as a pymatgen structure graph. For example, generated using the
            CrystalNN.get_bonded_structure() method.
        inc_orientation (bool, optional): Whether to include the orientation
            of the structure component. For surfaces, the miller index is given,
            for one-dimensional structures, the direction of the chain is given.
        inc_site_ids (bool, optional): Whether to include the site indices
            of the sites in the structure component.
        inc_molecule_graph (bool, optional): Whether to include MoleculeGraph
            objects for zero-dimensional components.

    Returns:
        (list of dict): Information on the components in a structure as a list
        of dictionaries with the keys:

        - "structure_graph": A pymatgen StructureGraph object for the
            component.
        - "dimensionality": The dimensionality of the structure component as an
            int.
        - "orientation": If inc_orientation is `True`, the orientation of the
            component as a tuple. E.g. (1, 1, 1)
        - "site_ids": If inc_site_ids is `True`, the site indices of the
            sites in the component as a tuple.
        - "molecule_graph": If inc_molecule_graph is `True`, the site a
            MoleculeGraph object for zero-dimensional components.
    """
    import networkx as nx  # optional dependency therefore not top level import

    comp_graphs = (bonded_structure.graph.subgraph(c) for c in
                   nx.weakly_connected_components(bonded_structure.graph))

    components = []
    for graph in comp_graphs:
        dimensionality, vertices = calculate_dimensionality_of_site(
            bonded_structure, list(graph.nodes())[0], inc_vertices=True)

        component = {'dimensionality': dimensionality}

        if inc_orientation:
            if dimensionality in [1, 2]:
                vertices = np.array(vertices)

                g = vertices.sum(axis=0) / vertices.shape[0]

                # run singular value decomposition
                _, _, vh = np.linalg.svd(vertices - g)

                # get direction (first column is best fit line,
                # 3rd column is unitary norm)
                index = 2 if dimensionality == 2 else 0
                orientation = get_integer_index(vh[index, :])
            else:
                orientation = None

            component['orientation'] = orientation

        if inc_site_ids:
            component['site_ids'] = tuple(graph.nodes())

        if inc_molecule_graph and dimensionality == 0:
            component['molecule_graph'] = zero_d_graph_to_molecule_graph(
                bonded_structure, graph)

        component_structure = Structure.from_sites(
            [bonded_structure.structure[n] for n in sorted(graph.nodes())])

        sorted_graph = nx.convert_node_labels_to_integers(
            graph, ordering="sorted")
        component_graph = StructureGraph(
            component_structure,
            graph_data=json_graph.adjacency_data(sorted_graph))
        component['structure_graph'] = component_graph

        components.append(component)
    return components


def calculate_dimensionality_of_site(bonded_structure, site_index,
                                     inc_vertices=False):
    """
    Calculates the dimensionality of the component containing the given site.

    Implements directly the modified breadth-first-search algorithm described in
    Algorithm 1 of:

    P. Larsem, M. Pandey, M. Strange, K. W. Jacobsen, 2018, arXiv:1808.02114

    Args:
        bonded_structure (StructureGraph): A structure with bonds, represented
            as a pymatgen structure graph. For example, generated using the
            CrystalNN.get_bonded_structure() method.
        site_index (int): The index of a site in the component of interest.
        inc_vertices (bool, optional): Whether to return the vertices (site
            images) of the component.

    Returns:
        (int or tuple): If inc_vertices is False, the dimensionality of the
        component will be returned as an int. If inc_vertices is true, the
        function will return a tuple of (dimensionality, vertices), where
        vertices is a list of tuples. E.g. [(0, 0, 0), (1, 1, 1)].
    """
    def neighbours(comp_index):
        return [(s.index, s.jimage) for s
                in bonded_structure.get_connected_sites(comp_index)]

    def rank(vertices):
        if len(vertices) == 0:
            return -1
        elif len(vertices) == 1:
            return 0
        else:
            vertices = np.array(list(vertices))
            return np.linalg.matrix_rank(vertices[1:] - vertices[0])

    connected_sites = {i: neighbours(i) for i in
                       range(bonded_structure.structure.num_sites)}

    seen_vertices = set()
    seen_comp_vertices = defaultdict(set)

    queue = [(site_index, (0, 0, 0))]
    while len(queue) > 0:
        comp_i, image_i = queue.pop(0)

        if (comp_i, image_i) in seen_vertices:
            continue
        seen_vertices.add((comp_i, image_i))

        if (rank(seen_comp_vertices[comp_i].union({image_i})) >
                rank(seen_comp_vertices[comp_i])):
            seen_comp_vertices[comp_i].add(image_i)

        for comp_j, image_j in connected_sites[comp_i]:

            image_j = tuple(np.add(image_j, image_i))

            if (comp_j, image_j) in seen_vertices:
                continue

            if (rank(seen_comp_vertices[comp_j].union({image_j})) >
                    rank(seen_comp_vertices[comp_j])):
                queue.append((comp_j, image_j))

    if inc_vertices:
        return (rank(seen_comp_vertices[site_index]),
                list(seen_comp_vertices[site_index]))
    else:
        return rank(seen_comp_vertices[site_index])


def zero_d_graph_to_molecule_graph(bonded_structure, graph):
    """
    Converts a zero-dimensional networkx Graph object into a MoleculeGraph.

    Implements a similar breadth-first search to that in
    calculate_dimensionality_of_site().

    Args:
        bonded_structure (StructureGraph): A structure with bonds, represented
            as a pymatgen structure graph. For example, generated using the
            CrystalNN.get_bonded_structure() method.
        graph (nx.Graph): A networkx `Graph` object for the component of
            interest.

    Returns:
        (MoleculeGraph): A MoleculeGraph object of the component.
    """
    import networkx as nx
    seen_indices = []
    sites = []

    start_index = list(graph.nodes())[0]
    queue = [(start_index, (0, 0, 0),
              bonded_structure.structure[start_index])]
    while len(queue) > 0:
        comp_i, image_i, site_i = queue.pop(0)

        if comp_i in [x[0] for x in seen_indices]:
            raise ValueError("Graph component is not 0D")

        seen_indices.append((comp_i, image_i))
        sites.append(site_i)

        for site_j in bonded_structure.get_connected_sites(
                comp_i, jimage=image_i):

            if ((site_j.index, site_j.jimage) not in seen_indices and
                    (site_j.index, site_j.jimage, site_j.site) not in queue):
                queue.append((site_j.index, site_j.jimage, site_j.site))

    # sort the list of indices and the graph by index to make consistent
    indices_ordering = np.argsort([x[0] for x in seen_indices])
    sorted_sites = np.array(sites, dtype=object)[indices_ordering]
    sorted_graph = nx.convert_node_labels_to_integers(graph, ordering="sorted")
    mol = Molecule([s.specie for s in sorted_sites],
                   [s.coords for s in sorted_sites])
    mol_graph = MoleculeGraph.with_edges(mol, nx.Graph(sorted_graph).edges())

    return mol_graph


def get_dimensionality_cheon(structure_raw, tolerance=0.45,
                             ldict=JmolNN().el_radius, standardize=True):
    """
    Algorithm for finding the dimensions of connected subunits in a structure.
    This method finds the dimensionality of the material even when the material
    is not layered along low-index planes, or does not have flat
    layers/molecular wires.

    Author: "Gowoon Cheon"
    Email: "gcheon@stanford.edu"

    See details at :

    Cheon, G.; Duerloo, K.-A. N.; Sendek, A. D.; Porter, C.; Chen, Y.; Reed,
    E. J. Data Mining for New Two- and One-Dimensional Weakly Bonded Solids and
    Lattice-Commensurate Heterostructures. Nano Lett. 2017.

    Args:
        structure_raw (Structure): A pymatgen Structure object.
        tolerance (float): length in angstroms used in finding bonded atoms.
            Two atoms are considered bonded if (radius of atom 1) + (radius of
            atom 2) + (tolerance) < (distance between atoms 1 and 2). Default
            value = 0.45, the value used by JMol and Cheon et al.
        ldict (dict): dictionary of bond lengths used in finding bonded atoms.
            Values from JMol are used as default
        standardize: works with conventional standard structures if True. It is
            recommended to keep this as True.

    Returns:
        (str): dimension of the largest cluster as a string. If there are ions
        or molecules it returns 'intercalated ion/molecule'
    """
    if standardize:
        structure = SpacegroupAnalyzer(structure_raw).\
            get_conventional_standard_structure()
    else:
        structure = structure_raw

    structure_save = copy.copy(structure_raw)
    connected_list1 = find_connected_atoms(
        structure, tolerance=tolerance, ldict=ldict)
    max1, min1, _ = find_clusters(structure, connected_list1)
    structure.make_supercell([[2, 0, 0], [0, 2, 0], [0, 0, 2]])
    connected_list2 = find_connected_atoms(
        structure, tolerance=tolerance, ldict=ldict)
    max2, min2, _ = find_clusters(structure, connected_list2)

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
            structure = copy.copy(structure_save)
            structure.make_supercell([[3, 0, 0], [0, 3, 0], [0, 0, 3]])
            connected_list3 = find_connected_atoms(
                structure, tolerance=tolerance, ldict=ldict)
            max3, min3, _ = find_clusters(structure, connected_list3)
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


def find_connected_atoms(struct, tolerance=0.45, ldict=JmolNN().el_radius):
    """
    Finds the list of bonded atoms.

    Author: "Gowoon Cheon"
    Email: "gcheon@stanford.edu"

    Args:
        struct (Structure): Input structure
        tolerance: length in angstroms used in finding bonded atoms. Two atoms
            are considered bonded if (radius of atom 1) + (radius of atom 2) +
            (tolerance) < (distance between atoms 1 and 2). Default
            value = 0.45, the value used by JMol and Cheon et al.
        ldict: dictionary of bond lengths used in finding bonded atoms. Values
            from JMol are used as default

    Returns:
        (np.ndarray): A numpy array of shape (number of bonded pairs, 2); each
        row of is of the form [atomi, atomj]. atomi and atomj are the indices of
        the atoms in the input structure. If any image of atomj is bonded to
        atomi with periodic boundary conditions, [atomi, atomj] is included in
        the list. If atomi is bonded to multiple images of atomj, it is only
        counted once.
    """
    n_atoms = len(struct.species)
    fc = np.array(struct.frac_coords)
    species = list(map(str, struct.species))

    # in case of charged species
    for i, item in enumerate(species):
        if item not in ldict.keys():
            species[i] = str(Specie.from_string(item).element)
    latmat = struct.lattice.matrix
    connected_list = []

    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            max_bond_length = ldict[species[i]] + ldict[species[j]] + tolerance
            add_ij = False
            for move_cell in itertools.product(
                    [0, 1, -1], [0, 1, -1], [0, 1, -1]):
                if not add_ij:
                    frac_diff = fc[j] + move_cell - fc[i]
                    distance_ij = np.dot(latmat.T, frac_diff)
                    if np.linalg.norm(distance_ij) < max_bond_length:
                        add_ij = True
            if add_ij:
                connected_list.append([i, j])
    return np.array(connected_list)


def find_clusters(struct, connected_list):
    """
    Finds bonded clusters of atoms in the structure with periodic boundary
    conditions.

    If there are atoms that are not bonded to anything, returns [0,1,0]. (For
    faster computation time in FindDimension())

    Author: "Gowoon Cheon"
    Email: "gcheon@stanford.edu"

    Args:
        struct (Structure): Input structure
        connected_list: Must be made from the same structure with
            FindConnected() function. An array of shape (number of bonded
            pairs, 2); each row of is of the form [atomi, atomj].

    Returns:
        max_cluster: the size of the largest cluster in the crystal structure
        min_cluster: the size of the smallest cluster in the crystal structure
        clusters: list of bonded clusters found here, clusters are formatted as
        sets of indices of atoms
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
            new_cluster_sizes = (clusters_noatom_sizes +
                                 [len(clusters_w_atom[0])])
        clusters = list(new_clusters)
        cluster_sizes = list(new_cluster_sizes)
        if n_atoms in cluster_sizes:
            break
    max_cluster = max(cluster_sizes)
    min_cluster = min(cluster_sizes)
    return [max_cluster, min_cluster, clusters]


def get_dimensionality_gorai(structure, max_hkl=2, el_radius_updates=None,
                             min_slab_size=5, min_vacuum_size=5,
                             standardize=True, bonds=None):
    """
    This method returns whether a structure is 3D, 2D (layered), or 1D (linear
    chains or molecules) according to the algorithm published in Gorai, P.,
    Toberer, E. & Stevanovic, V. Computational Identification of Promising
    Thermoelectric Materials Among Known Quasi-2D Binary Compounds. J. Mater.
    Chem. A 2, 4136 (2016).

    Note that a 1D structure detection might indicate problems in the bonding
    algorithm, particularly for ionic crystals (e.g., NaCl)

    Users can change the behavior of bonds detection by passing either
    el_radius_updates to update atomic radii for auto-detection of max bond
    distances, or bonds to explicitly specify max bond distances for atom pairs.
    Note that if you pass both, el_radius_updates are ignored.

    Args:
        structure: (Structure) structure to analyze dimensionality for
        max_hkl: (int) max index of planes to look for layers
        el_radius_updates: (dict) symbol->float to update atomic radii
        min_slab_size: (float) internal surface construction parameter
        min_vacuum_size: (float) internal surface construction parameter
        standardize (bool): whether to standardize the structure before
            analysis. Set to False only if you already have the structure in a
            convention where layers / chains will be along low <hkl> indexes.
        bonds ({(specie1, specie2): max_bond_dist}: bonds are
                specified as a dict of tuples: float of specie1, specie2
                and the max bonding distance. For example, PO4 groups may be
                defined as {("P", "O"): 3}.

    Returns: (int) the dimensionality of the structure - 1 (molecules/chains),
        2 (layered), or 3 (3D)

    """
    if standardize:
        structure = SpacegroupAnalyzer(structure). \
            get_conventional_standard_structure()

    if not bonds:
        bonds = get_max_bond_lengths(structure, el_radius_updates)

    num_surfaces = 0
    for h in range(max_hkl):
        for k in range(max_hkl):
            for l in range(max_hkl):
                if max([h, k, l]) > 0 and num_surfaces < 2:
                    sg = SlabGenerator(structure, (h, k, l),
                                       min_slab_size=min_slab_size,
                                       min_vacuum_size=min_vacuum_size)
                    slabs = sg.get_slabs(bonds)
                    for _ in slabs:
                        num_surfaces += 1

    return 3 - min(num_surfaces, 2)
