---
layout: default
title: pymatgen.analysis.dimensionality.md
nav_exclude: true
---

# pymatgen.analysis.dimensionality module

This module provides functions to get the dimensionality of a structure.

A number of different algorithms are implemented. These are based on the
following publications:

get_dimensionality_larsen:


    * P. M. Larsen, M. Pandey, M. Strange, K. W. Jacobsen. Definition of a
    scoring parameter to identify low-dimensional materials components.
    Phys. Rev. Materials 3, 034003 (2019).

get_dimensionality_cheon:


    * Cheon, G.; Duerloo, K.-A. N.; Sendek, A. D.; Porter, C.; Chen, Y.; Reed,
    E. J. Data Mining for New Two- and One-Dimensional Weakly Bonded Solids
    and Lattice-Commensurate Heterostructures. Nano Lett. 2017.

get_dimensionality_gorai:


    * Gorai, P., Toberer, E. & Stevanovic, V. Computational Identification of
    Promising Thermoelectric Materials Among Known Quasi-2D Binary Compounds.
    J. Mater. Chem. A 2, 4136 (2016).


### pymatgen.analysis.dimensionality.calculate_dimensionality_of_site(bonded_structure, site_index, inc_vertices=False)
Calculates the dimensionality of the component containing the given site.

Implements directly the modified breadth-first-search algorithm described in
Algorithm 1 of:

P. M. Larsen, M. Pandey, M. Strange, K. W. Jacobsen. Definition of a
scoring parameter to identify low-dimensional materials components.
Phys. Rev. Materials 3, 034003 (2019).


* **Parameters**


    * **bonded_structure** ([*StructureGraph*](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.StructureGraph)) – A structure with bonds, represented
    as a pymatgen structure graph. For example, generated using the
    CrystalNN.get_bonded_structure() method.


    * **site_index** (*int*) – The index of a site in the component of interest.


    * **inc_vertices** (*bool**, **optional*) – Whether to return the vertices (site
    images) of the component.



* **Returns**

    If inc_vertices is False, the dimensionality of the
    component will be returned as an int. If inc_vertices is true, the
    function will return a tuple of (dimensionality, vertices), where
    vertices is a list of tuples. E.g. [(0, 0, 0), (1, 1, 1)].



* **Return type**

    (int or tuple)



### pymatgen.analysis.dimensionality.find_clusters(struct, connected_matrix)
Finds bonded clusters of atoms in the structure with periodic boundary
conditions.

If there are atoms that are not bonded to anything, returns [0,1,0]. (For
faster computation time)

Author: Gowoon Cheon
Email: [gcheon@stanford.edu](mailto:gcheon@stanford.edu)


* **Parameters**


    * **struct** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input structure


    * **connected_matrix** – Must be made from the same structure with
    find_connected_atoms() function.



* **Returns**

    the size of the largest cluster in the crystal structure
    min_cluster: the size of the smallest cluster in the crystal structure
    clusters: list of bonded clusters found here, clusters are formatted as
    sets of indices of atoms



* **Return type**

    max_cluster



### pymatgen.analysis.dimensionality.find_connected_atoms(struct, tolerance=0.45, ldict=None)
Finds bonded atoms and returns a adjacency matrix of bonded atoms.

Author: “Gowoon Cheon”
Email: “[gcheon@stanford.edu](mailto:gcheon@stanford.edu)”


* **Parameters**


    * **struct** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input structure


    * **tolerance** – length in angstroms used in finding bonded atoms. Two atoms
    are considered bonded if (radius of atom 1) + (radius of atom 2) +
    (tolerance) < (distance between atoms 1 and 2). Default
    value = 0.45, the value used by JMol and Cheon et al.


    * **ldict** – dictionary of bond lengths used in finding bonded atoms. Values
    from JMol are used as default



* **Returns**

    A numpy array of shape (number of atoms, number of atoms);
    If any image of atom j is bonded to atom i with periodic boundary
    conditions, the matrix element [atom i, atom j] is 1.



* **Return type**

    (np.ndarray)



### pymatgen.analysis.dimensionality.get_dimensionality_cheon(structure_raw, tolerance=0.45, ldict=None, standardize=True, larger_cell=False)
Algorithm for finding the dimensions of connected subunits in a structure.
This method finds the dimensionality of the material even when the material
is not layered along low-index planes, or does not have flat
layers/molecular wires.

Author: “Gowoon Cheon”
Email: “[gcheon@stanford.edu](mailto:gcheon@stanford.edu)”

See details at :

Cheon, G.; Duerloo, K.-A. N.; Sendek, A. D.; Porter, C.; Chen, Y.; Reed,
E. J. Data Mining for New Two- and One-Dimensional Weakly Bonded Solids and
Lattice-Commensurate Heterostructures. Nano Lett. 2017.


* **Parameters**


    * **structure_raw** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – A pymatgen Structure object.


    * **tolerance** (*float*) – length in angstroms used in finding bonded atoms.
    Two atoms are considered bonded if (radius of atom 1) + (radius of
    atom 2) + (tolerance) < (distance between atoms 1 and 2). Default
    value = 0.45, the value used by JMol and Cheon et al.


    * **ldict** (*dict*) – dictionary of bond lengths used in finding bonded atoms.
    Values from JMol are used as default


    * **standardize** – works with conventional standard structures if True. It is
    recommended to keep this as True.


    * **larger_cell** – tests with 3x3x3 supercell instead of 2x2x2. Testing with
    2x2x2 supercell is faster but misclassifies rare interpenetrated 3D

    > structures. Testing with a larger cell circumvents this problem




* **Returns**

    dimension of the largest cluster as a string. If there are ions
    or molecules it returns ‘intercalated ion/molecule’



* **Return type**

    (str)



### pymatgen.analysis.dimensionality.get_dimensionality_gorai(structure, max_hkl=2, el_radius_updates=None, min_slab_size=5, min_vacuum_size=5, standardize=True, bonds=None)
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


* **Parameters**


    * **structure** – (Structure) structure to analyze dimensionality for


    * **max_hkl** – (int) max index of planes to look for layers


    * **el_radius_updates** – (dict) symbol->float to update atomic radii


    * **min_slab_size** – (float) internal surface construction parameter


    * **min_vacuum_size** – (float) internal surface construction parameter


    * **standardize** (*bool*) – whether to standardize the structure before
    analysis. Set to False only if you already have the structure in a
    convention where layers / chains will be along low <hkl> indexes.


    * **bonds** (*{**(**specie1**, **specie2*) – max_bond_dist}: bonds are
    specified as a dict of tuples: float of specie1, specie2
    and the max bonding distance. For example, PO4 groups may be
    defined as {(“P”, “O”): 3}.


Returns: (int) the dimensionality of the structure - 1 (molecules/chains),

    2 (layered), or 3 (3D)


### pymatgen.analysis.dimensionality.get_dimensionality_larsen(bonded_structure)
Gets the dimensionality of a bonded structure.

The dimensionality of the structure is the highest dimensionality of all
structure components. This method is very robust and can handle
many tricky structures, regardless of structure type or improper connections
due to periodic boundary conditions.

Requires a StructureGraph object as input. This can be generated using one
of the NearNeighbor classes. For example, using the CrystalNN class:

```default
bonded_structure = CrystalNN().get_bonded_structure(structure)
```

Based on the modified breadth-first-search algorithm described in:

P. M. Larsen, M. Pandey, M. Strange, K. W. Jacobsen. Definition of a
scoring parameter to identify low-dimensional materials components.
Phys. Rev. Materials 3, 034003 (2019).


* **Parameters**

    **bonded_structure** ([*StructureGraph*](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.StructureGraph)) – A structure with bonds, represented
    as a pymatgen structure graph. For example, generated using the
    CrystalNN.get_bonded_structure() method.



* **Returns**

    The dimensionality of the structure.



* **Return type**

    (int)



### pymatgen.analysis.dimensionality.get_structure_components(bonded_structure, inc_orientation=False, inc_site_ids=False, inc_molecule_graph=False)
Gets information on the components in a bonded structure.

Correctly determines the dimensionality of all structures, regardless of
structure type or improper connections due to periodic boundary conditions.

Requires a StructureGraph object as input. This can be generated using one
of the NearNeighbor classes. For example, using the CrystalNN class:

```default
bonded_structure = CrystalNN().get_bonded_structure(structure)
```

Based on the modified breadth-first-search algorithm described in:

P. M. Larsen, M. Pandey, M. Strange, K. W. Jacobsen. Definition of a
scoring parameter to identify low-dimensional materials components.
Phys. Rev. Materials 3, 034003 (2019).


* **Parameters**


    * **bonded_structure** ([*StructureGraph*](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.StructureGraph)) – A structure with bonds, represented
    as a pymatgen structure graph. For example, generated using the
    CrystalNN.get_bonded_structure() method.


    * **inc_orientation** (*bool**, **optional*) – Whether to include the orientation
    of the structure component. For surfaces, the miller index is given,
    for one-dimensional structures, the direction of the chain is given.


    * **inc_site_ids** (*bool**, **optional*) – Whether to include the site indices
    of the sites in the structure component.


    * **inc_molecule_graph** (*bool**, **optional*) – Whether to include MoleculeGraph
    objects for zero-dimensional components.



* **Returns**

    Information on the components in a structure as a list
    of dictionaries with the keys:


    * ”structure_graph”: A pymatgen StructureGraph object for the

        component.


    * ”dimensionality”: The dimensionality of the structure component as an

        int.


    * ”orientation”: If inc_orientation is True, the orientation of the

        component as a tuple. E.g. (1, 1, 1)


    * ”site_ids”: If inc_site_ids is True, the site indices of the

        sites in the component as a tuple.


    * ”molecule_graph”: If inc_molecule_graph is True, the site a

        MoleculeGraph object for zero-dimensional components.




* **Return type**

    (list of dict)



### pymatgen.analysis.dimensionality.zero_d_graph_to_molecule_graph(bonded_structure, graph)
Converts a zero-dimensional networkx Graph object into a MoleculeGraph.

Implements a similar breadth-first search to that in
calculate_dimensionality_of_site().


* **Parameters**


    * **bonded_structure** ([*StructureGraph*](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.StructureGraph)) – A structure with bonds, represented
    as a pymatgen structure graph. For example, generated using the
    CrystalNN.get_bonded_structure() method.


    * **graph** (*nx.Graph*) – A networkx Graph object for the component of
    interest.



* **Returns**

    A MoleculeGraph object of the component.



* **Return type**

    ([MoleculeGraph](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.MoleculeGraph))