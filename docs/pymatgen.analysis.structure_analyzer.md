---
layout: default
title: pymatgen.analysis.structure_analyzer.md
nav_exclude: true
---

# pymatgen.analysis.structure_analyzer module

This module provides classes to perform topological analyses of structures.


### _class_ pymatgen.analysis.structure_analyzer.OxideType(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), relative_cutoff=1.1)
Bases: `object`

Separate class for determining oxide type.


* **Parameters**


    * **structure** – Input structure.


    * **relative_cutoff** – Relative_cutoff \* act. cutoff stipulates the max.
    distance two O atoms must be from each other. Default value is
    1.1. At most 1.1 is recommended, nothing larger, otherwise the
    script cannot distinguish between superoxides and peroxides.



#### parse_oxide()
Determines if an oxide is a peroxide/superoxide/ozonide/normal oxide.


* **Returns**

    Type of oxide
    ozonide/peroxide/superoxide/hydroxide/None.
    nbonds (int): Number of peroxide/superoxide/hydroxide bonds in structure.



* **Return type**

    oxide_type (str)



### _class_ pymatgen.analysis.structure_analyzer.RelaxationAnalyzer(initial_structure, final_structure)
Bases: `object`

This class analyzes the relaxation in a calculation.

Please note that the input and final structures should have the same
ordering of sites. This is typically the case for most computational
codes.


* **Parameters**


    * **initial_structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Initial input structure to
    calculation.


    * **final_structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Final output structure from
    calculation.



#### get_percentage_bond_dist_changes(max_radius=3.0)
Returns the percentage bond distance changes for each site up to a
maximum radius for nearest neighbors.


* **Parameters**

    **max_radius** (*float*) – Maximum radius to search for nearest
    neighbors. This radius is applied to the initial structure,
    not the final structure.



* **Returns**

    Bond distance changes as a dict of dicts. E.g.,
    {index1: {index2: 0.011, …}}. For economy of representation, the
    index1 is always less than index2, i.e., since bonding between
    site1 and site_n is the same as bonding between site_n and site1,
    there is no reason to duplicate the information or computation.



#### get_percentage_lattice_parameter_changes()
Returns the percentage lattice parameter changes.


* **Returns**

    A dict of the percentage change in lattice parameter, e.g.,
    {‘a’: 0.012, ‘b’: 0.021, ‘c’: -0.031} implies a change of 1.2%,
    2.1% and -3.1% in the a, b and c lattice parameters respectively.



#### get_percentage_volume_change()
Returns the percentage volume change.


* **Returns**

    Volume change in percentage, e.g., 0.055 implies a 5.5% increase.



### _class_ pymatgen.analysis.structure_analyzer.VoronoiAnalyzer(cutoff=5.0, qhull_options='Qbb Qc Qz')
Bases: `object`

Performs a statistical analysis of Voronoi polyhedra around each site.
Each Voronoi polyhedron is described using Schaefli notation.
That is a set of indices {c_i} where c_i is the number of faces with i
number of vertices. E.g. for a bcc crystal, there is only one polyhedron
notation of which is [0,6,0,8,0,0,…].
In perfect crystals, these also corresponds to the Wigner-Seitz cells.
For distorted-crystals, liquids or amorphous structures, rather than one-type,
there is a statistical distribution of polyhedra.
See ref: Microstructure and its relaxation in Fe-B amorphous system
simulated by molecular dynamics,

> Stepanyuk et al., J. Non-cryst. Solids (1993), 159, 80-87.


* **Parameters**


    * **cutoff** (*float*) – cutoff distance to search for neighbors of a given atom
    (default = 5.0)


    * **qhull_options** (*str*) – options to pass to qhull (optional).



#### analyze(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), n=0)
Performs Voronoi analysis and returns the polyhedra around atom n
in Schlaefli notation.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – structure to analyze


    * **n** (*int*) – index of the center atom in structure



* **Returns**

    <c3,c4,c6,c6,c7,c8,c9,c10>

        where c_i denotes number of facets with i vertices.




* **Return type**

    voronoi index of n



#### analyze_structures(structures, step_freq=10, most_frequent_polyhedra=15)
Perform Voronoi analysis on a list of Structures.
Note that this might take a significant amount of time depending on the
size and number of structures.


* **Parameters**


    * **structures** (*list*) – list of Structures


    * **(****float** (*cutoff*) – cutoff distance around an atom to search for
    neighbors


    * **step_freq** (*int*) – perform analysis every step_freq steps


    * **qhull_options** (*str*) – options to pass to qhull


    * **most_frequent_polyhedra** (*int*) – this many unique polyhedra with
    highest frequencies is stored.



* **Returns**

    A list of tuples in the form (voronoi_index,frequency)



#### _static_ plot_vor_analysis(voronoi_ensemble)
Plot the Voronoi analysis.


* **Parameters**

    **voronoi_ensemble** –



* **Returns**

    matplotlib.pyplot



### _class_ pymatgen.analysis.structure_analyzer.VoronoiConnectivity(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), cutoff=10)
Bases: `object`

Computes the solid angles swept out by the shared face of the voronoi
polyhedron between two sites.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input structure


    * **cutoff** (*float*) –



#### _property_ connectivity_array()
Provides connectivity array.


* **Returns**

    An array of shape [atom_i, atom_j, image_j]. atom_i is
    the index of the atom in the input structure. Since the second
    atom can be outside of the unit cell, it must be described
    by both an atom index and an image index. Array data is the
    solid angle of polygon between atom_i and image_j of atom_j



* **Return type**

    connectivity



#### get_connections()
Returns a list of site pairs that are Voronoi Neighbors, along
with their real-space distances.


#### get_sitej(site_index, image_index)
Assuming there is some value in the connectivity array at indices
(1, 3, 12). sitei can be obtained directly from the input structure
(structure[1]). sitej can be obtained by passing 3, 12 to this function.


* **Parameters**


    * **site_index** (*int*) – index of the site (3 in the example)


    * **image_index** (*int*) – index of the image (12 in the example)



#### _property_ max_connectivity()
Returns the 2d array [site_i, site_j] that represents the maximum connectivity of
site i to any periodic image of site j.


### pymatgen.analysis.structure_analyzer.average_coordination_number(structures, freq=10)
Calculates the ensemble averaged Voronoi coordination numbers
of a list of Structures using VoronoiNN.
Typically used for analyzing the output of a Molecular Dynamics run.


* **Parameters**


    * **structures** (*list*) – list of Structures.


    * **freq** (*int*) – sampling frequency of coordination number [every freq steps].



* **Returns**

    Dictionary of elements as keys and average coordination numbers as values.



### pymatgen.analysis.structure_analyzer.contains_peroxide(structure, relative_cutoff=1.1)
Determines if a structure contains peroxide anions.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input structure.


    * **relative_cutoff** – The peroxide bond distance is 1.49 Angstrom.
    Relative_cutoff \* 1.49 stipulates the maximum distance two O
    atoms must be to each other to be considered a peroxide.



* **Returns**

    Boolean indicating if structure contains a peroxide anion.



### pymatgen.analysis.structure_analyzer.get_max_bond_lengths(structure, el_radius_updates=None)
Provides max bond length estimates for a structure based on the JMol
table and algorithms.


* **Parameters**


    * **structure** – (structure)


    * **el_radius_updates** – (dict) symbol->float to update atom_ic radii


Returns: (dict) - (Element1, Element2) -> float. The two elements are

    ordered by Z.


### pymatgen.analysis.structure_analyzer.oxide_type(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), relative_cutoff: float = 1.1, return_nbonds: bool = False)
Determines if an oxide is a peroxide/superoxide/ozonide/normal oxide.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input structure.


    * **relative_cutoff** (*float*) – Relative_cutoff \* act. cutoff stipulates the
    max distance two O atoms must be from each other.


    * **return_nbonds** (*bool*) – Should number of bonds be requested?



### pymatgen.analysis.structure_analyzer.solid_angle(center, coords)
Helper method to calculate the solid angle of a set of coords from the
center.


* **Parameters**


    * **center** (*3x1 array*) – Center to measure solid angle from.


    * **coords** (*Nx3 array*) – List of coords to determine solid angle.



* **Returns**

    The solid angle.



### pymatgen.analysis.structure_analyzer.sulfide_type(structure)
Determines if a structure is a sulfide/polysulfide/sulfate.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input structure.



* **Returns**

    (str) sulfide/polysulfide or None if structure is a sulfate.