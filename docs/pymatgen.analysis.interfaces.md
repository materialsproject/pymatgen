---
layout: default
title: pymatgen.analysis.interfaces.md
nav_exclude: true
---

1. TOC
{:toc}

# pymatgen.analysis.interfaces package

Module that implements various algorithms related to interface construction and analysis.


## pymatgen.analysis.interfaces.coherent_interfaces module

This module provides classes to store, generate, and manipulate material interfaces.


### _class_ CoherentInterfaceBuilder(substrate_structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), film_structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), film_miller: tuple[int, int, int], substrate_miller: tuple[int, int, int], zslgen: ZSLGenerator | None = None)
Bases: `object`

This class constructs the coherent interfaces between two crystalline slabs
Coherency is defined by matching lattices not sub-planes.


* **Parameters**


    * **substrate_structure** – structure of substrate


    * **film_structure** – structure of film


    * **film_miller** – miller index of the film layer


    * **substrate_miller** – miller index for the substrate layer


    * **zslgen** – BiDirectionalZSL if you want custom lattice matching tolerances for coherency.



#### _find_matches()
Finds and stores the ZSL matches.


#### _find_terminations()
Finds all terminations.


#### get_interfaces(termination: tuple[str, str], gap: float = 2.0, vacuum_over_film: float = 20.0, film_thickness: float = 1, substrate_thickness: float = 1, in_layers: bool = True)
Generates interface structures given the film and substrate structure
as well as the desired terminations.


* **Parameters**


    * **termination** (*tuple**[**str**, **str**]*) – termination from self.termination list


    * **gap** (*float**, **optional*) – gap between film and substrate. Defaults to 2.0.


    * **vacuum_over_film** (*float**, **optional*) – vacuum over the top of the film. Defaults to 20.0.


    * **film_thickness** (*float**, **optional*) – the film thickness. Defaults to 1.


    * **substrate_thickness** (*float**, **optional*) – substrate thickness. Defaults to 1.


    * **in_layers** (*bool**, **optional*) – set the thickness in layer units. Defaults to True.



* **Yields**

    *Iterator[Interface]* – interfaces from slabs



### from_2d_to_3d(mat: ndarray)
Converts a 2D matrix to a 3D matrix.


### get_2d_transform(start: Sequence, end: Sequence)
Gets a 2d transformation matrix
that converts start to end.


### get_rot_3d_for_2d(film_matrix, sub_matrix)
Find transformation matrix that will rotate and strain the film to the substrate while preserving the c-axis.

## pymatgen.analysis.interfaces.substrate_analyzer module

This module provides classes to identify optimal substrates for film growth.


### _class_ SubstrateAnalyzer(film_max_miller=1, substrate_max_miller=1, \*\*kwargs)
Bases: `ZSLGenerator`

This class applies a set of search criteria to identify suitable
substrates for film growth. It first uses a topological search by Zur
and McGill to identify matching super-lattices on various faces of the
two materials. Additional criteria can then be used to identify the most
suitable substrate. Currently, the only additional criteria is the
elastic strain energy of the super-lattices.

Initializes the substrate analyzer


* **Parameters**


    * **zslgen** (*ZSLGenerator*) – Defaults to a ZSLGenerator with standard
    tolerances, but can be fed one with custom tolerances


    * **film_max_miller** (*int*) – maximum miller index to generate for film
    surfaces


    * **substrate_max_miller** (*int*) – maximum miller index to generate for
    substrate surfaces.



#### calculate(film, substrate, elasticity_tensor=None, film_millers=None, substrate_millers=None, ground_state_energy=0, lowest=False)
Finds all topological matches for the substrate and calculates elastic
strain energy and total energy for the film if elasticity tensor and
ground state energy are provided:


* **Parameters**


    * **film** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – conventional standard structure for the film


    * **substrate** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – conventional standard structure for the
    substrate


    * **elasticity_tensor** ([*ElasticTensor*](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensor)) – elasticity tensor for the film
    in the IEEE orientation


    * **film_millers** (*array*) – film facets to consider in search as defined by
    miller indices


    * **substrate_millers** (*array*) – substrate facets to consider in search as
    defined by miller indices


    * **ground_state_energy** (*float*) – ground state energy for the film


    * **lowest** (*bool*) – only consider lowest matching area for each surface



#### generate_surface_vectors(film_millers, substrate_millers)
Generates the film/substrate slab combinations for a set of given
miller indices.


* **Parameters**


    * **film_millers** (*array*) – all miller indices to generate slabs for
    film


    * **substrate_millers** (*array*) – all miller indices to generate slabs
    for substrate



### _class_ SubstrateMatch(film_sl_vectors: list, substrate_sl_vectors: list, film_vectors: list, substrate_vectors: list, film_transformation: list, substrate_transformation: list, film_miller: tuple[int, int, int], substrate_miller: tuple[int, int, int], strain: [Strain](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.strain.Strain), von_mises_strain: float, ground_state_energy: float, elastic_energy: float)
Bases: `ZSLMatch`

A substrate match building on the Zur and McGill algorithm. This match class includes the miller
planes of the film and substrate the full strain tensor, the Von Mises strain, the ground state
energy if provided, and the elastic energy.


#### elastic_energy(_: floa_ )

#### film_miller(_: tuple[int, int, int_ )

#### _classmethod_ from_zsl(match: ZSLMatch, film: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), film_miller, substrate_miller, elasticity_tensor=None, ground_state_energy=0)
Generate a substrate match from a ZSL match plus metadata.


#### ground_state_energy(_: floa_ )

#### strain(_: [Strain](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.strain.Strain_ )

#### substrate_miller(_: tuple[int, int, int_ )

#### _property_ total_energy()
Total energy of this match.


#### von_mises_strain(_: floa_ )
## pymatgen.analysis.interfaces.zsl module

This module implements the Zur and McGill lattice matching algorithm.


### _class_ ZSLGenerator(max_area_ratio_tol=0.09, max_area=400, max_length_tol=0.03, max_angle_tol=0.01, bidirectional=False)
Bases: `MSONable`

This class generate matching interface super lattices based on the methodology
of lattice vector matching for heterostructural interfaces proposed by
Zur and McGill:
Journal of Applied Physics 55 (1984), 378 ; doi: 10.1063/1.333084
The process of generating all possible matching super lattices is:
1.) Reduce the surface lattice vectors and calculate area for the surfaces
2.) Generate all super lattice transformations within a maximum allowed area

> limit that give nearly equal area super-lattices for the two
> surfaces - generate_sl_transformation_sets

3.) For each superlattice set:

    1.) Reduce super lattice vectors
    2.) Check length and angle between film and substrate super lattice

    > vectors to determine if the super lattices are the nearly same
    > and therefore coincident - get_equiv_transformations.

Initialize a Zur Super Lattice Generator for a specific film and

    substrate


* **Parameters**


    * **max_area_ratio_tol** (*float*) – Max tolerance on ratio of
    super-lattices to consider equal


    * **max_area** (*float*) – max super lattice area to generate in search


    * **max_length_tol** – maximum length tolerance in checking if two
    vectors are of nearly the same length


    * **max_angle_tol** – maximum angle tolerance in checking of two sets
    of vectors have nearly the same angle between them.



#### generate_sl_transformation_sets(film_area, substrate_area)
Generates transformation sets for film/substrate pair given the
area of the unit cell area for the film and substrate. The
transformation sets map the film and substrate unit cells to super
lattices with a maximum area


* **Parameters**


    * **film_area** (*int*) – the unit cell area for the film


    * **substrate_area** (*int*) – the unit cell area for the substrate



* **Returns**

    a set of transformation_sets defined as:

        1.) the transformation matrices for the film to create a
        super lattice of area i\*film area
        2.) the transformation matrices for the substrate to create
        a super lattice of area j\*film area.




* **Return type**

    transformation_sets



#### get_equiv_transformations(transformation_sets, film_vectors, substrate_vectors)
Applies the transformation_sets to the film and substrate vectors
to generate super-lattices and checks if they matches.
Returns all matching vectors sets.


* **Parameters**


    * **transformation_sets** (*array*) – an array of transformation sets:
    each transformation set is an array with the (i,j)
    indicating the area multiples of the film and substrate it
    corresponds to, an array with all possible transformations
    for the film area multiple i and another array for the
    substrate area multiple j.


    * **film_vectors** (*array*) – film vectors to generate super lattices


    * **substrate_vectors** (*array*) – substrate vectors to generate super
    lattices



### _class_ ZSLMatch(film_sl_vectors: list, substrate_sl_vectors: list, film_vectors: list, substrate_vectors: list, film_transformation: list, substrate_transformation: list)
Bases: `MSONable`

A match from the Zur and McGill Algorithm. The super_lattice vectors are listed
as _sl_vectors. These are reduced according to the algorithm in the paper which
effectively a rotation in 3D space. Use the match_transformation property to get
the appropriate transformation matrix.


#### film_sl_vectors(_: lis_ )

#### film_transformation(_: lis_ )

#### film_vectors(_: lis_ )

#### _property_ match_area()
The area of the match between the substrate and film super lattice vectors.


#### _property_ match_transformation()
The transformation matrix to convert the film super lattice vectors to the substrate.


#### substrate_sl_vectors(_: lis_ )

#### substrate_transformation(_: lis_ )

#### substrate_vectors(_: lis_ )

### _bidirectional_same_vectors(vec_set1, vec_set2, max_length_tol, max_angle_tol)
Bidirectional version of above matching constraint check.


### _unidirectional_is_same_vectors(vec_set1, vec_set2, max_length_tol, max_angle_tol)
Determine if two sets of vectors are the same within length and angle
tolerances
:param vec_set1: an array of two vectors
:type vec_set1: array[array]
:param vec_set2: second array of two vectors.
:type vec_set2: array[array]


### fast_norm(a)
Much faster variant of numpy linalg norm.

Note that if numba is installed, this cannot be provided a list of ints;
please ensure input a is an np.array of floats.


### gen_sl_transform_matrices(area_multiple)
Generates the transformation matrices that convert a set of 2D
vectors into a super lattice of integer area multiple as proven
in Cassels:

Cassels, John William Scott. An introduction to the geometry of
numbers. Springer Science & Business Media, 2012.


* **Parameters**


    * **area_multiple** (*int*) – integer multiple of unit cell area for super


    * **area** (*lattice*) –



* **Returns**

    transformation matrices to convert unit vectors to
    super lattice vectors



* **Return type**

    matrix_list



### get_factors(n)
Generate all factors of n.


### is_same_vectors(vec_set1, vec_set2, bidirectional=False, max_length_tol=0.03, max_angle_tol=0.01)
Determine if two sets of vectors are the same within length and angle
tolerances
:param vec_set1: an array of two vectors
:type vec_set1: array[array]
:param vec_set2: second array of two vectors.
:type vec_set2: array[array]


### reduce_vectors(a, b)
Generate independent and unique basis vectors based on the
methodology of Zur and McGill.


### rel_angle(vec_set1, vec_set2)
Calculate the relative angle between two vector sets.


* **Parameters**


    * **vec_set1** (*array**[**array**]*) – an array of two vectors


    * **vec_set2** (*array**[**array**]*) – second array of two vectors



### rel_strain(vec1, vec2)
Calculate relative strain between two vectors.


### vec_angle(a, b)
Calculate angle between two vectors.


### vec_area(a, b)
Area of lattice plane defined by two vectors.