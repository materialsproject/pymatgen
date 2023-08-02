---
layout: default
title: pymatgen.analysis.interfaces.zsl.md
nav_exclude: true
---

# pymatgen.analysis.interfaces.zsl module

This module implements the Zur and McGill lattice matching algorithm.


### _class_ pymatgen.analysis.interfaces.zsl.ZSLGenerator(max_area_ratio_tol=0.09, max_area=400, max_length_tol=0.03, max_angle_tol=0.01, bidirectional=False)
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
:param film_area: the unit cell area for the film
:type film_area: int
:param substrate_area: the unit cell area for the substrate
:type substrate_area: int


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



### _class_ pymatgen.analysis.interfaces.zsl.ZSLMatch(film_sl_vectors: list, substrate_sl_vectors: list, film_vectors: list, substrate_vectors: list, film_transformation: list, substrate_transformation: list)
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

### pymatgen.analysis.interfaces.zsl.fast_norm(a)
Much faster variant of numpy linalg norm.

Note that if numba is installed, this cannot be provided a list of ints;
please ensure input a is an np.array of floats.


### pymatgen.analysis.interfaces.zsl.gen_sl_transform_matrices(area_multiple)
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



### pymatgen.analysis.interfaces.zsl.get_factors(n)
Generate all factors of n.


### pymatgen.analysis.interfaces.zsl.is_same_vectors(vec_set1, vec_set2, bidirectional=False, max_length_tol=0.03, max_angle_tol=0.01)
Determine if two sets of vectors are the same within length and angle
tolerances
:param vec_set1: an array of two vectors
:type vec_set1: array[array]
:param vec_set2: second array of two vectors.
:type vec_set2: array[array]


### pymatgen.analysis.interfaces.zsl.reduce_vectors(a, b)
Generate independent and unique basis vectors based on the
methodology of Zur and McGill.


### pymatgen.analysis.interfaces.zsl.rel_angle(vec_set1, vec_set2)
Calculate the relative angle between two vector sets.


* **Parameters**


    * **vec_set1** (*array**[**array**]*) – an array of two vectors


    * **vec_set2** (*array**[**array**]*) – second array of two vectors



### pymatgen.analysis.interfaces.zsl.rel_strain(vec1, vec2)
Calculate relative strain between two vectors.


### pymatgen.analysis.interfaces.zsl.vec_angle(a, b)
Calculate angle between two vectors.


### pymatgen.analysis.interfaces.zsl.vec_area(a, b)
Area of lattice plane defined by two vectors.