---
layout: default
title: pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md
nav_exclude: true
---

# pymatgen.analysis.chemenv.coordination_environments.coordination_geometries module

This module contains the class describing the coordination geometries that can exist in a given structure. These
“model” coordination geometries are described in the following articles :

>
> * Pure Appl. Chem., Vol. 79, No. 10, pp. 1779–1799, 2007.


> * Acta Cryst. A, Vol. 46, No. 1, pp. 1–11, 1990.

The module also contains descriptors of part of these geometries (plane of separation, …) that are used in the
identification algorithms.


### _class_ pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.AbstractChemenvAlgorithm(algorithm_type)
Bases: `MSONable`

Base class used to define a Chemenv algorithm used to identify the correct permutation for the computation
of the Continuous Symmetry Measure.

Base constructor for ChemenvAlgorithm.


* **Parameters**

    **algorithm_type** (*str*) – Type of algorithm.



#### _property_ algorithm_type()
Return the type of algorithm.

Returns: Type of the algorithm


#### _abstract_ as_dict()
A JSON-serializable dict representation of the algorithm.


### _class_ pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.AllCoordinationGeometries(permutations_safe_override=False, only_symbols=None)
Bases: `dict`

Class used to store all the reference “coordination geometries” (list with instances of the CoordinationGeometry
classes).

Initializes the list of Coordination Geometries.


* **Parameters**


    * **permutations_safe_override** – Whether to use safe permutations.


    * **only_symbols** – Whether to restrict the list of environments to be identified.



#### get_geometries(coordination=None, returned='cg')
Returns a list of coordination geometries with the given coordination number.


* **Parameters**


    * **coordination** – The coordination number of which the list of coordination geometries are returned.


    * **returned** – Type of objects in the list.



#### get_geometry_from_IUCr_symbol(IUCr_symbol)
Returns the coordination geometry of the given IUCr symbol.


* **Parameters**

    **IUCr_symbol** – The IUCr symbol of the coordination geometry.



#### get_geometry_from_IUPAC_symbol(IUPAC_symbol)
Returns the coordination geometry of the given IUPAC symbol.


* **Parameters**

    **IUPAC_symbol** – The IUPAC symbol of the coordination geometry.



#### get_geometry_from_mp_symbol(mp_symbol)
Returns the coordination geometry of the given mp_symbol.


* **Parameters**

    **mp_symbol** – The mp_symbol of the coordination geometry.



#### get_geometry_from_name(name)
Returns the coordination geometry of the given name.


* **Parameters**

    **name** – The name of the coordination geometry.



#### get_implemented_geometries(coordination=None, returned='cg', include_deactivated=False)
Returns a list of the implemented coordination geometries with the given coordination number.


* **Parameters**


    * **coordination** – The coordination number of which the list of implemented coordination geometries
    are returned.


    * **returned** – Type of objects in the list.


    * **include_deactivated** – Whether to include CoordinationGeometry that are deactivated.



#### get_not_implemented_geometries(coordination=None, returned='mp_symbol')
Returns a list of the implemented coordination geometries with the given coordination number.


* **Parameters**


    * **coordination** – The coordination number of which the list of implemented coordination geometries
    are returned.


    * **returned** – Type of objects in the list.



#### get_symbol_cn_mapping(coordination=None)
Return a dictionary mapping the symbol of a CoordinationGeometry to its coordination.


* **Parameters**

    **coordination** – Whether to restrict the dictionary to a given coordination.


Returns: Dictionary mapping the symbol of a CoordinationGeometry to its coordination.


#### get_symbol_name_mapping(coordination=None)
Return a dictionary mapping the symbol of a CoordinationGeometry to its name.


* **Parameters**

    **coordination** – Whether to restrict the dictionary to a given coordination.


Returns: Dictionary mapping the symbol of a CoordinationGeometry to its name.


#### is_a_valid_coordination_geometry(mp_symbol=None, IUPAC_symbol=None, IUCr_symbol=None, name=None, cn=None)
Checks whether a given coordination geometry is valid (exists) and whether the parameters are coherent with
each other.


* **Parameters**


    * **mp_symbol** – The mp_symbol of the coordination geometry.


    * **IUPAC_symbol** – The IUPAC_symbol of the coordination geometry.


    * **IUCr_symbol** – The IUCr_symbol of the coordination geometry.


    * **name** – The name of the coordination geometry.


    * **cn** – The coordination of the coordination geometry.



#### pretty_print(type='implemented_geometries', maxcn=8, additional_info=None)
Return a string with a list of the Coordination Geometries.


* **Parameters**


    * **type** – Type of string to be returned (all_geometries, all_geometries_latex_images, all_geometries_latex,
    implemented_geometries).


    * **maxcn** – Maximum coordination.


    * **additional_info** – Whether to add some additional info for each coordination geometry.


Returns: String describing the list of coordination geometries.


### _class_ pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry(mp_symbol, name, alternative_names=None, IUPAC_symbol=None, IUCr_symbol=None, coordination=None, central_site=None, points=None, solid_angles=None, permutations_safe_override=False, deactivate=False, faces=None, edges=None, algorithms=None, equivalent_indices=None, neighbors_sets_hints=None)
Bases: `object`

Class used to store the ideal representation of a chemical environment or “coordination geometry”.

Initializes one “coordination geometry” according to [Pure Appl. Chem., Vol. 79, No. 10, pp. 1779–1799, 2007]
and [Acta Cryst. A, Vol. 46, No. 1, pp. 1–11, 1990].


* **Parameters**


    * **mp_symbol** – Symbol used internally for the coordination geometry.


    * **name** – Name of the coordination geometry.


    * **alternative_names** – Alternative names for this coordination geometry.


    * **IUPAC_symbol** – The IUPAC symbol of this coordination geometry.


    * **IUCr_symbol** – The IUCr symbol of this coordination geometry.


    * **coordination** – The coordination number of this coordination geometry (number of neighboring atoms).


    * **central_site** – The coordinates of the central site of this coordination geometry.


    * **points** – The list of the coordinates of all the points of this coordination geometry.


    * **solid_angles** – The list of solid angles for each neighbor in this coordination geometry.


    * **permutations_safe_override** – Computes all the permutations if set to True (overrides the plane separation
    algorithms or any other algorithm, for testing purposes)


    * **deactivate** – Whether to deactivate this coordination geometry


    * **faces** – List of the faces with their vertices given in a clockwise or anticlockwise order, for drawing
    purposes.


    * **edges** – List of edges, for drawing purposes.


    * **algorithms** – Algorithms used to identify this coordination geometry.


    * **equivalent_indices** – The equivalent sets of indices in this coordination geometry (can be used to skip
    equivalent permutations that have already been performed).


    * **neighbors_sets_hints** – Neighbors sets hints for this coordination geometry.



#### CSM_SKIP_SEPARATION_PLANE_ALGO(_ = 10._ )

#### _property_ IUCr_symbol()
Returns the IUCr symbol of this coordination geometry.


#### _property_ IUCr_symbol_str()
Returns a string representation of the IUCr symbol of this coordination geometry.


#### _property_ IUPAC_symbol()
Returns the IUPAC symbol of this coordination geometry.


#### _property_ IUPAC_symbol_str()
Returns a string representation of the IUPAC symbol of this coordination geometry.


#### _class_ NeighborsSetsHints(hints_type, options)
Bases: `object`

Class used to describe neighbors sets hints.

This allows to possibly get a lower coordination from a capped-like model polyhedron.

Constructor for this NeighborsSetsHints.


* **Parameters**


    * **hints_type** – type of hint (single, double or triple cap)


    * **options** – options for the “hinting”, e.g. the maximum csm value beyond which no additional
    neighbors set could be found from a “cap hint”.



#### ALLOWED_HINTS_TYPES(_ = ('single_cap', 'double_cap', 'triple_cap'_ )

#### as_dict()
A JSON-serializable dict representation of this NeighborsSetsHints.


#### double_cap_hints(hints_info)
Return hints for an additional neighbors set, i.e. the voronoi indices that constitute this new
neighbors set, in case of a “Double cap” hint.


* **Parameters**

    **hints_info** – Info needed to build new “hinted” neighbors set.


Returns: Voronoi indices of the new “hinted” neighbors set.


#### _classmethod_ from_dict(dd)
Reconstructs the NeighborsSetsHints from its JSON-serializable dict representation.


* **Parameters**

    **dd** – a JSON-serializable dict representation of a NeighborsSetsHints.


Returns: a NeighborsSetsHints.


#### hints(hints_info)
Return hints for an additional neighbors set, i.e. the voronoi indices that constitute this new
neighbors set.


* **Parameters**

    **hints_info** – Info needed to build new “hinted” neighbors set.


Returns: Voronoi indices of the new “hinted” neighbors set.


#### single_cap_hints(hints_info)
Return hints for an additional neighbors set, i.e. the voronoi indices that constitute this new
neighbors set, in case of a “Single cap” hint.


* **Parameters**

    **hints_info** – Info needed to build new “hinted” neighbors set.


Returns: Voronoi indices of the new “hinted” neighbors set.


#### triple_cap_hints(hints_info)
Return hints for an additional neighbors set, i.e. the voronoi indices that constitute this new
neighbors set, in case of a “Triple cap” hint.


* **Parameters**

    **hints_info** – Info needed to build new “hinted” neighbors set.


Returns: Voronoi indices of the new “hinted” neighbors set.


#### _property_ algorithms()
Returns the list of algorithms that are used to identify this coordination geometry.


#### as_dict()
A JSON-serializable dict representation of this CoordinationGeometry.


#### _property_ ce_symbol()
Returns the symbol of this coordination geometry.


#### _property_ coordination_number()
Returns the coordination number of this coordination geometry.


#### _property_ distfactor_max()
The maximum distfactor for the perfect CoordinationGeometry.

Returns: Maximum distfactor for the perfect CoordinationGeometry (usually 1.0 for symmetric polyhedrons).


#### edges(sites, permutation=None, input='sites')
Returns the list of edges of this coordination geometry. Each edge is given as a
list of its end vertices coordinates.


#### faces(sites, permutation=None)
Returns the list of faces of this coordination geometry. Each face is given as a
list of its vertices coordinates.


#### _classmethod_ from_dict(dct)
Reconstructs the CoordinationGeometry from its JSON-serializable dict representation.


* **Parameters**

    **dct** – a JSON-serializable dict representation of a CoordinationGeometry.


Returns: a CoordinationGeometry.


#### get_central_site()
Returns the central site of this coordination geometry.


#### get_coordination_number()
Returns the coordination number of this coordination geometry.


#### get_name()
Returns the name of this coordination geometry.


#### get_pmeshes(sites, permutation=None)
Returns the pmesh strings used for jmol to show this geometry.


#### is_implemented()
Returns True if this coordination geometry is implemented.


#### _property_ mp_symbol()
Returns the MP symbol of this coordination geometry.


#### _property_ number_of_permutations()
Returns the number of permutations of this coordination geometry.


#### _property_ pauling_stability_ratio()
Returns the theoretical Pauling stability ratio (rC/rA) for this environment.


#### ref_permutation(permutation)
Returns the reference permutation for a set of equivalent permutations.

Can be useful to skip permutations that have already been performed.


* **Parameters**

    **permutation** – Current permutation


Returns: Reference permutation of the perfect CoordinationGeometry.


#### set_permutations_safe_override(permutations_safe_override)
Setup ChemEnv so that a safe set of permutations are used.


* **Parameters**

    **permutations_safe_override** – Whether to use safe permutations.



#### solid_angles(permutation=None)
Returns the list of “perfect” solid angles Each edge is given as a
list of its end vertices coordinates.


### _class_ pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.ExplicitPermutationsAlgorithm(permutations)
Bases: `AbstractChemenvAlgorithm`

Class representing the algorithm doing the explicit permutations for the calculation of
the Continuous Symmetry Measure.

> Initializes a separation plane for a given perfect coordination geometry.


* **Parameters**

    **permutations** – Permutations used for this algorithm.



#### _property_ as_dict()
Return the JSON-serializable dict representation of this ExplicitPermutationsAlgorithm algorithm.

Returns: a JSON-serializable dict representation of this ExplicitPermutationsAlgorithm algorithm.


#### _classmethod_ from_dict(dd)
Reconstructs the ExplicitPermutationsAlgorithm algorithm from its JSON-serializable dict representation.


* **Parameters**

    **dd** – a JSON-serializable dict representation of an ExplicitPermutationsAlgorithm algorithm.


Returns: an ExplicitPermutationsAlgorithm algorithm.


#### _property_ permutations()
Return the permutations to be performed for this algorithm.

Returns: Permutations to be performed.


### _class_ pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.SeparationPlane(plane_points, mirror_plane=False, ordered_plane=False, point_groups=None, ordered_point_groups=None, explicit_permutations=None, minimum_number_of_points=None, explicit_optimized_permutations=None, multiplicity=None, other_plane_points=None)
Bases: `AbstractChemenvAlgorithm`

Class representing the algorithm using separation planes for the calculation of
the Continuous Symmetry Measure.

> Initializes a separation plane for a given perfect coordination geometry.


* **Parameters**


    * **plane_points** – Indices of the points that are in the plane in the perfect structure (and should be
    found in the defective one as well).


    * **mirror_plane** – True if the separation plane is a mirror plane, in which case there is a correspondence
    of the points in each point_group (can reduce the number of permutations).


    * **ordered_plane** – True if the order of the points in the plane can be taken into account to reduce the
    number of permutations.


    * **point_groups** – Indices of the points in the two groups of points separated by the plane.


    * **ordered_point_groups** – Whether the order of the points in each group of points can be taken into account to
    reduce the number of permutations.


    * **explicit_permutations** – Explicit permutations to be performed in this separation plane algorithm.


    * **minimum_number_of_points** – Minimum number of points needed to initialize a separation plane
    for this algorithm.


    * **explicit_optimized_permutations** – Optimized set of explicit permutations to be performed in this
    separation plane algorithm.


    * **multiplicity** – Number of such planes in the model geometry.


    * **other_plane_points** – Indices of the points that are in the plane in the perfect structure for the other
    planes. The multiplicity should be equal to the length of this list + 1 (“main” separation plane +
    the other ones).



#### _property_ argsorted_ref_separation_perm()
“Arg sorted” ordered indices of the separation plane.

This is used in the identification of the final permutation to be used.

Returns: list of the “arg sorted” ordered indices of the separation plane.


#### _property_ as_dict()
Return the JSON-serializable dict representation of this SeparationPlane algorithm.

Returns: a JSON-serializable dict representation of this SeparationPlane algorithm.


#### _classmethod_ from_dict(dct)
Reconstructs the SeparationPlane algorithm from its JSON-serializable dict representation.


* **Parameters**

    **dd** – a JSON-serializable dict representation of an SeparationPlane algorithm.


Returns: a SeparationPlane algorithm.


#### _property_ permutations()
Permutations used for this separation plane algorithm.

Returns: List of permutations to be performed.


#### _property_ ref_separation_perm()
Ordered indices of the separation plane.

### Examples

For a separation plane of type 2|4|3, with plane_points indices [0, 3, 5, 8] and
point_groups indices [1, 4] and [2, 7, 6], the list of ordered indices is :
[0, 3, 5, 8, 1, 4, 2, 7, 6].

Returns: list of ordered indices of this separation plane.


#### safe_separation_permutations(ordered_plane=False, ordered_point_groups=None, add_opposite=False)
Simple and safe permutations for this separation plane.

This is not meant to be used in production. Default configuration for ChemEnv does not use this method.


* **Parameters**


    * **ordered_plane** – Whether the order of the points in the plane can be used to reduce the
    number of permutations.


    * **ordered_point_groups** – Whether the order of the points in each point group can be used to reduce the
    number of permutations.


    * **add_opposite** – Whether to add the permutations from the second group before the first group as well.


Returns: List of safe permutations.