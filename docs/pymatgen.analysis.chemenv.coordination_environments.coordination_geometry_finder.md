---
layout: default
title: pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md
nav_exclude: true
---

# pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder module

This module contains the main object used to identify the coordination environments in a given structure.
If you use this module, please cite:
David Waroquiers, Xavier Gonze, Gian-Marco Rignanese, Cathrin Welker-Nieuwoudt, Frank Rosowski,
Michael Goebel, Stephan Schenk, Peter Degelmann, Rute Andre, Robert Glaum, and Geoffroy Hautier,
“Statistical analysis of coordination environments in oxides”,
Chem. Mater., 2017, 29 (19), pp 8346-8360,
DOI: 10.1021/acs.chemmater.7b02766
D. Waroquiers, J. George, M. Horton, S. Schenk, K. A. Persson, G.-M. Rignanese, X. Gonze, G. Hautier
“ChemEnv: a fast and robust coordination environment identification tool”,
Acta Cryst. B 2020, 76, pp 683-695,
DOI: 10.1107/S2052520620007994.


### _class_ pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.AbstractGeometry(central_site=None, bare_coords=None, centering_type='standard', include_central_site_in_centroid=False, optimization=None)
Bases: `object`

Class used to describe a geometry (perfect or distorted).

Constructor for the abstract geometry
:param central_site: Coordinates of the central site
:param bare_coords: Coordinates of the neighbors of the central site
:param centering_type: How to center the abstract geometry
:param include_central_site_in_centroid: When the centering is on the centroid, the central site is included

> if this parameter is set to True.


* **Raise**

    ValueError if the parameters are not consistent.



#### _property_ cn()
Coordination number


* **Type**

    return



#### _property_ coordination_number()
Coordination number


* **Type**

    return



#### _classmethod_ from_cg(cg, centering_type='standard', include_central_site_in_centroid=False)

* **Parameters**


    * **cg** –


    * **centering_type** –


    * **include_central_site_in_centroid** –



* **Returns**




#### points_wcs_csc(permutation=None)

* **Parameters**

    **permutation** –



* **Returns**




#### points_wcs_ctwcc(permutation=None)

* **Parameters**

    **permutation** –



* **Returns**




#### points_wcs_ctwocc(permutation=None)

* **Parameters**

    **permutation** –



* **Returns**




#### points_wocs_csc(permutation=None)

* **Parameters**

    **permutation** –



* **Returns**




#### points_wocs_ctwcc(permutation=None)

* **Parameters**

    **permutation** –



* **Returns**




#### points_wocs_ctwocc(permutation=None)

* **Parameters**

    **permutation** –



* **Returns**




### _class_ pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder(permutations_safe_override: bool = False, plane_ordering_override: bool = True, plane_safe_permutations: bool = False, only_symbols=None)
Bases: `object`

Main class used to find the local environments in a structure.


* **Parameters**


    * **permutations_safe_override** – If set to True, all permutations are tested (very time-consuming for large


    * **numbers!****)** (*coordination*) –


    * **plane_ordering_override** – If set to False, the ordering of the points in the plane is disabled


    * **plane_safe_permutations** – Whether to use safe permutations.


    * **only_symbols** – Whether to restrict the list of environments to be identified.



#### BVA_DISTANCE_SCALE_FACTORS(_ = {'GGA_relaxed': 1.015, 'LDA_relaxed': 0.995, 'experimental': 1.0_ )

#### DEFAULT_BVA_DISTANCE_SCALE_FACTOR(_ = 1._ )

#### DEFAULT_SPG_ANALYZER_OPTIONS(_ = {'angle_tolerance': 5, 'symprec': 0.001_ )

#### DEFAULT_STRATEGY(_ = <pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.MultiWeightsChemenvStrategy object_ )

#### PRESETS(_ = {'DEFAULT': {'maximum_distance_factor': 2.0, 'minimum_angle_factor': 0.05, 'optimization': 2, 'voronoi_normalized_angle_tolerance': 0.03, 'voronoi_normalized_distance_tolerance': 0.05}_ )

#### STRUCTURE_REFINEMENT_NONE(_ = 'none_ )

#### STRUCTURE_REFINEMENT_REFINED(_ = 'refined_ )

#### STRUCTURE_REFINEMENT_SYMMETRIZED(_ = 'symmetrized_ )

#### compute_coordination_environments(structure, indices=None, only_cations=True, strategy=<pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.MultiWeightsChemenvStrategy object>, valences='bond-valence-analysis', initial_structure_environments=None)

* **Parameters**


    * **structure** –


    * **indices** –


    * **only_cations** –


    * **strategy** –


    * **valences** –


    * **initial_structure_environments** –



* **Returns**




#### compute_structure_environments(excluded_atoms=None, only_atoms=None, only_cations=True, only_indices=None, maximum_distance_factor=2.0, minimum_angle_factor=0.05, max_cn=None, min_cn=None, only_symbols=None, valences='undefined', additional_conditions=None, info=None, timelimit=None, initial_structure_environments=None, get_from_hints=False, voronoi_normalized_distance_tolerance=0.05, voronoi_normalized_angle_tolerance=0.03, voronoi_distance_cutoff=None, recompute=None, optimization=2)
Computes and returns the StructureEnvironments object containing all the information about the coordination
environments in the structure
:param excluded_atoms: Atoms for which the coordination geometries does not have to be identified
:param only_atoms: If not set to None, atoms for which the coordination geometries have to be identified
:param only_cations: If set to True, will only compute environments for cations
:param only_indices: If not set to None, will only compute environments the atoms of the given indices
:param maximum_distance_factor: If not set to None, neighbors beyond

> maximum_distance_factor\*closest_neighbor_distance are not considered


* **Parameters**


    * **minimum_angle_factor** – If not set to None, neighbors for which the angle is lower than
    minimum_angle_factor\*largest_angle_neighbor are not considered


    * **max_cn** – maximum coordination number to be considered


    * **min_cn** – minimum coordination number to be considered


    * **only_symbols** – if not set to None, consider only coordination environments with the given symbols


    * **valences** – valences of the atoms


    * **additional_conditions** – additional conditions to be considered in the bonds (example : only bonds
    between cation and anion


    * **info** – additional info about the calculation


    * **timelimit** – time limit (in secs) after which the calculation of the StructureEnvironments object stops


    * **initial_structure_environments** – initial StructureEnvironments object (most probably incomplete)


    * **get_from_hints** – whether to add neighbors sets from “hints” (e.g. capped environment => test the
    neighbors without the cap)


    * **voronoi_normalized_distance_tolerance** – tolerance for the normalized distance used to distinguish
    neighbors sets


    * **voronoi_normalized_angle_tolerance** – tolerance for the normalized angle used to distinguish
    neighbors sets


    * **voronoi_distance_cutoff** – determines distance of considered neighbors. Especially important to increase it
    for molecules in a box.


    * **recompute** – whether to recompute the sites already computed (when initial_structure_environments
    is not None)


    * **optimization** – optimization algorithm



* **Returns**

    The StructureEnvironments object containing all the information about the coordination
    environments in the structure.



#### coordination_geometry_symmetry_measures(coordination_geometry, tested_permutations=False, points_perfect=None, optimization=None)
Returns the symmetry measures of a given coordination_geometry for a set of permutations depending on
the permutation setup. Depending on the parameters of the LocalGeometryFinder and on the coordination

> geometry, different methods are called.


* **Parameters**

    **coordination_geometry** – Coordination geometry for which the symmetry measures are looked for



* **Returns**

    the symmetry measures of a given coordination_geometry for a set of permutations



* **Raise**

    NotImplementedError if the permutation_setup does not exists.



#### coordination_geometry_symmetry_measures_fallback_random(coordination_geometry, NRANDOM=10, points_perfect=None)
Returns the symmetry measures for a random set of permutations for the coordination geometry
“coordination_geometry”. Fallback implementation for the plane separation algorithms measures
of each permutation
:param coordination_geometry: The coordination geometry to be investigated
:param NRANDOM: Number of random permutations to be tested
:return: The symmetry measures for the given coordination geometry for each permutation investigated.


#### coordination_geometry_symmetry_measures_separation_plane(coordination_geometry, separation_plane_algo, testing=False, tested_permutations=False, points_perfect=None)
Returns the symmetry measures of the given coordination geometry “coordination_geometry” using separation
facets to reduce the complexity of the system. Caller to the refined 2POINTS, 3POINTS and other …
:param coordination_geometry: The coordination geometry to be investigated
:return: The symmetry measures for the given coordination geometry for each plane and permutation investigated.


#### coordination_geometry_symmetry_measures_separation_plane_optim(coordination_geometry, separation_plane_algo, points_perfect=None, nb_set=None, optimization=None)
Returns the symmetry measures of the given coordination geometry “coordination_geometry” using separation
facets to reduce the complexity of the system. Caller to the refined 2POINTS, 3POINTS and other …


* **Parameters**


    * **coordination_geometry** – The coordination geometry to be investigated.


    * **separation_plane_algo** – Separation Plane algorithm used.


    * **points_perfect** – Points corresponding to the perfect geometry.


    * **nb_set** – Neighbor set for this set of points. (used to store already computed separation planes)


    * **optimization** – Optimization level (1 or 2).



* **Returns**

    Continuous symmetry measures for the given coordination geometry for each plane and permutation

        investigated, corresponding permutations, corresponding algorithms,
        corresponding mappings from local to perfect environment and corresponding mappings
        from perfect to local environment.




* **Return type**

    tuple



#### coordination_geometry_symmetry_measures_sepplane_optim(coordination_geometry, points_perfect=None, nb_set=None, optimization=None)
Returns the symmetry measures of a given coordination_geometry for a set of permutations depending on
the permutation setup. Depending on the parameters of the LocalGeometryFinder and on the coordination

> geometry, different methods are called.


* **Parameters**

    **coordination_geometry** – Coordination geometry for which the symmetry measures are looked for



* **Returns**

    the symmetry measures of a given coordination_geometry for a set of permutations



* **Raise**

    NotImplementedError if the permutation_setup does not exists.



#### coordination_geometry_symmetry_measures_standard(coordination_geometry, algo, points_perfect=None, optimization=None)
Returns the symmetry measures for a set of permutations (whose setup depends on the coordination geometry)
for the coordination geometry “coordination_geometry”. Standard implementation looking for the symmetry
measures of each permutation
:param coordination_geometry: The coordination geometry to be investigated
:return: The symmetry measures for the given coordination geometry for each permutation investigated.


#### get_coordination_symmetry_measures(only_minimum=True, all_csms=True, optimization=None)
Returns the continuous symmetry measures of the current local geometry in a dictionary.
:return: the continuous symmetry measures of the current local geometry in a dictionary.


#### get_coordination_symmetry_measures_optim(only_minimum=True, all_csms=True, nb_set=None, optimization=None)
Returns the continuous symmetry measures of the current local geometry in a dictionary.
:return: the continuous symmetry measures of the current local geometry in a dictionary.


#### get_structure()
Returns the pymatgen Structure that has been setup for the identification of geometries (the initial one
might have been refined/symmetrized using the SpaceGroupAnalyzer).
:return: The pymatgen Structure that has been setup for the identification of geometries (the initial one
might have been refined/symmetrized using the SpaceGroupAnalyzer).


#### set_structure(lattice: [Lattice](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice), species, coords, coords_are_cartesian)
Sets up the pymatgen structure for which the coordination geometries have to be identified starting from the
lattice, the species and the coordinates
:param lattice: The lattice of the structure
:param species: The species on the sites
:param coords: The coordinates of the sites
:param coords_are_cartesian: If set to True, the coordinates are given in Cartesian coordinates.


#### setup_explicit_indices_local_geometry(explicit_indices)
Sets up explicit indices for the local geometry, for testing purposes
:param explicit_indices: explicit indices for the neighbors (set of numbers
from 0 to CN-1 in a given order).


#### setup_local_geometry(isite, coords, optimization=None)
Sets up the AbstractGeometry for the local geometry of site with index isite.
:param isite: Index of the site for which the local geometry has to be set up
:param coords: The coordinates of the (local) neighbors.


#### setup_ordered_indices_local_geometry(coordination)
Sets up ordered indices for the local geometry, for testing purposes
:param coordination: coordination of the local geometry.


#### setup_parameter(parameter, value)
Setup of one specific parameter to the given value. The other parameters are unchanged. See setup_parameters
method for the list of possible parameters
:param parameter: Parameter to setup/update
:param value: Value of the parameter.


#### setup_parameters(centering_type='standard', include_central_site_in_centroid=False, bva_distance_scale_factor=None, structure_refinement='refined', spg_analyzer_options=None)
Setup of the parameters for the coordination geometry finder. A reference point for the geometries has to be
chosen. This can be the centroid of the structure (including or excluding the atom for which the coordination
geometry is looked for) or the atom itself. In the ‘standard’ centering_type, the reference point is the central
atom for coordination numbers 1, 2, 3 and 4 and the centroid for coordination numbers > 4.
:param centering_type: Type of the reference point (centering) ‘standard’, ‘centroid’ or ‘central_site’
:param include_central_site_in_centroid: In case centering_type is ‘centroid’, the central site is included if

> this value is set to True.


* **Parameters**


    * **bva_distance_scale_factor** – Scaling factor for the bond valence analyzer (this might be different whether
    the structure is an experimental one, an LDA or a GGA relaxed one, or any other relaxation scheme (where
    under- or over-estimation of bond lengths is known).


    * **structure_refinement** – Refinement of the structure. Can be “none”, “refined” or “symmetrized”.


    * **spg_analyzer_options** – Options for the SpaceGroupAnalyzer (dictionary specifying “symprec”
    and “angle_tolerance”. See pymatgen’s SpaceGroupAnalyzer for more information.



#### setup_random_indices_local_geometry(coordination)
Sets up random indices for the local geometry, for testing purposes
:param coordination: coordination of the local geometry.


#### setup_random_structure(coordination)
Sets up a purely random structure with a given coordination.
:param coordination: coordination number for the random structure.


#### setup_structure(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Sets up the structure for which the coordination geometries have to be identified. The structure is analyzed
with the space group analyzer and a refined structure is used
:param structure: A pymatgen Structure.


#### setup_test_perfect_environment(symbol, randomness=False, max_random_dist=0.1, symbol_type='mp_symbol', indices='RANDOM', random_translation='NONE', random_rotation='NONE', random_scale='NONE', points=None)

* **Parameters**


    * **symbol** –


    * **randomness** –


    * **max_random_dist** –


    * **symbol_type** –


    * **indices** –


    * **random_translation** –


    * **random_rotation** –


    * **random_scale** –


    * **points** –



* **Returns**




#### update_nb_set_environments(se, isite, cn, inb_set, nb_set, recompute=False, optimization=None)

* **Parameters**


    * **se** –


    * **isite** –


    * **cn** –


    * **inb_set** –


    * **nb_set** –


    * **recompute** –


    * **optimization** –



* **Returns**




### pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.find_rotation(points_distorted, points_perfect)
This finds the rotation matrix that aligns the (distorted) set of points “points_distorted” with respect to the
(perfect) set of points “points_perfect” in a least-square sense.
:param points_distorted: List of points describing a given (distorted) polyhedron for which the rotation that

> aligns these points in a least-square sense to the set of perfect points “points_perfect”


* **Parameters**

    **points_perfect** – List of “perfect” points describing a given model polyhedron.



* **Returns**

    The rotation matrix.



### pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.find_scaling_factor(points_distorted, points_perfect, rot)
This finds the scaling factor between the (distorted) set of points “points_distorted” and the
(perfect) set of points “points_perfect” in a least-square sense.
:param points_distorted: List of points describing a given (distorted) polyhedron for which the scaling factor has

> to be obtained.


* **Parameters**


    * **points_perfect** – List of “perfect” points describing a given model polyhedron.


    * **rot** – The rotation matrix



* **Returns**

    The scaling factor between the two structures and the rotated set of (distorted) points.



### pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.symmetry_measure(points_distorted, points_perfect)
Computes the continuous symmetry measure of the (distorted) set of points “points_distorted” with respect to the
(perfect) set of points “points_perfect”.
:param points_distorted: List of points describing a given (distorted) polyhedron for which the symmetry measure

> has to be computed with respect to the model polyhedron described by the list of points
> “points_perfect”.


* **Parameters**

    **points_perfect** – List of “perfect” points describing a given model polyhedron.



* **Returns**

    The continuous symmetry measure of the distorted polyhedron with respect to the perfect polyhedron.