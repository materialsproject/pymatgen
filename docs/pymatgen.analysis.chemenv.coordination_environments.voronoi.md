---
layout: default
title: pymatgen.analysis.chemenv.coordination_environments.voronoi.md
nav_exclude: true
---

# pymatgen.analysis.chemenv.coordination_environments.voronoi module

This module contains the object used to describe the possible bonded atoms based on a Voronoi analysis.


### _class_ pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer(structure=None, voronoi_list2=None, voronoi_cutoff=10.0, isites=None, normalized_distance_tolerance=1e-05, normalized_angle_tolerance=0.001, additional_conditions=None, valences=None, maximum_distance_factor=None, minimum_angle_factor=None)
Bases: `MSONable`

Class used to store the full Voronoi of a given structure.

Constructor for the VoronoiContainer object. Either a structure is given, in which case the Voronoi is
computed, or the different components of the VoronoiContainer are given (used in the from_dict method).


* **Parameters**


    * **structure** – Structure for which the Voronoi is computed.


    * **voronoi_list2** – List of voronoi polyhedrons for each site.


    * **voronoi_cutoff** – cutoff used for the voronoi.


    * **isites** – indices of sites for which the Voronoi has to be computed.


    * **normalized_distance_tolerance** – Tolerance for two normalized distances to be considered equal.


    * **normalized_angle_tolerance** – Tolerance for two normalized angles to be considered equal.


    * **additional_conditions** – Additional conditions to be used.


    * **valences** – Valences of all the sites in the structure (used when additional conditions require it).


    * **maximum_distance_factor** – The maximum distance factor to be considered.


    * **minimum_angle_factor** – The minimum angle factor to be considered.



* **Raises**

    **RuntimeError if the Voronoi cannot be constructed.** –



#### AC(_ = <pymatgen.analysis.chemenv.utils.defs_utils.AdditionalConditions object_ )

#### as_dict()
Bson-serializable dict representation of the VoronoiContainer.


* **Returns**

    dictionary that is BSON-encodable.



#### default_normalized_angle_tolerance(_ = 0.00_ )

#### default_normalized_distance_tolerance(_ = 1e-0_ )

#### default_voronoi_cutoff(_ = 10._ )

#### _classmethod_ from_dict(dct)
Reconstructs the VoronoiContainer object from a dict representation of the VoronoiContainer created using
the as_dict method.


* **Parameters**

    **dct** – dict representation of the VoronoiContainer object.



* **Returns**

    VoronoiContainer object.



#### get_rdf_figure(isite, normalized=True, figsize=None, step_function=None)
Get the Radial Distribution Figure for a given site.


* **Parameters**


    * **isite** – Index of the site.


    * **normalized** – Whether to normalize distances.


    * **figsize** – Size of the figure.


    * **step_function** – Type of step function to be used for the RDF.



* **Returns**

    Matplotlib figure.



#### get_sadf_figure(isite, normalized=True, figsize=None, step_function=None)
Get the Solid Angle Distribution Figure for a given site.


* **Parameters**


    * **isite** – Index of the site.


    * **normalized** – Whether to normalize angles.


    * **figsize** – Size of the figure.


    * **step_function** – Type of step function to be used for the SADF.



* **Returns**

    Matplotlib figure.



#### is_close_to(other, rtol=0.0, atol=1e-08)
Whether two DetailedVoronoiContainer objects are close to each other.


* **Parameters**


    * **other** – Another DetailedVoronoiContainer to be compared with.


    * **rtol** – Relative tolerance to compare values.


    * **atol** – Absolute tolerance to compare values.



* **Returns**

    True if the two DetailedVoronoiContainer are close to each other.



#### maps_and_surfaces(isite, surface_calculation_type=None, max_dist=2.0, additional_conditions=None)
Get the different surfaces and their cn_map corresponding to the different distance-angle cutoffs
for a given site.


* **Parameters**


    * **isite** – Index of the site


    * **surface_calculation_type** – How to compute the surface.


    * **max_dist** – The maximum distance factor to be considered.


    * **additional_conditions** – If additional conditions have to be considered.



* **Returns**

    Surfaces and cn_map’s for each distance-angle cutoff.



#### maps_and_surfaces_bounded(isite, surface_calculation_options=None, additional_conditions=None)
Get the different surfaces (using boundaries) and their cn_map corresponding to the different
distance-angle cutoffs for a given site.


* **Parameters**


    * **isite** – Index of the site


    * **surface_calculation_options** – Options for the boundaries.


    * **additional_conditions** – If additional conditions have to be considered.



* **Returns**

    Surfaces and cn_map’s for each distance-angle cutoff.



#### neighbors(isite, distfactor, angfactor, additional_condition=None)
Get the neighbors of a given site corresponding to a given distance and angle factor.


* **Parameters**


    * **isite** – Index of the site.


    * **distfactor** – Distance factor.


    * **angfactor** – Angle factor.


    * **additional_condition** – Additional condition to be used (currently not implemented).



* **Returns**

    List of neighbors of the given site for the given distance and angle factors.



#### neighbors_surfaces(isite, surface_calculation_type=None, max_dist=2.0)
Get the different surfaces corresponding to the different distance-angle cutoffs for a given site.


* **Parameters**


    * **isite** – Index of the site


    * **surface_calculation_type** – How to compute the surface.


    * **max_dist** – The maximum distance factor to be considered.



* **Returns**

    Surfaces for each distance-angle cutoff.



#### neighbors_surfaces_bounded(isite, surface_calculation_options=None)
Get the different surfaces (using boundaries) corresponding to the different distance-angle cutoffs
for a given site.


* **Parameters**


    * **isite** – Index of the site.


    * **surface_calculation_options** – Options for the boundaries.



* **Returns**

    Surfaces for each distance-angle cutoff.



#### setup_neighbors_distances_and_angles(indices)
Initializes the angle and distance separations.


* **Parameters**

    **indices** – Indices of the sites for which the Voronoi is needed.



#### setup_voronoi_list(indices, voronoi_cutoff)
Set up of the voronoi list of neighbors by calling qhull.


* **Parameters**


    * **indices** – indices of the sites for which the Voronoi is needed.


    * **voronoi_cutoff** – Voronoi cutoff for the search of neighbors.



* **Raises**

    **RuntimeError** – If an infinite vertex is found in the voronoi construction.



#### to_bson_voronoi_list2()
Transforms the voronoi_list into a vlist + bson_nb_voro_list, that are BSON-encodable.


* **Returns**

    [vlist, bson_nb_voro_list], to be used in the as_dict method.



#### voronoi_parameters_bounds_and_limits(isite, plot_type, max_dist)
Get the different boundaries and limits of the distance and angle factors for the given site.


* **Parameters**


    * **isite** – Index of the site.


    * **plot_type** – Types of distance/angle parameters to get.


    * **max_dist** – Maximum distance factor.



* **Returns**

    Distance and angle bounds and limits.



### pymatgen.analysis.chemenv.coordination_environments.voronoi.from_bson_voronoi_list2(bson_nb_voro_list2, structure)
Returns the voronoi_list needed for the VoronoiContainer object from a bson-encoded voronoi_list.


* **Parameters**


    * **bson_nb_voro_list2** – List of periodic sites involved in the Voronoi.


    * **structure** – Structure object.



* **Returns**

    The voronoi_list needed for the VoronoiContainer (with PeriodicSites as keys of the dictionary - not
    allowed in the BSON format).