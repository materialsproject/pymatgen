---
layout: default
title: pymatgen.analysis.chemenv.coordination_environments.structure_environments.md
nav_exclude: true
---

# pymatgen.analysis.chemenv.coordination_environments.structure_environments module

This module contains objects that are used to describe the environments in a structure. The most detailed object
(StructureEnvironments) contains a very thorough analysis of the environments of a given atom but is difficult to
used as such. The LightStructureEnvironments object is a lighter version that is obtained by applying a “strategy”
on the StructureEnvironments object. Basically, the LightStructureEnvironments provides the coordination environment(s)
and possibly some fraction corresponding to these.


### _class_ pymatgen.analysis.chemenv.coordination_environments.structure_environments.ChemicalEnvironments(coord_geoms=None)
Bases: `MSONable`

Class used to store all the information about the chemical environment of a given site for a given list of
coordinated neighbors (internally called “cn_map”).

Initializes the ChemicalEnvironments object containing all the information about the chemical
environment of a given site.


* **Parameters**

    **coord_geoms** – coordination geometries to be added to the chemical environment.



#### add_coord_geom(mp_symbol, symmetry_measure, algo='UNKNOWN', permutation=None, override=False, local2perfect_map=None, perfect2local_map=None, detailed_voronoi_index=None, other_symmetry_measures=None, rotation_matrix=None, scaling_factor=None)
Adds a coordination geometry to the ChemicalEnvironments object.


* **Parameters**


    * **mp_symbol** – Symbol of the coordination geometry added.


    * **symmetry_measure** – Symmetry measure of the coordination geometry added.


    * **algo** – Algorithm used for the search of the coordination geometry added.


    * **permutation** – Permutation of the neighbors that leads to the csm stored.


    * **override** – If set to True, the coordination geometry will override the existent one if present.


    * **local2perfect_map** – Mapping of the local indices to the perfect indices.


    * **perfect2local_map** – Mapping of the perfect indices to the local indices.


    * **detailed_voronoi_index** – Index in the voronoi containing the neighbors set.


    * **other_symmetry_measures** – Other symmetry measure of the coordination geometry added (with/without the
    central atom, centered on the central atom or on the centroid with/without the central atom).


    * **rotation_matrix** – Rotation matrix mapping the local geometry to the perfect geometry.


    * **scaling_factor** – Scaling factor mapping the local geometry to the perfect geometry.



* **Raises**

    **ChemenvError if the coordination geometry is already added and override is set to False** –



#### as_dict()
Returns a dictionary representation of the ChemicalEnvironments object.


* **Returns**

    A dictionary representation of the ChemicalEnvironments object.



#### _classmethod_ from_dict(d)
Reconstructs the ChemicalEnvironments object from a dict representation of the ChemicalEnvironments created
using the as_dict method.


* **Parameters**

    **d** – dict representation of the ChemicalEnvironments object.



* **Returns**

    ChemicalEnvironments object.



#### is_close_to(other, rtol=0.0, atol=1e-08)
Whether this ChemicalEnvironments object is close to another one.


* **Parameters**


    * **other** – Another ChemicalEnvironments object.


    * **rtol** – Relative tolerance for the comparison of Continuous Symmetry Measures.


    * **atol** – Absolute tolerance for the comparison of Continuous Symmetry Measures.



* **Returns**

    True if the two ChemicalEnvironments objects are close to each other.



#### minimum_geometries(n=None, symmetry_measure_type=None, max_csm=None)
Returns a list of geometries with increasing continuous symmetry measure in this ChemicalEnvironments object.


* **Parameters**

    **n** – Number of geometries to be included in the list.



* **Returns**

    List of geometries with increasing continuous symmetry measure in this ChemicalEnvironments object.



* **Raises**

    **ValueError if no coordination geometry is found in this ChemicalEnvironments object.** –



#### minimum_geometry(symmetry_measure_type=None, max_csm=None)
Returns the geometry with the minimum continuous symmetry measure of this ChemicalEnvironments.


* **Returns**

    tuple (symbol, csm) with symbol being the geometry with the minimum continuous symmetry measure and
    csm being the continuous symmetry measure associated to it.



* **Raises**

    **ValueError if no coordination geometry is found in this ChemicalEnvironments object.** –



### _class_ pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments(strategy, coordination_environments=None, all_nbs_sites=None, neighbors_sets=None, structure=None, valences=None, valences_origin=None)
Bases: `MSONable`

Class used to store the chemical environments of a given structure obtained from a given ChemenvStrategy. Currently,
only strategies leading to the determination of a unique environment for each site is allowed
This class does not store all the information contained in the StructureEnvironments object, only the coordination
environment found.

Constructor for the LightStructureEnvironments object.


* **Parameters**


    * **strategy** – ChemEnv strategy used to get the environments.


    * **coordination_environments** – The coordination environments identified.


    * **all_nbs_sites** – All the possible neighbors for each site in the structure.


    * **neighbors_sets** – The neighbors sets of each site in the structure.


    * **structure** – The structure.


    * **valences** – The valences used to get the environments (if needed).


    * **valences_origin** – How the valences were obtained (e.g. from the Bond-valence analysis or from the original
    structure).



#### DEFAULT_STATISTICS_FIELDS(_ = ('anion_list', 'anion_atom_list', 'cation_list', 'cation_atom_list', 'neutral_list', 'neutral_atom_list', 'atom_coordination_environments_present', 'ion_coordination_environments_present', 'fraction_atom_coordination_environments_present', 'fraction_ion_coordination_environments_present', 'coordination_environments_atom_present', 'coordination_environments_ion_present'_ )

#### DELTA_MAX_OXIDATION_STATE(_ = 0._ )

#### _class_ NeighborsSet(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), isite, all_nbs_sites, all_nbs_sites_indices)
Bases: `object`

Class used to store a given set of neighbors of a given site (based on a list of sites, the voronoi
container is not part of the LightStructureEnvironments object).

Constructor for NeighborsSet.


* **Parameters**


    * **structure** – Structure object.


    * **isite** – Index of the site for which neighbors are stored in this NeighborsSet.


    * **all_nbs_sites** – All the possible neighbors for this site.


    * **all_nbs_sites_indices** – Indices of the sites in all_nbs_sites that make up this NeighborsSet.



#### as_dict()
A JSON-serializable dict representation of the NeighborsSet.


#### _classmethod_ from_dict(dd, structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), all_nbs_sites)
Reconstructs the NeighborsSet algorithm from its JSON-serializable dict representation, together with
the structure and all the possible neighbors sites.

As an inner (nested) class, the NeighborsSet is not supposed to be used anywhere else that inside the
LightStructureEnvironments. The from_dict method is thus using the structure and all_nbs_sites when
reconstructing itself. These two are both in the LightStructureEnvironments object.


* **Parameters**


    * **dd** – a JSON-serializable dict representation of a NeighborsSet.


    * **structure** – The structure.


    * **all_nbs_sites** – The list of all the possible neighbors for a given site.


Returns: a NeighborsSet.


#### _property_ neighb_coords()
Coordinates of neighbors for this NeighborsSet.


#### _property_ neighb_indices_and_images(_: list[dict[str, int]_ )
List of indices and images with respect to the original unit cell sites for this NeighborsSet.


#### _property_ neighb_sites()
Neighbors for this NeighborsSet as pymatgen Sites.


#### _property_ neighb_sites_and_indices()
List of neighbors for this NeighborsSet as pymatgen Sites and their index in the original structure.


#### as_dict()

* **Returns**

    Bson-serializable representation of the LightStructureEnvironments object.



* **Return type**

    dict



#### clear_environments(conditions=None)
Get the clear environments in the structure.


* **Parameters**

    **conditions** – Conditions to be checked for an environment to be “clear”.


Returns: Set of clear environments in this structure.


#### contains_only_one_anion(anion)
Whether this LightStructureEnvironments concerns a structure with only one given anion type.


* **Parameters**

    **anion** – Anion (e.g. O2-, …).


Returns: True if this LightStructureEnvironments concerns a structure with only one given anion.


#### contains_only_one_anion_atom(anion_atom)
Whether this LightStructureEnvironments concerns a structure with only one given anion atom type.


* **Parameters**

    **anion_atom** – Anion (e.g. O, …). The structure could contain O2- and O- though.


Returns: True if this LightStructureEnvironments concerns a structure with only one given anion_atom.


#### environments_identified()
Return the set of environments identified in this structure.

Returns: Set of environments identified in this structure.


#### _classmethod_ from_dict(d)
Reconstructs the LightStructureEnvironments object from a dict representation of the
LightStructureEnvironments created using the as_dict method.


* **Parameters**

    **d** – dict representation of the LightStructureEnvironments object.



* **Returns**

    LightStructureEnvironments object.



#### _classmethod_ from_structure_environments(strategy, structure_environments, valences=None, valences_origin=None)
Construct a LightStructureEnvironments object from a strategy and a StructureEnvironments object.


* **Parameters**


    * **strategy** – ChemEnv strategy used.


    * **structure_environments** – StructureEnvironments object from which to construct the LightStructureEnvironments.


    * **valences** – The valences of each site in the structure.


    * **valences_origin** – How the valences were obtained (e.g. from the Bond-valence analysis or from the original
    structure).


Returns: a LightStructureEnvironments object.


#### get_site_info_for_specie_allces(specie, min_fraction=0)
Get list of indices that have the given specie.


* **Parameters**


    * **specie** – Species to get.


    * **min_fraction** – Minimum fraction of the coordination environment.


Returns: Dictionary with the list of coordination environments for the given species, the indices of the sites

    in which they appear, their fractions and continuous symmetry measures.


#### get_site_info_for_specie_ce(specie, ce_symbol)
Get list of indices that have the given specie with a given Coordination environment.


* **Parameters**


    * **specie** – Species to get.


    * **ce_symbol** – Symbol of the coordination environment to get.


Returns: Dictionary with the list of indices in the structure that have the given specie in the given

    environment, their fraction and continuous symmetry measures.


#### get_statistics(statistics_fields=('anion_list', 'anion_atom_list', 'cation_list', 'cation_atom_list', 'neutral_list', 'neutral_atom_list', 'atom_coordination_environments_present', 'ion_coordination_environments_present', 'fraction_atom_coordination_environments_present', 'fraction_ion_coordination_environments_present', 'coordination_environments_atom_present', 'coordination_environments_ion_present'), bson_compatible=False)
Get the statistics of environments for this structure.


* **Parameters**


    * **statistics_fields** – Which statistics to get.


    * **bson_compatible** – Whether to make the dictionary BSON-compatible.



* **Returns**

    A dictionary with the requested statistics.



#### setup_statistic_lists()
Set up the statistics of environments for this LightStructureEnvironments.


#### site_contains_environment(isite, ce_symbol)
Whether a given site contains a given coordination environment.


* **Parameters**


    * **isite** – Index of the site.


    * **ce_symbol** – Symbol of the coordination environment.


Returns: True if the site contains the given coordination environment.


#### site_has_clear_environment(isite, conditions=None)
Whether a given site has a “clear” environments.

A “clear” environment is somewhat arbitrary. You can pass (multiple) conditions, e.g. the environment should
have a continuous symmetry measure lower than this, a fraction higher than that, …


* **Parameters**


    * **isite** – Index of the site.


    * **conditions** – Conditions to be checked for an environment to be “clear”.


Returns: True if the site has a clear environment.


#### structure_contains_atom_environment(atom_symbol, ce_symbol)
Checks whether the structure contains a given atom in a given environment.


* **Parameters**


    * **atom_symbol** – Symbol of the atom.


    * **ce_symbol** – Symbol of the coordination environment.



* **Returns**

    True if the coordination environment is found, False otherwise



#### structure_has_clear_environments(conditions=None, skip_none=True, skip_empty=False)
Whether all sites in a structure have “clear” environments.


* **Parameters**


    * **conditions** – Conditions to be checked for an environment to be “clear”.


    * **skip_none** – Whether to skip sites for which no environments have been computed.


    * **skip_empty** – Whether to skip sites for which no environments could be found.



* **Returns**

    True if all the sites in the structure have clear environments.



* **Return type**

    bool



#### _property_ uniquely_determines_coordination_environments()
True if the coordination environments are uniquely determined.


### _class_ pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments(voronoi, valences, sites_map, equivalent_sites, ce_list, structure, neighbors_sets=None, info=None)
Bases: `MSONable`

Class used to store the chemical environments of a given structure.

Constructor for the StructureEnvironments object.


* **Parameters**


    * **voronoi** – VoronoiContainer object for the structure.


    * **valences** – Valences provided.


    * **sites_map** – Mapping of equivalent sites to the unequivalent sites that have been computed.


    * **equivalent_sites** – List of list of equivalent sites of the structure.


    * **ce_list** – List of chemical environments.


    * **structure** – Structure object.


    * **neighbors_sets** – List of neighbors sets.


    * **info** – Additional information for this StructureEnvironments object.



#### AC(_ = <pymatgen.analysis.chemenv.utils.defs_utils.AdditionalConditions object_ )

#### _class_ NeighborsSet(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), isite, detailed_voronoi, site_voronoi_indices, sources=None)
Bases: `object`

Class used to store a given set of neighbors of a given site (based on the detailed_voronoi).

Constructor for NeighborsSet.


* **Parameters**


    * **structure** – Structure object.


    * **isite** – Index of the site for which neighbors are stored in this NeighborsSet.


    * **detailed_voronoi** – Corresponding DetailedVoronoiContainer object containing all the possible
    neighbors of the give site.


    * **site_voronoi_indices** – Indices of the voronoi sites in the DetailedVoronoiContainer object that
    make up this NeighborsSet.


    * **sources** – Sources for this NeighborsSet, i.e. how this NeighborsSet was generated.



#### add_source(source)
Add a source to this NeighborsSet.


* **Parameters**

    **source** – Information about the generation of this NeighborsSet.



#### angle_plateau()
Returns the angles plateau’s for this NeighborsSet.


#### _property_ angles()
Angles for each neighbor in this NeighborsSet.


#### as_dict()
A JSON-serializable dict representation of the NeighborsSet.


#### _property_ coords()
Coordinates of the current central atom and its neighbors for this NeighborsSet.


#### distance_plateau()
Returns the distances plateau’s for this NeighborsSet.


#### _property_ distances()
Distances to each neighbor in this NeighborsSet.


#### _classmethod_ from_dict(dd, structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), detailed_voronoi)
Reconstructs the NeighborsSet algorithm from its JSON-serializable dict representation, together with
the structure and the DetailedVoronoiContainer.

As an inner (nested) class, the NeighborsSet is not supposed to be used anywhere else that inside the
StructureEnvironments. The from_dict method is thus using the structure and  detailed_voronoi when
reconstructing itself. These two are both in the StructureEnvironments object.


* **Parameters**


    * **dd** – a JSON-serializable dict representation of a NeighborsSet.


    * **structure** – The structure.


    * **detailed_voronoi** – The Voronoi object containing all the neighboring atoms from which the subset of
    neighbors for this NeighborsSet is extracted.


Returns: a NeighborsSet.


#### get_neighb_voronoi_indices(permutation)
Return the indices in the detailed_voronoi corresponding to the current permutation.


* **Parameters**

    **permutation** – Current permutation for which the indices in the detailed_voronoi are needed.


Returns: List of indices in the detailed_voronoi.


#### _property_ info()
Summarized information about this NeighborsSet.


#### _property_ neighb_coords()
Coordinates of neighbors for this NeighborsSet.


#### _property_ neighb_coordsOpt()
Optimized access to the coordinates of neighbors for this NeighborsSet.


#### _property_ neighb_sites()
Neighbors for this NeighborsSet as pymatgen Sites.


#### _property_ neighb_sites_and_indices()
List of neighbors for this NeighborsSet as pymatgen Sites and their index in the original structure.


#### _property_ normalized_angles()
Normalized angles for each neighbor in this NeighborsSet.


#### _property_ normalized_distances()
Normalized distances to each neighbor in this NeighborsSet.


#### _property_ source()
Returns the source of this NeighborsSet (how it was generated, e.g. from which Voronoi cut-offs, or from
hints).


#### voronoi_grid_surface_points(additional_condition=1, other_origins='DO_NOTHING')
Get the surface points in the Voronoi grid for this neighbor from the sources.
The general shape of the points should look like a staircase such as in the following figure :

> ^

0.0|

    B—-C|    ||    |a  |      k    D——-E
n  |      |            |
g  |      |            |
l  |      |            |
e  |      j            F—-n———G

> |                           ||                           |A—-g——-h—-i———H1.0+————————————————->

    1.0              distance              2.0   ->+Inf


* **Parameters**


    * **additional_condition** – Additional condition for the neighbors.


    * **other_origins** – What to do with sources that do not come from the Voronoi grid (e.g. “from hints”).



#### add_neighbors_set(isite, nb_set)
Adds a neighbor set to the list of neighbors sets for this site.


* **Parameters**


    * **isite** – Index of the site under consideration.


    * **nb_set** – NeighborsSet to be added.



#### as_dict()
Bson-serializable dict representation of the StructureEnvironments object.


* **Returns**

    Bson-serializable dict representation of the StructureEnvironments object.



#### differences_wrt(other)
Return differences found in the current StructureEnvironments with respect to another StructureEnvironments.


* **Parameters**

    **other** – A StructureEnvironments object.



* **Returns**

    List of differences between the two StructureEnvironments objects.



#### _classmethod_ from_dict(d)
Reconstructs the StructureEnvironments object from a dict representation of the StructureEnvironments created
using the as_dict method.


* **Parameters**

    **d** – dict representation of the StructureEnvironments object.



* **Returns**

    StructureEnvironments object.



#### get_coordination_environments(isite, cn, nb_set)
Get the ChemicalEnvironments for a given site, coordination and neighbors set.


* **Parameters**


    * **isite** – Index of the site for which the ChemicalEnvironments is looked for.


    * **cn** – Coordination for which the ChemicalEnvironments is looked for.


    * **nb_set** – Neighbors set for which the ChemicalEnvironments is looked for.


Returns: a ChemicalEnvironments object.


#### get_csm(isite, mp_symbol)
Get the continuous symmetry measure for a given site in the given coordination environment.


* **Parameters**


    * **isite** – Index of the site.


    * **mp_symbol** – Symbol of the coordination environment for which we want the continuous symmetry measure.


Returns: Continuous symmetry measure of the given site in the given environment.


#### get_csm_and_maps(isite, max_csm=8.0, figsize=None, symmetry_measure_type=None)
Plotting of the coordination numbers of a given site for all the distfactor/angfactor parameters. If the
chemical environments are given, a color map is added to the plot, with the lowest continuous symmetry measure
as the value for the color of that distfactor/angfactor set.


* **Parameters**


    * **isite** – Index of the site for which the plot has to be done.


    * **max_csm** – Maximum continuous symmetry measure to be shown.


    * **figsize** – Size of the figure.


    * **symmetry_measure_type** – Type of continuous symmetry measure to be used.



* **Returns**

    Matplotlib figure and axes representing the csm and maps.



#### get_csms(isite, mp_symbol)
Returns the continuous symmetry measure(s) of site with index isite with respect to the

    perfect coordination environment with mp_symbol. For some environments, a given mp_symbol might not
    be available (if there is no voronoi parameters leading to a number of neighbors corresponding to
    the coordination number of environment mp_symbol). For some environments, a given mp_symbol might
    lead to more than one csm (when two or more different voronoi parameters lead to different neighbors
    but with same number of neighbors).


* **Parameters**


    * **isite** – Index of the site.


    * **mp_symbol** – MP symbol of the perfect environment for which the csm has to be given.



* **Returns**

    List of csms for site isite with respect to geometry mp_symbol



#### get_environments_figure(isite, plot_type=None, title='Coordination numbers', max_dist=2.0, colormap=None, figsize=None, strategy=None)
Plotting of the coordination environments of a given site for all the distfactor/angfactor regions. The
chemical environments with the lowest continuous symmetry measure is shown for each distfactor/angfactor
region as the value for the color of that distfactor/angfactor region (using a colormap).


* **Parameters**


    * **isite** – Index of the site for which the plot has to be done.


    * **plot_type** – How to plot the coordinations.


    * **title** – Title for the figure.


    * **max_dist** – Maximum distance to be plotted when the plotting of the distance is set to ‘initial_normalized’
    or ‘initial_real’ (Warning: this is not the same meaning in both cases! In the first case, the
    closest atom lies at a “normalized” distance of 1.0 so that 2.0 means refers to this normalized
    distance while in the second case, the real distance is used).


    * **colormap** – Color map to be used for the continuous symmetry measure.


    * **figsize** – Size of the figure.


    * **strategy** – Whether to plot information about one of the Chemenv Strategies.



* **Returns**

    Matplotlib figure and axes representing the environments.



#### init_neighbors_sets(isite, additional_conditions=None, valences=None)
Initialize the list of neighbors sets for the current site.


* **Parameters**


    * **isite** – Index of the site under consideration.


    * **additional_conditions** – Additional conditions to be used for the initialization of the list of
    neighbors sets, e.g. “Only anion-cation bonds”, …


    * **valences** – List of valences for each site in the structure (needed if an additional condition based on the
    valence is used, e.g. only anion-cation bonds).



#### plot_csm_and_maps(isite, max_csm=8.0)
Plotting of the coordination numbers of a given site for all the distfactor/angfactor parameters. If the
chemical environments are given, a color map is added to the plot, with the lowest continuous symmetry measure
as the value for the color of that distfactor/angfactor set.


* **Parameters**


    * **isite** – Index of the site for which the plot has to be done


    * **max_csm** – Maximum continuous symmetry measure to be shown.



#### plot_environments(isite, plot_type=None, title='Coordination numbers', max_dist=2.0, figsize=None, strategy=None)
Plotting of the coordination numbers of a given site for all the distfactor/angfactor parameters. If the
chemical environments are given, a color map is added to the plot, with the lowest continuous symmetry measure
as the value for the color of that distfactor/angfactor set.


* **Parameters**


    * **isite** – Index of the site for which the plot has to be done.


    * **plot_type** – How to plot the coordinations.


    * **title** – Title for the figure.


    * **max_dist** – Maximum distance to be plotted when the plotting of the distance is set to ‘initial_normalized’
    or ‘initial_real’ (Warning: this is not the same meaning in both cases! In the first case, the
    closest atom lies at a “normalized” distance of 1.0 so that 2.0 means refers to this normalized
    distance while in the second case, the real distance is used).


    * **figsize** – Size of the figure.


    * **strategy** – Whether to plot information about one of the Chemenv Strategies.



#### save_environments_figure(isite, imagename='image.png', plot_type=None, title='Coordination numbers', max_dist=2.0, figsize=None)
Saves the environments figure to a given file.


* **Parameters**


    * **isite** – Index of the site for which the plot has to be done.


    * **imagename** – Name of the file to which the figure has to be saved.


    * **plot_type** – How to plot the coordinations.


    * **title** – Title for the figure.


    * **max_dist** – Maximum distance to be plotted when the plotting of the distance is set to ‘initial_normalized’
    or ‘initial_real’ (Warning: this is not the same meaning in both cases! In the first case, the
    closest atom lies at a “normalized” distance of 1.0 so that 2.0 means refers to this normalized
    distance while in the second case, the real distance is used).


    * **figsize** – Size of the figure.



#### update_coordination_environments(isite, cn, nb_set, ce)
Updates the coordination environment for this site, coordination and neighbor set.


* **Parameters**


    * **isite** – Index of the site to be updated.


    * **cn** – Coordination to be updated.


    * **nb_set** – Neighbors set to be updated.


    * **ce** – ChemicalEnvironments object for this neighbors set.



#### update_site_info(isite, info_dict)
Update information about this site.


* **Parameters**


    * **isite** – Index of the site for which info has to be updated.


    * **info_dict** – Dictionary of information to be added for this site.