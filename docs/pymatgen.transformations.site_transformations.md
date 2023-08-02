---
layout: default
title: pymatgen.transformations.site_transformations.md
nav_exclude: true
---

# pymatgen.transformations.site_transformations module

This module defines site transformations which transforms a structure into
another structure. Site transformations differ from standard transformations
in that they operate in a site-specific manner.
All transformations should inherit the AbstractTransformation ABC.


### _class_ pymatgen.transformations.site_transformations.AddSitePropertyTransformation(site_properties)
Bases: [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)

Simple transformation to add site properties to a given structure.


* **Parameters**

    **site_properties** (*dict*) – site properties to be added to a structure.



#### apply_transformation(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Apply the transformation.

Arg:

    structure (Structure): A structurally similar structure in

        regards to crystal and site positions.


* **Returns**

    Returns a copy of structure with sites properties added.



#### _property_ inverse()
None.


* **Type**

    Return



#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Return



### _class_ pymatgen.transformations.site_transformations.InsertSitesTransformation(species, coords, coords_are_cartesian=False, validate_proximity=True)
Bases: [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)

This transformation substitutes certain sites with certain species.


* **Parameters**


    * **species** – A list of species. e.g., [“Li”, “Fe”]


    * **coords** – A list of coords corresponding to those species. e.g.,
    [[0,0,0],[0.5,0.5,0.5]].


    * **coords_are_cartesian** (*bool*) – Set to True if coords are given in
    Cartesian coords. Defaults to False.


    * **validate_proximity** (*bool*) – Set to False if you do not wish to ensure
    that added sites are not too close to other sites. Defaults to True.



#### apply_transformation(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Apply the transformation.

Arg:

    structure (Structure): A structurally similar structure in

        regards to crystal and site positions.


* **Returns**

    Returns a copy of structure with sites inserted.



#### _property_ inverse()
None.


* **Type**

    Return



#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Return



### _class_ pymatgen.transformations.site_transformations.PartialRemoveSitesTransformation(indices, fractions, algo=1)
Bases: [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)

Remove fraction of specie from a structure.
Requires an oxidation state decorated structure for Ewald sum to be
computed.

Given that the solution to selecting the right removals is NP-hard, there
are several algorithms provided with varying degrees of accuracy and speed.
The options are as follows:

ALGO_FAST:

    This is a highly optimized algorithm to quickly go through the search
    tree. It is guaranteed to find the optimal solution, but will return
    only a single lowest energy structure. Typically, you will want to use
    this.

ALGO_COMPLETE:

    The complete algo ensures that you get all symmetrically distinct
    orderings, ranked by the estimated Ewald energy. But this can be an
    extremely time-consuming process if the number of possible orderings is
    very large. Use this if you really want all possible orderings. If you
    want just the lowest energy ordering, ALGO_FAST is accurate and faster.

ALGO_BEST_FIRST:

    This algorithm is for ordering the really large cells that defeats even
    ALGO_FAST. For example, if you have 48 sites of which you want to
    remove 16 of them, the number of possible orderings is around
    2 x 10^12. ALGO_BEST_FIRST shortcircuits the entire search tree by
    removing the highest energy site first, then followed by the next
    highest energy site, and so on. It is guaranteed to find a solution
    in a reasonable time, but it is also likely to be highly inaccurate.

ALGO_ENUMERATE:

    This algorithm uses the EnumerateStructureTransformation to perform
    ordering. This algo returns *complete* orderings up to a single unit
    cell size. It is more robust than the ALGO_COMPLETE, but requires
    Gus Hart’s enumlib to be installed.


* **Parameters**


    * **indices** – A list of list of indices.
    e.g. [[0, 1], [2, 3, 4, 5]]


    * **fractions** – The corresponding fractions to remove. Must be same length as
    indices. e.g., [0.5, 0.25]


    * **algo** – This parameter allows you to choose the algorithm to perform
    ordering. Use one of PartialRemoveSpecieTransformation.ALGO_\*
    variables to set the algo.



#### ALGO_BEST_FIRST(_ = _ )

#### ALGO_COMPLETE(_ = _ )

#### ALGO_ENUMERATE(_ = _ )

#### ALGO_FAST(_ = _ )

#### apply_transformation(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), return_ranked_list: bool | int = False)
Apply the transformation.


* **Parameters**


    * **structure** – input structure


    * **return_ranked_list** (*bool** | **int*) – Whether or not multiple structures are returned.
    If return_ranked_list is int, that number of structures is returned.



* **Returns**

    Depending on returned_ranked list, either a transformed structure
    or a list of dictionaries, where each dictionary is of the form
    {“structure” = …. , “other_arguments”}
    the key “transformation” is reserved for the transformation that
    was actually applied to the structure.
    This transformation is parsed by the alchemy classes for generating
    a more specific transformation history. Any other information will
    be stored in the transformation_parameters dictionary in the
    transmuted structure class.



#### _property_ inverse()
None.


* **Type**

    Return



#### _property_ is_one_to_many(_: boo_ )
True.


* **Type**

    Return



### _class_ pymatgen.transformations.site_transformations.RadialSiteDistortionTransformation(site_index, displacement=0.1, nn_only=False)
Bases: [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)

Radially perturbs atoms around a site. Can be used to create spherical distortion due to a
point defect.


* **Parameters**


    * **site_index** (*int*) – index of the site in structure to place at the center of the distortion (will
    not be distorted). This index must be provided before the structure is provided in
    apply_transformation in order to keep in line with the base class.


    * **displacement** (*float*) – distance to perturb the atoms around the objective site


    * **nn_only** (*bool*) – Whether or not to perturb beyond the nearest neighbors. If True, then only the
    nearest neighbors will be perturbed, leaving the other sites undisturbed. If False, then
    the nearest neighbors will receive the full displacement, and then subsequent sites will receive
    a displacement=0.1 / r, where r is the distance each site to the origin site. For small displacements,
    atoms beyond the NN environment will receive very small displacements, and these are almost equal.
    For large displacements, this difference is noticeable.



#### apply_transformation(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Apply the transformation.


* **Parameters**

    **structure** – Structure or Molecule to apply the transformation to



* **Returns**

    the transformed structure



#### _property_ inverse()
Returns the inverse transformation if available.
Otherwise, should return None.


#### _property_ is_one_to_many(_: boo_ )
Determines if a Transformation is a one-to-many transformation. If a
Transformation is a one-to-many transformation, the
apply_transformation method should have a keyword arg
“return_ranked_list” which allows for the transformed structures to be
returned as a ranked list.


#### _property_ use_multiprocessing()
Indicates whether the transformation can be applied by a
subprocessing pool. This should be overridden to return True for
transformations that the transmuter can parallelize.


### _class_ pymatgen.transformations.site_transformations.RemoveSitesTransformation(indices_to_remove)
Bases: [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)

Remove certain sites in a structure.


* **Parameters**

    **indices_to_remove** – List of indices to remove. E.g., [0, 1, 2].



#### apply_transformation(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Apply the transformation.

Arg:

    structure (Structure): A structurally similar structure in

        regards to crystal and site positions.


* **Returns**

    Returns a copy of structure with sites removed.



#### _property_ inverse()
None.


* **Type**

    Return



#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Return



### _class_ pymatgen.transformations.site_transformations.ReplaceSiteSpeciesTransformation(indices_species_map)
Bases: [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)

This transformation substitutes certain sites with certain species.


* **Parameters**

    **indices_species_map** – A dict containing the species mapping in
    int-string pairs. E.g., { 1:”Na”} or {2:”Mn2+”}. Multiple
    substitutions can be done. Overloaded to accept sp_and_occu
    dictionary. E.g. {1: {“Ge”:0.75, “C”:0.25} }, which
    substitutes a single species with multiple species to generate a
    disordered structure.



#### apply_transformation(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Apply the transformation.

Arg:

    structure (Structure): A structurally similar structure in

        regards to crystal and site positions.


* **Returns**

    Returns a copy of structure with sites replaced.



#### _property_ inverse()
None.


* **Type**

    Return



#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Return



### _class_ pymatgen.transformations.site_transformations.TranslateSitesTransformation(indices_to_move, translation_vector, vector_in_frac_coords=True)
Bases: [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)

This class translates a set of sites by a certain vector.


* **Parameters**


    * **indices_to_move** – The indices of the sites to move


    * **translation_vector** – Vector to move the sites. If a list of list or numpy
    array of shape, (len(indices_to_move), 3), is provided then each
    translation vector is applied to the corresponding site in the
    indices_to_move.


    * **vector_in_frac_coords** – Set to True if the translation vector is in
    fractional coordinates, and False if it is in cartesian
    coordinations. Defaults to True.



#### apply_transformation(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Apply the transformation.

Arg:

    structure (Structure): A structurally similar structure in

        regards to crystal and site positions.


* **Returns**

    Returns a copy of structure with sites translated.



#### as_dict()
JSON-serializable dict representation.


#### _property_ inverse()
Returns:
TranslateSitesTransformation with the reverse translation.


#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Return