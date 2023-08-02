---
layout: default
title: pymatgen.transformations.standard_transformations.md
nav_exclude: true
---

# pymatgen.transformations.standard_transformations module

This module defines standard transformations which transforms a structure into
another structure. Standard transformations operate in a structure-wide manner,
rather than site-specific manner.
All transformations should inherit the AbstractTransformation ABC.


### _class_ pymatgen.transformations.standard_transformations.AutoOxiStateDecorationTransformation(symm_tol=0.1, max_radius=4, max_permutations=100000, distance_scale_factor=1.015)
Bases: [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)

This transformation automatically decorates a structure with oxidation
states using a bond valence approach.


* **Parameters**


    * **symm_tol** (*float*) – Symmetry tolerance used to determine which sites are
    symmetrically equivalent. Set to 0 to turn off symmetry.


    * **max_radius** (*float*) – Maximum radius in Angstrom used to find nearest
    neighbors.


    * **max_permutations** (*int*) – Maximum number of permutations of oxidation
    states to test.


    * **distance_scale_factor** (*float*) – A scale factor to be applied. This is
    useful for scaling distances, esp in the case of
    calculation-relaxed structures, which may tend to under (GGA) or
    over bind (LDA). The default of 1.015 works for GGA. For
    experimental structure, set this to 1.



#### apply_transformation(structure)
Apply the transformation.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input Structure



* **Returns**

    Oxidation state decorated Structure.



#### _property_ inverse()
None.


* **Type**

    Returns



#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Returns



### _class_ pymatgen.transformations.standard_transformations.ChargedCellTransformation(charge=0)
Bases: [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)

The ChargedCellTransformation applies a charge to a structure (or defect
object).


* **Parameters**

    **charge** – A integer charge to apply to the structure.
    Defaults to zero. Has to be a single integer. e.g. 2.



#### apply_transformation(structure)
Apply the transformation.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input Structure



* **Returns**

    Charged Structure.



#### _property_ inverse()
NotImplementedError.


* **Type**

    Raises



#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Returns



### _class_ pymatgen.transformations.standard_transformations.ConventionalCellTransformation(symprec: float = 0.01, angle_tolerance=5, international_monoclinic=True)
Bases: [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)

This class finds the conventional cell of the input structure.


* **Parameters**


    * **symprec** (*float*) – tolerance as in SpacegroupAnalyzer


    * **angle_tolerance** (*float*) – angle tolerance as in SpacegroupAnalyzer


    * **international_monoclinic** (*bool*) – whether to use beta (True) or alpha (False)


as the non-right-angle in the unit cell.


#### apply_transformation(structure)
Returns most primitive cell for structure.


* **Parameters**

    **structure** – A structure



* **Returns**

    The same structure in a conventional standard setting



#### _property_ inverse()
None.


* **Type**

    Returns



#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Returns



### _class_ pymatgen.transformations.standard_transformations.DeformStructureTransformation(deformation=((1, 0, 0), (0, 1, 0), (0, 0, 1)))
Bases: [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)

This transformation deforms a structure by a deformation gradient matrix.


* **Parameters**

    **deformation** (*array*) – deformation gradient for the transformation.



#### apply_transformation(structure)
Apply the transformation.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input Structure



* **Returns**

    Deformed Structure.



#### _property_ inverse()
Returns:
Inverse Transformation.


#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Returns



### _class_ pymatgen.transformations.standard_transformations.DiscretizeOccupanciesTransformation(max_denominator=5, tol: float | None = None, fix_denominator=False)
Bases: [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)

Discretizes the site occupancies in a disordered structure; useful for
grouping similar structures or as a pre-processing step for order-disorder
transformations.


* **Parameters**


    * **max_denominator** – An integer maximum denominator for discretization. A higher
    denominator allows for finer resolution in the site occupancies.


    * **tol** – A float that sets the maximum difference between the original and
    discretized occupancies before throwing an error. If None, it is
    set to 1 / (4 \* max_denominator).


    * **fix_denominator** (*bool*) – If True, will enforce a common denominator for all species.
    This prevents a mix of denominators (for example, 1/3, 1/4)
    that might require large cell sizes to perform an enumeration.
    ‘tol’ needs to be > 1.0 in some cases.



#### apply_transformation(structure)
Discretizes the site occupancies in the structure.


* **Parameters**

    **structure** – disordered Structure to discretize occupancies



* **Returns**

    A new disordered Structure with occupancies discretized



#### _property_ inverse()
None.


* **Type**

    Returns



#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Returns



### _class_ pymatgen.transformations.standard_transformations.OrderDisorderedStructureTransformation(algo=0, symmetrized_structures=False, no_oxi_states=False)
Bases: [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)

Order a disordered structure. The disordered structure must be oxidation
state decorated for Ewald sum to be computed. No attempt is made to perform
symmetry determination to reduce the number of combinations.

Hence, attempting to performing ordering on a large number of disordered
sites may be extremely expensive. The time scales approximately with the
number of possible combinations. The algorithm can currently compute
approximately 5,000,000 permutations per minute.

Also, simple rounding of the occupancies are performed, with no attempt
made to achieve a target composition. This is usually not a problem for
most ordering problems, but there can be times where rounding errors may
result in structures that do not have the desired composition.
This second step will be implemented in the next iteration of the code.

If multiple fractions for a single species are found for different sites,
these will be treated separately if the difference is above a threshold
tolerance. currently this is .1

For example, if a fraction of .25 Li is on sites 0,1,2,3  and .5 on sites
4, 5, 6, 7 then 1 site from [0,1,2,3] will be filled and 2 sites from [4,5,6,7]
will be filled, even though a lower energy combination might be found by
putting all lithium in sites [4,5,6,7].

USE WITH CARE.


* **Parameters**


    * **algo** (*int*) – Algorithm to use.


    * **symmetrized_structures** (*bool*) – Whether the input structures are
    instances of SymmetrizedStructure, and that their symmetry
    should be used for the grouping of sites.


    * **no_oxi_states** (*bool*) – Whether to remove oxidation states prior to
    ordering.



#### ALGO_BEST_FIRST(_ = _ )

#### ALGO_COMPLETE(_ = _ )

#### ALGO_FAST(_ = _ )

#### apply_transformation(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), return_ranked_list: bool | int = False)
For this transformation, the apply_transformation method will return
only the ordered structure with the lowest Ewald energy, to be
consistent with the method signature of the other transformations.
However, all structures are stored in the  all_structures attribute in
the transformation object for easy access.


* **Parameters**


    * **structure** – Oxidation state decorated disordered structure to order


    * **return_ranked_list** (*bool** | **int**, **optional*) – If return_ranked_list is int, that number of structures

    is returned. If False, only the single lowest energy structure is returned. Defaults to False.




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

    Returns



#### _property_ is_one_to_many(_: boo_ )
True.


* **Type**

    Returns



#### _property_ lowest_energy_structure()
Lowest energy structure found.


* **Type**

    return



### _class_ pymatgen.transformations.standard_transformations.OxidationStateDecorationTransformation(oxidation_states)
Bases: [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)

This transformation decorates a structure with oxidation states.


* **Parameters**


    * **oxidation_states** (*dict*) – Oxidation states supplied as a dict,


    * **e.g.** – 1, “O”:-2}.


    * **{"Li"** – 1, “O”:-2}.



#### apply_transformation(structure)
Apply the transformation.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input Structure



* **Returns**

    Oxidation state decorated Structure.



#### _property_ inverse()
None.


* **Type**

    Returns



#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Returns



### _class_ pymatgen.transformations.standard_transformations.OxidationStateRemovalTransformation()
Bases: [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)

This transformation removes oxidation states from a structure.

No arg needed.


#### apply_transformation(structure)
Apply the transformation.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input Structure



* **Returns**

    Non-oxidation state decorated Structure.



#### _property_ inverse()
None.


* **Type**

    Returns



#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Returns



### _class_ pymatgen.transformations.standard_transformations.PartialRemoveSpecieTransformation(specie_to_remove, fraction_to_remove, algo=0)
Bases: [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)

Remove fraction of specie from a structure.

Requires an oxidation state decorated structure for Ewald sum to be
computed.

Given that the solution to selecting the right removals is NP-hard, there
are several algorithms provided with varying degrees of accuracy and speed.
Please see
[`pymatgen.transformations.site_transformations.PartialRemoveSitesTransformation`](pymatgen.transformations.site_transformations.md#pymatgen.transformations.site_transformations.PartialRemoveSitesTransformation).


* **Parameters**


    * **specie_to_remove** – Species to remove. Must have oxidation state E.g.,
    “Li+”


    * **fraction_to_remove** – Fraction of specie to remove. E.g., 0.5


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


    * **return_ranked_list** (*bool** | **int**, **optional*) – If return_ranked_list is int, that number of structures

    is returned. If False, only the single lowest energy structure is returned. Defaults to False.




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

    Returns



#### _property_ is_one_to_many(_: boo_ )
True.


* **Type**

    Returns



### _class_ pymatgen.transformations.standard_transformations.PerturbStructureTransformation(distance: float = 0.01, min_distance: int | float | None = None)
Bases: [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)

This transformation perturbs a structure by a specified distance in random
directions. Used for breaking symmetries.


* **Parameters**


    * **distance** – Distance of perturbation in angstroms. All sites
    will be perturbed by exactly that distance in a random
    direction.


    * **min_distance** – if None, all displacements will be equidistant. If int
    or float, perturb each site a distance drawn from the uniform
    distribution between ‘min_distance’ and ‘distance’.



#### apply_transformation(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Apply the transformation.


* **Parameters**

    **structure** – Input Structure



* **Returns**

    Structure with sites perturbed.



#### _property_ inverse()
None.


* **Type**

    Returns



#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Returns



### _class_ pymatgen.transformations.standard_transformations.PrimitiveCellTransformation(tolerance=0.5)
Bases: [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)

This class finds the primitive cell of the input structure.
It returns a structure that is not necessarily orthogonalized
Author: Will Richards.


* **Parameters**

    **tolerance** (*float*) – Tolerance for each coordinate of a particular
    site. For example, [0.5, 0, 0.5] in Cartesian coordinates will be
    considered to be on the same coordinates as [0, 0, 0] for a
    tolerance of 0.5. Defaults to 0.5.



#### apply_transformation(structure)
Returns most primitive cell for structure.


* **Parameters**

    **structure** – A structure



* **Returns**

    The most primitive structure found. The returned structure is
    guaranteed to have len(new structure) <= len(structure).



#### _property_ inverse()
None.


* **Type**

    Returns



#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Returns



### _class_ pymatgen.transformations.standard_transformations.RemoveSpeciesTransformation(species_to_remove)
Bases: [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)

Remove all occurrences of some species from a structure.


* **Parameters**

    **species_to_remove** – List of species to remove. E.g., [“Li”, “Mn”].



#### apply_transformation(structure)
Apply the transformation.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input Structure



* **Returns**

    Structure with species removed.



#### _property_ inverse()
None.


* **Type**

    Returns



#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Returns



### _class_ pymatgen.transformations.standard_transformations.RotationTransformation(axis, angle, angle_in_radians=False)
Bases: [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)

The RotationTransformation applies a rotation to a structure.


* **Parameters**


    * **axis** (*3x1 array*) – Axis of rotation, e.g., [1, 0, 0]


    * **angle** (*float*) – Angle to rotate


    * **angle_in_radians** (*bool*) – Set to True if angle is supplied in radians.
    Else degrees are assumed.



#### apply_transformation(structure)
Apply the transformation.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input Structure



* **Returns**

    Rotated Structure.



#### _property_ inverse()
Returns:
Inverse Transformation.


#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Returns



### _class_ pymatgen.transformations.standard_transformations.ScaleToRelaxedTransformation(unrelaxed_structure, relaxed_structure, species_map=None)
Bases: [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)

Takes the unrelaxed and relaxed structure and applies its site and volume
relaxation to a structurally similar structures (e.g. bulk: NaCl and PbTe
(rock-salt), slab: Sc(10-10) and Mg(10-10) (hcp), GB: Mo(001) sigma 5 GB,
Fe(001) sigma 5). Useful for finding an initial guess of a set of similar
structures closer to its most relaxed state.


* **Parameters**


    * **unrelaxed_structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Initial, unrelaxed structure


    * **relaxed_structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Relaxed structure


    * **species_map** (*dict*) – A dict or list of tuples containing the species mapping in
    string-string pairs. The first species corresponds to the relaxed
    structure while the second corresponds to the species in the
    structure to be scaled. E.g., {“Li”:”Na”} or [(“Fe2+”,”Mn2+”)].
    Multiple substitutions can be done. Overloaded to accept
    sp_and_occu dictionary E.g. {“Si: {“Ge”:0.75, “C”:0.25}},
    which substitutes a single species with multiple species to
    generate a disordered structure.



#### apply_transformation(structure)
Returns a copy of structure with lattice parameters
and sites scaled to the same degree as the relaxed_structure.

Arg:

    structure (Structure): A structurally similar structure in

        regards to crystal and site positions.


#### _property_ inverse()
None.


* **Type**

    Returns



#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Returns



### _class_ pymatgen.transformations.standard_transformations.SubstitutionTransformation(species_map: dict[SpeciesLike, SpeciesLike | dict[SpeciesLike, float]] | list[tuple[SpeciesLike, SpeciesLike]])
Bases: [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)

This transformation substitutes species for one another.


* **Parameters**

    **species_map** – A dict or list of tuples containing the species mapping in
    string-string pairs. E.g., {“Li”: “Na”} or [(“Fe2+”,”Mn2+”)].
    Multiple substitutions can be done. Overloaded to accept
    sp_and_occu dictionary E.g. {“Si: {“Ge”:0.75, “C”:0.25}},
    which substitutes a single species with multiple species to
    generate a disordered structure.



#### apply_transformation(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Apply the transformation.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input Structure



* **Returns**

    Substituted Structure.



#### _property_ inverse()
Returns:
Inverse Transformation.


#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Returns



### _class_ pymatgen.transformations.standard_transformations.SupercellTransformation(scaling_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)))
Bases: [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)

The RotationTransformation applies a rotation to a structure.


* **Parameters**

    **scaling_matrix** – A matrix of transforming the lattice vectors.
    Defaults to the identity matrix. Has to be all integers. e.g.,
    [[2,1,0],[0,3,0],[0,0,1]] generates a new structure with
    lattice vectors a” = 2a + b, b” = 3b, c” = c where a, b, and c
    are the lattice vectors of the original structure.



#### apply_transformation(structure)
Apply the transformation.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input Structure



* **Returns**

    Supercell Structure.



#### _static_ from_scaling_factors(scale_a=1, scale_b=1, scale_c=1)
Convenience method to get a SupercellTransformation from a simple
series of three numbers for scaling each lattice vector. Equivalent to
calling the normal with [[scale_a, 0, 0], [0, scale_b, 0],
[0, 0, scale_c]].


* **Parameters**


    * **scale_a** – Scaling factor for lattice direction a. Defaults to 1.


    * **scale_b** – Scaling factor for lattice direction b. Defaults to 1.


    * **scale_c** – Scaling factor for lattice direction c. Defaults to 1.



* **Returns**

    SupercellTransformation.



#### _property_ inverse()
NotImplementedError.


* **Type**

    Raises



#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Returns