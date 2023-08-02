---
layout: default
title: pymatgen.transformations.transformation_abc.md
nav_exclude: true
---

# pymatgen.transformations.transformation_abc module

Defines an abstract base class contract for Transformation object.


### _class_ pymatgen.transformations.transformation_abc.AbstractTransformation()
Bases: `MSONable`

Abstract transformation class.


#### _abstract_ apply_transformation(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Applies the transformation to a structure. Depending on whether a
transformation is one-to-many, there may be an option to return a
ranked list of structures.


* **Parameters**


    * **structure** – input structure


    * **return_ranked_list** (*bool** | **int**, **optional*) – If return_ranked_list is int, that number of structures

    is returned. If False, only the single lowest energy structure is returned. Defaults to False.




* **Returns**

    depending on returned_ranked list, either a transformed structure
    or
    a list of dictionaries, where each dictionary is of the form
    {‘structure’ = …. , ‘other_arguments’}
    the key ‘transformation’ is reserved for the transformation that
    was actually applied to the structure.
    This transformation is parsed by the alchemy classes for generating
    a more specific transformation history. Any other information will
    be stored in the transformation_parameters dictionary in the
    transmuted structure class.



#### _abstract property_ inverse(_: AbstractTransformation | Non_ )
Returns the inverse transformation if available.
Otherwise, should return None.


#### _abstract property_ is_one_to_many(_: boo_ )
Determines if a Transformation is a one-to-many transformation. If a
Transformation is a one-to-many transformation, the
apply_transformation method should have a keyword arg
“return_ranked_list” which allows for the transformed structures to be
returned as a ranked list.


#### _property_ use_multiprocessing(_: boo_ )
Indicates whether the transformation can be applied by a
subprocessing pool. This should be overridden to return True for
transformations that the transmuter can parallelize.