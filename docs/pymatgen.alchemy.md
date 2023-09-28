---
layout: default
title: pymatgen.alchemy.md
nav_exclude: true
---

1. TOC
{:toc}

# pymatgen.alchemy package

This package provides the modules for performing large scale transformations on
a large number of structures.


## pymatgen.alchemy.filters module

This module defines filters for Transmuter object.


### _class_ AbstractStructureFilter()
Bases: `MSONable`

AbstractStructureFilter that defines an API to perform testing of
Structures. Structures that return True to a test are retained during
transmutation while those that return False are removed.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### _abstract_ test(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Method to execute the test.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input structure to test



* **Returns**

    (bool) Structures that return true are kept in the Transmuter
    object during filtering.



### _class_ ChargeBalanceFilter()
Bases: `AbstractStructureFilter`

This filter removes structures that are not charge balanced from the
transmuter. This only works if the structure is oxidation state
decorated, as structures with only elemental sites are automatically
assumed to have net charge of 0.

No args required.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### test(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Method to execute the test.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input structure to test



* **Returns**

    True if structure is neutral.



* **Return type**

    bool



### _class_ ContainsSpecieFilter(species, strict_compare=False, AND=True, exclude=False)
Bases: `AbstractStructureFilter`

Filter for structures containing certain elements or species.
By default compares by atomic number.


* **Parameters**


    * **species** (*[**Species/Element**]*) – list of species to look for


    * **AND** – whether all species must be present to pass (or fail) filter.


    * **strict_compare** – if true, compares objects by specie or element
    object if false, compares atomic number


    * **exclude** – If true, returns false for any structures with the specie
    (excludes them from the Transmuter).



#### _abc_impl(_ = <_abc._abc_data object_ )

#### as_dict()
Returns: MSONable dict.


#### _classmethod_ from_dict(dct)

* **Parameters**

    **dct** (*dict*) – Dict representation.



* **Returns**

    Filter



#### test(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Method to execute the test.


* **Returns**

    True if structure does not contain specified species.



* **Return type**

    bool



### _class_ RemoveDuplicatesFilter(structure_matcher: dict | [StructureMatcher](pymatgen.analysis.md#pymatgen.analysis.structure_matcher.StructureMatcher) | None = None, symprec: float | None = None)
Bases: `AbstractStructureFilter`

This filter removes exact duplicate structures from the transmuter.

Remove duplicate structures based on the structure matcher
and symmetry (if symprec is given).


* **Parameters**


    * **structure_matcher** (*dict** | *[*StructureMatcher*](pymatgen.analysis.md#pymatgen.analysis.structure_matcher.StructureMatcher)*, **optional*) – Provides a structure matcher to be used for
    structure comparison.


    * **symprec** (*float**, **optional*) – The precision in the symmetry finder algorithm if None (
    default value), no symmetry check is performed and only the
    structure matcher is used. A recommended value is 1e-5.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### test(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))

* **Parameters**

    **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input structure to test.



* **Returns**

    True if structure is not in list.



* **Return type**

    bool



### _class_ RemoveExistingFilter(existing_structures, structure_matcher=None, symprec=None)
Bases: `AbstractStructureFilter`

This filter removes structures existing in a given list from the transmuter.

Remove existing structures based on the structure matcher
and symmetry (if symprec is given).


* **Parameters**


    * **existing_structures** – List of existing structures to compare with


    * **structure_matcher** – Provides a structure matcher to be used for
    structure comparison.


    * **symprec** – The precision in the symmetry finder algorithm if None (
    default value), no symmetry check is performed and only the
    structure matcher is used. A recommended value is 1e-5.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### as_dict()
Returns: MSONable dict.


#### test(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Method to execute the test.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input structure to test



* **Returns**

    True if structure is not in existing list.



* **Return type**

    bool



### _class_ SpecieProximityFilter(specie_and_min_dist_dict)
Bases: `AbstractStructureFilter`

This filter removes structures that have certain species that are too close
together.


* **Parameters**

    **specie_and_min_dist_dict** (*dict*) – A species string to float mapping. For
    example, {“Na+”: 1} means that all Na+ ions must be at least 1
    Angstrom away from each other. Multiple species criteria can be
    applied. Note that the testing is done based on the actual object
    . If you have a structure with Element, you must use {“Na”:1}
    instead to filter based on Element and not Species.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### as_dict()
Returns: MSONable dict.


#### _classmethod_ from_dict(dct)

* **Parameters**

    **dct** (*dict*) – Dict representation.



* **Returns**

    Filter



#### test(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Method to execute the test.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input structure to test



* **Returns**

    True if structure does not contain species within specified distances.



* **Return type**

    bool



### _class_ SpeciesMaxDistFilter(sp1, sp2, max_dist)
Bases: `AbstractStructureFilter`

This filter removes structures that do have two particular species that are
not nearest neighbors by a predefined max_dist. For instance, if you are
analyzing Li battery materials, you would expect that each Li+ would be
nearest neighbor to lower oxidation state transition metal for
electrostatic reasons. This only works if the structure is oxidation state
decorated, as structures with only elemental sites are automatically
assumed to have net charge of 0.


* **Parameters**


    * **sp1** ([*Species*](pymatgen.core.md#pymatgen.core.periodic_table.Species)) – First specie


    * **sp2** ([*Species*](pymatgen.core.md#pymatgen.core.periodic_table.Species)) – Second specie


    * **max_dist** (*float*) – Maximum distance between species.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### test(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Method to execute the test.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input structure to test



* **Returns**

    True if structure does not contain the two species are distances

        greater than max_dist.




* **Return type**

    bool


## pymatgen.alchemy.materials module

This module provides various representations of transformed structures. A
TransformedStructure is a structure that has been modified by undergoing a
series of transformations.


### _class_ TransformedStructure(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), transformations: list[[AbstractTransformation](pymatgen.transformations.md#pymatgen.transformations.transformation_abc.AbstractTransformation)] | None = None, history: list[[AbstractTransformation](pymatgen.transformations.md#pymatgen.transformations.transformation_abc.AbstractTransformation) | dict[str, Any]] | None = None, other_parameters: dict[str, Any] | None = None)
Bases: `MSONable`

Container object for new structures that include history of
transformations.

Each transformed structure is made up of a sequence of structures with
associated transformation history.

Initializes a transformed structure from a structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input structure


    * **transformations** (*list**[**Transformation**]*) – List of transformations to
    apply.


    * **history** (*list**[**Transformation**]*) – Previous history.


    * **other_parameters** (*dict*) – Additional parameters to be added.



#### append_filter(structure_filter: AbstractStructureFilter)
Adds a filter.


* **Parameters**

    **structure_filter** (*StructureFilter*) – A filter implementing the
    AbstractStructureFilter API. Tells transmuter what structures to retain.



#### append_transformation(transformation, return_alternatives: bool = False, clear_redo: bool = True)
Appends a transformation to the TransformedStructure.


* **Parameters**


    * **transformation** – Transformation to append


    * **return_alternatives** – Whether to return alternative
    TransformedStructures for one-to-many transformations.
    return_alternatives can be a number, which stipulates the
    total number of structures to return.


    * **clear_redo** – Boolean indicating whether to clear the redo list.
    By default, this is True, meaning any appends clears the
    history of undoing. However, when using append_transformation
    to do a redo, the redo list should not be cleared to allow
    multiple redos.



#### as_dict()
Dict representation of the TransformedStructure.


#### extend_transformations(transformations: list[[AbstractTransformation](pymatgen.transformations.md#pymatgen.transformations.transformation_abc.AbstractTransformation)], return_alternatives: bool = False)
Extends a sequence of transformations to the TransformedStructure.


* **Parameters**


    * **transformations** – Sequence of Transformations


    * **return_alternatives** – Whether to return alternative
    TransformedStructures for one-to-many transformations.
    return_alternatives can be a number, which stipulates the
    total number of structures to return.



#### _classmethod_ from_cif_str(cif_string: str, transformations: list[[AbstractTransformation](pymatgen.transformations.md#pymatgen.transformations.transformation_abc.AbstractTransformation)] | None = None, primitive: bool = True, occupancy_tolerance: float = 1.0)
Generates TransformedStructure from a cif string.


* **Parameters**


    * **cif_string** (*str*) – Input cif string. Should contain only one
    structure. For CIFs containing multiple structures, please use
    CifTransmuter.


    * **transformations** (*list**[**Transformation**]*) – Sequence of transformations
    to be applied to the input structure.


    * **primitive** (*bool*) – Option to set if the primitive cell should be
    extracted. Defaults to True. However, there are certain
    instances where you might want to use a non-primitive cell,
    e.g., if you are trying to generate all possible orderings of
    partial removals or order a disordered structure.


    * **occupancy_tolerance** (*float*) – If total occupancy of a site is
    between 1 and occupancy_tolerance, the occupancies will be
    scaled down to 1.



* **Returns**

    TransformedStructure



#### _classmethod_ from_cif_string(\*args, \*\*kwds)
from_cif_string is deprecated!
Use from_cif_str instead


#### _classmethod_ from_dict(d)
Creates a TransformedStructure from a dict.


#### _classmethod_ from_poscar_str(poscar_string: str, transformations: list[[AbstractTransformation](pymatgen.transformations.md#pymatgen.transformations.transformation_abc.AbstractTransformation)] | None = None)
Generates TransformedStructure from a poscar string.


* **Parameters**


    * **poscar_string** (*str*) – Input POSCAR string.


    * **transformations** (*list**[**Transformation**]*) – Sequence of transformations
    to be applied to the input structure.



#### _classmethod_ from_poscar_string(\*args, \*\*kwds)
from_poscar_string is deprecated!
Use from_poscar_str instead


#### _classmethod_ from_snl(snl: [StructureNL](pymatgen.util.md#pymatgen.util.provenance.StructureNL))
Create TransformedStructure from SNL.


* **Parameters**

    **snl** ([*StructureNL*](pymatgen.util.md#pymatgen.util.provenance.StructureNL)) – Starting snl



* **Returns**

    TransformedStructure



#### get_vasp_input(vasp_input_set: type[~pymatgen.io.vasp.sets.VaspInputSet] = <class 'pymatgen.io.vasp.sets.MPRelaxSet'>, \*\*kwargs)
Returns VASP input as a dict of VASP objects.


* **Parameters**


    * **vasp_input_set** ([*pymatgen.io.vasp.sets.VaspInputSet*](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.VaspInputSet)) – input set
    to create VASP input files from structures


    * **\*\*kwargs** – All keyword args supported by the VASP input set.



#### redo_next_change()
Redo the last undone change in the TransformedStructure.


* **Raises**

    **IndexError** – If already at the latest change.



#### set_parameter(key: str, value: Any)
Set a parameter.


* **Parameters**


    * **key** – The string key


    * **value** – The value.



#### _property_ structures(_: list[[Structure](pymatgen.core.md#pymatgen.core.structure.Structure)_ )
Copy of all structures in the TransformedStructure. A
structure is stored after every single transformation.


#### to_snl(authors, \*\*kwargs)
Generate SNL from TransformedStructure.


* **Parameters**


    * **authors** – List of authors


    * **\*\*kwargs** – All kwargs supported by StructureNL.




* **Returns**

    StructureNL



#### undo_last_change()
Undo the last change in the TransformedStructure.


* **Raises**

    **IndexError** – If already at the oldest change.



#### _property_ was_modified(_: boo_ )
Boolean describing whether the last transformation on the structure
made any alterations to it one example of when this would return false
is in the case of performing a substitution transformation on the
structure when the specie to replace isn’t in the structure.


#### write_vasp_input(vasp_input_set: type[~pymatgen.io.vasp.sets.VaspInputSet] = <class 'pymatgen.io.vasp.sets.MPRelaxSet'>, output_dir: str = '.', create_directory: bool = True, \*\*kwargs)
Writes VASP input to an output_dir.


* **Parameters**


    * **vasp_input_set** – pymatgen.io.vasp.sets.VaspInputSet like object that creates vasp input files from
    structures.


    * **output_dir** – Directory to output files


    * **create_directory** – Create the directory if not present. Defaults to
    True.


    * **\*\*kwargs** – All keyword args supported by the VASP input set.


## pymatgen.alchemy.transmuters module

This module implements various transmuter classes.
Transmuters are essentially classes that generate TransformedStructures from
various data sources. They enable the high-throughput generation of new
structures and input files.

It also includes the helper function, batch_write_vasp_input to generate an
entire directory of vasp input files for running.


### _class_ CifTransmuter(cif_string, transformations=None, primitive=True, extend_collection=False)
Bases: `StandardTransmuter`

Generates a Transmuter from a cif string, possibly containing multiple
structures.

Generates a Transmuter from a cif string, possibly
containing multiple structures.


* **Parameters**


    * **cif_string** – A string containing a cif or a series of CIFs


    * **transformations** – New transformations to be applied to all
    structures


    * **primitive** – Whether to generate the primitive cell from the cif.


    * **extend_collection** – Whether to use more than one output structure
    from one-to-many transformations. extend_collection can be a
    number, which determines the maximum branching for each
    transformation.



#### _classmethod_ from_filenames(filenames, transformations=None, primitive=True, extend_collection=False)
Generates a TransformedStructureCollection from a cif, possibly
containing multiple structures.


* **Parameters**


    * **filenames** – List of strings of the cif files


    * **transformations** – New transformations to be applied to all
    structures


    * **primitive** – Same meaning as in __init__.


    * **extend_collection** – Same meaning as in __init__.



### _class_ PoscarTransmuter(poscar_string, transformations=None, extend_collection=False)
Bases: `StandardTransmuter`

Generates a transmuter from a sequence of POSCARs.


* **Parameters**


    * **poscar_string** – List of POSCAR strings


    * **transformations** – New transformations to be applied to all
    structures.


    * **extend_collection** – Whether to use more than one output structure
    from one-to-many transformations.



#### _static_ from_filenames(poscar_filenames, transformations=None, extend_collection=False)
Convenient constructor to generates a POSCAR transmuter from a list of
POSCAR filenames.


* **Parameters**


    * **poscar_filenames** – List of POSCAR filenames


    * **transformations** – New transformations to be applied to all
    structures.


    * **extend_collection** – Same meaning as in __init__.



### _class_ StandardTransmuter(transformed_structures, transformations=None, extend_collection: int = 0, ncores: int | None = None)
Bases: `object`

An example of a Transmuter object, which performs a sequence of
transformations on many structures to generate TransformedStructures.


#### transformed_structures()
List of all transformed structures.


* **Type**

    list[[Structure](pymatgen.core.md#pymatgen.core.structure.Structure)]


Initializes a transmuter from an initial list of
pymatgen.alchemy.materials.TransformedStructure.


* **Parameters**


    * **transformed_structures** (*[**TransformedStructure**]*) – Input transformed
    structures


    * **transformations** (*[**Transformations**]*) – New transformations to be
    applied to all structures.


    * **extend_collection** (*int*) – Whether to use more than one output
    structure from one-to-many transformations. extend_collection
    can be an int, which determines the maximum branching for each
    transformation.


    * **ncores** (*int*) – Number of cores to use for applying transformations.
    Uses multiprocessing.Pool. Default is None, which implies
    serial.



#### add_tags(tags)
Add tags for the structures generated by the transmuter.


* **Parameters**

    **tags** – A sequence of tags. Note that this should be a sequence of
    strings, e.g., [“My awesome structures”, “Project X”].



#### append_transformation(transformation, extend_collection=False, clear_redo=True)
Appends a transformation to all TransformedStructures.


* **Parameters**


    * **transformation** – Transformation to append


    * **extend_collection** – Whether to use more than one output structure
    from one-to-many transformations. extend_collection can be a
    number, which determines the maximum branching for each
    transformation.


    * **clear_redo** (*bool*) – Whether to clear the redo list. By default,
    this is True, meaning any appends clears the history of
    undoing. However, when using append_transformation to do a
    redo, the redo list should not be cleared to allow multiple
    redos.



* **Returns**

    List of booleans corresponding to initial transformed structures
    each boolean describes whether the transformation altered the
    structure



#### append_transformed_structures(trafo_structs_or_transmuter)
Method is overloaded to accept either a list of transformed structures
or transmuter, it which case it appends the second transmuter”s
structures.


* **Parameters**

    **trafo_structs_or_transmuter** – A list of transformed structures or a
    transmuter.



#### apply_filter(structure_filter)
Applies a structure_filter to the list of TransformedStructures
in the transmuter.


* **Parameters**

    **structure_filter** – StructureFilter to apply.



#### extend_transformations(transformations)
Extends a sequence of transformations to the TransformedStructure.


* **Parameters**

    **transformations** – Sequence of Transformations



#### _classmethod_ from_structures(structures, transformations=None, extend_collection=0)
Alternative constructor from structures rather than
TransformedStructures.


* **Parameters**


    * **structures** – Sequence of structures


    * **transformations** – New transformations to be applied to all
    structures


    * **extend_collection** – Whether to use more than one output structure
    from one-to-many transformations. extend_collection can be a
    number, which determines the maximum branching for each
    transformation.



* **Returns**

    StandardTransmuter



#### redo_next_change()
Redo the last undone transformation in the TransformedStructure.


* **Raises**

    **IndexError if already at the latest change.** –



#### set_parameter(key, value)
Add parameters to the transmuter. Additional parameters are stored in
the as_dict() output.


* **Parameters**


    * **key** – The key for the parameter.


    * **value** – The value for the parameter.



#### undo_last_change()
Undo the last transformation in the TransformedStructure.


* **Raises**

    **IndexError if already at the oldest change.** –



#### write_vasp_input(\*\*kwargs)
Batch write vasp input for a sequence of transformed structures to
output_dir, following the format output_dir/{formula}_{number}.


* **Parameters**

    **kwargs** – All kwargs supported by batch_write_vasp_input.



### _apply_transformation(inputs)
Helper method for multiprocessing of apply_transformation. Must not be
in the class so that it can be pickled.


* **Parameters**

    **inputs** – Tuple containing the transformed structure, the transformation
    to be applied, a boolean indicating whether to extend the
    collection, and a boolean indicating whether to clear the redo



* **Returns**

    List of output structures (the modified initial structure, plus
    any new structures created by a one-to-many transformation)



### batch_write_vasp_input(transformed_structures: Sequence[TransformedStructure], vasp_input_set: type[VaspInputSet] = <class 'pymatgen.io.vasp.sets.MPRelaxSet'>, output_dir: str = '.', create_directory: bool = True, subfolder: Callable[[TransformedStructure], str] | None = None, include_cif: bool = False, \*\*kwargs)
Batch write vasp input for a sequence of transformed structures to
output_dir, following the format output_dir/{group}/{formula}_{number}.


* **Parameters**


    * **transformed_structures** – Sequence of TransformedStructures.


    * **vasp_input_set** – pymatgen.io.vasp.sets.VaspInputSet to creates
    vasp input files from structures.


    * **output_dir** – Directory to output files


    * **create_directory** (*bool*) – Create the directory if not present.
    Defaults to True.


    * **subfolder** – Function to create subdirectory name from
    transformed_structure.
    e.g., lambda x: x.other_parameters[“tags”][0] to use the first
    tag.


    * **include_cif** (*bool*) – Boolean indication whether to output a CIF as
    well. CIF files are generally better supported in visualization
    programs.


    * **\*\*kwargs** – Any kwargs supported by vasp_input_set.