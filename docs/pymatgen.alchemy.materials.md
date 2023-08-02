---
layout: default
title: pymatgen.alchemy.materials.md
nav_exclude: true
---

# pymatgen.alchemy.materials module

This module provides various representations of transformed structures. A
TransformedStructure is a structure that has been modified by undergoing a
series of transformations.


### _class_ pymatgen.alchemy.materials.TransformedStructure(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), transformations: list[[AbstractTransformation](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)] | None = None, history: list[[AbstractTransformation](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation) | dict[str, Any]] | None = None, other_parameters: dict[str, Any] | None = None)
Bases: `MSONable`

Container object for new structures that include history of
transformations.

Each transformed structure is made up of a sequence of structures with
associated transformation history.

Initializes a transformed structure from a structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input structure


    * **transformations** (*list**[**Transformation**]*) – List of transformations to
    apply.


    * **history** (*list**[**Transformation**]*) – Previous history.


    * **other_parameters** (*dict*) – Additional parameters to be added.



#### append_filter(structure_filter: [AbstractStructureFilter](pymatgen.alchemy.filters.md#pymatgen.alchemy.filters.AbstractStructureFilter))
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


#### extend_transformations(transformations: list[[AbstractTransformation](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)], return_alternatives: bool = False)
Extends a sequence of transformations to the TransformedStructure.


* **Parameters**


    * **transformations** – Sequence of Transformations


    * **return_alternatives** – Whether to return alternative
    TransformedStructures for one-to-many transformations.
    return_alternatives can be a number, which stipulates the
    total number of structures to return.



#### _static_ from_cif_string(cif_string: str, transformations: list[[AbstractTransformation](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)] | None = None, primitive: bool = True, occupancy_tolerance: float = 1.0)
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



#### _classmethod_ from_dict(d)
Creates a TransformedStructure from a dict.


#### _static_ from_poscar_string(poscar_string: str, transformations: list[[AbstractTransformation](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)] | None = None)
Generates TransformedStructure from a poscar string.


* **Parameters**


    * **poscar_string** (*str*) – Input POSCAR string.


    * **transformations** (*list**[**Transformation**]*) – Sequence of transformations
    to be applied to the input structure.



#### _classmethod_ from_snl(snl: [StructureNL](pymatgen.util.provenance.md#pymatgen.util.provenance.StructureNL))
Create TransformedStructure from SNL.


* **Parameters**

    **snl** ([*StructureNL*](pymatgen.util.provenance.md#pymatgen.util.provenance.StructureNL)) – Starting snl



* **Returns**

    TransformedStructure



#### get_vasp_input(vasp_input_set: type[pymatgen.io.vasp.sets.VaspInputSet] = <class 'pymatgen.io.vasp.sets.MPRelaxSet'>, \*\*kwargs)
Returns VASP input as a dict of VASP objects.


* **Parameters**


    * **vasp_input_set** ([*pymatgen.io.vasp.sets.VaspInputSet*](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.VaspInputSet)) – input set
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



#### _property_ structures(_: list[[pymatgen.core.structure.Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure)_ )
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


#### write_vasp_input(vasp_input_set: type[pymatgen.io.vasp.sets.VaspInputSet] = <class 'pymatgen.io.vasp.sets.MPRelaxSet'>, output_dir: str = '.', create_directory: bool = True, \*\*kwargs)
Writes VASP input to an output_dir.


* **Parameters**


    * **vasp_input_set** – pymatgen.io.vasp.sets.VaspInputSet like object that creates vasp input files from
    structures.


    * **output_dir** – Directory to output files


    * **create_directory** – Create the directory if not present. Defaults to
    True.


    * **\*\*kwargs** – All keyword args supported by the VASP input set.