---
layout: default
title: pymatgen.alchemy.filters.md
nav_exclude: true
---

# pymatgen.alchemy.filters module

This module defines filters for Transmuter object.


### _class_ pymatgen.alchemy.filters.AbstractStructureFilter()
Bases: `MSONable`

AbstractStructureFilter that defines an API to perform testing of
Structures. Structures that return True to a test are retained during
transmutation while those that return False are removed.


#### _abstract_ test(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Method to execute the test.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input structure to test



* **Returns**

    (bool) Structures that return true are kept in the Transmuter
    object during filtering.



### _class_ pymatgen.alchemy.filters.ChargeBalanceFilter()
Bases: `AbstractStructureFilter`

This filter removes structures that are not charge balanced from the
transmuter. This only works if the structure is oxidation state
decorated, as structures with only elemental sites are automatically
assumed to have net charge of 0.

No args required.


#### test(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Method to execute the test.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input structure to test


Returns: True if structure is neutral.


### _class_ pymatgen.alchemy.filters.ContainsSpecieFilter(species, strict_compare=False, AND=True, exclude=False)
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



#### as_dict()
Returns: MSONable dict.


#### _classmethod_ from_dict(dct)

* **Parameters**

    **dct** (*dict*) – Dict representation.



* **Returns**

    Filter



#### test(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Method to execute the test.

Returns: True if structure do not contain specified species.


### _class_ pymatgen.alchemy.filters.RemoveDuplicatesFilter(structure_matcher: dict | [StructureMatcher](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.StructureMatcher) | None = None, symprec: float | None = None)
Bases: `AbstractStructureFilter`

This filter removes exact duplicate structures from the transmuter.

Remove duplicate structures based on the structure matcher
and symmetry (if symprec is given).


* **Parameters**


    * **structure_matcher** (*dict** | *[*StructureMatcher*](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.StructureMatcher)*, **optional*) – Provides a structure matcher to be used for
    structure comparison.


    * **symprec** (*float**, **optional*) – The precision in the symmetry finder algorithm if None (
    default value), no symmetry check is performed and only the
    structure matcher is used. A recommended value is 1e-5.



#### test(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))

* **Parameters**

    **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input structure to test.


Returns: True if structure is not in list.


### _class_ pymatgen.alchemy.filters.RemoveExistingFilter(existing_structures, structure_matcher=None, symprec=None)
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



#### as_dict()
Returns: MSONable dict.


#### test(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Method to execute the test.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input structure to test


Returns: True if structure is not in existing list.


### _class_ pymatgen.alchemy.filters.SpecieProximityFilter(specie_and_min_dist_dict)
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



#### as_dict()
Returns: MSONable dict.


#### _classmethod_ from_dict(dct)

* **Parameters**

    **dct** (*dict*) – Dict representation.



* **Returns**

    Filter



#### test(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Method to execute the test.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input structure to test


Returns: True if structure does not contain species within specified

    distances.


### _class_ pymatgen.alchemy.filters.SpeciesMaxDistFilter(sp1, sp2, max_dist)
Bases: `AbstractStructureFilter`

This filter removes structures that do have two particular species that are
not nearest neighbors by a predefined max_dist. For instance, if you are
analyzing Li battery materials, you would expect that each Li+ would be
nearest neighbor to lower oxidation state transition metal for
electrostatic reasons. This only works if the structure is oxidation state
decorated, as structures with only elemental sites are automatically
assumed to have net charge of 0.


* **Parameters**


    * **sp1** ([*Species*](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species)) – First specie


    * **sp2** ([*Species*](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species)) – Second specie


    * **max_dist** (*float*) – Maximum distance between species.



#### test(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Method to execute the test.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input structure to test


Returns: True if structure does not contain the two species are distances

    greater than max_dist.