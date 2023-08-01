---
layout: default
title: pymatgen.alchemy.tests.md
nav_exclude: true
---

# pymatgen.alchemy.tests package


## pymatgen.alchemy.tests.test_filters module


### _class_ pymatgen.alchemy.tests.test_filters.ContainsSpecieFilterTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_filtering()

#### test_to_from_dict()

### _class_ pymatgen.alchemy.tests.test_filters.RemoveDuplicatesFilterTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_filter()

#### test_to_from_dict()

### _class_ pymatgen.alchemy.tests.test_filters.RemoveExistingFilterTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_filter()

### _class_ pymatgen.alchemy.tests.test_filters.SpecieProximityFilterTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_filter()

#### test_to_from_dict()
## pymatgen.alchemy.tests.test_materials module


### _class_ pymatgen.alchemy.tests.test_materials.TransformedStructureTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_append_filter()

#### test_append_transformation()

#### test_as_dict()

#### test_final_structure()

#### test_from_dict()

#### test_get_vasp_input()

#### test_snl()

#### test_undo_and_redo_last_change()
## pymatgen.alchemy.tests.test_transmuters module


### _class_ pymatgen.alchemy.tests.test_transmuters.CifTransmuterTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_init()

### _class_ pymatgen.alchemy.tests.test_transmuters.PoscarTransmuterTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_init()

#### test_transmuter()