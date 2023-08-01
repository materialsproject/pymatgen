---
layout: default
title: pymatgen.io.feff.tests.md
nav_exclude: true
---

# pymatgen.io.feff.tests package


## pymatgen.io.feff.tests.test_inputs module


### _class_ pymatgen.io.feff.tests.test_inputs.FeffAtomsTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### _classmethod_ setUpClass()
Hook method for setting up class fixture before running tests in the class.


#### test_absorber_line()

#### test_absorbing_atom()

#### test_as_dict_and_from_dict()

#### test_atoms_from_file()

#### test_cluster_from_file()

#### test_distances()

#### test_get_string()

#### test_single_absorbing_atom()
When there is only one absorbing atom in the structure, it should not appear
in the pot_dict to avoid an error


### _class_ pymatgen.io.feff.tests.test_inputs.FeffPotTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_as_dict_and_from_dict()

#### test_init()

#### test_single_absorbing_atom()
When there is only one absorbing atom in the structure, it should not appear
in the pot_dict to avoid an error


### _class_ pymatgen.io.feff.tests.test_inputs.FeffTagsTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_as_dict_and_from_dict()

#### test_diff()

#### test_eels_tags()

#### test_init()

### _class_ pymatgen.io.feff.tests.test_inputs.HeaderTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_as_dict_and_from_dict()

#### test_from_string()

#### test_get_string()

#### test_init()

### _class_ pymatgen.io.feff.tests.test_inputs.PathsTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_paths_string()
## pymatgen.io.feff.tests.test_outputs module


### _class_ pymatgen.io.feff.tests.test_outputs.FeffLdosTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### filepath1(_ = '/Users/shyue/repos/pymatgen/pymatgen/util/../../test_files/feff.inp_ )

#### filepath2(_ = '/Users/shyue/repos/pymatgen/pymatgen/util/../../test_files/ldos_ )

#### ldos(_ = <pymatgen.io.feff.outputs.LDos object_ )

#### reci_dos(_ = <pymatgen.io.feff.outputs.LDos object_ )

#### reci_feffinp(_ = '/Users/shyue/repos/pymatgen/pymatgen/util/../../test_files/feff_reci_dos/feff.inp_ )

#### reci_ldos(_ = '/Users/shyue/repos/pymatgen/pymatgen/util/../../test_files/feff_reci_dos/ldos_ )

#### test_as_dict_and_from_dict()

#### test_complete_dos()

#### test_init()

#### test_reci_charge()

#### test_reci_complete_dos()

#### test_reci_init()

### _class_ pymatgen.io.feff.tests.test_outputs.XmuTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_as_dict_and_from_dict()

#### test_init()
## pymatgen.io.feff.tests.test_sets module


### _class_ pymatgen.io.feff.tests.test_sets.FeffInputSetTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### _classmethod_ setUpClass()
Hook method for setting up class fixture before running tests in the class.


#### test_charged_structure()

#### test_eels_tags_set()

#### test_eels_to_from_dict()

#### test_get_feffPot()

#### test_get_feff_atoms()

#### test_get_header()

#### test_getfefftags()

#### test_large_systems()

#### test_number_of_kpoints()

#### test_post_distdiff()

#### test_postfeffset()

#### test_reciprocal_tags_and_input()

#### test_small_system_EXAFS()

#### test_to_and_from_dict()

#### test_user_tag_settings()