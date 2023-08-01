---
layout: default
title: pymatgen.command_line.tests.md
nav_exclude: true
---

# pymatgen.command_line.tests package

Command line tests.


## pymatgen.command_line.tests.test_bader_caller module


### _class_ pymatgen.command_line.tests.test_bader_caller.BaderAnalysisTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_atom_parsing()

#### test_automatic_runner()

#### test_from_path()

#### test_init()

#### test_missing_file_bader_exe_path()
## pymatgen.command_line.tests.test_chargemol_caller module


### _class_ pymatgen.command_line.tests.test_chargemol_caller.ChargemolAnalysisTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_parse_chargemol()

#### test_parse_chargemol2()
## pymatgen.command_line.tests.test_critic2_caller module


### _class_ pymatgen.command_line.tests.test_critic2_caller.Critic2AnalysisTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_graph_output()

#### test_properties_to_from_dict()

### _class_ pymatgen.command_line.tests.test_critic2_caller.Critic2CallerTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_from_path()

#### test_from_structure()
## pymatgen.command_line.tests.test_enumlib_caller module


### _class_ pymatgen.command_line.tests.test_enumlib_caller.EnumlibAdaptorTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_init()

#### test_partial_disorder()

#### test_rounding_errors()

#### test_timeout()
## pymatgen.command_line.tests.test_gulp_caller module

Created on Jan 22, 2013.

@author: Bharat Medasani


### _class_ pymatgen.command_line.tests.test_gulp_caller.BuckinghamPotentialBushTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_element_different_valence()

#### test_existing_element()

#### test_non_exisitng_element()

#### test_spring()

### _class_ pymatgen.command_line.tests.test_gulp_caller.BuckinghamPotentialLewisTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_element_different_valence()

#### test_existing_element()

#### test_non_exisitng_element()

#### test_spring()

#### test_values()

### _class_ pymatgen.command_line.tests.test_gulp_caller.GlobalFunctionsTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_get_energy_buckingham()

#### test_get_energy_relax_structure_buckingham()

#### test_get_energy_tersoff()

### _class_ pymatgen.command_line.tests.test_gulp_caller.GulpCallerTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_decimal()

#### test_run()

### _class_ pymatgen.command_line.tests.test_gulp_caller.GulpIOTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_buckingham_input()

#### test_buckingham_potential()

#### test_get_energy()

#### test_get_relaxed_structure()

#### test_keyword_line_with_correct_keywords()

#### test_library_line_explicit_path()

#### test_library_line_wrong_file()

#### test_specie_potential()

#### test_structure_lines_default_options()

#### test_structure_lines_no_frac_coords()

#### test_structure_lines_no_unitcell()

#### test_tersoff_inpt()

#### test_tersoff_potential()
## pymatgen.command_line.tests.test_mcsqs_caller module


### _class_ pymatgen.command_line.tests.test_mcsqs_caller.McsqsCallerTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_mcsqs_caller_parallel()

#### test_mcsqs_caller_runtime_error()

#### test_mcsqs_caller_supercell()

#### test_mcsqs_caller_total_atoms()

#### test_mcsqs_caller_total_atoms_auto_instances()

#### test_mcsqs_perfect_match_error()

#### test_mcsqs_perfect_match_error_parallel()
## pymatgen.command_line.tests.test_vampire_caller module


### _class_ pymatgen.command_line.tests.test_vampire_caller.VampireCallerTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### _classmethod_ setUpClass()
Hook method for setting up class fixture before running tests in the class.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_vampire()