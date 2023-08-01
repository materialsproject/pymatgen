---
layout: default
title: pymatgen.io.lammps.tests.md
nav_exclude: true
---

# pymatgen.io.lammps.tests package


## pymatgen.io.lammps.tests.test_data module


### _class_ pymatgen.io.lammps.tests.test_data.CombinedDataTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### _classmethod_ setUpClass()
Hook method for setting up class fixture before running tests in the class.


#### test_as_lammpsdata()

#### test_disassemble()

#### test_from_ff_and_topologies()

#### test_from_files()

#### test_from_lammpsdata()

#### test_from_structure()

#### test_get_string()

#### test_json_dict()

#### test_structure()

### _class_ pymatgen.io.lammps.tests.test_data.ForceFieldTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### _classmethod_ setUpClass()
Hook method for setting up class fixture before running tests in the class.


#### _classmethod_ tearDownClass()
Hook method for deconstructing the class fixture after running all tests in the class.


#### test_from_dict()

#### test_from_file()

#### test_init()

#### test_to_file()

### _class_ pymatgen.io.lammps.tests.test_data.FuncTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_lattice_2_lmpbox()

### _class_ pymatgen.io.lammps.tests.test_data.LammpsBoxTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### _classmethod_ setUpClass()
Hook method for setting up class fixture before running tests in the class.


#### test_get_box_shift()

#### test_get_string()

#### test_to_lattice()

#### test_volume()

### _class_ pymatgen.io.lammps.tests.test_data.LammpsDataTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### _classmethod_ setUpClass()
Hook method for setting up class fixture before running tests in the class.


#### _classmethod_ tearDownClass()
Hook method for deconstructing the class fixture after running all tests in the class.


#### test_disassemble()

#### test_from_ff_and_topologies()

#### test_from_file()

#### test_from_structure()

#### test_get_string()

#### test_json_dict()

#### test_set_charge_atom()

#### test_set_charge_atom_type()

#### test_sort_structure()

#### test_structure()

#### test_write_file()

### _class_ pymatgen.io.lammps.tests.test_data.TopologyTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_from_bonding()

#### test_init()
## pymatgen.io.lammps.tests.test_generators module


### _class_ pymatgen.io.lammps.tests.test_generators.LammpsMinimizationTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### _classmethod_ setUpClass()
Hook method for setting up class fixture before running tests in the class.


#### test_get_input_set()
## pymatgen.io.lammps.tests.test_inputs module


### _class_ pymatgen.io.lammps.tests.test_inputs.FuncTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### _classmethod_ tearDownClass()
Hook method for deconstructing the class fixture after running all tests in the class.


#### test_write_lammps_inputs()

### _class_ pymatgen.io.lammps.tests.test_inputs.LammpsInputFileTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### _classmethod_ setUpClass()
Hook method for setting up class fixture before running tests in the class.


#### test_add_command()

#### test_add_commands()

#### test_add_comment()

#### test_add_stage()

#### test_append()

#### test_contains_command()

#### test_from_file()

#### test_from_string()

#### test_get_args()

#### test_get_string()

#### test_merge_stages()

#### test_ncomments()

#### test_nstages()

#### test_remove_command()

#### test_remove_stage()

#### test_rename_stage()

#### test_set_args()

#### test_stages_names()

### _class_ pymatgen.io.lammps.tests.test_inputs.LammpsRunTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### maxDiff(_ = Non_ )

#### _classmethod_ tearDownClass()
Hook method for deconstructing the class fixture after running all tests in the class.


#### test_md()

### _class_ pymatgen.io.lammps.tests.test_inputs.LammpsTemplateGenTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_write_inputs()
## pymatgen.io.lammps.tests.test_outputs module


### _class_ pymatgen.io.lammps.tests.test_outputs.FuncTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_parse_lammps_dumps()

#### test_parse_lammps_log()

### _class_ pymatgen.io.lammps.tests.test_outputs.LammpsDumpTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### _classmethod_ setUpClass()
Hook method for setting up class fixture before running tests in the class.


#### test_from_string()

#### test_json_dict()
## pymatgen.io.lammps.tests.test_utils module


### _class_ pymatgen.io.lammps.tests.test_utils.TestPackmolOutput(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### _classmethod_ setUpClass()
Hook method for setting up class fixture before running tests in the class.


#### test_packed_molecule()

### _class_ pymatgen.io.lammps.tests.test_utils.TestPolymer(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### _classmethod_ setUpClass()
Hook method for setting up class fixture before running tests in the class.


#### test_polymer_chain_lengths()

#### test_polymer_chain_topologies()