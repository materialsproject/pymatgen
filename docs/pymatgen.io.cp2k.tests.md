---
layout: default
title: pymatgen.io.cp2k.tests.md
nav_exclude: true
---

# pymatgen.io.cp2k.tests package


## pymatgen.io.cp2k.tests.test_inputs module


### _class_ pymatgen.io.cp2k.tests.test_inputs.BasisAndPotentialTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_basis()

#### test_basis_info()

#### test_potential_info()

#### test_potentials()

### _class_ pymatgen.io.cp2k.tests.test_inputs.InputTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_basic_keywords()

#### test_basic_sections()

#### test_ci_file()

#### test_coords()

#### test_kind()

#### test_mongo()

#### test_odd_file()

#### test_preprocessor()

#### test_sectionlist()
## pymatgen.io.cp2k.tests.test_outputs module


### _class_ pymatgen.io.cp2k.tests.test_outputs.SetTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### energy_force()
Can get energy and forces


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_band()
Can parse bandstructure files


#### test_chi()

#### test_dos()
Can parse dos files


#### test_files()
Can find files successfully


#### test_gtensor()

#### test_hyperfine()

#### test_run_info()
Can extract run info from out file

## pymatgen.io.cp2k.tests.test_sets module


### _class_ pymatgen.io.cp2k.tests.test_sets.SetTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_dft_set()