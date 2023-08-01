---
layout: default
title: pymatgen.io.xtb.tests.md
nav_exclude: true
---

# pymatgen.io.xtb.tests package

This package implements tests for CREST outputs.


## pymatgen.io.xtb.tests.test_inputs module


### _class_ pymatgen.io.xtb.tests.test_inputs.TestCRESTInput(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Checks that all attributes of CRESTInput match the expected values for
sample inputs

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_constraints_file()

#### test_coordinates_file()
## pymatgen.io.xtb.tests.test_outputs module


### _class_ pymatgen.io.xtb.tests.test_outputs.TestCRESTOutput(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Checks that all attributes of CRESTOutput match the expected values for a
sample CREST output directory.

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_all()