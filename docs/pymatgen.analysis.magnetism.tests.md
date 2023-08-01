---
layout: default
title: pymatgen.analysis.magnetism.tests.md
nav_exclude: true
---

# pymatgen.analysis.magnetism.tests namespace


## pymatgen.analysis.magnetism.tests.test_analyzer module


### _class_ pymatgen.analysis.magnetism.tests.test_analyzer.CollinearMagneticStructureAnalyzerTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_get_ferromagnetic_structure()

#### test_get_representations()

#### test_magnetic_properties()

#### test_matches()

#### test_missing_spin()

#### test_modes()

#### test_net_positive()

#### test_round_magmoms()

#### test_str()

### _class_ pymatgen.analysis.magnetism.tests.test_analyzer.MagneticDeformationTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_magnetic_deformation()

### _class_ pymatgen.analysis.magnetism.tests.test_analyzer.MagneticStructureEnumeratorTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_ordering_enumeration()
## pymatgen.analysis.magnetism.tests.test_heisenberg module


### _class_ pymatgen.analysis.magnetism.tests.test_heisenberg.HeisenbergMapperTest(methodName='runTest')
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


#### test_exchange_params()

#### test_get_igraph()

#### test_graphs()

#### test_heisenberg_model()

#### test_mean_field()

#### test_nn_interactions()

#### test_sites()
## pymatgen.analysis.magnetism.tests.test_jahnteller module


### _class_ pymatgen.analysis.magnetism.tests.test_jahnteller.JahnTellerTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_jahn_teller_species_analysis()

#### test_jahn_teller_structure_analysis()

#### test_mu_so()