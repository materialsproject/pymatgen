---
layout: default
title: pymatgen.analysis.elasticity.tests.md
nav_exclude: true
---

# pymatgen.analysis.elasticity.tests namespace


## pymatgen.analysis.elasticity.tests.test_elastic module


### _class_ pymatgen.analysis.elasticity.tests.test_elastic.DiffFitTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Tests various functions related to diff fitting

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_find_eq_stress()

#### test_fit()

#### test_generate_pseudo()

#### test_get_diff_coeff()

#### test_get_strain_state_dict()

### _class_ pymatgen.analysis.elasticity.tests.test_elastic.ElasticTensorExpansionTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_calculate_stress()

#### test_energy_density()

#### test_from_diff_fit()

#### test_get_compliance_expansion()

#### test_get_effective_ecs()

#### test_get_strain_from_stress()

#### test_get_yield_stress()

#### test_gruneisen()

#### test_init()

#### test_thermal_expansion_coeff()

### _class_ pymatgen.analysis.elasticity.tests.test_elastic.ElasticTensorTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_compliance_tensor()

#### test_directional_elastic_mod()

#### test_directional_poisson_ratio()

#### test_energy_density()

#### test_from_independent_strains()

#### test_from_pseudoinverse()

#### test_new()

#### test_properties()

#### test_structure_based_methods()

### _class_ pymatgen.analysis.elasticity.tests.test_elastic.NthOrderElasticTensorTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_calculate_stress()

#### test_energy_density()

#### test_from_diff_fit()

#### test_init()
## pymatgen.analysis.elasticity.tests.test_strain module


### _class_ pymatgen.analysis.elasticity.tests.test_strain.DeformationTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_apply_to_structure()

#### test_independence()

#### test_properties()

### _class_ pymatgen.analysis.elasticity.tests.test_strain.DeformedStructureSetTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_init()

### _class_ pymatgen.analysis.elasticity.tests.test_strain.StrainTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_convert_strain_to_deformation()

#### test_from_deformation()

#### test_from_index_amount()

#### test_new()

#### test_properties()
## pymatgen.analysis.elasticity.tests.test_stress module


### _class_ pymatgen.analysis.elasticity.tests.test_stress.StressTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_properties()