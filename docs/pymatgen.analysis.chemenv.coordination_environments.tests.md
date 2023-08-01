---
layout: default
title: pymatgen.analysis.chemenv.coordination_environments.tests.md
nav_exclude: true
---

# pymatgen.analysis.chemenv.coordination_environments.tests package

Test package.


## pymatgen.analysis.chemenv.coordination_environments.tests.test_chemenv_strategies module


### _class_ pymatgen.analysis.chemenv.coordination_environments.tests.test_chemenv_strategies.StrategyOptionsTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_options()

#### test_strategies()
## pymatgen.analysis.chemenv.coordination_environments.tests.test_coordination_geometries module


### _class_ pymatgen.analysis.chemenv.coordination_environments.tests.test_coordination_geometries.CoordinationGeometriesTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_algorithms()

#### test_coordination_geometry()

#### test_hints()

### _class_ pymatgen.analysis.chemenv.coordination_environments.tests.test_coordination_geometries.FakeSite(coords)
Bases: `object`

## pymatgen.analysis.chemenv.coordination_environments.tests.test_coordination_geometry_finder module


### _class_ pymatgen.analysis.chemenv.coordination_environments.tests.test_coordination_geometry_finder.CoordinationGeometryFinderTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_abstract_geometry()

#### test_disable_hints()

#### test_perfect_environments()
## pymatgen.analysis.chemenv.coordination_environments.tests.test_read_write module


### _class_ pymatgen.analysis.chemenv.coordination_environments.tests.test_read_write.ReadWriteChemenvTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### _classmethod_ setUpClass()
Hook method for setting up class fixture before running tests in the class.


#### _classmethod_ tearDownClass()
Hook method for deconstructing the class fixture after running all tests in the class.


#### test_read_write_structure_environments()

#### test_read_write_voronoi()

#### test_strategies()

#### test_structure_environments_neighbors_sets()
## pymatgen.analysis.chemenv.coordination_environments.tests.test_structure_environments module


### _class_ pymatgen.analysis.chemenv.coordination_environments.tests.test_structure_environments.StructureEnvironmentsTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_from_structure_environments()

#### test_light_structure_environments()

#### test_structure_environments()
## pymatgen.analysis.chemenv.coordination_environments.tests.test_voronoi module


### _class_ pymatgen.analysis.chemenv.coordination_environments.tests.test_voronoi.VoronoiContainerTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_get_vertices_dist_ang_indices()

#### test_voronoi()
## pymatgen.analysis.chemenv.coordination_environments.tests.test_weights module


### _class_ pymatgen.analysis.chemenv.coordination_environments.tests.test_weights.DummyStructureEnvironments()
Bases: `object`


### _class_ pymatgen.analysis.chemenv.coordination_environments.tests.test_weights.DummyVoronoiContainer()
Bases: `object`


### _class_ pymatgen.analysis.chemenv.coordination_environments.tests.test_weights.FakeNbSet(cn=None)
Bases: `object`


### _class_ pymatgen.analysis.chemenv.coordination_environments.tests.test_weights.StrategyWeightsTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_CN_bias_weight()

#### test_angle_weight()

#### test_delta_csms_weight()

#### test_dist_angle_area_weight()

#### test_dist_nb_set_weight()

#### test_normalized_angle_distance_weight()

#### test_self_csms_weight()