---
layout: default
title: pymatgen.analysis.chemenv.utils.tests.md
nav_exclude: true
---

# pymatgen.analysis.chemenv.utils.tests package

Test package.


## pymatgen.analysis.chemenv.utils.tests.test_chemenv_config module


### _class_ pymatgen.analysis.chemenv.utils.tests.test_chemenv_config.ChemenvConfigTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_chemenv_config()
## pymatgen.analysis.chemenv.utils.tests.test_coordination_geometry_utils module


### _class_ pymatgen.analysis.chemenv.utils.tests.test_coordination_geometry_utils.PlanesUtilsTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_distances()

#### test_factors_abcd_normal_vector()

#### test_from_npoints_plane()

#### test_indices_separate()

#### test_is_in_plane()

#### test_normal_vector_is_normed()

#### test_orthonormal_vectors()

#### test_plane_2_coefficients()

#### test_plane_3_coefficients()

#### test_plane_comparison()

#### test_plane_is_in_list_of_planes()

#### test_projections()
## pymatgen.analysis.chemenv.utils.tests.test_func_utils module


### _class_ pymatgen.analysis.chemenv.utils.tests.test_func_utils.FuncUtilsTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_CSMFiniteRatioFunction()

#### test_CSMInfiniteRatioFunction()

#### test_DeltaCSMRatioFunction()
## pymatgen.analysis.chemenv.utils.tests.test_graph_utils module


### _class_ pymatgen.analysis.chemenv.utils.tests.test_graph_utils.EnvironmentNodesGraphUtilsTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_cycle()

### _class_ pymatgen.analysis.chemenv.utils.tests.test_graph_utils.FakeNode(isite)
Bases: `object`


### _class_ pymatgen.analysis.chemenv.utils.tests.test_graph_utils.FakeNodeWithEqLtMethods(isite)
Bases: `object`


### _class_ pymatgen.analysis.chemenv.utils.tests.test_graph_utils.FakeNodeWithEqLtMethodsBis(isite)
Bases: `FakeNodeWithEqLtMethods`


### _class_ pymatgen.analysis.chemenv.utils.tests.test_graph_utils.FakeNodeWithEqMethod(isite)
Bases: `object`


### _class_ pymatgen.analysis.chemenv.utils.tests.test_graph_utils.FakeNodeWithEqMethodWrongSortable(isite)
Bases: `object`


### _class_ pymatgen.analysis.chemenv.utils.tests.test_graph_utils.GraphUtilsTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_get_delta()

#### test_multigraph_cycle()

#### test_simple_graph_cycle()
## pymatgen.analysis.chemenv.utils.tests.test_math_utils module


### _class_ pymatgen.analysis.chemenv.utils.tests.test_math_utils.MathUtilsTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_cosinus_step()

#### test_linearly_independent_vectors()

#### test_list_cartesian_product()

#### test_math_utils()

#### test_power3_step()

#### test_powern_parts_step()

#### test_scale_and_clamp()

#### test_smootherstep()

#### test_smoothstep()