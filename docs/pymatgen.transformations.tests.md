---
layout: default
title: pymatgen.transformations.tests.md
nav_exclude: true
---

# pymatgen.transformations.tests package


## pymatgen.transformations.tests.test_advanced_transformations module


### _class_ pymatgen.transformations.tests.test_advanced_transformations.AddAdsorbateTransformationTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_apply_transformation()

### _class_ pymatgen.transformations.tests.test_advanced_transformations.ChargeBalanceTransformationTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_apply_transformation()

### _class_ pymatgen.transformations.tests.test_advanced_transformations.CubicSupercellTransformationTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_apply_transformation()

### _class_ pymatgen.transformations.tests.test_advanced_transformations.DisorderedOrderedTransformationTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_apply_transformation()

### _class_ pymatgen.transformations.tests.test_advanced_transformations.DopingTransformationTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_apply_transformation()

#### test_as_from_dict()

#### test_find_codopant()

### _class_ pymatgen.transformations.tests.test_advanced_transformations.EnumerateStructureTransformationTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_apply_transformation()

#### test_callable_sort_criteria()

#### test_m3gnet()

#### test_max_disordered_sites()

#### test_to_from_dict()

### _class_ pymatgen.transformations.tests.test_advanced_transformations.GrainBoundaryTransformationTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_apply_transformation()

### _class_ pymatgen.transformations.tests.test_advanced_transformations.MagOrderingTransformationTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_advanced_usage()

#### test_apply_transformation()

#### test_as_from_dict()

#### test_ferrimagnetic()

#### test_zero_spin_case()

### _class_ pymatgen.transformations.tests.test_advanced_transformations.MonteCarloRattleTransformationTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_apply_transformation()

### _class_ pymatgen.transformations.tests.test_advanced_transformations.MultipleSubstitutionTransformationTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_apply_transformation()

### _class_ pymatgen.transformations.tests.test_advanced_transformations.SQSTransformationTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_apply_transformation()

#### test_return_ranked_list()

#### test_spin()

### _class_ pymatgen.transformations.tests.test_advanced_transformations.SlabTransformationTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_apply_transformation()

### _class_ pymatgen.transformations.tests.test_advanced_transformations.SubstituteSurfaceSiteTransformationTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_apply_transformation()

### _class_ pymatgen.transformations.tests.test_advanced_transformations.SubstitutionPredictorTransformationTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_apply_transformation()

#### test_as_dict()

### _class_ pymatgen.transformations.tests.test_advanced_transformations.SuperTransformationTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_apply_transformation()

#### test_apply_transformation_mult()

### pymatgen.transformations.tests.test_advanced_transformations.get_table()
Loads a lightweight lambda table for use in unit tests to reduce
initialization time, and make unit tests insensitive to changes in the
default lambda table.

## pymatgen.transformations.tests.test_site_transformations module


### _class_ pymatgen.transformations.tests.test_site_transformations.AddSitePropertyTransformationTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_apply_transformation()

### _class_ pymatgen.transformations.tests.test_site_transformations.InsertSitesTransformationTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_apply_transformation()

#### test_to_from_dict()

### _class_ pymatgen.transformations.tests.test_site_transformations.PartialRemoveSitesTransformationTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_apply_transformation_best_first()

#### test_apply_transformation_complete()

#### test_apply_transformation_enumerate()

#### test_apply_transformation_fast()

#### test_str()

#### test_to_from_dict()

### _class_ pymatgen.transformations.tests.test_site_transformations.RadialSiteDistortionTransformationTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test()

#### test_second_nn()

### _class_ pymatgen.transformations.tests.test_site_transformations.RemoveSitesTransformationTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_apply_transformation()

#### test_to_from_dict()

### _class_ pymatgen.transformations.tests.test_site_transformations.ReplaceSiteSpeciesTransformationTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_apply_transformation()

#### test_to_from_dict()

### _class_ pymatgen.transformations.tests.test_site_transformations.TranslateSitesTransformationTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_apply_transformation()

#### test_apply_transformation_site_by_site()

#### test_to_from_dict()
## pymatgen.transformations.tests.test_standard_transformations module


### _class_ pymatgen.transformations.tests.test_standard_transformations.AutoOxiStateDecorationTransformationTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_apply_transformation()

#### test_as_from_dict()

### _class_ pymatgen.transformations.tests.test_standard_transformations.ChargedCellTransformationTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_apply_transformation()

### _class_ pymatgen.transformations.tests.test_standard_transformations.ConventionalCellTransformationTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_apply_transformation()

### _class_ pymatgen.transformations.tests.test_standard_transformations.DeformStructureTransformationTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_apply_transformation()

### _class_ pymatgen.transformations.tests.test_standard_transformations.DiscretizeOccupanciesTransformationTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_apply_transformation()

### _class_ pymatgen.transformations.tests.test_standard_transformations.OrderDisorderedStructureTransformationTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_apply_transformation()

#### test_best_first()

#### test_no_oxidation()

#### test_symmetrized_structure()

#### test_too_small_cell()

### _class_ pymatgen.transformations.tests.test_standard_transformations.OxidationStateDecorationTransformationTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_apply_transformation()

### _class_ pymatgen.transformations.tests.test_standard_transformations.OxidationStateRemovalTransformationTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_apply_transformation()

### _class_ pymatgen.transformations.tests.test_standard_transformations.PartialRemoveSpecieTransformationTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_apply_transformation()

#### test_apply_transformation_fast()

#### test_apply_transformations_best_first()

#### test_apply_transformations_complete_ranking()

### _class_ pymatgen.transformations.tests.test_standard_transformations.PerturbStructureTransformationTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_apply_transformation()

### _class_ pymatgen.transformations.tests.test_standard_transformations.PrimitiveCellTransformationTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_apply_transformation()

### _class_ pymatgen.transformations.tests.test_standard_transformations.RemoveSpeciesTransformationTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_apply_transformation()

### _class_ pymatgen.transformations.tests.test_standard_transformations.RotationTransformationsTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_as_from_dict()

#### test_rotation_transformation()

### _class_ pymatgen.transformations.tests.test_standard_transformations.ScaleToRelaxedTransformationTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_apply_transformation()

### _class_ pymatgen.transformations.tests.test_standard_transformations.SubstitutionTransformationTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_apply_transformation()

#### test_fractional_substitution()

### _class_ pymatgen.transformations.tests.test_standard_transformations.SupercellTransformationTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_apply_transformation()

#### test_from_scaling_factors()