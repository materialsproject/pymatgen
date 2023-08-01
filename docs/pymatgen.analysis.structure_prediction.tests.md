---
layout: default
title: pymatgen.analysis.structure_prediction.tests.md
nav_exclude: true
---

# pymatgen.analysis.structure_prediction.tests package


## pymatgen.analysis.structure_prediction.tests.test_dopant_predictor module


### _class_ pymatgen.analysis.structure_prediction.tests.test_dopant_predictor.DopantPredictionTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_dopants_from_shannon_radii()

#### test_dopants_from_substitution_probabilities()
## pymatgen.analysis.structure_prediction.tests.test_substitution_probability module


### _class_ pymatgen.analysis.structure_prediction.tests.test_substitution_probability.SubstitutionPredictorTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_prediction()

### _class_ pymatgen.analysis.structure_prediction.tests.test_substitution_probability.SubstitutionProbabilityTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_full_lambda_table()
This test tests specific values in the data folder. If the
json is updated, these tests will have to be as well


#### test_mini_lambda_table()

### pymatgen.analysis.structure_prediction.tests.test_substitution_probability.get_table()
Loads a lightweight lambda table for use in unit tests to reduce
initialization time, and make unit tests insensitive to changes in the
default lambda table.

## pymatgen.analysis.structure_prediction.tests.test_substitutor module


### _class_ pymatgen.analysis.structure_prediction.tests.test_substitutor.SubstitutorTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_as_dict()

#### test_substitutor()

### pymatgen.analysis.structure_prediction.tests.test_substitutor.get_table()
Loads a lightweight lambda table for use in unit tests to reduce
initialization time, and make unit tests insensitive to changes in the
default lambda table.

## pymatgen.analysis.structure_prediction.tests.test_volume_predictor module


### _class_ pymatgen.analysis.structure_prediction.tests.test_volume_predictor.DLSVolumePredictorTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_predict()

### _class_ pymatgen.analysis.structure_prediction.tests.test_volume_predictor.RLSVolumePredictorTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_modes()

#### test_predict()