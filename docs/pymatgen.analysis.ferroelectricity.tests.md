---
layout: default
title: pymatgen.analysis.ferroelectricity.tests.md
nav_exclude: true
---

# pymatgen.analysis.ferroelectricity.tests namespace


## pymatgen.analysis.ferroelectricity.tests.test_polarization module


### _class_ pymatgen.analysis.ferroelectricity.tests.test_polarization.EnergyTrendTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_endpoints_minima()

#### test_max_spline_jump()

#### test_smoothness()

### _class_ pymatgen.analysis.ferroelectricity.tests.test_polarization.PolarizationTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_from_outcars_and_structures()

#### test_get_lattice_quanta()

#### test_get_polarization_change()

#### test_get_polarization_change_norm()

#### test_get_same_branch_polarization_data()

#### test_max_spline_jumps()

#### test_smoothness()

### _class_ pymatgen.analysis.ferroelectricity.tests.test_polarization.UtilsTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_get_total_ionic_dipole()

#### test_zval_dict_from_potcar()