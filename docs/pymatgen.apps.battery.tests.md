---
layout: default
title: pymatgen.apps.battery.tests.md
nav_exclude: true
---

# pymatgen.apps.battery.tests package

Tests for battery classes.


## pymatgen.apps.battery.tests.test_analyzer module


### _class_ pymatgen.apps.battery.tests.test_analyzer.BatteryAnalyzerTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### load_from_cif(filename, oxidations, working_ion='Li')

#### load_from_internal(name, oxidations, working_ion='Li')

#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_capacitygrav_calculations()

#### test_capacityvol_calculations()

#### test_ion_removal()

#### test_oxide_check()
## pymatgen.apps.battery.tests.test_conversion_battery module


### _class_ pymatgen.apps.battery.tests.test_conversion_battery.ConversionElectrodeTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_composite()

#### test_init()

#### test_summary()
## pymatgen.apps.battery.tests.test_insertion_battery module


### _class_ pymatgen.apps.battery.tests.test_insertion_battery.InsertionElectrodeTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_capacities()

#### test_entries()

#### test_get_all_entries()

#### test_get_instability()

#### test_get_muO2()

#### test_get_summary_dict()

#### test_init_no_structure()

#### test_to_from_dict()

#### test_voltage()

#### test_voltage_pair()
## pymatgen.apps.battery.tests.test_plotter module


### _class_ pymatgen.apps.battery.tests.test_plotter.VoltageProfilePlotterTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### testName()

#### testPlotly()