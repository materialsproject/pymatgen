---
layout: default
title: pymatgen.phonon.tests.md
nav_exclude: true
---

# pymatgen.phonon.tests package


## pymatgen.phonon.tests.test_bandstructure module


### _class_ pymatgen.phonon.tests.test_bandstructure.PhononBandStructureSymmLineTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_basic()

#### test_branches()

#### test_dict_methods()

#### test_nac()

#### test_write_methods()
## pymatgen.phonon.tests.test_dos module


### _class_ pymatgen.phonon.tests.test_dos.CompleteDosTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_dict_methods()

#### test_properties()

#### test_str()

### _class_ pymatgen.phonon.tests.test_dos.DosTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_dict_methods()

#### test_get_smeared_densities()

#### test_properties()

#### test_thermodynamic_functions()
## pymatgen.phonon.tests.test_gruneisen module


### _class_ pymatgen.phonon.tests.test_gruneisen.GruneisenParameterTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_acoustic_debye_temp()

#### test_average_gruneisen()

#### test_debye_temp_phonopy()

#### test_frequencies()

#### test_fromdict_asdict()

#### test_gruneisen()

#### test_multi()

#### test_phdos()

#### test_plot()

#### test_tdos()

#### test_thermal_conductivity_slack()

### _class_ pymatgen.phonon.tests.test_gruneisen.GruneisenPhononBandStructureSymmLineTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_as_dict_from_dict()

#### test_plot()
## pymatgen.phonon.tests.test_ir_spectra module


### _class_ pymatgen.phonon.tests.test_ir_spectra.IRDielectricTensorTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_basic()
## pymatgen.phonon.tests.test_plotter module


### _class_ pymatgen.phonon.tests.test_plotter.PhononBSPlotterTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_bs_plot_data()

#### test_plot()

#### test_plot_compare()

#### test_proj_plot()

### _class_ pymatgen.phonon.tests.test_plotter.PhononDosPlotterTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_add_dos_dict()

#### test_get_dos_dict()

#### test_plot()

### _class_ pymatgen.phonon.tests.test_plotter.ThermoPlotterTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_plot_functions()
## pymatgen.phonon.tests.test_thermal_displacements module


### _class_ pymatgen.phonon.tests.test_thermal_displacements.ThermalDisplacementTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Test data from J. George’s matlab code [https://github.com/JaGeo/MolecularToolbox](https://github.com/JaGeo/MolecularToolbox).

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_B()

#### test_U1U2U3()

#### test_Ucart()

#### test_Ucif()

#### test_Ustar()

#### test_angle()

#### test_beta()

#### test_compute_directionality_quality_criterion()

#### test_from_cif_P1()

#### test_from_ucif()

#### test_ratio_prolate()

#### test_to_structure_with_site_properties()

#### test_visualization_directionality_criterion()

#### test_write_file()