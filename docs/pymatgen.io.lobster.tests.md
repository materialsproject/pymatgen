---
layout: default
title: pymatgen.io.lobster.tests.md
nav_exclude: true
---

# pymatgen.io.lobster.tests package


## pymatgen.io.lobster.tests.test_lobster module


### _class_ pymatgen.io.lobster.tests.test_lobster.BandoverlapsTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_attributes()

#### test_has_good_quality()

### _class_ pymatgen.io.lobster.tests.test_lobster.ChargeTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_get_structure_with_charges()

#### testattributes()

### _class_ pymatgen.io.lobster.tests.test_lobster.CohpcarTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_attributes()

#### test_cohp_data()

#### test_energies()

#### test_orbital_resolved_cohp()

### _class_ pymatgen.io.lobster.tests.test_lobster.DoscarTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_complete_dos()

#### test_energies()

#### test_is_spin_polarized()

#### test_itdensities()

#### test_pdos()

#### test_tdensities()

#### test_tdos()

### _class_ pymatgen.io.lobster.tests.test_lobster.FatbandTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_attributes()

#### test_get_bandstructure()

#### test_raises()

### _class_ pymatgen.io.lobster.tests.test_lobster.GrosspopTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_structure_with_grosspop()

#### testattributes()

### _class_ pymatgen.io.lobster.tests.test_lobster.IcohplistTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_attributes()

#### test_values()

### _class_ pymatgen.io.lobster.tests.test_lobster.LobsterinTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### is_kpoint_in_list(kpoint, kpointlist, weight, weightlist)

#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_MSONable_implementation()

#### test_diff()

#### test_from_file()

#### test_get_all_possible_basis_functions()

#### test_get_basis()

#### test_get_potcar_symbols()

#### test_getitem()

#### test_initialize_from_dict()

#### test_setitem()

#### test_standard_settings()

#### test_standard_with_energy_range_from_vasprun()

#### test_write_INCAR()

#### test_write_KPOINTS()

#### test_write_lobsterin()

### _class_ pymatgen.io.lobster.tests.test_lobster.LobsteroutTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_get_doc()

#### testattributes()

### _class_ pymatgen.io.lobster.tests.test_lobster.MadelungEnergiesTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_attributes()

### _class_ pymatgen.io.lobster.tests.test_lobster.SitePotentialsTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_attributes()

#### test_get_structure()

### _class_ pymatgen.io.lobster.tests.test_lobster.TestUtils(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_get_all_possible_basis_combinations()

### _class_ pymatgen.io.lobster.tests.test_lobster.WavefunctionTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_get_volumetricdata_density()

#### test_get_volumetricdata_imaginary()

#### test_get_volumetricdata_real()

#### test_parse_file()

#### test_set_volumetric_data()

#### test_write_file()
## pymatgen.io.lobster.tests.test_lobsterenv module


### _class_ pymatgen.io.lobster.tests.test_lobsterenv.TestLobsterNeighbors(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_cation_anion_mode_without_ions()

#### test_extended_structure_graph()

#### test_get_anion_types()

#### test_get_info_cohps_to_neighbors()

#### test_get_info_icohps_neighbors()

#### test_get_nn_info()

#### test_get_plot_label()

#### test_get_structure_environments()

#### test_get_strucuture_environments_further_tests()

#### test_get_sum_icohps_between_neighbors_of_atom()

#### test_molecules_allowed()

#### test_order_parameter()

#### test_raises_extended_structure_graph()

#### test_set_limits()

#### test_structure_graph()

#### test_wrong_additional_correction()