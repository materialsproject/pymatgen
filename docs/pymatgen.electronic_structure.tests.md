---
layout: default
title: pymatgen.electronic_structure.tests.md
nav_exclude: true
---

# pymatgen.electronic_structure.tests package

Tests for pymatgen.electronic_structure.


## pymatgen.electronic_structure.tests.test_bandstructure module


### _class_ pymatgen.electronic_structure.tests.test_bandstructure.BandStructureSymmLineTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_apply_scissor_insulator()

#### test_apply_scissor_spin_polarized()

#### test_as_dict()

#### test_basic()

#### test_get_band_gap()

#### test_get_branch()

#### test_get_cbm()

#### test_get_direct_band_gap()

#### test_get_direct_band_gap_dict()

#### test_get_sym_eq_kpoints_and_degeneracy()

#### test_get_vbm()

#### test_is_metal()

#### test_old_format_load()

#### test_properties()

### _class_ pymatgen.electronic_structure.tests.test_bandstructure.KpointTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_as_dict()

#### test_from_dict()

#### test_properties()

### _class_ pymatgen.electronic_structure.tests.test_bandstructure.LobsterBandStructureSymmLineTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_as_dict()

#### test_basic()

#### test_get_band_gap()

#### test_get_branch()

#### test_get_cbm()

#### test_get_direct_band_gap()

#### test_get_direct_band_gap_dict()

#### test_get_sym_eq_kpoints_and_degeneracy()

#### test_get_vbm()

#### test_is_metal()

#### test_old_format_load()

#### test_proj_bandstructure_plot()

### _class_ pymatgen.electronic_structure.tests.test_bandstructure.ReconstructBandStructureTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_reconstruct_band_structure()

#### test_vasprun_bs()
## pymatgen.electronic_structure.tests.test_boltztrap module


### _class_ pymatgen.electronic_structure.tests.test_boltztrap.BoltztrapAnalyzerTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### _classmethod_ setUpClass()
Hook method for setting up class fixture before running tests in the class.


#### _classmethod_ tearDownClass()
Hook method for deconstructing the class fixture after running all tests in the class.


#### test_extreme()

#### test_get_average_eff_mass()

#### test_get_carrier_concentration()

#### test_get_complete_dos()

#### test_get_complexity_factor()

#### test_get_conductivity()

#### test_get_hall_carrier_concentration()

#### test_get_power_factor()

#### test_get_seebeck()

#### test_get_seebeck_eff_mass()

#### test_get_symm_bands()

#### test_get_thermal_conductivity()

#### test_get_zt()

#### test_properties()

#### test_to_from_dict()
## pymatgen.electronic_structure.tests.test_boltztrap2 module


### _class_ pymatgen.electronic_structure.tests.test_boltztrap2.BandstructureLoaderTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_get_volume()

#### test_properties()

### _class_ pymatgen.electronic_structure.tests.test_boltztrap2.BztInterpolatorTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_get_band_structure()

#### test_properties()

#### test_tot_dos()

#### test_tot_proj_dos()

### _class_ pymatgen.electronic_structure.tests.test_boltztrap2.BztPlotterTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_plot()

### _class_ pymatgen.electronic_structure.tests.test_boltztrap2.BztTransportPropertiesTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_compute_properties_doping()

#### test_properties()

### _class_ pymatgen.electronic_structure.tests.test_boltztrap2.VasprunBSLoaderTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_get_volume()

#### test_properties()

### _class_ pymatgen.electronic_structure.tests.test_boltztrap2.VasprunLoaderTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_from_file()

#### test_get_volume()

#### test_properties()
## pymatgen.electronic_structure.tests.test_cohp module


### _class_ pymatgen.electronic_structure.tests.test_cohp.CohpTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_antibnd_states_below_efermi()

#### test_as_from_dict()

#### test_attributes()

#### test_get_icohp()

#### test_get_interpolated_value()

#### test_str()

### _class_ pymatgen.electronic_structure.tests.test_cohp.CombinedIcohpTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_extremum_icohpvalue()

#### test_get_icohp_by_label()

#### test_get_icohp_dict_by_bondlengths()

#### test_get_icohp_dict_of_site()

#### test_get_summed_icohp_by_label_list()

### _class_ pymatgen.electronic_structure.tests.test_cohp.CompleteCohpTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_attributes()

#### test_dict()

#### test_get_cohp_by_label()

#### test_get_cohp_by_label_summed_spin()

#### test_get_summed_cohp_by_label_and_orbital_list()

#### test_get_summed_cohp_by_label_and_orbital_list_summed_spin_channels()

#### test_get_summed_cohp_by_label_list()

#### test_get_summed_cohp_by_label_list_summed_spin()

#### test_icohp_values()

#### test_orbital_resolved_cohp()

#### test_orbital_resolved_cohp_summed_spin_channels()

### _class_ pymatgen.electronic_structure.tests.test_cohp.IcohpValueTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_attributes()

#### test_icohpvalue()

#### test_str()

#### test_summed_icohp()

### _class_ pymatgen.electronic_structure.tests.test_cohp.MethodTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_get_integrated_cohp_in_energy_range_full()

#### test_get_integrated_cohp_in_energy_range_onefloat()

#### test_get_integrated_cohp_in_energy_range_whole_range()
## pymatgen.electronic_structure.tests.test_core module


### _class_ pymatgen.electronic_structure.tests.test_core.MagmomTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_equality()

#### test_get_consistent_set_and_saxis()

#### test_get_moments()

#### test_have_consistent_saxis()

#### test_init()

#### test_is_collinear()

#### test_negative()

#### test_relative_to_crystal_axes()

### _class_ pymatgen.electronic_structure.tests.test_core.OrbitalTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_cached()

#### test_init()

### _class_ pymatgen.electronic_structure.tests.test_core.SpinTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_cached()

#### test_from_int()

#### test_init()
## pymatgen.electronic_structure.tests.test_dos module


### _class_ pymatgen.electronic_structure.tests.test_dos.CompleteDosTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_as_dict()

#### test_band_filling()

#### test_band_width()

#### test_dos_fp_exceptions()

#### test_get_band_center()

#### test_get_dos_fp()

#### test_get_dos_fp_similarity()

#### test_get_gap()

#### test_get_n_moment()

#### test_get_upper_band_edge()

#### test_kurtosis()

#### test_skewness()

#### test_str()

#### test_to_from_dict()

### _class_ pymatgen.electronic_structure.tests.test_dos.DOSTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_get_gap()

### _class_ pymatgen.electronic_structure.tests.test_dos.DosTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_as_dict()

#### test_get_gap()

#### test_get_smeared_densities()

### _class_ pymatgen.electronic_structure.tests.test_dos.FermiDosTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_as_dict()

#### test_doping_fermi()

### _class_ pymatgen.electronic_structure.tests.test_dos.LobsterCompleteDosTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_get_element_spd_dos()

#### test_get_site_orbital_dos()

#### test_get_site_t2g_eg_resolved_dos()

#### test_get_spd_dos()

### _class_ pymatgen.electronic_structure.tests.test_dos.SpinPolarizationTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_spin_polarization()
## pymatgen.electronic_structure.tests.test_plotter module


### _class_ pymatgen.electronic_structure.tests.test_plotter.BSDOSPlotterTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_methods()

### _class_ pymatgen.electronic_structure.tests.test_plotter.BSPlotterProjectedTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_methods()

### _class_ pymatgen.electronic_structure.tests.test_plotter.BSPlotterTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_add_bs()

#### test_bs_plot_data()

#### test_get_branch_steps()

#### test_get_plot()

#### test_get_ticks()

#### test_interpolate_bands()

#### test_rescale_distances()

### _class_ pymatgen.electronic_structure.tests.test_plotter.BoltztrapPlotterTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_plot_carriers()

#### test_plot_complexity_factor_mu()

#### test_plot_conductivity_dop()

#### test_plot_conductivity_mu()

#### test_plot_conductivity_temp()

#### test_plot_dos()

#### test_plot_eff_mass_dop()

#### test_plot_eff_mass_temp()

#### test_plot_hall_carriers()

#### test_plot_power_factor_dop()

#### test_plot_power_factor_mu()

#### test_plot_power_factor_temp()

#### test_plot_seebeck_dop()

#### test_plot_seebeck_eff_mass_mu()

#### test_plot_seebeck_mu()

#### test_plot_seebeck_temp()

#### test_plot_zt_dop()

#### test_plot_zt_mu()

#### test_plot_zt_temp()

### _class_ pymatgen.electronic_structure.tests.test_plotter.CohpPlotterTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_add_cohp_dict()

#### test_attributes()

#### test_get_cohp_dict()

#### test_get_plot()

#### test_save_plot()

### _class_ pymatgen.electronic_structure.tests.test_plotter.DosPlotterTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### _static_ get_plot_attributes(plt)

#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_add_dos_dict()

#### test_get_dos_dict()

#### test_get_plot()

#### test_get_plot_limits()

### _class_ pymatgen.electronic_structure.tests.test_plotter.PlotBZTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_bz_plot()

#### test_fold_point()