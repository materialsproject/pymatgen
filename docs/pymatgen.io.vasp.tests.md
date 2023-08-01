---
layout: default
title: pymatgen.io.vasp.tests.md
nav_exclude: true
---

# pymatgen.io.vasp.tests package


## pymatgen.io.vasp.tests.test_inputs module


### _class_ pymatgen.io.vasp.tests.test_inputs.IncarTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_as_dict_and_from_dict()

#### test_check_params()

#### test_diff()

#### test_get_string()

#### test_init()

#### test_lsorbit_magmom()

#### test_proc_types()

#### test_quad_efg()

#### test_types()

#### test_write()

### _class_ pymatgen.io.vasp.tests.test_inputs.KpointsTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_as_dict_from_dict()

#### test_automatic_kpoint()

#### test_init()

#### test_kpt_bands_as_dict_from_dict()

#### test_pickle()

#### test_static_constructors()

#### test_style_setter()

### _class_ pymatgen.io.vasp.tests.test_inputs.PoscarTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_cart_scale()

#### test_from_file()

#### test_from_md_run()

#### test_init()

#### test_selective_dynamics()

#### test_setattr()

#### test_significant_figures()

#### test_str()

#### test_to_from_dict()

#### test_velocities()

#### test_write()

#### test_write_md_poscar()

### _class_ pymatgen.io.vasp.tests.test_inputs.PotcarSingleTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_attributes()

#### test_bad_value()

#### test_electron_config()

#### test_found_unknown_key()

#### test_functional_types()

#### test_hash()

#### test_identify_potcar()

#### test_keywords()

#### test_multi_potcar_with_and_without_hash()

#### test_nelectrons()

#### test_potcar_file_hash_warning()

#### test_potcar_hash_warning()

#### test_psctr()

#### test_verify_correct_potcar_with_hash()

#### test_verify_faulty_potcar_with_hash()

### _class_ pymatgen.io.vasp.tests.test_inputs.PotcarTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_from_file()

#### test_init()

#### test_pickle()

#### test_potcar_map()

#### test_set_symbol()

#### test_to_from_dict()

#### test_write()

### _class_ pymatgen.io.vasp.tests.test_inputs.VaspInputTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_from_directory()

#### test_run_vasp()

#### test_to_from_dict()

#### test_write()
## pymatgen.io.vasp.tests.test_optics module


### _class_ pymatgen.io.vasp.tests.test_optics.VasprunTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_optics()

### pymatgen.io.vasp.tests.test_optics.test_delta_func()

### pymatgen.io.vasp.tests.test_optics.test_step_func()
## pymatgen.io.vasp.tests.test_outputs module


### _class_ pymatgen.io.vasp.tests.test_outputs.BSVasprunTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_get_band_structure()

### _class_ pymatgen.io.vasp.tests.test_outputs.ChgcarTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### _classmethod_ setUpClass()
Hook method for setting up class fixture before running tests in the class.


#### test_add()

#### test_as_dict_and_from_dict()

#### test_hdf5()

#### test_init()

#### test_soc_chgcar()

#### test_spin_data()

#### test_write()

### _class_ pymatgen.io.vasp.tests.test_outputs.DynmatTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_init()

### _class_ pymatgen.io.vasp.tests.test_outputs.EigenvalTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_eigenvalue_band_properties()

#### test_eigenvalue_band_properties_separate_spins()

#### test_init()

#### test_ispin2()

### _class_ pymatgen.io.vasp.tests.test_outputs.ElfcarTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_alpha()

#### test_init()

#### test_interpolation()

### _class_ pymatgen.io.vasp.tests.test_outputs.LocpotTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_init()

### _class_ pymatgen.io.vasp.tests.test_outputs.OszicarTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_init()

### _class_ pymatgen.io.vasp.tests.test_outputs.OutcarTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_avg_core_poten()

#### test_chemical_shielding()

#### test_chemical_shielding_with_different_core_contribution()

#### test_core_state_eigen()

#### test_cs_core_contribution()

#### test_cs_g0_contribution()

#### test_cs_raw_tensors()

#### test_dielectric()

#### test_drift()

#### test_electrostatic_potential()

#### test_energies()

#### test_freq_dielectric()

#### test_freq_dielectric_vasp544()

#### test_init()

#### test_mag_electrostatic_error()

#### test_nmr_efg()

#### test_nplwvs()

#### test_onsite_density_matrix()

#### test_parse_sci_notation()

#### test_polarization()

#### test_pseudo_zval()

#### test_read_elastic_tensor()

#### test_read_fermi_contact_shift()

#### test_read_lcalcpol()

#### test_read_piezo_tensor()

#### test_read_table_pattern()

#### test_serial_compilation()

#### test_single_atom()

#### test_soc()

#### test_stopped()

#### test_stopped_old()

#### test_vasp620_format()

### _class_ pymatgen.io.vasp.tests.test_outputs.ProcarTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_init()

#### test_phase_factors()

### _class_ pymatgen.io.vasp.tests.test_outputs.VasprunTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_Xe()

#### test_as_dict()

#### test_bad_random_seed()

#### test_bad_vasprun()

#### test_charge_charge_dielectric()
VASP 5.4.4 writes out two dielectric functions to vasprun.xml
These are the “density-density” and “velocity-velocity” linear response functions.
See the comments in linear_optics.F for details.


#### test_charged_structure()

#### test_chi()

#### test_dfpt()

#### test_dfpt_ionic()

#### test_dfpt_unconverged()

#### test_dielectric()

#### test_dielectric_vasp608()

#### test_eigenvalue_band_properties_separate_spins()

#### test_energies()

#### test_force_constants()

#### test_get_band_structure()

#### test_indirect_vasprun()

#### test_invalid_element()

#### test_kpointset_electronvelocities()

#### test_multiple_dielectric()

#### test_no_projected()

#### test_nonlmn()

#### test_optical_absorption_coeff()

#### test_optical_vasprun()

#### test_parsing_chemical_shift_calculations()

#### test_parsing_efg_calcs()

#### test_potcar_not_found()

#### test_projected_magnetisation()

#### test_runtype()

#### test_sc_step_overflow()

#### test_search_for_potcar()

#### test_selective_dynamics()

#### test_smart_efermi()

#### test_standard()

#### test_unconverged()

#### test_uniform()

#### test_update_potcar()

#### test_vasprun_ml()

#### test_vasprun_with_more_than_two_unlabelled_dielectric_functions()

#### test_vdw()

### _class_ pymatgen.io.vasp.tests.test_outputs.WSWQTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_consistency()

### _class_ pymatgen.io.vasp.tests.test_outputs.WavecarTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test__generate_G_points()

#### test__generate_nbmax()

#### test_evaluate_wavefunc()

#### test_fft_mesh_advanced()

#### test_fft_mesh_basic()

#### test_get_parchg()

#### test_n2_45210()

#### test_n2_spin()

#### test_standard()

#### test_write_unks()

### _class_ pymatgen.io.vasp.tests.test_outputs.WavederTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_consistency()

### _class_ pymatgen.io.vasp.tests.test_outputs.XdatcarTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_init()
## pymatgen.io.vasp.tests.test_sets module


### _class_ pymatgen.io.vasp.tests.test_sets.DictSetTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### _classmethod_ setUpClass()
Hook method for setting up class fixture before running tests in the class.


#### test_as_dict()

### _class_ pymatgen.io.vasp.tests.test_sets.FuncTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_batch_write_input()

### _class_ pymatgen.io.vasp.tests.test_sets.LobsterSetTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_as_from_dict()

#### test_incar()

#### test_kpoints()

#### test_potcar()

### _class_ pymatgen.io.vasp.tests.test_sets.MITMDSetTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_as_from_dict()

#### test_params()

### _class_ pymatgen.io.vasp.tests.test_sets.MITMPRelaxSetTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### _classmethod_ setUpClass()
Hook method for setting up class fixture before running tests in the class.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_MPMetalRelaxSet()

#### test_as_from_dict()

#### test_estimate_nbands()

#### test_get_incar()

#### test_get_kpoints()

#### test_get_vasp_input()

#### test_hubbard_off_and_ediff_override()

#### test_incar_lmaxmix()

#### test_lda_potcar()

#### test_metal_check()

#### test_nelect()

#### test_poscar()

#### test_potcar_special_defaults()

#### test_potcar_symbols()

#### test_potcar_validation()

#### test_user_potcar_settings()

#### test_valid_magmom_struct()

#### test_write_input()

### _class_ pymatgen.io.vasp.tests.test_sets.MITNEBSetTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_as_from_dict()

#### test_incar()

#### test_kpoints()

#### test_potcar_symbols()

#### test_write_input()

### _class_ pymatgen.io.vasp.tests.test_sets.MPAbsorptionSetTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### pytestmark(_ = [Mark(name='skipif', args=(False,), kwargs={'reason': 'PMG_VASP_PSP_DIR is not set.'})_ )

#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_ipa()

#### test_rpa()

### _class_ pymatgen.io.vasp.tests.test_sets.MPHSEBSTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_init()

#### test_override_from_prev_calc()

### _class_ pymatgen.io.vasp.tests.test_sets.MPNMRSetTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_incar()

### _class_ pymatgen.io.vasp.tests.test_sets.MPNonSCFSetTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_init()

#### test_kpoints()

#### test_optics()

#### test_override_from_prev()

#### test_user_kpoint_override()

### _class_ pymatgen.io.vasp.tests.test_sets.MPSOCSetTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_from_prev_calc()

#### test_override_from_prev_calc()

### _class_ pymatgen.io.vasp.tests.test_sets.MPScanRelaxSetTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_as_from_dict()

#### test_incar()

#### test_incar_overrides()

#### test_kspacing_cap()

#### test_nonmetal()

#### test_other_vdw()

#### test_potcar()

#### test_rvv10()

#### test_scan_substitute()

#### test_write_input()

### _class_ pymatgen.io.vasp.tests.test_sets.MPScanStaticSetTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_init()

#### test_override_from_prev_calc()

### _class_ pymatgen.io.vasp.tests.test_sets.MPStaticSetTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_grid_size_from_struct()

#### test_init()

#### test_kspacing_override()

#### test_override_from_prev_calc()

#### test_standardize_structure()

#### test_user_incar_kspacing()

#### test_write_input_zipped()

### _class_ pymatgen.io.vasp.tests.test_sets.MVLElasticSetTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_incar()

### _class_ pymatgen.io.vasp.tests.test_sets.MVLGBSetTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### pytestmark(_ = [Mark(name='skipif', args=(False,), kwargs={'reason': 'PMG_VASP_PSP_DIR is not set.'})_ )

#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_bulk()

#### test_kpoints()

#### test_slab()

### _class_ pymatgen.io.vasp.tests.test_sets.MVLGWSetTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### pytestmark(_ = [Mark(name='skipif', args=(False,), kwargs={'reason': 'PMG_VASP_PSP_DIR is not set.'})_ )

#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_bse()

#### test_diag()

#### test_static()

### _class_ pymatgen.io.vasp.tests.test_sets.MVLNPTMDSetTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### pytestmark(_ = [Mark(name='skipif', args=(False,), kwargs={'reason': 'PMG_VASP_PSP_DIR is not set.'})_ )

#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_as_from_dict()

#### test_incar()

### _class_ pymatgen.io.vasp.tests.test_sets.MVLRelax52SetTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_as_from_dict()

#### test_incar()

#### test_potcar()

### _class_ pymatgen.io.vasp.tests.test_sets.MVLScanRelaxSetTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_as_from_dict()

#### test_incar()

#### test_potcar()

### _class_ pymatgen.io.vasp.tests.test_sets.MVLSlabSetTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### pytestmark(_ = [Mark(name='skipif', args=(False,), kwargs={'reason': 'PMG_VASP_PSP_DIR is not set.'})_ )

#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_as_dict()

#### test_bulk()

#### test_kpoints()

#### test_slab()

#### test_user_incar_settings()

### _class_ pymatgen.io.vasp.tests.test_sets.MagmomLdauTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_ln_magmom()

#### test_structure_from_prev_run()

### _class_ pymatgen.io.vasp.tests.test_sets.SetChangeCheckTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_sets_changed()

### pymatgen.io.vasp.tests.test_sets.test_Yb_2_warning(input_set: [VaspInputSet](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.VaspInputSet))