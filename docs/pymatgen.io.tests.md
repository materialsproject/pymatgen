---
layout: default
title: pymatgen.io.tests.md
nav_exclude: true
---

# pymatgen.io.tests package


## pymatgen.io.tests.test_adf module


### _class_ pymatgen.io.tests.test_adf.AdfInputTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_main()

### _class_ pymatgen.io.tests.test_adf.AdfKeyTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_atom_block_key()

#### test_end()

#### test_from_string()

#### test_option_operations()

#### test_options()

#### test_simple()

#### test_subkeys()

#### test_subkeys_subkeys()

### _class_ pymatgen.io.tests.test_adf.AdfOutputTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_analytical_freq()

#### test_numerical_freq()

#### test_single_point()

### _class_ pymatgen.io.tests.test_adf.AdfTaskTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_energy()

#### test_serialization()

### pymatgen.io.tests.test_adf.readfile(file_object)
Return the content of the file as a string.


* **Parameters**


    * **file_object** (*file** or **str*) – The file to read. This can be either a File object or a file path.


    * **Returns** –


    * **-------** –


    * **content** (*str*) – The content of the file.


## pymatgen.io.tests.test_ase module


### _class_ pymatgen.io.tests.test_ase.AseAtomsAdaptorTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_back_forth()

#### test_get_atoms_from_molecule()

#### test_get_atoms_from_molecule_dyn()

#### test_get_atoms_from_molecule_mags()

#### test_get_atoms_from_structure()

#### test_get_atoms_from_structure_charge()

#### test_get_atoms_from_structure_dyn()

#### test_get_atoms_from_structure_mags()

#### test_get_atoms_from_structure_oxi_states()

#### test_get_molecule()

#### test_get_structure()

#### test_get_structure_dyn()

#### test_get_structure_mag()
## pymatgen.io.tests.test_atat module


### _class_ pymatgen.io.tests.test_atat.AtatTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_mcsqs_cif_nacl()

#### test_mcsqs_cif_pzt()

#### test_mcsqs_export()

#### test_mcsqs_import()
## pymatgen.io.tests.test_babel module

## pymatgen.io.tests.test_cif module


### _class_ pymatgen.io.tests.test_cif.CifBlockTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_double_quoted_data()

#### test_double_quotes_and_underscore_data()

#### test_long_loop()

#### test_nested_fake_multiline_quotes()

#### test_to_string()

### _class_ pymatgen.io.tests.test_cif.CifIOTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_CifParser()

#### test_CifParserCod()
Parsing problematic cif files from the COD database


#### test_CifWriter()

#### test_bad_cif()

#### test_cif_parser_springer_pauling()

#### test_cif_writer_labeled()

#### test_cifwrite_without_refinement()

#### test_disordered()

#### test_dot_positions()

#### test_empty()

#### test_empty_deque()

#### test_get_lattice_from_lattice_type()

#### test_get_symmetrized_structure()

#### test_implicit_hydrogen()

#### test_missing_atom_site_type_with_oxi_states()

#### test_no_coords_or_species()

#### test_no_symmops()

#### test_one_line_symm()

#### test_parse_symbol()
Test the _parse_symbol function with several potentially
problematic examples of symbols and labels.


#### test_primes()

#### test_replacing_finite_precision_frac_coords()

#### test_site_labels()

#### test_site_symbol_preference()

#### test_specie_cifwriter()

#### test_symmetrized()

### _class_ pymatgen.io.tests.test_cif.MagCifTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_bibtex()

#### test_get_structures()

#### test_mcif_detection()

#### test_write()
## pymatgen.io.tests.test_common module


### pymatgen.io.tests.test_common.test_cube_io_faithful(tmp_path: Path)
## pymatgen.io.tests.test_core module


### _class_ pymatgen.io.tests.test_core.FakeClass(a, b)
Bases: `object`


#### write_file()

### _class_ pymatgen.io.tests.test_core.StructInputFile(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Bases: [`InputFile`](pymatgen.io.md#pymatgen.io.core.InputFile)

Test implementation of an InputFile object for CIF.


#### _classmethod_ from_str(contents: str)
Create an InputFile object from a string.


* **Parameters**

    **contents** – The contents of the file as a single string



* **Returns**

    InputFile



#### get_string()
Return a string representation of an entire input file.


### _class_ pymatgen.io.tests.test_core.TestInputFile(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_file_io()

#### test_msonable()

### _class_ pymatgen.io.tests.test_core.TestInputSet(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### _classmethod_ setUpClass()
Hook method for setting up class fixture before running tests in the class.


#### test_copy()

#### test_deepcopy()

#### test_equality()

#### test_mapping()

#### test_msonable()

#### test_write()

#### test_write_from_str()
## pymatgen.io.tests.test_cssr module

Created on Jan 24, 2012


### _class_ pymatgen.io.tests.test_cssr.CssrTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_from_file()

#### test_str()
## pymatgen.io.tests.test_fiesta module


### _class_ pymatgen.io.tests.test_fiesta.FiestaInputTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_init()

#### test_str_and_from_string()

### _class_ pymatgen.io.tests.test_fiesta.FiestaOutputTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_props()
## pymatgen.io.tests.test_gaussian module


### _class_ pymatgen.io.tests.test_gaussian.GaussianInputTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_from_cart_coords()

#### test_from_file()

#### test_from_string()

#### test_gen_basis()

#### test_init()

#### test_multiple_paramaters()
This test makes sure that input files with multi-parameter keywords
and route cards with multiple lines can be parsed accurately.


#### test_no_molecule()
Test that we can write input files without a geometry


#### test_no_molecule_func_bset_charge_mult()
Test that we can write inputs files without a geometry,
functional, basis set, charge or multiplicity
(mainly for postprocessing jobs where this info is read from .chk)


#### test_str_and_from_string()

### _class_ pymatgen.io.tests.test_gaussian.GaussianOutputTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_geo_opt()
Test an optimization where no “input orientation” is outputted


#### test_multiple_parameters()
This test makes sure that input files with multi-parameter keywords
and route cards with multiple lines can be parsed accurately.


#### test_multiple_parameters_with_multiple_completed_lines()
This test makes sure that input files with multi-parameter keywords
and route cards with multiple completed lines which are split by line break parse correctly.


#### test_pop()

#### test_props()

#### test_resume()

#### test_scan()

#### test_td()
## pymatgen.io.tests.test_jarvis module


### _class_ pymatgen.io.tests.test_jarvis.JarvisAtomsAdaptorTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_get_atoms_from_structure()

#### test_get_structure()
## pymatgen.io.tests.test_lmto module


### _class_ pymatgen.io.tests.test_lmto.CoplTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_attributes()

#### test_cohp_data()

#### test_energies()

### _class_ pymatgen.io.tests.test_lmto.CtrlTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_dict()

#### test_read_write()

#### test_structure()
## pymatgen.io.tests.test_nwchem module


### _class_ pymatgen.io.tests.test_nwchem.NwInputTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_from_string_and_file()

#### test_str()

#### test_to_from_dict()

### _class_ pymatgen.io.tests.test_nwchem.NwOutputTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_get_excitation_spectrum()

#### test_parse_tddft()

#### test_read()

### _class_ pymatgen.io.tests.test_nwchem.NwTaskTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_dft_cosmo_task()

#### test_dft_task()

#### test_esp_task()

#### test_init()

#### test_multi_bset()

#### test_str_and_from_string()

#### test_to_from_dict()
## pymatgen.io.tests.test_packmol module

## pymatgen.io.tests.test_phonopy module


### _class_ pymatgen.io.tests.test_phonopy.GetDisplacedStructuresTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_get_displaced_structures()

### _class_ pymatgen.io.tests.test_phonopy.PhonopyParserTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_get_complete_dos()

#### test_get_ph_bs()

#### test_get_ph_dos()

### _class_ pymatgen.io.tests.test_phonopy.StructureConversionTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_structure_conversion()

### _class_ pymatgen.io.tests.test_phonopy.TestGruneisen(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_gruneisen_parameter()

#### test_ph_bs_symm_line()

### _class_ pymatgen.io.tests.test_phonopy.TestPhonopyFromForceConstants(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_get_phonon_band_structure_from_fc()

#### test_get_phonon_band_structure_symm_line_from_fc()

#### test_get_phonon_dos_from_fc()

### _class_ pymatgen.io.tests.test_phonopy.TestThermalDisplacementMatrices(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_get_thermal_displacement_matrix()
## pymatgen.io.tests.test_prismatic module


### _class_ pymatgen.io.tests.test_prismatic.PrismaticTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_to_string()
## pymatgen.io.tests.test_pwscf module


### _class_ pymatgen.io.tests.test_pwscf.PWInputTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_init()

#### test_proc_val()

#### test_read_str()

#### test_str_mixed_oxidation()

#### test_str_with_oxidation()

#### test_str_without_oxidation()

#### test_write_str_with_kpoints()

### _class_ pymatgen.io.tests.test_pwscf.PWOuputTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_get_celldm()

#### test_properties()
## pymatgen.io.tests.test_res module


### _class_ pymatgen.io.tests.test_res.TestAirssProvider()
Bases: `object`


#### pytestmark(_ = [Mark(name='parametrize', args=('provider', [<pymatgen.io.res.AirssProvider object>]), kwargs={})_ )

#### test_as_dict(provider: [AirssProvider](pymatgen.io.md#pymatgen.io.res.AirssProvider))

#### test_cggf(provider: [AirssProvider](pymatgen.io.md#pymatgen.io.res.AirssProvider))

#### test_entry(provider: [AirssProvider](pymatgen.io.md#pymatgen.io.res.AirssProvider))

#### test_lattice(provider: [AirssProvider](pymatgen.io.md#pymatgen.io.res.AirssProvider))

#### test_misc(provider: [AirssProvider](pymatgen.io.md#pymatgen.io.res.AirssProvider))

#### test_moks(provider: [AirssProvider](pymatgen.io.md#pymatgen.io.res.AirssProvider))

#### test_pspots(provider: [AirssProvider](pymatgen.io.md#pymatgen.io.res.AirssProvider))

#### test_raise(provider: [AirssProvider](pymatgen.io.md#pymatgen.io.res.AirssProvider))

#### test_titl(provider: [AirssProvider](pymatgen.io.md#pymatgen.io.res.AirssProvider))

### _class_ pymatgen.io.tests.test_res.TestSpin()
Bases: `object`


#### test_gh_2938_example()

#### test_read_spin()

### _class_ pymatgen.io.tests.test_res.TestStructureModule()
Bases: `object`


#### test_structure_from_file()
## pymatgen.io.tests.test_shengbte module


### _class_ pymatgen.io.tests.test_shengbte.TestShengBTE(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_MSONable_implementation()

#### test_from_dict()

#### test_from_file()
## pymatgen.io.tests.test_template_input module


### _class_ pymatgen.io.tests.test_template_input.TestTemplateInputGen(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_write_inputs()
## pymatgen.io.tests.test_wannier90 module

Tests for pymatgen.io.wannier90


### _class_ pymatgen.io.tests.test_wannier90.UnkTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_eq()

#### test_from_file()

#### test_init()

#### test_read_write()

#### test_repr()

#### test_write_file()
## pymatgen.io.tests.test_xcrysden module


### _class_ pymatgen.io.tests.test_xcrysden.XSFTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_to_string()

#### test_xsf()

#### test_xsf_symbolparse()
Ensure that the same structure is parsed
even if the atomic symbol / number convention
is different.

## pymatgen.io.tests.test_xr module


### _class_ pymatgen.io.tests.test_xr.XrTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_from_file()

#### test_str()
## pymatgen.io.tests.test_xyz module


### _class_ pymatgen.io.tests.test_xyz.XYZTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_as_dataframe()

#### test_from_file()

#### test_from_string()

#### test_init_from_structure()

#### test_str()
## pymatgen.io.tests.test_zeopp module


### _class_ pymatgen.io.tests.test_zeopp.GetFreeSphereParamsTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_get_free_sphere_params()

### _class_ pymatgen.io.tests.test_zeopp.GetHighAccuracyVoronoiNodesTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_get_voronoi_nodes()

### _class_ pymatgen.io.tests.test_zeopp.GetVoronoiNodesMultiOxiTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_get_voronoi_nodes()

### _class_ pymatgen.io.tests.test_zeopp.GetVoronoiNodesTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_get_voronoi_nodes()

### _class_ pymatgen.io.tests.test_zeopp.ZeoCssrOxiTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_from_file()

#### test_str()

### _class_ pymatgen.io.tests.test_zeopp.ZeoCssrTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_from_file()

#### test_str()

### _class_ pymatgen.io.tests.test_zeopp.ZeoVoronoiXYZTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_from_file()

#### test_str()