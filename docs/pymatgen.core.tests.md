---
layout: default
title: pymatgen.core.tests.md
nav_exclude: true
---

# pymatgen.core.tests package

Unittests for pymatgen/core.


## pymatgen.core.tests.test_bonds module


### _class_ pymatgen.core.tests.test_bonds.CovalentBondTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_get_bond_order()

#### test_is_bonded()

#### test_length()

#### test_str()

### _class_ pymatgen.core.tests.test_bonds.FuncTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_get_bond_length()

#### test_get_bond_order()

#### test_obtain_all_bond_lengths()
## pymatgen.core.tests.test_composition module

Created on Nov 10, 2012.

@author: Shyue Ping Ong


### _class_ pymatgen.core.tests.test_composition.ChemicalPotentialTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_init()

#### test_math()

### _class_ pymatgen.core.tests.test_composition.CompositionTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_add()

#### test_almost_equals()

#### test_alphabetical_formula()

#### test_anonymized_formula()

#### test_as_dict()

#### test_average_electroneg()

#### test_chemical_system()

#### test_comparisons()

#### test_contains_element_type()

#### test_div()

#### test_equality()

#### test_equals()

#### test_formula()

#### test_fractional_composition()

#### test_from_dict()

#### test_from_weight_dict()

#### test_get_atomic_fraction()

#### test_get_wt_fraction()

#### test_hash_robustness()

#### test_hill_formula()

#### test_immutable()

#### test_in()

#### test_indeterminate_formula()

#### test_init()

#### test_init_numerical_tolerance()

#### test_integer_formula()

#### test_is_valid()

#### test_iupac_formula()

#### test_metallofullerene()

#### test_mixed_valence()

#### test_mul()

#### test_negative_compositions()

#### test_num_atoms()

#### test_oxi_state_decoration()

#### test_oxi_state_guesses()

#### test_pickle()

#### test_reduced_composition()

#### test_reduced_formula()

#### test_remove_charges()

#### test_replace()

#### test_special_formulas()

#### test_str_and_repr()

#### test_sub()

#### test_to_data_dict()

#### test_to_latex_html_unicode()

#### test_tofrom_weight_dict()

#### test_total_electrons()

#### test_weight()
## pymatgen.core.tests.test_interface module


### _class_ pymatgen.core.tests.test_interface.InterfaceTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_basic_props()

#### test_from_slabs()

#### test_gap_setter()

#### test_get_shifts_based_on_adsorbate_sites()

#### test_in_plane_offset_setter()

#### test_vacuum_over_film_setter()
## pymatgen.core.tests.test_ion module


### _class_ pymatgen.core.tests.test_ion.IonTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_alphabetical_formula()

#### test_anonymized_formula()

#### test_as_dict()

#### test_charge_from_formula()

#### test_equality()

#### test_equals()

#### test_formula()

#### test_from_dict()

#### test_init()

#### test_len()

#### test_mixed_valence()

#### test_mul()

#### test_num_atoms()

#### test_oxi_state_guesses()

#### test_special_formulas()

#### test_to_latex_string()
## pymatgen.core.tests.test_lattice module


### _class_ pymatgen.core.tests.test_lattice.LatticeTestCase(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_attributes()
Docstring for test_attributes.


#### test_copy()

#### test_d_hkl()

#### test_dot_and_norm()

#### test_equal()

#### test_find_all_mappings()

#### test_find_mapping()

#### test_format()

#### test_get_all_distances()

#### test_get_cartesian_or_frac_coord()

#### test_get_distance_and_image()

#### test_get_distance_and_image_strict()

#### test_get_lll_reduced_lattice()

#### test_get_miller_index_from_sites()

#### test_get_niggli_reduced_lattice()

#### test_get_points_in_sphere()

#### test_get_vector_along_lattice_directions()

#### test_get_wigner_seitz_cell()

#### test_init()

#### test_is_3d_periodic()

#### test_is_hexagonal()

#### test_lattice_matrices()
If alpha == 90 and beta == 90, two matrices are identical.


#### test_lll_basis()

#### test_mapping_symmetry()

#### test_monoclinic()

#### test_points_in_spheres()

#### test_reciprocal_lattice()

#### test_scale()

#### test_selling_dist()

#### test_selling_vector()

#### test_static_methods()

#### test_to_from_dict()
## pymatgen.core.tests.test_libxcfunc module


### _class_ pymatgen.core.tests.test_libxcfunc.LibxcFuncTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_libxcfunc_api()
Testing libxcfunc_api.

## pymatgen.core.tests.test_molecular_orbitals module


### _class_ pymatgen.core.tests.test_molecular_orbitals.MolecularOrbitalTestCase(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_aos_as_list()

#### test_fractional_compositions()

#### test_max_electronegativity()

#### test_obtain_band_edges()
## pymatgen.core.tests.test_operations module


### _class_ pymatgen.core.tests.test_operations.MagSymmOpTestCase(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_operate_magmom()

#### test_to_from_dict()

#### test_xyzt_string()

### _class_ pymatgen.core.tests.test_operations.SymmOpTestCase(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_apply_rotation_only()

#### test_are_symmetrically_related()

#### test_are_symmetrically_related_vectors()

#### test_inverse()

#### test_inversion()

#### test_operate()

#### test_operate_multi()

#### test_properties()

#### test_reflection()

#### test_to_from_dict()

#### test_transform_tensor()

#### test_xyz()
## pymatgen.core.tests.test_periodic_table module


### _class_ pymatgen.core.tests.test_periodic_table.DummySpeciesTestCase(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_attr()

#### test_eq()

#### test_from_string()

#### test_immutable()

#### test_init()

#### test_pickle()

#### test_sort()

### _class_ pymatgen.core.tests.test_periodic_table.ElementTestCase(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_attributes()

#### test_block()

#### test_data()

#### test_deepcopy()

#### test_dict()

#### test_from_name()

#### test_from_row_and_group()

#### test_full_electronic_structure()

#### test_ground_state_term_symbol()

#### test_group()

#### test_ie_ea()

#### test_init()

#### test_is()

#### test_is_metal()

#### test_nan_X()

#### test_oxidation_states()

#### test_pickle()

#### test_print_periodic_table()

#### test_radii()

#### test_row()

#### test_sort()

#### test_term_symbols()

#### test_valence()

### _class_ pymatgen.core.tests.test_periodic_table.FuncTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_get_el_sp()

### _class_ pymatgen.core.tests.test_periodic_table.SpeciesTestCase(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_attr()

#### test_cmp()

#### test_deepcopy()

#### test_eq()

#### test_get_crystal_field_spin()

#### test_get_nmr_mom()

#### test_get_shannon_radius()

#### test_init()

#### test_ionic_radius()

#### test_no_oxidation_state()

#### test_pickle()

#### test_sort()

#### test_stringify()

#### test_to_from_string()

### pymatgen.core.tests.test_periodic_table.test_symbol_oxi_state_str(symbol_oxi, expected_element, expected_oxi_state)
## pymatgen.core.tests.test_settings module


### pymatgen.core.tests.test_settings.test_load_settings(tmp_path: Path, monkeypatch: MonkeyPatch)
Test .pmgrc.yaml file is loaded correctly and env vars take precedence.

## pymatgen.core.tests.test_sites module


### _class_ pymatgen.core.tests.test_sites.PeriodicSiteTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_as_from_dict()

#### test_distance()

#### test_distance_and_image()

#### test_distance_from_point()

#### test_equality()

#### test_equality_with_label()

#### test_is_periodic_image()

#### test_properties()
Test the properties for a site.


#### test_repr()

#### test_setters()

#### test_to_unit_cell()

### _class_ pymatgen.core.tests.test_sites.SiteTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_distance()

#### test_gt_lt()

#### test_hash()

#### test_pickle()

#### test_properties()

#### test_setters()

#### test_to_from_dict()

### pymatgen.core.tests.test_sites.get_distance_and_image_old(site1, site2, jimage=None)
Gets distance between two sites assuming periodic boundary conditions.
If the index jimage of two sites atom j is not specified it selects the
j image nearest to the i atom and returns the distance and jimage
indices in terms of lattice vector translations. If the index jimage of
atom j is specified it returns the distance between the i atom and the
specified jimage atom, the given jimage is also returned.


* **Parameters**


    * **other** – other site to get distance from.


    * **jimage** – specific periodic image in terms of lattice translations,
    e.g., [1,0,0] implies to take periodic image that is one
    a-lattice vector away. If jimage is None, the image that is
    nearest to the site is found.



* **Returns**

    distance and periodic lattice translations of the other site
    for which the distance applies.



* **Return type**

    (distance, jimage)


**NOTE**: Assumes the primitive cell vectors are sufficiently not skewed such
that the condition |a|cos(ab_angle) < |b| for all possible cell
vector pairs. \*\* this method does not check this condition \*\*

## pymatgen.core.tests.test_spectrum module


### _class_ pymatgen.core.tests.test_spectrum.SpectrumTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_copy()

#### test_normalize()

#### test_operators()

#### test_smear()

#### test_str()
## pymatgen.core.tests.test_structure module


### _class_ pymatgen.core.tests.test_structure.IMoleculeTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_bad_molecule()

#### test_break_bond()

#### test_default_dict_attrs()

#### test_equal()

#### test_get_angle_dihedral()

#### test_get_boxed_structure()

#### test_get_centered_molecule()

#### test_get_covalent_bonds()

#### test_get_dist_matrix()

#### test_get_distance()

#### test_get_neighbors()

#### test_get_neighbors_in_shell()

#### test_get_zmatrix()

#### test_no_spin_check()

#### test_prop()

#### test_properties()

#### test_repr_str()

#### test_set_item()

#### test_site_properties()

#### test_to_from_dict()

#### test_to_from_file_string()

### _class_ pymatgen.core.tests.test_structure.IStructureTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_as_dataframe()

#### test_as_dict()

#### test_bad_structure()

#### test_coincide_sites()

#### test_copy()

#### test_equal()

#### test_fractional_occupations()

#### test_from_dict()

#### test_get_all_neighbors_and_get_neighbors()

#### test_get_all_neighbors_equal()

#### test_get_all_neighbors_outside_cell()

#### test_get_all_neighbors_small_cutoff()

#### test_get_dist_matrix()

#### test_get_distance()

#### test_get_miller_index()
Test for get miller index convenience method.


#### test_get_neighbor_list()

#### test_get_orderings()

#### test_get_primitive_structure()

#### test_get_sorted_structure()

#### test_get_space_group_data()

#### test_get_symmetric_neighbor_list()

#### test_interpolate()

#### test_interpolate_lattice()

#### test_interpolate_lattice_rotation()

#### test_labeled_structure()

#### test_matches()

#### test_pbc()

#### test_primitive_cell_site_merging()

#### test_primitive_on_large_supercell()

#### test_primitive_positions()

#### test_primitive_structure_volume_check()

#### test_site_properties()

#### test_sites_setter()

#### test_specie_init()

#### test_to_from_file_string()

#### test_volume_and_density()

### _class_ pymatgen.core.tests.test_structure.MoleculeTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_add_site_property()

#### test_apply_operation()

#### test_calculate_ase_mol()

#### test_calculate_gfnxtb()

#### test_extract_cluster()

#### test_insert_remove_append()

#### test_mutable_sequence_methods()

#### test_no_spin_check()

#### test_relax_ase_mol()

#### test_relax_ase_mol_return_traj()

#### test_relax_gfnxtb()

#### test_replace()

#### test_rotate_sites()

#### test_substitute()

#### test_to_from_dict()

#### test_to_from_file_string()

#### test_translate_sites()

### _class_ pymatgen.core.tests.test_structure.NeighborTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_msonable()

#### test_neighbor_labels()

### _class_ pymatgen.core.tests.test_structure.StructureTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_add_oxidation_states_by_element()

#### test_add_oxidation_states_by_guess()

#### test_add_oxidation_states_by_site()

#### test_add_remove_site_property()

#### test_add_remove_spin_states()

#### test_another_supercell()

#### test_append_insert_remove_replace_substitute()

#### test_apply_operation()

#### test_apply_strain()

#### test_calculate_ase()

#### test_calculate_chgnet()

#### test_calculate_m3gnet()

#### test_charge()

#### test_default_dict_attrs()

#### test_disordered_supercell_primitive_cell()

#### test_extract_cluster()

#### test_from_magnetic_spacegroup()

#### test_from_prototype()

#### test_from_sites()

#### test_from_spacegroup()

#### test_init_error()

#### test_make_supercell()

#### test_make_supercell_labeled()

#### test_merge_sites()

#### test_mul()

#### test_mutable_sequence_methods()

#### test_not_hashable()

#### test_perturb()

#### test_propertied_structure()

#### test_properties()

#### test_relax_ase()

#### test_relax_ase_opt_kwargs()

#### test_relax_ase_return_traj()

#### test_relax_chgnet()

#### test_relax_m3gnet()

#### test_relax_m3gnet_fixed_lattice()

#### test_relax_m3gnet_with_traj()

#### test_remove_oxidation_states()

#### test_rotate_sites()

#### test_scale_lattice()

#### test_set_item()

#### test_sort()

#### test_species()

#### test_to_from_abivars()
Test as_dict, from_dict with fmt == abivars.


#### test_to_from_dict()

#### test_to_from_file_string()

#### test_translate_sites()

#### test_vesta_lattice_matrix()
## pymatgen.core.tests.test_surface module


### _class_ pymatgen.core.tests.test_surface.MillerIndexFinderTests(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_generate_all_slabs()

#### test_get_symmetrically_distinct_miller_indices()

#### test_get_symmetrically_equivalent_miller_indices()

#### test_miller_index_from_sites()
Test surface miller index convenience function.


### _class_ pymatgen.core.tests.test_surface.ReconstructionGeneratorTests(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_build_slab()

#### test_get_d()

#### test_previous_reconstructions()

### _class_ pymatgen.core.tests.test_surface.SlabGeneratorTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_bonds_broken()

#### test_get_orthogonal_c_slab()

#### test_get_orthogonal_c_slab_site_props()

#### test_get_slab()

#### test_get_slabs()

#### test_get_tasker2_slabs()

#### test_move_to_other_side()

#### test_nonstoichiometric_symmetrized_slab()

#### test_normal_search()

#### test_triclinic_TeI()

### _class_ pymatgen.core.tests.test_surface.SlabTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_add_adsorbate_atom()

#### test_as_dict()

#### test_as_from_dict()

#### test_dipole_and_is_polar()

#### test_get_slab_regions()

#### test_get_sorted_structure()

#### test_get_symmetric_sites()

#### test_init()

#### test_methods()

#### test_oriented_unit_cell()

#### test_surface_sites_and_symmetry()

#### test_symmetrization()

### pymatgen.core.tests.test_surface.get_path(path_str)
## pymatgen.core.tests.test_tensors module


### _class_ pymatgen.core.tests.test_tensors.SquareTensorTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_get_scaled()

#### test_is_rotation()

#### test_new()

#### test_polar_decomposition()

#### test_properties()

#### test_refine_rotation()

#### test_serialization()

### _class_ pymatgen.core.tests.test_tensors.TensorCollectionTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### list_based_function_check(attribute, coll, \*args, \*\*kwargs)
This function allows for more efficient testing of list-based
functions in a “collection”-style class like TensorCollection.

It ensures that the test function


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_list_based_functions()

#### test_serialization()

### _class_ pymatgen.core.tests.test_tensors.TensorTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_convert_to_ieee()

#### test_einsum_sequence()

#### test_fit_to_structure()

#### test_from_values_indices()

#### test_from_voigt()

#### test_is_fit_to_structure()

#### test_is_symmetric()

#### test_new()

#### test_populate()

#### test_projection_methods()

#### test_rotate()

#### test_round()

#### test_serialization()

#### test_structure_transform()

#### test_summary_methods()

#### test_symmetrized()

#### test_symmetry_reduce()

#### test_tensor_mapping()

#### test_transform()

#### test_zeroed()
## pymatgen.core.tests.test_trajectory module


### _class_ pymatgen.core.tests.test_trajectory.TrajectoryTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_conversion()

#### test_displacements()

#### test_extend()

#### test_extend_frame_props()

#### test_extend_site_props()

#### test_frame_properties()

#### test_length()

#### test_list_slice()

#### test_single_index_slice()

#### test_site_properties()

#### test_slice()

#### test_to_from_dict()

#### test_variable_lattice()

#### test_xdatcar_write()
## pymatgen.core.tests.test_units module


### _class_ pymatgen.core.tests.test_units.ArrayWithFloatWithUnitTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_array_algebra()

#### test_as_base_units()

#### test_energy()
Similar to FloatWithUnitTest.test_energy.
Check whether EnergyArray and FloatWithUnit have same behavior.

# TODO
One can merge the two tests easily:

for obj in [Energy, EnergyArray]:

    a = obj(…)
    self.assert(…)


#### test_factors()

#### test_length()
Similar to FloatWithUnitTest.test_time.
Check whether EnergyArray and FloatWithUnit have same behavior.


#### test_time()
Similar to FloatWithUnitTest.test_time.
Check whether EnergyArray and FloatWithUnit have same behavior.


### _class_ pymatgen.core.tests.test_units.DataPersistenceTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_pickle()
Test whether FloatWithUnit and ArrayWithUnit support pickle.


### _class_ pymatgen.core.tests.test_units.FloatWithUnitTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_as_base_units()

#### test_compound_operations()

#### test_energy()

#### test_length()

#### test_memory()

#### test_time()

#### test_unitized()

### _class_ pymatgen.core.tests.test_units.UnitTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_init()
## pymatgen.core.tests.test_xcfunc module


### _class_ pymatgen.core.tests.test_xcfunc.LibxcFuncTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_xcfunc_api()
Testing XcFunc API.