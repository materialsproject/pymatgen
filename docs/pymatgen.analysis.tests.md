---
layout: default
title: pymatgen.analysis.tests.md
nav_exclude: true
---

# pymatgen.analysis.tests namespace


## pymatgen.analysis.tests.test_adsorption module


### _class_ pymatgen.analysis.tests.test_adsorption.AdsorbateSiteFinderTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_adsorb_both_surfaces()

#### test_find_adsorption_sites()

#### test_from_bulk_and_miller()

#### test_functions()

#### test_generate_adsorption_structures()

#### test_generate_substitution_structures()

#### test_init()
## pymatgen.analysis.tests.test_bond_dissociation module


### _class_ pymatgen.analysis.tests.test_bond_dissociation.BondDissociationTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_ec_neg_pcm_40()

#### test_pc_neutral_pcm_65()

#### test_tfsi_neg_no_pcm()
## pymatgen.analysis.tests.test_bond_valence module


### _class_ pymatgen.analysis.tests.test_bond_valence.BVAnalyzerTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_get_oxi_state_structure()

#### test_get_valence()

### _class_ pymatgen.analysis.tests.test_bond_valence.BondValenceSumTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_calculate_bv_sum()

#### test_calculate_bv_sum_unordered()
## pymatgen.analysis.tests.test_chempot_diagram module


### _class_ pymatgen.analysis.tests.test_chempot_diagram.ChemicalPotentialDiagramTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_border_hyperplanes()

#### test_centroid()

#### test_dim()

#### test_domains()

#### test_el_refs()

#### test_el_refs_formal()

#### test_get_2d_orthonormal_vector()

#### test_get_plot()

#### test_lims()

#### test_pca()
## pymatgen.analysis.tests.test_cost module


### _class_ pymatgen.analysis.tests.test_cost.CostAnalyzerTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_cost_per_kg()

#### test_cost_per_mol()

#### test_sanity()

### _class_ pymatgen.analysis.tests.test_cost.CostDBTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_sanity()
## pymatgen.analysis.tests.test_dimensionality module


### _class_ pymatgen.analysis.tests.test_dimensionality.CheonDimensionalityTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_get_dimensionality()

#### test_get_dimensionality_with_bonds()

#### test_tricky_structure()

### _class_ pymatgen.analysis.tests.test_dimensionality.GoraiDimensionalityTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_get_dimensionality()

#### test_get_dimensionality_with_bonds()

### _class_ pymatgen.analysis.tests.test_dimensionality.LarsenDimensionalityTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_calculate_dimensionality_of_site()

#### test_get_dimensionality()

#### test_get_structure_components()

#### test_tricky_structure()
Test for a tricky structure that other dimensionality finders say is
2D but is actually an interpenetrated 3D structure.


#### test_zero_d_to_molecule_graph()
## pymatgen.analysis.tests.test_disorder module


### _class_ pymatgen.analysis.tests.test_disorder.OrderParameterTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_compute_warren_cowley_parameters()
## pymatgen.analysis.tests.test_energy_models module


### _class_ pymatgen.analysis.tests.test_energy_models.EwaldElectrostaticModelTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_get_energy()

#### test_to_from_dict()

### _class_ pymatgen.analysis.tests.test_energy_models.IsingModelTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_get_energy()

#### test_to_from_dict()

### _class_ pymatgen.analysis.tests.test_energy_models.SymmetryModelTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_get_energy()

#### test_to_from_dict()
## pymatgen.analysis.tests.test_eos module


### _class_ pymatgen.analysis.tests.test_eos.EOSTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_eos_func()

#### test_eos_func_call()

#### test_fitting()

#### test_numerical_eos_values()

#### test_numerical_eoswrapper()

#### test_run_all_models()

#### test_summary_dict()
## pymatgen.analysis.tests.test_ewald module


### _class_ pymatgen.analysis.tests.test_ewald.EwaldMinimizerTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_init()

#### test_site()
Test that uses an uncharged structure


### _class_ pymatgen.analysis.tests.test_ewald.EwaldSummationTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_as_dict()

#### test_from_dict()

#### test_init()
## pymatgen.analysis.tests.test_fragmenter module


### _class_ pymatgen.analysis.tests.test_fragmenter.TestFragmentMolecule(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### _classmethod_ setUpClass()
Hook method for setting up class fixture before running tests in the class.


#### test_babel_pc_defaults()

#### test_babel_pc_frag1()

#### test_babel_pc_old_defaults()

#### test_babel_pc_with_ro_depth_0_vs_depth_10()

#### test_babel_tfsi()

#### test_edges_given_pc_frag1()

#### test_edges_given_pc_not_defaults()

#### test_edges_given_tfsi()

#### test_pc_depth_0_vs_depth_10()

#### test_pc_frag1_then_pc()

#### test_pc_then_ec_depth_10()
## pymatgen.analysis.tests.test_functional_groups module

## pymatgen.analysis.tests.test_graphs module


### _class_ pymatgen.analysis.tests.test_graphs.MoleculeGraphTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_as_from_dict()

#### test_build_unique_fragments()

#### test_construction()

#### test_coordination()

#### test_edge_editing()

#### test_find_rings()

#### test_get_disconnected()

#### test_insert_remove()

#### test_isomorphic()

#### test_properties()

#### test_replace()

#### test_set_node_attributes()

#### test_sort()

#### test_split()

#### test_substitute()

### _class_ pymatgen.analysis.tests.test_graphs.StructureGraphTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_auto_image_detection()

#### test_draw()

#### test_edge_editing()

#### test_extract_molecules()

#### test_from_edges()

#### test_from_local_env_and_equality_and_diff()

#### test_inappropriate_construction()

#### test_insert_remove()

#### test_mul()

#### test_no_duplicate_hops()

#### test_properties()

#### test_set_node_attributes()

#### test_sort()

#### test_str()

#### test_substitute()

#### test_to_from_dict()

#### test_types_and_weights_of_connections()

#### test_types_of_coordination_environments()

#### test_weight_statistics()
## pymatgen.analysis.tests.test_hhi module


### _class_ pymatgen.analysis.tests.test_hhi.HHIModelTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_hhi()
## pymatgen.analysis.tests.test_interface_reactions module


### _class_ pymatgen.analysis.tests.test_interface_reactions.InterfaceReactionTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_convert()

#### test_convexity()

#### test_get_chempot_correction()

#### test_get_critical_original_kink_ratio()

#### test_get_dataframe()

#### test_get_energy()

#### test_get_entry_energy()

#### test_get_get_elmt_amt_in_rxt()

#### test_get_grand_potential()

#### test_get_kinks()

#### test_get_no_mixing_energy()

#### test_get_original_composition_ratio()

#### test_get_reaction()

#### test_labels()

#### test_minimum()

#### test_plot()

#### test_products_property()

#### test_reverse_convert()
## pymatgen.analysis.tests.test_local_env module


### _class_ pymatgen.analysis.tests.test_local_env.CovalentBondNNTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_bonded_structure()

#### test_nn_length()

#### test_nn_orders()

### _class_ pymatgen.analysis.tests.test_local_env.Critic2NNTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_cn()

### _class_ pymatgen.analysis.tests.test_local_env.CrystalNNTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_cation_anion()

#### test_discrete_cn()

#### test_fixed_length()

#### test_get_bonded_structure()

#### test_get_cn()

#### test_noble_gas_material()

#### test_sanity()

#### test_shifted_sites()

#### test_weighted_cn()

#### test_weighted_cn_no_oxid()

#### test_x_diff_weight()

### _class_ pymatgen.analysis.tests.test_local_env.CutOffDictNNTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_cn()

#### test_from_preset()

### _class_ pymatgen.analysis.tests.test_local_env.JmolNNTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_get_nn()

### _class_ pymatgen.analysis.tests.test_local_env.LocalStructOrderParamsTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_get_order_parameters()

#### test_init()

### _class_ pymatgen.analysis.tests.test_local_env.MetalEdgeExtenderTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_custom_coordinators()

#### test_custom_cutoff()

#### test_custom_metals()

#### test_metal_edge_extender()

#### test_oxygen_edge_extender()

### _class_ pymatgen.analysis.tests.test_local_env.MiniDistNNTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_all_nn_classes()

#### test_get_local_order_params()

### _class_ pymatgen.analysis.tests.test_local_env.MotifIdentificationTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_get_neighbors_of_site_with_index()

#### test_site_is_of_motif_type()

### _class_ pymatgen.analysis.tests.test_local_env.NearNeighborTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### set_nn_info()

#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_on_disorder_options()

### _class_ pymatgen.analysis.tests.test_local_env.OpenBabelNNTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_nn_length()

#### test_nn_orders()

### _class_ pymatgen.analysis.tests.test_local_env.TestIsayevNN(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_get_nn()

### _class_ pymatgen.analysis.tests.test_local_env.ValenceIonicRadiusEvaluatorTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Setup MgO rocksalt structure for testing Vacancy


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_radii_ionic_structure()

#### test_valences_ionic_structure()

### _class_ pymatgen.analysis.tests.test_local_env.VoronoiNNTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_Cs2O()
A problematic structure in the Materials Project


#### test_adj_neighbors()

#### test_all_at_once()

#### test_filtered()

#### test_get_cn()

#### test_get_coordinated_sites()

#### test_get_voronoi_polyhedra()

#### test_nn_shell()

#### test_solid_angle()

#### test_volume()
## pymatgen.analysis.tests.test_molecule_matcher module


### _class_ pymatgen.analysis.tests.test_molecule_matcher.BruteForceOrderMatcherSi2OTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### _classmethod_ setUpClass()
Hook method for setting up class fixture before running tests in the class.


#### test_mismatched_atoms()

#### test_permuted_atoms_order()

#### test_perturbed_atom_position()

#### test_random_match()

#### test_rotated_molecule()

### _class_ pymatgen.analysis.tests.test_molecule_matcher.BruteForceOrderMatcherSiTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### _classmethod_ setUpClass()
Hook method for setting up class fixture before running tests in the class.


#### test_random_match()

#### test_to_and_from_dict()

### _class_ pymatgen.analysis.tests.test_molecule_matcher.GeneticOrderMatcherSi2OTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### _classmethod_ setUpClass()
Hook method for setting up class fixture before running tests in the class.


#### test_mismatched_atoms()

#### test_permuted_atoms_order()

#### test_perturbed_atom_position()

#### test_random_match()

#### test_rotated_molecule()

### _class_ pymatgen.analysis.tests.test_molecule_matcher.GeneticOrderMatcherSiTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### _classmethod_ setUpClass()
Hook method for setting up class fixture before running tests in the class.


#### test_mismatched_atoms()

#### test_permuted_atoms_order()

#### test_perturbed_atom_position()

#### test_random_match()

#### test_rotated_molecule()

#### test_to_and_from_dict()

### _class_ pymatgen.analysis.tests.test_molecule_matcher.GeneticOrderMatcherTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_fit()

#### test_get_rmsd()

#### test_mismatched_atom_composition()

#### test_rotated_molecule()

#### test_to_and_from_dict()

### _class_ pymatgen.analysis.tests.test_molecule_matcher.HungarianOrderMatcherSi2OTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### _classmethod_ setUpClass()
Hook method for setting up class fixture before running tests in the class.


#### test_mismatched_atoms()

#### test_permuted_atoms_order()

#### test_perturbed_atom_position()

#### test_random_match()

#### test_rotated_molecule()

### _class_ pymatgen.analysis.tests.test_molecule_matcher.HungarianOrderMatcherSiTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### _classmethod_ setUpClass()
Hook method for setting up class fixture before running tests in the class.


#### test_mismatched_atoms()

#### test_permuted_atoms_order()

#### test_perturbed_atom_position()

#### test_random_match()

#### test_rotated_molecule()

#### test_to_and_from_dict()

### _class_ pymatgen.analysis.tests.test_molecule_matcher.HungarianOrderMatcherTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_fit()

#### test_get_rmsd()

#### test_mismatched_atom_composition()

#### test_rotated_molecule()

#### test_to_and_from_dict()

### _class_ pymatgen.analysis.tests.test_molecule_matcher.KabschMatcherSi2OTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### _classmethod_ setUpClass()
Hook method for setting up class fixture before running tests in the class.


#### test_mismatched_atoms()

#### test_permuted_atoms_order()

#### test_perturbed_atom_position()

#### test_rotated_molecule()

### _class_ pymatgen.analysis.tests.test_molecule_matcher.KabschMatcherSiTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### _classmethod_ setUpClass()
Hook method for setting up class fixture before running tests in the class.


#### test_mismatched_atoms()

#### test_permuted_atoms_order()

#### test_perturbed_atom_position()

#### test_rotated_molecule()

#### test_to_and_from_dict()

### _class_ pymatgen.analysis.tests.test_molecule_matcher.KabschMatcherTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_fit()

#### test_get_rmsd()

#### test_mismatched_atom_composition()

#### test_mismatched_atom_order()

#### test_rotated_molecule()

#### test_to_and_from_dict()

### _class_ pymatgen.analysis.tests.test_molecule_matcher.MoleculeMatcherTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### fit_with_mapper(mapper)

#### test_cdi_23()

#### test_fit()

#### test_get_rmsd()

#### test_group_molecules()

#### test_strange_inchi()

#### test_thiane()

#### test_thiane_ethynyl()

#### test_to_and_from_dict()

### pymatgen.analysis.tests.test_molecule_matcher.generate_Si2O_cluster()

### pymatgen.analysis.tests.test_molecule_matcher.generate_Si_cluster()

### pymatgen.analysis.tests.test_molecule_matcher.permute(mol, seed)
Performs a random permutation of the sites in a structure.


* **Parameters**

    **seed** (*int*) – The seed value for the random generator.



### pymatgen.analysis.tests.test_molecule_matcher.perturb(mol, scale, seed)
Performs a random perturbation of the sites in a structure.


* **Parameters**


    * **scale** (*float*) – Distance in angstroms by which to perturb each site.


    * **seed** (*int*) – The seed value for the random generator.



### pymatgen.analysis.tests.test_molecule_matcher.rotate(mol, seed)
Performs a random rotation of the sites in a structure.


* **Parameters**


    * **mol** ([*Molecule*](pymatgen.core.md#pymatgen.core.structure.Molecule)) – The Molecule object which will be transformed.


    * **seed** (*int*) – The seed value for the random generator.


## pymatgen.analysis.tests.test_molecule_structure_comparator module


### _class_ pymatgen.analysis.tests.test_molecule_structure_comparator.TestMoleculeStructureComparator(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_are_equal()

#### test_get_13_bonds()

#### test_get_bonds()

#### test_to_and_from_dict()
## pymatgen.analysis.tests.test_nmr module


### _class_ pymatgen.analysis.tests.test_nmr.TestChemicalShieldingNotation(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_construction()

#### test_notations()

#### test_principal_axis_system()

### _class_ pymatgen.analysis.tests.test_nmr.TestElectricFieldGradient(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_Attributes()

#### test_construction()

#### test_principal_axis_system()
## pymatgen.analysis.tests.test_path_finder module


### _class_ pymatgen.analysis.tests.test_path_finder.PathFinderTest(methodName='runTest')
Bases: `TestCase`

Uses Li migration in LiFePO4

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_image_num()
## pymatgen.analysis.tests.test_phase_diagram module


### _class_ pymatgen.analysis.tests.test_phase_diagram.CompoundPhaseDiagramTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_get_formation_energy()

#### test_stable_entries()

#### test_str()

### _class_ pymatgen.analysis.tests.test_phase_diagram.GrandPotentialPhaseDiagramTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_get_formation_energy()

#### test_stable_entries()

#### test_str()

### _class_ pymatgen.analysis.tests.test_phase_diagram.PDEntryTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_get_chemical_energy()

#### test_get_composition()

#### test_get_energy()

#### test_get_energy_per_atom()

#### test_get_name()

#### test_is_element()

#### test_read_csv()

#### test_str()

#### test_to_from_dict()

### _class_ pymatgen.analysis.tests.test_phase_diagram.PDPlotterTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_mpl_plots()

#### test_pd_plot_data()

#### test_plot_pd_with_no_unstable()

#### test_plotly_plots()

### _class_ pymatgen.analysis.tests.test_phase_diagram.PatchedPhaseDiagramTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_contains()

#### test_dimensionality()

#### test_get_decomp_and_e_above_hull()

#### test_get_decomposition()

#### test_get_equilibrium_reaction_energy()

#### test_get_form_energy()

#### test_get_hull_energy()

#### test_get_pd_for_entry()

#### test_get_phase_separation_energy()

#### test_get_qhull_entries()

#### test_get_stable_entries()

#### test_getitem()

#### test_iter()

#### test_len()

#### test_raises_on_missing_terminal_entries()

#### test_repr()

#### test_setitem_and_delitem()

#### test_to_from_dict()

### _class_ pymatgen.analysis.tests.test_phase_diagram.PhaseDiagramTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_1d_pd()

#### test_all_entries_hulldata()

#### test_dim1()

#### test_downstream_methods_can_also_ignore_errors()

#### test_el_refs()

#### test_get_all_chempots()

#### test_get_composition_chempots()

#### test_get_critical_compositions()

#### test_get_critical_compositions_fractional()

#### test_get_decomp_and_e_above_hull_on_error()

#### test_get_decomposition()

#### test_get_e_above_hull()

#### test_get_element_profile()

#### test_get_equilibrium_reaction_energy()

#### test_get_formation_energy()

#### test_get_get_chempot_range_map()

#### test_get_hull_energy()

#### test_get_hull_energy_per_atom()

#### test_get_phase_separation_energy()

#### test_get_plot()

#### test_get_transition_chempots()

#### test_getmu_range_stability_phase()

#### test_getmu_vertices_stability_phase()

#### test_init()

#### test_ordering()

#### test_planar_inputs()

#### test_read_json()

#### test_repr()

#### test_stable_entries()

#### test_str()

#### test_to_from_dict()

#### test_val_err_on_no_entries()

### _class_ pymatgen.analysis.tests.test_phase_diagram.ReactionDiagramTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_formula()

#### test_get_compound_pd()

### _class_ pymatgen.analysis.tests.test_phase_diagram.TransformedPDEntryTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_get_composition()

#### test_get_energy()

#### test_get_energy_per_atom()

#### test_get_name()

#### test_is_element()

#### test_normalize()

#### test_str()

#### test_to_from_dict()

### _class_ pymatgen.analysis.tests.test_phase_diagram.UtilityFunctionTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_tet_coord()

#### test_triangular_coord()

#### test_unique_lines()
## pymatgen.analysis.tests.test_piezo module

Test for the piezo tensor class


### _class_ pymatgen.analysis.tests.test_piezo.PiezoTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_from_vasp_voigt()

#### test_from_voigt()

#### test_new()
## pymatgen.analysis.tests.test_piezo_sensitivity module

Test for the piezo tensor class


### _class_ pymatgen.analysis.tests.test_piezo_sensitivity.PiezoSensitivityTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_BornEffectiveChargeTensor()

#### test_ForceConstantMatrix()

#### test_InternalStrainTensor()

#### test_get_BEC_operations()

#### test_get_FCM_operations()

#### test_get_FCM_symmetry()

#### test_get_asum_FCM()

#### test_get_piezo()

#### test_get_rand_BEC()

#### test_get_rand_IST()

#### test_get_stable_FCM()

#### test_get_unstable_FCM()

#### test_rand_FCM()

#### test_rand_piezo()
## pymatgen.analysis.tests.test_pourbaix_diagram module


### _class_ pymatgen.analysis.tests.test_pourbaix_diagram.PourbaixDiagramTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### _classmethod_ setUpClass()
Hook method for setting up class fixture before running tests in the class.


#### test_get_decomposition()

#### test_get_pourbaix_domains()

#### test_get_stable_entry()

#### test_multicomponent()

#### test_multielement_parallel()

#### test_pourbaix_diagram()

#### test_properties()

#### test_serialization()

#### test_solid_filter()

### _class_ pymatgen.analysis.tests.test_pourbaix_diagram.PourbaixEntryTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_calc_coeff_terms()

#### test_energy_functions()

#### test_get_elt_fraction()

#### test_multi_entry()

#### test_pourbaix_entry()

#### test_to_from_dict()

### _class_ pymatgen.analysis.tests.test_pourbaix_diagram.PourbaixPlotterTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_plot_entry_stability()

#### test_plot_pourbaix()
## pymatgen.analysis.tests.test_prototypes module


### _class_ pymatgen.analysis.tests.test_prototypes.AflowPrototypeMatcherTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_prototype_matching()
## pymatgen.analysis.tests.test_quasiharmonic_debye_approx module


### _class_ pymatgen.analysis.tests.test_quasiharmonic_debye_approx.TestAnharmonicQuasiharmociDebyeApprox(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_debye_temperature()

#### test_gruneisen_parameter()

#### test_optimum_volume()

#### test_thermal_conductivity()

#### test_vibrational_free_energy()

#### test_vibrational_internal_energy()

### _class_ pymatgen.analysis.tests.test_quasiharmonic_debye_approx.TestQuasiharmociDebyeApprox(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_bulk_modulus()

#### test_debye_temperature()

#### test_gruneisen_parameter()

#### test_optimum_volume()

#### test_thermal_conductivity()

#### test_vibrational_free_energy()

#### test_vibrational_internal_energy()
## pymatgen.analysis.tests.test_reaction_calculator module


### _class_ pymatgen.analysis.tests.test_reaction_calculator.BalancedReactionTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_from_str()

#### test_init()

#### test_remove_spectator_species()

#### test_to_from_dict()

### _class_ pymatgen.analysis.tests.test_reaction_calculator.ComputedReactionTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_all_entries()

#### test_calculated_reaction_energy()

#### test_calculated_reaction_energy_uncertainty()

#### test_calculated_reaction_energy_uncertainty_for_nan()

#### test_calculated_reaction_energy_uncertainty_for_no_uncertainty()

#### test_init()

#### test_to_from_dict()

### _class_ pymatgen.analysis.tests.test_reaction_calculator.ReactionTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_as_entry()

#### test_calculate_energy()

#### test_equals()

#### test_init()

#### test_normalize_to()

#### test_overdetermined()

#### test_products_reactants()

#### test_rank()

#### test_scientific_notation()

#### test_singular_case()

#### test_to_from_dict()

#### test_underdetermined()

#### test_underdetermined_reactants()
## pymatgen.analysis.tests.test_structure_analyzer module


### _class_ pymatgen.analysis.tests.test_structure_analyzer.MiscFunctionTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_average_coordination_number()

#### test_contains_peroxide()

#### test_oxide_type()

#### test_solid_angle()

#### test_sulfide_type()

### _class_ pymatgen.analysis.tests.test_structure_analyzer.RelaxationAnalyzerTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_get_percentage_bond_dist_changes()

#### test_vol_and_para_changes()

### _class_ pymatgen.analysis.tests.test_structure_analyzer.VoronoiAnalyzerTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_analyze()

### _class_ pymatgen.analysis.tests.test_structure_analyzer.VoronoiConnectivityTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_connectivity_array()
## pymatgen.analysis.tests.test_structure_matcher module


### _class_ pymatgen.analysis.tests.test_structure_matcher.StructureMatcherTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_as_dict_and_from_dict()

#### test_cart_dists()

#### test_class()

#### test_cmp_fstruct()

#### test_disordered_primitive_to_ordered_supercell()

#### test_disordered_to_disordered()

#### test_electronegativity()

#### test_find_match1()

#### test_find_match2()

#### test_fit()
Take two known matched structures


    1. Ensure match


    2. Ensure match after translation and rotations


    3. Ensure no-match after large site translation


    4. Ensure match after site shuffling


#### test_get_lattices()

#### test_get_mapping()

#### test_get_mask()

#### test_get_s2_large_s2()

#### test_get_supercell_matrix()

#### test_get_supercell_size()

#### test_get_supercells()

#### test_ignore_species()

#### test_left_handed_lattice()
Ensure Left handed lattices are accepted


#### test_mix()

#### test_no_scaling()

#### test_occupancy_comparator()

#### test_ordered_primitive_to_disordered_supercell()

#### test_out_of_cell_s2_like_s1()

#### test_oxi()
Test oxidation state removal matching


#### test_primitive()
Test primitive cell reduction


#### test_rms_vs_minimax()

#### test_subset()

#### test_supercell_fit()

#### test_supercell_subsets()
## pymatgen.analysis.tests.test_surface_analysis module


### _class_ pymatgen.analysis.tests.test_surface_analysis.NanoscaleStabilityTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_scaled_wulff()

#### test_stability_at_r()

### _class_ pymatgen.analysis.tests.test_surface_analysis.SlabEntryTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_cleaned_up_slab()

#### test_create_slab_label()

#### test_properties()

#### test_surface_energy()

### _class_ pymatgen.analysis.tests.test_surface_analysis.SurfaceEnergyPlotterTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_color_palette_dict()

#### test_entry_dict_from_list()

#### test_get_stable_entry_at_u()

#### test_get_surface_equilibrium()

#### test_stable_u_range_dict()

#### test_wulff_from_chempot()

### _class_ pymatgen.analysis.tests.test_surface_analysis.WorkfunctionAnalyzerTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_is_converged()

#### test_shift()

### pymatgen.analysis.tests.test_surface_analysis.get_entry_dict(filename)

### pymatgen.analysis.tests.test_surface_analysis.get_path(path_str)

### pymatgen.analysis.tests.test_surface_analysis.load_O_adsorption()
## pymatgen.analysis.tests.test_transition_state module


### _class_ pymatgen.analysis.tests.test_transition_state.NEBAnalysisTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### runTest()

#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_combine_neb_plots()
## pymatgen.analysis.tests.test_wulff module


### _class_ pymatgen.analysis.tests.test_wulff.WulffShapeTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### consistency_tests()

#### setUp()
Hook method for setting up the test fixture before exercising it.


#### symm_check(ucell, wulff_vertices)
# Checks if the point group of the Wulff shape matches
# the point group of its conventional unit cell


* **Parameters**


    * **ucell** (*str*) – Unit cell that the Wulff shape is based on.


    * **wulff_vertices** (*list*) – List of all vertices on the Wulff
    shape. Use wulff.wulff_pt_list to obtain the list
    (see wulff_generator.py).


return (bool)


#### symmetry_test()

#### test_corner_and_edges()

#### test_get_azimuth_elev()

#### test_get_plot()

#### test_get_plotly()

#### test_properties()
## pymatgen.analysis.tests.test_xps module


### _class_ pymatgen.analysis.tests.test_xps.XPSTestCase(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_from_dos()