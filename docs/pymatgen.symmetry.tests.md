---
layout: default
title: pymatgen.symmetry.tests.md
nav_exclude: true
---

# pymatgen.symmetry.tests package


## pymatgen.symmetry.tests.test_analyzer module


### _class_ pymatgen.symmetry.tests.test_analyzer.FuncTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_cluster_sites()

### _class_ pymatgen.symmetry.tests.test_analyzer.PointGroupAnalyzerTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_asym_top()

#### test_cyclic()

#### test_dihedral()

#### test_get_kpoint_weights()

#### test_linear()

#### test_spherical()

#### test_symmetrize_molecule1()

#### test_symmetrize_molecule2()

#### test_tricky()

### _class_ pymatgen.symmetry.tests.test_analyzer.SpacegroupAnalyzerTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_find_primitive()
F m -3 m Li2O testing of converting to primitive cell.


#### test_get_conventional_standard_structure()

#### test_get_crystal_system()

#### test_get_hall()

#### test_get_ir_reciprocal_mesh()

#### test_get_ir_reciprocal_mesh_map()

#### test_get_point_group_operations()

#### test_get_point_group_operations_uniq()

#### test_get_pointgroup()

#### test_get_primitive_standard_structure()

#### test_get_refined_structure()

#### test_get_space_number()

#### test_get_space_symbol()

#### test_get_symmetry()

#### test_get_symmetry_dataset()

#### test_init_cell()

#### test_is_laue()

#### test_magnetic()

#### test_primitive()

#### test_symmetrized_structure()

#### test_tricky_structure()

### _class_ pymatgen.symmetry.tests.test_analyzer.SpacegroupTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_are_symmetrically_equivalent()
## pymatgen.symmetry.tests.test_groups module


### _class_ pymatgen.symmetry.tests.test_groups.PointGroupTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_get_orbit()

#### test_is_sub_super_group()

#### test_order()

### _class_ pymatgen.symmetry.tests.test_groups.SpaceGroupTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_abbrev_symbols()

#### test_attr()

#### test_crystal_system()

#### test_full_symbols()

#### test_get_orbit()

#### test_get_orbit_and_generators()

#### test_get_settings()

#### test_hexagonal()

#### test_is_compatible()

#### test_order_symm_ops()

#### test_other_settings()

#### test_point_group_is_set()

#### test_raises_on_bad_int_number()

#### test_renamed_e_symbols()

#### test_repr()

#### test_string()

#### test_subgroup_supergroup()

#### test_symmops()
## pymatgen.symmetry.tests.test_kpath_hin module


### _class_ pymatgen.symmetry.tests.test_kpath_hin.KPathSeekTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_kpath_acentered()

#### test_kpath_generation()
## pymatgen.symmetry.tests.test_kpath_lm module


### _class_ pymatgen.symmetry.tests.test_kpath_lm.KPathLatimerMunroTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_kpath_acentered()

#### test_kpath_generation()

#### test_magnetic_kpath_generation()
## pymatgen.symmetry.tests.test_kpath_sc module


### _class_ pymatgen.symmetry.tests.test_kpath_sc.BandStructureSCTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_kpath_acentered()

#### test_kpath_generation()
## pymatgen.symmetry.tests.test_kpaths module


### _class_ pymatgen.symmetry.tests.test_kpaths.HighSymmKpathTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_continuous_kpath()

#### test_kpath_generation()
## pymatgen.symmetry.tests.test_maggroups module


### _class_ pymatgen.symmetry.tests.test_maggroups.MagneticSpaceGroupTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_crystal_system()

#### test_equivalence_to_spacegroup()

#### test_init()

#### test_is_compatible()

#### test_sg_symbol()

#### test_str()

#### test_symmetry_ops()
## pymatgen.symmetry.tests.test_settings module


### _class_ pymatgen.symmetry.tests.test_settings.JonesFaithfulTransformationTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_init()

#### test_inverse()

#### test_transform_coords()

#### test_transform_lattice()

#### test_transform_symmops()
## pymatgen.symmetry.tests.test_site_symmetries module


### _class_ pymatgen.symmetry.tests.test_site_symmetries.SiteSymmetriesTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_get_shared_symmetries_operations()

#### test_get_site_symmetries()