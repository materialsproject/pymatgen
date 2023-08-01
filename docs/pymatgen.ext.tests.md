---
layout: default
title: pymatgen.ext.tests.md
nav_exclude: true
---

# pymatgen.ext.tests namespace


## pymatgen.ext.tests.test_cod module


### _class_ pymatgen.ext.tests.test_cod.CODTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_get_cod_ids()

#### test_get_structure_by_formula()

#### test_get_structure_by_id()
## pymatgen.ext.tests.test_matproj module


### _class_ pymatgen.ext.tests.test_matproj.MPResterOldTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### pytestmark(_ = [Mark(name='skipif', args=(False,), kwargs={'reason': 'PMG_MAPI_KEY environment variable not set or MP API is down.'})_ )

#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_api_key_is_none()

#### test_database_version()

#### test_download_info()

#### test_find_structure()

#### test_get_all_materials_ids_doc()

#### test_get_bandstructure_by_material_id()

#### test_get_cohesive_energy()

#### test_get_data()

#### test_get_dos_by_id()

#### test_get_entries()

#### test_get_entries_in_chemsys()

#### test_get_entry_by_material_id()

#### test_get_exp_entry()

#### test_get_exp_thermo_data()

#### test_get_gb_data()

#### test_get_interface_reactions()

#### test_get_materials_id_from_task_id()

#### test_get_materials_id_references()

#### test_get_phonon_data_by_material_id()

#### test_get_pourbaix_entries()

#### test_get_reaction()

#### test_get_stability()

#### test_get_structure_by_material_id()

#### test_get_structures()

#### test_get_substrates()

#### test_get_surface_data()

#### test_get_wulff_shape()

#### test_get_xas_data()

#### test_include_user_agent()

#### test_parse_criteria()

#### test_pourbaix_heavy()

#### test_pourbaix_mpr_pipeline()

#### test_query()

#### test_query_chunk_size()
## pymatgen.ext.tests.test_optimade module


### _class_ pymatgen.ext.tests.test_optimade.OptimadeTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_build_filter()

#### test_get_snls_mp()

#### test_get_structures_mp()