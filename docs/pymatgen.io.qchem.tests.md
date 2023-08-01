---
layout: default
title: pymatgen.io.qchem.tests.md
nav_exclude: true
---

# pymatgen.io.qchem.tests package

This package implements modules test the input and output modules for QChem


## pymatgen.io.qchem.tests.test_inputs module


### _class_ pymatgen.io.qchem.tests.test_inputs.TestQCInput(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test__str__()

#### test_almo_template()

#### test_cdft_template()

#### test_find_sections()

#### test_from_multi_jobs_file()

#### test_from_string()

#### test_molecule_template()

#### test_multi_job_string()

#### test_multi_molecule_template()

#### test_opt_template()

#### test_pcm_nonels_template()

#### test_pcm_template()

#### test_read_almo()

#### test_read_bad_pcm()

#### test_read_bad_scan()

#### test_read_bad_smx()

#### test_read_bad_solvent()

#### test_read_cdft()

#### test_read_molecule()

#### test_read_multi_molecule()

#### test_read_nbo()

#### test_read_negative()

#### test_read_only_rem()

#### test_read_opt()

#### test_read_pcm()

#### test_read_pcm_nonels()

#### test_read_plots()

#### test_read_rem()

#### test_read_scan()

#### test_read_smx()

#### test_read_solvent()

#### test_read_svp()

#### test_read_write_custom_smd()

#### test_read_write_nbo7()

#### test_read_write_nbo_e2pert()

#### test_rem_template()

#### test_scan_template()

#### test_smx_template()

#### test_solvent_template()

#### test_svp_template()

#### test_van_der_waals_template()

#### test_write_file_from_OptSet()

#### test_write_file_from_OptSet_with_vdw()
## pymatgen.io.qchem.tests.test_outputs module


### _class_ pymatgen.io.qchem.tests.test_outputs.TestQCOutput(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### _static_ generate_multi_job_dict()
Used to generate test dictionary for multiple jobs.


#### _static_ generate_single_job_dict()
Used to generate test dictionary for single jobs.


#### test_NBO5_vs_NBO7_hybridization_character()

#### test_NBO7_infinite_e2pert()

#### test_NBO7_parsing()

#### test_NBO_3C()

#### test_NBO_hyperbonds()

#### test_NBO_parsing()

#### test_all()

#### test_almo_msdft2_parsing()

#### test_cdft_dc_parsing()

#### test_cdft_parsing()

#### test_cmirs_benzene()

#### test_cmirs_dielst10()

#### test_cmirs_water()

#### test_fodft_parsing()

#### test_isosvp_dielst10()

#### test_isosvp_water()

#### test_pod_parsing()

#### test_structural_change()
## pymatgen.io.qchem.tests.test_sets module


### _class_ pymatgen.io.qchem.tests.test_sets.ForceSetTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_init()

#### test_pcm_init()

#### test_smd_init()

### _class_ pymatgen.io.qchem.tests.test_sets.FreqSetTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_init()

#### test_pcm_init()

#### test_smd_init()

### _class_ pymatgen.io.qchem.tests.test_sets.OptSetTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_init()

#### test_nbo_init()

#### test_overwrite_opt_input()

#### test_pcm_init()

#### test_smd_init()

#### test_v5_vs_v6()

### _class_ pymatgen.io.qchem.tests.test_sets.PESScanSetTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_init()

#### test_pcm_init()

#### test_smd_init()

### _class_ pymatgen.io.qchem.tests.test_sets.QChemDictSetTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_cmirs_write()
Also tests overwrite_inputs with a RHOISO value


#### test_custom_smd_write()

#### test_double_solvation()

#### test_full_init()

#### test_init()

#### test_isosvp_write()
Also tests overwrite_inputs with a RHOISO value


#### test_overwrite_input()

#### test_pcm_write()

#### test_smd_write()

#### test_solvation_warnings()
Tests warnings / errors resulting from nonsensical overwrite_inputs


### _class_ pymatgen.io.qchem.tests.test_sets.SinglePointSetTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_cmirs_init()

#### test_init()

#### test_isosvp_init()

#### test_pcm_init()

#### test_plots_init()

#### test_scf_extra_print()

#### test_smd_init()

### _class_ pymatgen.io.qchem.tests.test_sets.TransitionStateSetTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_init()

#### test_pcm_init()

#### test_smd_init()
## pymatgen.io.qchem.tests.test_utils module


### _class_ pymatgen.io.qchem.tests.test_utils.UtilTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

test utils

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_lower_and_check_unique()

#### test_process_parsed_HESS()