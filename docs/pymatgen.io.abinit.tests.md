---
layout: default
title: pymatgen.io.abinit.tests.md
nav_exclude: true
---

# pymatgen.io.abinit.tests package


## pymatgen.io.abinit.tests.test_abiobjects module


### _class_ pymatgen.io.abinit.tests.test_abiobjects.ElectronsAlgorithmTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_base()

### _class_ pymatgen.io.abinit.tests.test_abiobjects.ElectronsTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_base()

### _class_ pymatgen.io.abinit.tests.test_abiobjects.KSamplingTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_base()

### _class_ pymatgen.io.abinit.tests.test_abiobjects.LatticeFromAbivarsTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_rprim_acell()

#### test_znucl_typat()
Test the order of typat and znucl in the Abinit input and enforce_typat, enforce_znucl.


### _class_ pymatgen.io.abinit.tests.test_abiobjects.PPModelTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_base()

### _class_ pymatgen.io.abinit.tests.test_abiobjects.RelaxationTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_base()

### _class_ pymatgen.io.abinit.tests.test_abiobjects.SmearingTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_base()

### _class_ pymatgen.io.abinit.tests.test_abiobjects.SpinModeTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_base()
## pymatgen.io.abinit.tests.test_inputs module


### _class_ pymatgen.io.abinit.tests.test_inputs.AbinitInputTestCase(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Unit tests for BasicAbinitInput.

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_api()
Testing BasicAbinitInput API.


#### test_helper_functions()
Testing BasicAbinitInput helper functions.


#### test_input_errors()
Testing typical BasicAbinitInput Error


### _class_ pymatgen.io.abinit.tests.test_inputs.FactoryTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_ebands_input()
Testing ebands_input factory.


#### test_gs_input()
Testing gs_input factory.


#### test_ion_ioncell_relax_input()
Testing ion_ioncell_relax_input factory.


### _class_ pymatgen.io.abinit.tests.test_inputs.ShiftModeTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_shiftmode()
Testing shiftmode


### _class_ pymatgen.io.abinit.tests.test_inputs.TestMultiDataset(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Unit tests for BasicMultiDataset.

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_api()
Testing BasicMultiDataset API.


### pymatgen.io.abinit.tests.test_inputs.abiref_file(filename)
Return absolute path to filename in ~pymatgen/test_files/abinit


### pymatgen.io.abinit.tests.test_inputs.abiref_files(\*filenames)
Return list of absolute paths to filenames in ~pymatgen/test_files/abinit

## pymatgen.io.abinit.tests.test_netcdf module


### _class_ pymatgen.io.abinit.tests.test_netcdf.ETSF_Reader_TestCase(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_read_Si2()

### _class_ pymatgen.io.abinit.tests.test_netcdf.TestAbinitHeader(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_api()

### pymatgen.io.abinit.tests.test_netcdf.ref_file(filename)
## pymatgen.io.abinit.tests.test_pseudos module


### _class_ pymatgen.io.abinit.tests.test_pseudos.PseudoTableTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_methods()
Test PseudoTable methods


### _class_ pymatgen.io.abinit.tests.test_pseudos.PseudoTestCase(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_nc_pseudos()
Test norm-conserving pseudopotentials


#### test_oncvpsp_pseudo_fr()
Test the ONCVPSP Pb pseudo (relativistic version with SO).


#### test_oncvpsp_pseudo_sr()
Test the ONCVPSP Ge pseudo (scalar relativistic version).


#### test_pawxml_pseudos()
Test O.GGA_PBE-JTH-paw.xml.


### pymatgen.io.abinit.tests.test_pseudos.ref_file(filename)

### pymatgen.io.abinit.tests.test_pseudos.ref_files(\*filenames)