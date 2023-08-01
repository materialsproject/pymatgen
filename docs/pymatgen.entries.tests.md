---
layout: default
title: pymatgen.entries.tests.md
nav_exclude: true
---

# pymatgen.entries.tests package

tests for pymatgen.entries.


## pymatgen.entries.tests.test_compatibility module


### _class_ pymatgen.entries.tests.test_compatibility.AqueousCorrectionTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_compound_energy()

### _class_ pymatgen.entries.tests.test_compatibility.CorrectionErrors2020CompatibilityTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_errors()

### _class_ pymatgen.entries.tests.test_compatibility.CorrectionSpecificityTest(methodName='runTest')
Bases: `TestCase`

Make sure corrections are only applied to GGA or GGA+U entries.

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_correction_specificity()

### _class_ pymatgen.entries.tests.test_compatibility.DummyCompatibility()
Bases: [`Compatibility`](pymatgen.entries.md#pymatgen.entries.compatibility.Compatibility)

Dummy class to test abstract Compatibility interface.


#### get_adjustments(entry)
Get the energy adjustments for a ComputedEntry.

This method must generate a list of EnergyAdjustment objects
of the appropriate type (constant, composition-based, or temperature-based)
to be applied to the ComputedEntry, and must raise a CompatibilityError
if the entry is not compatible.


* **Parameters**

    **entry** – A ComputedEntry object.



* **Returns**

    A list of EnergyAdjustment to be applied to the

        Entry.




* **Return type**

    list[[EnergyAdjustment](pymatgen.entries.md#pymatgen.entries.computed_entries.EnergyAdjustment)]



* **Raises**

    **CompatibilityError if the entry is not compatible** –



### _class_ pymatgen.entries.tests.test_compatibility.MITAqueousCompatibilityTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_aqueous_compat()

#### test_dont_error_on_weird_elements()

#### test_msonable()

#### test_potcar_doenst_match_structure()

### _class_ pymatgen.entries.tests.test_compatibility.MITCompatibilityTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_U_value()

#### test_correction_value()

#### test_element_processing()

#### test_get_explanation_dict()

#### test_msonable()

#### test_potcar_doenst_match_structure()

#### test_potcar_spec_is_none()

#### test_process_entry()

#### test_revert_to_symbols()

#### test_same_potcar_symbol()

#### test_wrong_U_value()

#### test_wrong_psp()

### _class_ pymatgen.entries.tests.test_compatibility.MaterialsProjectCompatibility2020Test(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_U_values()

#### test_check_potcar()

#### test_config_file()

#### test_correction_values()

#### test_element_processing()

#### test_energy_adjustments()

#### test_get_explanation_dict()

#### test_msonable()

#### test_oxdiation()

#### test_oxdiation_by_electronegativity()

#### test_oxi_state_guess()

#### test_process_entries()

#### test_process_entry()

#### test_process_entry_with_oxidation_state()

#### test_processing_entries_inplace()

#### test_wrong_psp()

### _class_ pymatgen.entries.tests.test_compatibility.MaterialsProjectCompatibilityTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### tearDown()
Hook method for deconstructing the test fixture after testing it.


#### test_U_values()

#### test_correction_values()

#### test_element_processing()

#### test_get_corrections_dict()

#### test_get_explanation_dict()

#### test_msonable()

#### test_process_entries()

#### test_process_entry()

#### test_wrong_psp()

### _class_ pymatgen.entries.tests.test_compatibility.OxideTypeCorrectionNoPeroxideCorrTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_oxide_energy_corr()

#### test_ozonide()

#### test_peroxide_energy_corr()

### _class_ pymatgen.entries.tests.test_compatibility.OxideTypeCorrectionTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_no_struct_compat()

#### test_process_entry_oxide()

#### test_process_entry_ozonide()

#### test_process_entry_peroxide()

#### test_process_entry_superoxide()

### _class_ pymatgen.entries.tests.test_compatibility.SulfideTypeCorrection2020Test(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_struct_no_struct()

### _class_ pymatgen.entries.tests.test_compatibility.TestMaterialsProjectAqueousCompatibility()
Bases: `object`

Test MaterialsProjectAqueousCompatibility.

-x- formation energy of H2O should always be -2.458 eV/H2O
-x- H2 energy should always be the same value
-x- H2O energy should always be the same value
-x- Should get warnings if you init without all energy args
-x- Should get CompatibilityError if you get_entry without all energy args
-x- energy args should auto-populate from entries passed to process_entries
-x- check compound entropies appropriately added
-x- check hydrate adjustment appropriately applied

### Notes

Argument values from MaterialsProjectCompatibility as of April 2020:

    corrected DFT energy of H2O = -15.5875 eV/H2O (mp-697111) or -5.195 eV/atom
    corrected DFT energy of O2 = -4.9276 eV/atom (mp-12957)
    total energy corrections applied to H2O (eV/H2O) -0.70229 eV/H2O or -0.234 eV/atom


#### test_compound_entropy()

#### test_h_h2o_energy_no_args()

#### test_h_h2o_energy_with_args_multi()

#### test_h_h2o_energy_with_args_single()

#### test_hydrate_adjustment()

#### test_processing_entries_inplace()

### pymatgen.entries.tests.test_compatibility.test_clean_arg()
clean=False should preserve existing corrections, clean=True should delete
them before processing.


### pymatgen.entries.tests.test_compatibility.test_energy_adjustment_normalize()
Both manual and automatically generated energy adjustments should be scaled
by the normalize method.


### pymatgen.entries.tests.test_compatibility.test_no_duplicate_corrections()
Compatibility should never apply the same correction twice.


### pymatgen.entries.tests.test_compatibility.test_overlapping_adjustments()
Compatibility should raise a CompatibilityError if there is already a
correction with the same name, but a different value, and process_entries
should skip that entry.


### pymatgen.entries.tests.test_compatibility.test_process_entries_return_type()
process_entries should accept single entries or a list, and always return a list.

## pymatgen.entries.tests.test_computed_entries module


### _class_ pymatgen.entries.tests.test_computed_entries.ComputedEntryTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_composition()

#### test_conflicting_correction_adjustment()
Should raise a ValueError if a user tries to manually set both the correction
and energy_adjustment, even if the values match.


#### test_copy()

#### test_energy()

#### test_entry_id()

#### test_is_element()

#### test_normalize()

#### test_normalize_energy_adjustments()

#### test_per_atom_props()

#### test_str()

#### test_sulfide_energy()

#### test_to_from_dict()

#### test_to_from_dict_with_adjustment()
Legacy case where adjustment was provided manually.


#### test_to_from_dict_with_adjustment_2()
Modern case where correction was provided manually.


#### test_to_from_dict_with_adjustment_3()
Legacy case where the entry was serialized before the energy_adjustment
attribute was part of ComputedEntry.


### _class_ pymatgen.entries.tests.test_computed_entries.ComputedStructureEntryTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_composition()

#### test_copy()

#### test_energy()

#### test_eq()

#### test_str()

#### test_to_from_dict()

#### test_to_from_dict_structure_with_adjustment_3()
Legacy case where the structure entry was serialized before the energy_adjustment
attribute was part of ComputedEntry.


### _class_ pymatgen.entries.tests.test_computed_entries.GibbsComputedStructureEntryTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_expt_gas_entry()

#### test_from_entries()

#### test_from_pd()

#### test_gf_sisso()

#### test_interpolation()

#### test_normalize()

#### test_str()

#### test_to_from_dict()

### pymatgen.entries.tests.test_computed_entries.test_composition_energy_adjustment()

### pymatgen.entries.tests.test_computed_entries.test_constant_energy_adjustment()

### pymatgen.entries.tests.test_computed_entries.test_energy_adjustment()

### pymatgen.entries.tests.test_computed_entries.test_manual_energy_adjustment()

### pymatgen.entries.tests.test_computed_entries.test_temp_energy_adjustment()
## pymatgen.entries.tests.test_correction_calculator module


### _class_ pymatgen.entries.tests.test_correction_calculator.CorrectionCalculatorTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_missing_entry_response()
Test that correct error is raised (ValueError) if the input is missing a computed entry.


#### test_no_uncertainties()
Test that corrections can be calculated with no uncertainties.


#### test_normal_corrections()
Test that the values in MPCompatibility.yaml are reproduced correctly.


#### test_warnings_options()
Test that compounds can be included/excluded using the allow_{warning} optional parameters.

## pymatgen.entries.tests.test_entry_tools module


### _class_ pymatgen.entries.tests.test_entry_tools.EntrySetTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_as_dict()

#### test_chemsys()

#### test_get_subset()

#### test_remove_non_ground_states()

### _class_ pymatgen.entries.tests.test_entry_tools.FuncTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_group_entries_by_composition()

#### test_group_entries_by_structure()
## pymatgen.entries.tests.test_exp_entries module


### _class_ pymatgen.entries.tests.test_exp_entries.ExpEntryTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_energy()

#### test_str()

#### test_to_from_dict()
## pymatgen.entries.tests.test_mixing_scheme module

Tests for the Materials Project DFT mixing scheme.

**NOTE FOR FUTURE DEVELOPERS**

PLEASE DO NOT CHANGE THESE TESTS WITHOUT READING THIS ENTIRE DOCUMENTATION!

The MP DFT mixing scheme is conceptually simple but highly intricate to implement.
The tests in this file were crafted carefully over a period of many months
in order to capture as many nuances of the desired behavior as possible, and were
integral to developing the mixing code itself.

### Test Structure

The majority of the tests use “mixing states” to check behavior. Mixing states
are merely combinations of different ComputedStructureEntry with different run_type
- e.g. all ground states present for both run_types, only one run_type, etc. Mixing
states are defined using the MixingState utility class, which has attributes that
return 1) the GGA entries, 2) the R2SCAN entries, 3) all entries, and 4) the pandas
DataFrame that represents the mixing state. Note that ‘GGA’ and ‘R2SCAN’ are used
throughout this test file to represent run_type_1 and run_type_2, respectively,
but the mixing scheme is design to be able to mix any two functionals.

Most mixing states are subsets of the ms_complete mixing state. ms_complete
was crafted to capture most of the scenarios that may be encountered when mixing
ComputedStructureEntry from different functionals. It comprises a complete binary
phase diagram, with all entries present as both GGA and R2SCAN calculations.

Note that these entries are inspired by, but NOT equivalent to, the real SnBr2 phase
diagram. Rather than use real energies or structures, arbitrary energies and structures
have been used to keep this test file cleaner and easier to understand. The Bromine
structures are the one exception to this. These structures are taken from real calculations
and will fail to structure match with default arguments (they are included here to
test “fuzzy matching” behavior). Entry-id’s are assigned numerically, with the
same number corresponding to equivalent materials. So “gga-1” should have the
same structure as “r2scan-1”.

### Description of the ms_complete mixing state

Future developers are HIGHLY encouraged to plot the PhaseDiagram associated with both
the R2SCAN and the GGA entries in ms_complete before attempting to modify any of these
tests.

**GGA entries**


* gga-1: stable ground state of Sn


* gga-2: unstable polymorph of Br


* gga-3: stable polymorph of Br


* gga-4: stable polymorph of SnBr2


* gga-5: unstable polymorph of SnBr2, 1 eV/atom above hull


* gga-6: unstable polymorph of SnBr2, 2 eV/atom above hull


* gga-7: unstable composition (SnBr4), 0.6 eV/atom above hull

**R2SCAN entries**


* All the same entries exist as in GGA


* Entries with corresponding numbers have matching structures
(i.e. gga-1 and r2scan-1)


* Unless otherwise listed below, energies are 1 eV/atom lower than those in GGA


* for Br, the GGA ground state structure gga-3 does not match r2scan-3 unless you
use fuzzy matching


* for SnBr2, a different polymorph is stabilized than in GGA (r2scan-5 whereas
r2scan-4 was the GGA ground state)


* entry r2scan-6 (unstable SnBr2 polymorph) is scaled to 25% the size of gga-6.
This will match with default StructureMatcher settings but not with customized
settings.


* SnBr4 (r2scan-7) appears on the hull whereas it is not stable in GGA.


* A few mixing states (but not ms_complete) also include an unstable polymorph
of SnBr4, entry r2scan-8, that does not exist in GGA.

### Types of Tests

In general there are 3 types of tests. Most tests follow this pattern.


1. Unit tests of get_mixing_state_data, to verify that it generates
the correct DataFrame when provided with a particular list of entries


2. Unit tests that verify get_adjustments behaves correctly, when provided with
a ComputedEntry and a pre-generated mixing_scheme_state_data DataFrame


3. Functional tests of the behavior of the process_entries method, which is the
primary user-facing interface

### Implementation Notes


* Tests are organized into two classes. One class collects tests of the different
args / kwargs that can be passed to the mixing scheme. The other class collects
tests of different mixing states.


* Tests are written in pure pytest format, so the class organization is merely
a convenience to facilitate code folding, etc.


* pytest fixtures are used to define the various mixing states. Using pytest
fixtures is helpful here since it ensures every test receives fresh, unmodified
copies of the respective entries. process_entries modifies entry energies
in place, so tests could cross-contaminate one another if a fixture were not used.


### _class_ pymatgen.entries.tests.test_mixing_scheme.DummyCompatibility()
Bases: [`Compatibility`](pymatgen.entries.md#pymatgen.entries.compatibility.Compatibility)

Dummy class to test compat1 and compat2 kwargs.


#### get_adjustments(entry)
Get the energy adjustments for a ComputedEntry.

This method must generate a list of EnergyAdjustment objects
of the appropriate type (constant, composition-based, or temperature-based)
to be applied to the ComputedEntry, and must raise a CompatibilityError
if the entry is not compatible.


* **Parameters**

    **entry** – A ComputedEntry object.



* **Returns**

    A list of EnergyAdjustment to be applied to the

        Entry.




* **Return type**

    list[[EnergyAdjustment](pymatgen.entries.md#pymatgen.entries.computed_entries.EnergyAdjustment)]



* **Raises**

    **CompatibilityError if the entry is not compatible** –



### _class_ pymatgen.entries.tests.test_mixing_scheme.MixingState(gga_entries, scan_entries, dataframe)
Bases: `object`

Lightweight container class for a mixing state, used by the tests below.


* **Parameters**


    * **gga_entries** – list of run_type_1 (usually GGA or GGA+U) ComputedEntry


    * **scan_entries** – list of run_type_2 (usually R2SCAN) ComputedEntry


    * **DataFrame** – pandas DataFrame representing the mixing state, of the
    format returned by get_mixing_state_data.



#### _property_ all_entries()

### _class_ pymatgen.entries.tests.test_mixing_scheme.TestMaterialsProjectDFTMixingSchemeArgs()
Bases: `object`

Test behavior of args / kwargs passed to the mixing scheme.
In general, these tests are for the interface rather than the actual
mixing.


#### test_alternate_run_types()
Test behavior with alternate run_types.


#### test_alternate_structure_matcher(ms_complete)
Test alternate structure matcher kwargs. By setting scale to False, entries
gga-6 and r2scan-6 will no longer match and will be listed as separate rows
in the mixing state DataFrame.


#### test_check_potcar(ms_complete)
Entries with invalid or missing POTCAR raise error by default but should be ignored if
check_potcar=False in MaterialsProjectDFTMixingScheme.


#### test_clean(mixing_scheme_no_compat)

#### test_compat_args(ms_complete)
Test the behavior of compat1 and compat2 kwargs
The DummyCompatibility class defined in this test file should lower the
energies of all entries by 10 eV/atom.


#### test_empty_entries(mixing_scheme_no_compat)

#### test_fuzzy_matching(ms_complete)
Test fuzzy diatomic matching. If fuzzy matching is disabled, then entry r2scan-3 will not
match GGA ground state gga-3, preventing the mixing scheme from using the R2SCAN energy
scale.

In this situation, the mixing scheme should adjust the energies of all R2SCAN entries that
match GGA ones to be identical to the GGA’s, and discard the corresponding GGA entries. Entry
r2scan-3 should be discarded because there is no reference ground state energy for it.


#### test_incompatible_run_type(mixing_scheme_no_compat)
If entry.parameters.run_type is not “GGA”, “GGA+U”, or “R2SCAN”, we should get a
warning from process_entries and a CompatibilityError from get_adjustments.


#### test_multiple_matching_structures(mixing_scheme_no_compat, ms_complete_duplicate_structs)
Test behavior when the entries contain many structures that match to
the same material. For this test, entry gga-4 (SnBr2 ground state) and
its matching r2scan-4 are each duplicated into new entries gga-10 and
r2scan-10, respectively.

In this situation, the mixing scheme should only keep the lowest energy entry
within each run_type, so gga-10 and r2scan-10 should be discarded, and the
results should otherwise be identical to ms_scan_complete.


#### test_no_entry_ids(mixing_scheme_no_compat)
Unique entry_ids are required.


#### test_no_foreign_entries(mixing_scheme_no_compat, ms_complete)
If process_entries or get_adjustments is called with a populated mixing_state_data
kwarg and one or more of the entry_ids is not present in the mixing_state_data,
raise CompatbilityError.


#### test_no_mixing_data(ms_complete)
Test the behavior of get_adjustments when mixing_state_data is None.


#### test_no_run_type(mixing_scheme_no_compat)
If one of the entries doesn’t have a run_type attribute, we should get a warning
from process_entries and a CompatibilityError from get_adjustments.


#### test_no_single_entry(mixing_scheme_no_compat)
Raise CompatibilityError if process_entries is called on a single Entry
without any state_data.


#### test_no_structure(mixing_scheme_no_compat)

#### test_processing_entries_inplace()

#### test_same_run_type()
Test behavior when run_type_1 and run_type_2 are the same
or overlap.


### _class_ pymatgen.entries.tests.test_mixing_scheme.TestMaterialsProjectDFTMixingSchemeStates()
Bases: `object`

Test the behavior of the mixing scheme under different mixing states (i.e., different
combinations of GGA and R2SCAN entries).


#### test_chemsys_mismatch(mixing_scheme_no_compat, ms_scan_chemsys_superset)
Test what happens if the entries aren’t in the same chemsys. run_type_2 entries
that are outside the run_type_1 chemsys should be discarded.


#### test_state_all_gga_scan_gs(mixing_scheme_no_compat, ms_all_gga_scan_gs)
Mixing state in which we have a complete GGA PhaseDiagram and all GGA
ground states present as R2SCAN entries.

In this situation, we should build the hull using R2SCAN energies, discard
the GGA ground state materials, and correct the remaining GGA energies onto
the R2SCAN hull such that they maintain their e_above_hull.


#### test_state_complete_entries(mixing_scheme_no_compat, ms_complete)
Mixing state in which every material is present in both GGA and R2SCAN.

In this state, the mixing scheme should return only the R2SCAN entries, unmodified


#### test_state_energy_modified(mixing_scheme_no_compat, ms_complete)
Mixing state in which we try to process an entry whose energy has been
changed and is no longer consistent with the mixing state data.

This is a corner case that should never occur, because if the entry_id
of this entry is not already in the state_data, a CompatibilityError
should be raised. If the Entry was passed to get_mixing_state_data, then
its energy should already be there. Hence, this is testing a case
in which either 1) get_mixing_state_data fails to work properly or 2)
the energy of the Entry is somehow modified in between calling
get_mixing_state_data and process_entries. Such a situation could
potentially arise in e.g. the build pipeline if one is calling
process_entries with a separately-calculated state_data DataFrame.


#### test_state_gga_1_scan(mixing_scheme_no_compat, ms_gga_1_scan)
Mixing state in which we have a complete GGA PhaseDiagram and 1 R2SCAN entry
The R2SCAN entry chosen is the GGA ground state for SnBr2 (r2scan-4).

In this state, the mixing scheme should adjust the entry of r2scan-4 to
match the GGA energy and discard entry gga-4.


#### test_state_gga_1_scan_plus_novel(mixing_scheme_no_compat, ms_gga_1_scan_novel)
Mixing state in which we have a complete GGA PhaseDiagram and 1 R2SCAN entry
at a composition not in the GGA phase diagram.

In this state, the mixing scheme should discard the R2SCAN entry


#### test_state_gga_2_scan_diff_match(mixing_scheme_no_compat, ms_gga_2_scan_diff_match)
Mixing state in which we have a complete GGA PhaseDiagram and 2 R2SCAN entries
at different compositions, where both R2SCAN materials match GGA materials but
only one matches a GGA ground state.

In this state, the energies of both R2SCAN entries should be set equal to the
corresponding GGA energies, and the GGA entries discarded.


#### test_state_gga_2_scan_diff_nomatch(mixing_scheme_no_compat, ms_gga_2_scan_diff_no_match)
Mixing state in which we have a complete GGA PhaseDiagram and 2 R2SCAN entries
at different compositions, where one of the R2SCAN materials does not match
any GGA material.

In this state, the energy of the matching R2SCAN entry should be adjusted
to the GGA value, the corresponding GGA entry should be discarded, and the
novel R2SCAN material that doesn’t match anything should be discarded


#### test_state_gga_2_scan_same(mixing_scheme_no_compat, ms_gga_2_scan_same)
Mixing state in which we have a complete GGA PhaseDiagram and 2 R2SCAN entries
at a single composition, one of which is the GGA ground state.

In this state, the mixing scheme should correct the energy of unstable polymorph
r2scan-6 to maintain the same e_above_hull (r2scan-6 minus r2scan-4 which is the ground
state). Entry r2scan-4 (the GGA ground state) should be corrected to the GGA energy.
Entries gga-4 and gga-6 should be discarded.


#### test_state_gga_only(mixing_scheme_no_compat, ms_gga_only)
Mixing state in which we only have GGA entries, forming a complete PhaseDiagram.

In this state, the mixing scheme should not do anything


#### test_state_incomplete_gga_all_scan(mixing_scheme_no_compat, ms_incomplete_gga_all_scan)
Mixing state in which we have an incomplete GGA PhaseDiagram and all entries
present in R2SCAN.

This case should fail, because a complete run_type_1 phase diagram is required by the
mixing scheme


#### test_state_novel_scan_comp(mixing_scheme_no_compat, ms_all_gga_scan_gs_plus_novel)
Mixing state in which we have all GGA ground states in R2SCAN and then we try to
process a R2SCAN entry at a composition that is not in the GGA PhaseDiagram.

In this case, the mixing scheme should preserve the energy of the novel R2SCAN
entry and discard the 3 GGA ground states


#### test_state_scan_only(mixing_scheme_no_compat, ms_scan_only)
Mixing state in which we only have R2SCAN entries, forming a complete PhaseDiagram.

In this case, the mixing scheme should not do anything


### pymatgen.entries.tests.test_mixing_scheme.mixing_scheme_no_compat()
Return an instance of MaterialsProjectDFTMixingScheme with no additional
compatibility schemes (e.g., compat_1=None). Used by most of the tests where
we are manually supplying energies.


### pymatgen.entries.tests.test_mixing_scheme.ms_all_gga_scan_gs(ms_complete)
Mixing state with all GGA entries and R2SCAN entries corresponding to all GGA
ground states, but no others.


### pymatgen.entries.tests.test_mixing_scheme.ms_all_gga_scan_gs_plus_novel(ms_all_gga_scan_gs)
Mixing state with all GGA entries and R2SCAN entries corresponding to all GGA
ground states, plus one R2SCAN entry at a novel composition not in the GGA
phase diagram.


### pymatgen.entries.tests.test_mixing_scheme.ms_all_scan_novel(ms_complete)
Mixing state with all GGA entries and all R2SCAN, with an additional unstable
polymorphs of SnBr4 (r2scan-8) only in R2SCAN.


### pymatgen.entries.tests.test_mixing_scheme.ms_complete()
Mixing state where we have R2SCAN for all GGA.


### pymatgen.entries.tests.test_mixing_scheme.ms_complete_duplicate_structs(ms_complete)
Mixing state where we have R2SCAN for all GGA, plus extra entries that duplicate
the structures of gga-4 and r2scan-4, and have slightly higher energies.


### pymatgen.entries.tests.test_mixing_scheme.ms_gga_1_scan(ms_complete)
Mixing state with all GGA entries and one R2SCAN, corresponding to the GGA
ground state of SnBr2 (r2scan-4).


### pymatgen.entries.tests.test_mixing_scheme.ms_gga_1_scan_novel(ms_complete)
Mixing state with all GGA entries and 1 R2SCAN, corresponding to a composition
(SnBr) that is not present in the GGA entries.


### pymatgen.entries.tests.test_mixing_scheme.ms_gga_2_scan_diff_match(ms_complete)
Mixing state with all GGA entries and 2 R2SCAN entries corresponding to
different compositions, where both R2SCAN materials match GGA materials, but
only one matches a GGA ground state.
r2scan-4 and r2scan-7.


### pymatgen.entries.tests.test_mixing_scheme.ms_gga_2_scan_diff_no_match(ms_complete)
Mixing state with all GGA entries and 2 R2SCAN, corresponding to the GGA
ground state of SnBr2 (r2scan-4) and one unstable polymorph of SnBr4
that does not match any GGA material (r2scan-8).


### pymatgen.entries.tests.test_mixing_scheme.ms_gga_2_scan_same(ms_complete)
Mixing state with all GGA entries and 2 R2SCAN, corresponding to the GGA
ground state and one unstable polymorph of SnBr2 (r2scan-4 and r2scan-6).


### pymatgen.entries.tests.test_mixing_scheme.ms_gga_only(ms_complete)
Mixing state with only GGA entries.


### pymatgen.entries.tests.test_mixing_scheme.ms_incomplete_gga_all_scan(ms_complete)
Mixing state with an incomplete GGA phase diagram.


### pymatgen.entries.tests.test_mixing_scheme.ms_scan_chemsys_superset(ms_complete)
Mixing state where we have R2SCAN for all GGA, and there is an additional R2SCAN
entry outside the GGA chemsys.


### pymatgen.entries.tests.test_mixing_scheme.ms_scan_only(ms_complete)
Mixing state with only R2SCAN entries.


### pymatgen.entries.tests.test_mixing_scheme.test_data_ms_complete(ms_complete)
Verify that the test chemical system
ComputedStructureEntry match (or don’t match) as intended.