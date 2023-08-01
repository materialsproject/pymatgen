---
layout: default
title: pymatgen.entries.md
nav_exclude: true
---

# pymatgen.entries package

Entries are containers for calculated information, which is used in
many analyses. This module contains entry related tools and implements
the base Entry class, which is the basic entity that can be used to
store calculated information. Other Entry classes such as ComputedEntry
and PDEntry inherit from this class.


### _class_ pymatgen.entries.Entry(composition: [Composition](pymatgen.core.md#pymatgen.core.composition.Composition) | str | dict[str, float], energy: float)
Bases: `MSONable`

A lightweight object containing the energy associated with
a specific chemical composition. This base class is not
intended to be instantiated directly. Note that classes
which inherit from Entry must define a .energy property.

Initializes an Entry.


* **Parameters**


    * **composition** ([*Composition*](pymatgen.core.md#pymatgen.core.composition.Composition)) – Composition of the entry. For
    flexibility, this can take the form of all the typical input
    taken by a Composition, including a {symbol: amt} dict,
    a string formula, and others.


    * **energy** (*float*) – Energy of the entry.



#### as_dict()
MSONable dict.


#### _property_ composition(_: [Composition](pymatgen.core.md#pymatgen.core.composition.Composition_ )
the composition of the entry.


* **Type**

    return



#### _abstract property_ energy(_: floa_ )
the energy of the entry.


* **Type**

    return



#### _property_ energy_per_atom(_: floa_ )
the energy per atom of the entry.


* **Type**

    return



#### _property_ is_element(_: boo_ )
Whether composition of entry is an element.


* **Type**

    return



#### normalize(mode: Literal['formula_unit', 'atom'] = 'formula_unit')
Normalize the entry’s composition and energy.


* **Parameters**

    **mode** (*"formula_unit"** | **"atom"*) – “formula_unit” (the default) normalizes to composition.reduced_formula.
    “atom” normalizes such that the composition amounts sum to 1.


## Subpackages


* [pymatgen.entries.tests package](pymatgen.entries.tests.md)




    * [pymatgen.entries.tests.test_compatibility module](pymatgen.entries.tests.md#module-pymatgen.entries.tests.test_compatibility)


        * [`AqueousCorrectionTest`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.AqueousCorrectionTest)


            * [`AqueousCorrectionTest.setUp()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.AqueousCorrectionTest.setUp)


            * [`AqueousCorrectionTest.test_compound_energy()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.AqueousCorrectionTest.test_compound_energy)


        * [`CorrectionErrors2020CompatibilityTest`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.CorrectionErrors2020CompatibilityTest)


            * [`CorrectionErrors2020CompatibilityTest.setUp()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.CorrectionErrors2020CompatibilityTest.setUp)


            * [`CorrectionErrors2020CompatibilityTest.tearDown()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.CorrectionErrors2020CompatibilityTest.tearDown)


            * [`CorrectionErrors2020CompatibilityTest.test_errors()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.CorrectionErrors2020CompatibilityTest.test_errors)


        * [`CorrectionSpecificityTest`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.CorrectionSpecificityTest)


            * [`CorrectionSpecificityTest.setUp()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.CorrectionSpecificityTest.setUp)


            * [`CorrectionSpecificityTest.test_correction_specificity()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.CorrectionSpecificityTest.test_correction_specificity)


        * [`DummyCompatibility`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.DummyCompatibility)


            * [`DummyCompatibility.get_adjustments()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.DummyCompatibility.get_adjustments)


        * [`MITAqueousCompatibilityTest`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MITAqueousCompatibilityTest)


            * [`MITAqueousCompatibilityTest.setUp()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MITAqueousCompatibilityTest.setUp)


            * [`MITAqueousCompatibilityTest.test_aqueous_compat()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MITAqueousCompatibilityTest.test_aqueous_compat)


            * [`MITAqueousCompatibilityTest.test_dont_error_on_weird_elements()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MITAqueousCompatibilityTest.test_dont_error_on_weird_elements)


            * [`MITAqueousCompatibilityTest.test_msonable()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MITAqueousCompatibilityTest.test_msonable)


            * [`MITAqueousCompatibilityTest.test_potcar_doenst_match_structure()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MITAqueousCompatibilityTest.test_potcar_doenst_match_structure)


        * [`MITCompatibilityTest`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MITCompatibilityTest)


            * [`MITCompatibilityTest.setUp()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MITCompatibilityTest.setUp)


            * [`MITCompatibilityTest.tearDown()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MITCompatibilityTest.tearDown)


            * [`MITCompatibilityTest.test_U_value()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MITCompatibilityTest.test_U_value)


            * [`MITCompatibilityTest.test_correction_value()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MITCompatibilityTest.test_correction_value)


            * [`MITCompatibilityTest.test_element_processing()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MITCompatibilityTest.test_element_processing)


            * [`MITCompatibilityTest.test_get_explanation_dict()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MITCompatibilityTest.test_get_explanation_dict)


            * [`MITCompatibilityTest.test_msonable()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MITCompatibilityTest.test_msonable)


            * [`MITCompatibilityTest.test_potcar_doenst_match_structure()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MITCompatibilityTest.test_potcar_doenst_match_structure)


            * [`MITCompatibilityTest.test_potcar_spec_is_none()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MITCompatibilityTest.test_potcar_spec_is_none)


            * [`MITCompatibilityTest.test_process_entry()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MITCompatibilityTest.test_process_entry)


            * [`MITCompatibilityTest.test_revert_to_symbols()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MITCompatibilityTest.test_revert_to_symbols)


            * [`MITCompatibilityTest.test_same_potcar_symbol()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MITCompatibilityTest.test_same_potcar_symbol)


            * [`MITCompatibilityTest.test_wrong_U_value()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MITCompatibilityTest.test_wrong_U_value)


            * [`MITCompatibilityTest.test_wrong_psp()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MITCompatibilityTest.test_wrong_psp)


        * [`MaterialsProjectCompatibility2020Test`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MaterialsProjectCompatibility2020Test)


            * [`MaterialsProjectCompatibility2020Test.setUp()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MaterialsProjectCompatibility2020Test.setUp)


            * [`MaterialsProjectCompatibility2020Test.tearDown()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MaterialsProjectCompatibility2020Test.tearDown)


            * [`MaterialsProjectCompatibility2020Test.test_U_values()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MaterialsProjectCompatibility2020Test.test_U_values)


            * [`MaterialsProjectCompatibility2020Test.test_check_potcar()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MaterialsProjectCompatibility2020Test.test_check_potcar)


            * [`MaterialsProjectCompatibility2020Test.test_config_file()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MaterialsProjectCompatibility2020Test.test_config_file)


            * [`MaterialsProjectCompatibility2020Test.test_correction_values()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MaterialsProjectCompatibility2020Test.test_correction_values)


            * [`MaterialsProjectCompatibility2020Test.test_element_processing()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MaterialsProjectCompatibility2020Test.test_element_processing)


            * [`MaterialsProjectCompatibility2020Test.test_energy_adjustments()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MaterialsProjectCompatibility2020Test.test_energy_adjustments)


            * [`MaterialsProjectCompatibility2020Test.test_get_explanation_dict()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MaterialsProjectCompatibility2020Test.test_get_explanation_dict)


            * [`MaterialsProjectCompatibility2020Test.test_msonable()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MaterialsProjectCompatibility2020Test.test_msonable)


            * [`MaterialsProjectCompatibility2020Test.test_oxdiation()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MaterialsProjectCompatibility2020Test.test_oxdiation)


            * [`MaterialsProjectCompatibility2020Test.test_oxdiation_by_electronegativity()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MaterialsProjectCompatibility2020Test.test_oxdiation_by_electronegativity)


            * [`MaterialsProjectCompatibility2020Test.test_oxi_state_guess()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MaterialsProjectCompatibility2020Test.test_oxi_state_guess)


            * [`MaterialsProjectCompatibility2020Test.test_process_entries()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MaterialsProjectCompatibility2020Test.test_process_entries)


            * [`MaterialsProjectCompatibility2020Test.test_process_entry()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MaterialsProjectCompatibility2020Test.test_process_entry)


            * [`MaterialsProjectCompatibility2020Test.test_process_entry_with_oxidation_state()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MaterialsProjectCompatibility2020Test.test_process_entry_with_oxidation_state)


            * [`MaterialsProjectCompatibility2020Test.test_processing_entries_inplace()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MaterialsProjectCompatibility2020Test.test_processing_entries_inplace)


            * [`MaterialsProjectCompatibility2020Test.test_wrong_psp()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MaterialsProjectCompatibility2020Test.test_wrong_psp)


        * [`MaterialsProjectCompatibilityTest`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MaterialsProjectCompatibilityTest)


            * [`MaterialsProjectCompatibilityTest.setUp()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MaterialsProjectCompatibilityTest.setUp)


            * [`MaterialsProjectCompatibilityTest.tearDown()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MaterialsProjectCompatibilityTest.tearDown)


            * [`MaterialsProjectCompatibilityTest.test_U_values()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MaterialsProjectCompatibilityTest.test_U_values)


            * [`MaterialsProjectCompatibilityTest.test_correction_values()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MaterialsProjectCompatibilityTest.test_correction_values)


            * [`MaterialsProjectCompatibilityTest.test_element_processing()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MaterialsProjectCompatibilityTest.test_element_processing)


            * [`MaterialsProjectCompatibilityTest.test_get_corrections_dict()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MaterialsProjectCompatibilityTest.test_get_corrections_dict)


            * [`MaterialsProjectCompatibilityTest.test_get_explanation_dict()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MaterialsProjectCompatibilityTest.test_get_explanation_dict)


            * [`MaterialsProjectCompatibilityTest.test_msonable()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MaterialsProjectCompatibilityTest.test_msonable)


            * [`MaterialsProjectCompatibilityTest.test_process_entries()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MaterialsProjectCompatibilityTest.test_process_entries)


            * [`MaterialsProjectCompatibilityTest.test_process_entry()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MaterialsProjectCompatibilityTest.test_process_entry)


            * [`MaterialsProjectCompatibilityTest.test_wrong_psp()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.MaterialsProjectCompatibilityTest.test_wrong_psp)


        * [`OxideTypeCorrectionNoPeroxideCorrTest`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.OxideTypeCorrectionNoPeroxideCorrTest)


            * [`OxideTypeCorrectionNoPeroxideCorrTest.setUp()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.OxideTypeCorrectionNoPeroxideCorrTest.setUp)


            * [`OxideTypeCorrectionNoPeroxideCorrTest.test_oxide_energy_corr()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.OxideTypeCorrectionNoPeroxideCorrTest.test_oxide_energy_corr)


            * [`OxideTypeCorrectionNoPeroxideCorrTest.test_ozonide()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.OxideTypeCorrectionNoPeroxideCorrTest.test_ozonide)


            * [`OxideTypeCorrectionNoPeroxideCorrTest.test_peroxide_energy_corr()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.OxideTypeCorrectionNoPeroxideCorrTest.test_peroxide_energy_corr)


        * [`OxideTypeCorrectionTest`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.OxideTypeCorrectionTest)


            * [`OxideTypeCorrectionTest.setUp()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.OxideTypeCorrectionTest.setUp)


            * [`OxideTypeCorrectionTest.test_no_struct_compat()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.OxideTypeCorrectionTest.test_no_struct_compat)


            * [`OxideTypeCorrectionTest.test_process_entry_oxide()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.OxideTypeCorrectionTest.test_process_entry_oxide)


            * [`OxideTypeCorrectionTest.test_process_entry_ozonide()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.OxideTypeCorrectionTest.test_process_entry_ozonide)


            * [`OxideTypeCorrectionTest.test_process_entry_peroxide()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.OxideTypeCorrectionTest.test_process_entry_peroxide)


            * [`OxideTypeCorrectionTest.test_process_entry_superoxide()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.OxideTypeCorrectionTest.test_process_entry_superoxide)


        * [`SulfideTypeCorrection2020Test`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.SulfideTypeCorrection2020Test)


            * [`SulfideTypeCorrection2020Test.setUp()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.SulfideTypeCorrection2020Test.setUp)


            * [`SulfideTypeCorrection2020Test.test_struct_no_struct()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.SulfideTypeCorrection2020Test.test_struct_no_struct)


        * [`TestMaterialsProjectAqueousCompatibility`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.TestMaterialsProjectAqueousCompatibility)


            * [`TestMaterialsProjectAqueousCompatibility.test_compound_entropy()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.TestMaterialsProjectAqueousCompatibility.test_compound_entropy)


            * [`TestMaterialsProjectAqueousCompatibility.test_h_h2o_energy_no_args()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.TestMaterialsProjectAqueousCompatibility.test_h_h2o_energy_no_args)


            * [`TestMaterialsProjectAqueousCompatibility.test_h_h2o_energy_with_args_multi()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.TestMaterialsProjectAqueousCompatibility.test_h_h2o_energy_with_args_multi)


            * [`TestMaterialsProjectAqueousCompatibility.test_h_h2o_energy_with_args_single()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.TestMaterialsProjectAqueousCompatibility.test_h_h2o_energy_with_args_single)


            * [`TestMaterialsProjectAqueousCompatibility.test_hydrate_adjustment()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.TestMaterialsProjectAqueousCompatibility.test_hydrate_adjustment)


            * [`TestMaterialsProjectAqueousCompatibility.test_processing_entries_inplace()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.TestMaterialsProjectAqueousCompatibility.test_processing_entries_inplace)


        * [`test_clean_arg()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.test_clean_arg)


        * [`test_energy_adjustment_normalize()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.test_energy_adjustment_normalize)


        * [`test_no_duplicate_corrections()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.test_no_duplicate_corrections)


        * [`test_overlapping_adjustments()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.test_overlapping_adjustments)


        * [`test_process_entries_return_type()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_compatibility.test_process_entries_return_type)


    * [pymatgen.entries.tests.test_computed_entries module](pymatgen.entries.tests.md#module-pymatgen.entries.tests.test_computed_entries)


        * [`ComputedEntryTest`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.ComputedEntryTest)


            * [`ComputedEntryTest.setUp()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.ComputedEntryTest.setUp)


            * [`ComputedEntryTest.test_composition()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.ComputedEntryTest.test_composition)


            * [`ComputedEntryTest.test_conflicting_correction_adjustment()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.ComputedEntryTest.test_conflicting_correction_adjustment)


            * [`ComputedEntryTest.test_copy()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.ComputedEntryTest.test_copy)


            * [`ComputedEntryTest.test_energy()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.ComputedEntryTest.test_energy)


            * [`ComputedEntryTest.test_entry_id()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.ComputedEntryTest.test_entry_id)


            * [`ComputedEntryTest.test_is_element()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.ComputedEntryTest.test_is_element)


            * [`ComputedEntryTest.test_normalize()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.ComputedEntryTest.test_normalize)


            * [`ComputedEntryTest.test_normalize_energy_adjustments()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.ComputedEntryTest.test_normalize_energy_adjustments)


            * [`ComputedEntryTest.test_per_atom_props()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.ComputedEntryTest.test_per_atom_props)


            * [`ComputedEntryTest.test_str()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.ComputedEntryTest.test_str)


            * [`ComputedEntryTest.test_sulfide_energy()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.ComputedEntryTest.test_sulfide_energy)


            * [`ComputedEntryTest.test_to_from_dict()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.ComputedEntryTest.test_to_from_dict)


            * [`ComputedEntryTest.test_to_from_dict_with_adjustment()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.ComputedEntryTest.test_to_from_dict_with_adjustment)


            * [`ComputedEntryTest.test_to_from_dict_with_adjustment_2()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.ComputedEntryTest.test_to_from_dict_with_adjustment_2)


            * [`ComputedEntryTest.test_to_from_dict_with_adjustment_3()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.ComputedEntryTest.test_to_from_dict_with_adjustment_3)


        * [`ComputedStructureEntryTest`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.ComputedStructureEntryTest)


            * [`ComputedStructureEntryTest.setUp()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.ComputedStructureEntryTest.setUp)


            * [`ComputedStructureEntryTest.test_composition()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.ComputedStructureEntryTest.test_composition)


            * [`ComputedStructureEntryTest.test_copy()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.ComputedStructureEntryTest.test_copy)


            * [`ComputedStructureEntryTest.test_energy()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.ComputedStructureEntryTest.test_energy)


            * [`ComputedStructureEntryTest.test_eq()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.ComputedStructureEntryTest.test_eq)


            * [`ComputedStructureEntryTest.test_str()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.ComputedStructureEntryTest.test_str)


            * [`ComputedStructureEntryTest.test_to_from_dict()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.ComputedStructureEntryTest.test_to_from_dict)


            * [`ComputedStructureEntryTest.test_to_from_dict_structure_with_adjustment_3()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.ComputedStructureEntryTest.test_to_from_dict_structure_with_adjustment_3)


        * [`GibbsComputedStructureEntryTest`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.GibbsComputedStructureEntryTest)


            * [`GibbsComputedStructureEntryTest.setUp()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.GibbsComputedStructureEntryTest.setUp)


            * [`GibbsComputedStructureEntryTest.test_expt_gas_entry()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.GibbsComputedStructureEntryTest.test_expt_gas_entry)


            * [`GibbsComputedStructureEntryTest.test_from_entries()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.GibbsComputedStructureEntryTest.test_from_entries)


            * [`GibbsComputedStructureEntryTest.test_from_pd()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.GibbsComputedStructureEntryTest.test_from_pd)


            * [`GibbsComputedStructureEntryTest.test_gf_sisso()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.GibbsComputedStructureEntryTest.test_gf_sisso)


            * [`GibbsComputedStructureEntryTest.test_interpolation()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.GibbsComputedStructureEntryTest.test_interpolation)


            * [`GibbsComputedStructureEntryTest.test_normalize()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.GibbsComputedStructureEntryTest.test_normalize)


            * [`GibbsComputedStructureEntryTest.test_str()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.GibbsComputedStructureEntryTest.test_str)


            * [`GibbsComputedStructureEntryTest.test_to_from_dict()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.GibbsComputedStructureEntryTest.test_to_from_dict)


        * [`test_composition_energy_adjustment()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.test_composition_energy_adjustment)


        * [`test_constant_energy_adjustment()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.test_constant_energy_adjustment)


        * [`test_energy_adjustment()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.test_energy_adjustment)


        * [`test_manual_energy_adjustment()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.test_manual_energy_adjustment)


        * [`test_temp_energy_adjustment()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_computed_entries.test_temp_energy_adjustment)


    * [pymatgen.entries.tests.test_correction_calculator module](pymatgen.entries.tests.md#module-pymatgen.entries.tests.test_correction_calculator)


        * [`CorrectionCalculatorTest`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_correction_calculator.CorrectionCalculatorTest)


            * [`CorrectionCalculatorTest.setUp()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_correction_calculator.CorrectionCalculatorTest.setUp)


            * [`CorrectionCalculatorTest.test_missing_entry_response()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_correction_calculator.CorrectionCalculatorTest.test_missing_entry_response)


            * [`CorrectionCalculatorTest.test_no_uncertainties()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_correction_calculator.CorrectionCalculatorTest.test_no_uncertainties)


            * [`CorrectionCalculatorTest.test_normal_corrections()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_correction_calculator.CorrectionCalculatorTest.test_normal_corrections)


            * [`CorrectionCalculatorTest.test_warnings_options()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_correction_calculator.CorrectionCalculatorTest.test_warnings_options)


    * [pymatgen.entries.tests.test_entry_tools module](pymatgen.entries.tests.md#module-pymatgen.entries.tests.test_entry_tools)


        * [`EntrySetTest`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_entry_tools.EntrySetTest)


            * [`EntrySetTest.setUp()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_entry_tools.EntrySetTest.setUp)


            * [`EntrySetTest.test_as_dict()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_entry_tools.EntrySetTest.test_as_dict)


            * [`EntrySetTest.test_chemsys()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_entry_tools.EntrySetTest.test_chemsys)


            * [`EntrySetTest.test_get_subset()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_entry_tools.EntrySetTest.test_get_subset)


            * [`EntrySetTest.test_remove_non_ground_states()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_entry_tools.EntrySetTest.test_remove_non_ground_states)


        * [`FuncTest`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_entry_tools.FuncTest)


            * [`FuncTest.test_group_entries_by_composition()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_entry_tools.FuncTest.test_group_entries_by_composition)


            * [`FuncTest.test_group_entries_by_structure()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_entry_tools.FuncTest.test_group_entries_by_structure)


    * [pymatgen.entries.tests.test_exp_entries module](pymatgen.entries.tests.md#module-pymatgen.entries.tests.test_exp_entries)


        * [`ExpEntryTest`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_exp_entries.ExpEntryTest)


            * [`ExpEntryTest.setUp()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_exp_entries.ExpEntryTest.setUp)


            * [`ExpEntryTest.test_energy()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_exp_entries.ExpEntryTest.test_energy)


            * [`ExpEntryTest.test_str()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_exp_entries.ExpEntryTest.test_str)


            * [`ExpEntryTest.test_to_from_dict()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_exp_entries.ExpEntryTest.test_to_from_dict)


    * [pymatgen.entries.tests.test_mixing_scheme module](pymatgen.entries.tests.md#module-pymatgen.entries.tests.test_mixing_scheme)


        * [Test Structure](pymatgen.entries.tests.md#test-structure)


        * [Description of the ms_complete mixing state](pymatgen.entries.tests.md#description-of-the-ms-complete-mixing-state)


        * [Types of Tests](pymatgen.entries.tests.md#types-of-tests)


        * [Implementation Notes](pymatgen.entries.tests.md#implementation-notes)


        * [`DummyCompatibility`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.DummyCompatibility)


            * [`DummyCompatibility.get_adjustments()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.DummyCompatibility.get_adjustments)


        * [`MixingState`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.MixingState)


            * [`MixingState.all_entries`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.MixingState.all_entries)


        * [`TestMaterialsProjectDFTMixingSchemeArgs`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.TestMaterialsProjectDFTMixingSchemeArgs)


            * [`TestMaterialsProjectDFTMixingSchemeArgs.test_alternate_run_types()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.TestMaterialsProjectDFTMixingSchemeArgs.test_alternate_run_types)


            * [`TestMaterialsProjectDFTMixingSchemeArgs.test_alternate_structure_matcher()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.TestMaterialsProjectDFTMixingSchemeArgs.test_alternate_structure_matcher)


            * [`TestMaterialsProjectDFTMixingSchemeArgs.test_check_potcar()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.TestMaterialsProjectDFTMixingSchemeArgs.test_check_potcar)


            * [`TestMaterialsProjectDFTMixingSchemeArgs.test_clean()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.TestMaterialsProjectDFTMixingSchemeArgs.test_clean)


            * [`TestMaterialsProjectDFTMixingSchemeArgs.test_compat_args()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.TestMaterialsProjectDFTMixingSchemeArgs.test_compat_args)


            * [`TestMaterialsProjectDFTMixingSchemeArgs.test_empty_entries()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.TestMaterialsProjectDFTMixingSchemeArgs.test_empty_entries)


            * [`TestMaterialsProjectDFTMixingSchemeArgs.test_fuzzy_matching()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.TestMaterialsProjectDFTMixingSchemeArgs.test_fuzzy_matching)


            * [`TestMaterialsProjectDFTMixingSchemeArgs.test_incompatible_run_type()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.TestMaterialsProjectDFTMixingSchemeArgs.test_incompatible_run_type)


            * [`TestMaterialsProjectDFTMixingSchemeArgs.test_multiple_matching_structures()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.TestMaterialsProjectDFTMixingSchemeArgs.test_multiple_matching_structures)


            * [`TestMaterialsProjectDFTMixingSchemeArgs.test_no_entry_ids()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.TestMaterialsProjectDFTMixingSchemeArgs.test_no_entry_ids)


            * [`TestMaterialsProjectDFTMixingSchemeArgs.test_no_foreign_entries()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.TestMaterialsProjectDFTMixingSchemeArgs.test_no_foreign_entries)


            * [`TestMaterialsProjectDFTMixingSchemeArgs.test_no_mixing_data()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.TestMaterialsProjectDFTMixingSchemeArgs.test_no_mixing_data)


            * [`TestMaterialsProjectDFTMixingSchemeArgs.test_no_run_type()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.TestMaterialsProjectDFTMixingSchemeArgs.test_no_run_type)


            * [`TestMaterialsProjectDFTMixingSchemeArgs.test_no_single_entry()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.TestMaterialsProjectDFTMixingSchemeArgs.test_no_single_entry)


            * [`TestMaterialsProjectDFTMixingSchemeArgs.test_no_structure()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.TestMaterialsProjectDFTMixingSchemeArgs.test_no_structure)


            * [`TestMaterialsProjectDFTMixingSchemeArgs.test_processing_entries_inplace()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.TestMaterialsProjectDFTMixingSchemeArgs.test_processing_entries_inplace)


            * [`TestMaterialsProjectDFTMixingSchemeArgs.test_same_run_type()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.TestMaterialsProjectDFTMixingSchemeArgs.test_same_run_type)


        * [`TestMaterialsProjectDFTMixingSchemeStates`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.TestMaterialsProjectDFTMixingSchemeStates)


            * [`TestMaterialsProjectDFTMixingSchemeStates.test_chemsys_mismatch()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.TestMaterialsProjectDFTMixingSchemeStates.test_chemsys_mismatch)


            * [`TestMaterialsProjectDFTMixingSchemeStates.test_state_all_gga_scan_gs()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.TestMaterialsProjectDFTMixingSchemeStates.test_state_all_gga_scan_gs)


            * [`TestMaterialsProjectDFTMixingSchemeStates.test_state_complete_entries()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.TestMaterialsProjectDFTMixingSchemeStates.test_state_complete_entries)


            * [`TestMaterialsProjectDFTMixingSchemeStates.test_state_energy_modified()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.TestMaterialsProjectDFTMixingSchemeStates.test_state_energy_modified)


            * [`TestMaterialsProjectDFTMixingSchemeStates.test_state_gga_1_scan()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.TestMaterialsProjectDFTMixingSchemeStates.test_state_gga_1_scan)


            * [`TestMaterialsProjectDFTMixingSchemeStates.test_state_gga_1_scan_plus_novel()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.TestMaterialsProjectDFTMixingSchemeStates.test_state_gga_1_scan_plus_novel)


            * [`TestMaterialsProjectDFTMixingSchemeStates.test_state_gga_2_scan_diff_match()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.TestMaterialsProjectDFTMixingSchemeStates.test_state_gga_2_scan_diff_match)


            * [`TestMaterialsProjectDFTMixingSchemeStates.test_state_gga_2_scan_diff_nomatch()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.TestMaterialsProjectDFTMixingSchemeStates.test_state_gga_2_scan_diff_nomatch)


            * [`TestMaterialsProjectDFTMixingSchemeStates.test_state_gga_2_scan_same()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.TestMaterialsProjectDFTMixingSchemeStates.test_state_gga_2_scan_same)


            * [`TestMaterialsProjectDFTMixingSchemeStates.test_state_gga_only()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.TestMaterialsProjectDFTMixingSchemeStates.test_state_gga_only)


            * [`TestMaterialsProjectDFTMixingSchemeStates.test_state_incomplete_gga_all_scan()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.TestMaterialsProjectDFTMixingSchemeStates.test_state_incomplete_gga_all_scan)


            * [`TestMaterialsProjectDFTMixingSchemeStates.test_state_novel_scan_comp()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.TestMaterialsProjectDFTMixingSchemeStates.test_state_novel_scan_comp)


            * [`TestMaterialsProjectDFTMixingSchemeStates.test_state_scan_only()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.TestMaterialsProjectDFTMixingSchemeStates.test_state_scan_only)


        * [`mixing_scheme_no_compat()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.mixing_scheme_no_compat)


        * [`ms_all_gga_scan_gs()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.ms_all_gga_scan_gs)


        * [`ms_all_gga_scan_gs_plus_novel()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.ms_all_gga_scan_gs_plus_novel)


        * [`ms_all_scan_novel()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.ms_all_scan_novel)


        * [`ms_complete()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.ms_complete)


        * [`ms_complete_duplicate_structs()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.ms_complete_duplicate_structs)


        * [`ms_gga_1_scan()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.ms_gga_1_scan)


        * [`ms_gga_1_scan_novel()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.ms_gga_1_scan_novel)


        * [`ms_gga_2_scan_diff_match()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.ms_gga_2_scan_diff_match)


        * [`ms_gga_2_scan_diff_no_match()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.ms_gga_2_scan_diff_no_match)


        * [`ms_gga_2_scan_same()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.ms_gga_2_scan_same)


        * [`ms_gga_only()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.ms_gga_only)


        * [`ms_incomplete_gga_all_scan()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.ms_incomplete_gga_all_scan)


        * [`ms_scan_chemsys_superset()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.ms_scan_chemsys_superset)


        * [`ms_scan_only()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.ms_scan_only)


        * [`test_data_ms_complete()`](pymatgen.entries.tests.md#pymatgen.entries.tests.test_mixing_scheme.test_data_ms_complete)



## pymatgen.entries.compatibility module

This module implements Compatibility corrections for mixing runs of different
functionals.


### _class_ pymatgen.entries.compatibility.AnionCorrection(\*args, \*\*kwargs)
Bases: `AnionCorrection`

Correct anion energies to obtain the right formation energies. Note that
this depends on calculations being run within the same input set.

Used by legacy MaterialsProjectCompatibility and MITCompatibility.


* **Parameters**


    * **config_file** – Path to the selected compatibility.yaml config file.


    * **correct_peroxide** – Specify whether peroxide/superoxide/ozonide
    corrections are to be applied or not.



### _class_ pymatgen.entries.compatibility.AqueousCorrection(\*args, \*\*kwargs)
Bases: `AqueousCorrection`

This class implements aqueous phase compound corrections for elements
and H2O.

Used only by MITAqueousCompatibility.


* **Parameters**


    * **config_file** – Path to the selected compatibility.yaml config file.


    * **error_file** – Path to the selected compatibilityErrors.yaml config file.



### _class_ pymatgen.entries.compatibility.Compatibility()
Bases: `MSONable`

Abstract Compatibility class, not intended for direct use.
Compatibility classes are used to correct the energies of an entry or a set
of entries. All Compatibility classes must implement get_adjustments() method.


#### _static_ explain(entry)
Prints an explanation of the energy adjustments applied by the
Compatibility class. Inspired by the “explain” methods in many database
methodologies.


* **Parameters**

    **entry** – A ComputedEntry.



#### _abstract_ get_adjustments(entry: ComputedEntry | ComputedStructureEntry)
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

    list[EnergyAdjustment]



* **Raises**

    **CompatibilityError if the entry is not compatible** –



#### process_entries(entries: AnyComputedEntry | list[AnyComputedEntry], clean: bool = True, verbose: bool = False, inplace: bool = True, on_error: Literal['ignore', 'warn', 'raise'] = 'ignore')
Process a sequence of entries with the chosen Compatibility scheme.

Warning: This method changes entries in place! All changes can be undone and original entries
restored by setting entry.energy_adjustments = [].


* **Parameters**


    * **entries** (*AnyComputedEntry** | **list**[**AnyComputedEntry**]*) – A sequence of
    Computed(Structure)Entry objects.


    * **clean** (*bool*) – Whether to remove any previously-applied energy adjustments.
    If True, all EnergyAdjustment are removed prior to processing the Entry.
    Defaults to True.


    * **verbose** (*bool*) – Whether to display progress bar for processing multiple entries.
    Defaults to False.


    * **inplace** (*bool*) – Whether to adjust input entries in place. Defaults to True.


    * **on_error** (*'ignore'** | **'warn'** | **'raise'*) – What to do when get_adjustments(entry)
    raises CompatibilityError. Defaults to ‘ignore’.



* **Returns**

    Adjusted entries. Entries in the original list incompatible with

        chosen correction scheme are excluded from the returned list.




* **Return type**

    list[AnyComputedEntry]



#### process_entry(entry: ComputedEntry, \*\*kwargs)
Process a single entry with the chosen Corrections. Note
that this method will change the data of the original entry.


* **Parameters**


    * **entry** – A ComputedEntry object.


    * **\*\*kwargs** – Will be passed to process_entries().



* **Returns**

    An adjusted entry if entry is compatible, else None.



### _exception_ pymatgen.entries.compatibility.CompatibilityError()
Bases: `Exception`

Exception class for Compatibility. Raised by attempting correction
on incompatible calculation.


### _class_ pymatgen.entries.compatibility.Correction()
Bases: `object`

A Correction class is a pre-defined scheme for correction a computed
entry based on the type and chemistry of the structure and the
calculation parameters. All Correction classes must implement a
correct_entry method.


#### correct_entry(entry)
Corrects a single entry.


* **Parameters**

    **entry** – A ComputedEntry object.



* **Returns**

    An processed entry.



* **Raises**

    **CompatibilityError if entry is not compatible.** –



#### _abstract_ get_correction(entry: ComputedEntry | ComputedStructureEntry)
Returns correction and uncertainty for a single entry.


* **Parameters**

    **entry** – A ComputedEntry object.



* **Returns**

    The energy correction to be applied and the uncertainty of the correction.



* **Raises**

    **CompatibilityError if entry is not compatible.** –



### _class_ pymatgen.entries.compatibility.CorrectionsList(corrections: Sequence[Correction])
Bases: `Compatibility`

The CorrectionsList class combines a list of corrections to be applied to
an entry or a set of entries. Note that some of the Corrections have
interdependencies. For example, PotcarCorrection must always be used
before any other compatibility. Also, AnionCorrection(“MP”) must be used
with PotcarCorrection(“MP”) (similarly with “MIT”). Typically,
you should use the specific MaterialsProjectCompatibility and
MITCompatibility subclasses instead.


* **Parameters**

    **corrections** (*list**[**Correction**]*) – Correction objects to apply.



#### explain(entry)
Prints an explanation of the corrections that are being applied for a
given compatibility scheme. Inspired by the “explain” methods in many
database methodologies.


* **Parameters**

    **entry** – A ComputedEntry.



#### get_adjustments(entry: ComputedEntry | ComputedStructureEntry)
Get the list of energy adjustments to be applied to an entry.


#### get_corrections_dict(entry: ComputedEntry | ComputedStructureEntry)
Returns the correction values and uncertainties applied to a particular entry.


* **Parameters**

    **entry** – A ComputedEntry object.



* **Returns**

    Map from correction names to values

        (1st) and uncertainties (2nd).




* **Return type**

    tuple[dict[str, float], dict[str, float]]



#### get_explanation_dict(entry)
Provides an explanation dict of the corrections that are being applied
for a given compatibility scheme. Inspired by the “explain” methods
in many database methodologies.


* **Parameters**

    **entry** – A ComputedEntry.



* **Returns**

    (dict) of the form
    {“Compatibility”: “string”,
    “Uncorrected_energy”: float,
    “Corrected_energy”: float,
    “correction_uncertainty:” float,
    “Corrections”: [{“Name of Correction”: {
    “Value”: float, “Explanation”: “string”, “Uncertainty”: float}]}



### _class_ pymatgen.entries.compatibility.GasCorrection(\*args, \*\*kwargs)
Bases: `GasCorrection`

Correct gas energies to obtain the right formation energies. Note that
this depends on calculations being run within the same input set.
Used by legacy MaterialsProjectCompatibility and MITCompatibility.


* **Parameters**

    **config_file** – Path to the selected compatibility.yaml config file.



### _class_ pymatgen.entries.compatibility.MITAqueousCompatibility(compat_type: Literal['GGA', 'Advanced'] = 'Advanced', correct_peroxide: bool = True, check_potcar_hash: bool = False)
Bases: `CorrectionsList`

This class implements the GGA/GGA+U mixing scheme, which allows mixing of
entries. Note that this should only be used for VASP calculations using the
MIT parameters (see pymatgen.io.vasp.sets MITVaspInputSet). Using
this compatibility scheme on runs with different parameters is not valid.


* **Parameters**


    * **compat_type** – Two options, GGA or Advanced. GGA means all GGA+U
    entries are excluded. Advanced means mixing scheme is
    implemented to make entries compatible with each other,
    but entries which are supposed to be done in GGA+U will have the
    equivalent GGA entries excluded. For example, Fe oxides should
    have a U value under the Advanced scheme. A GGA Fe oxide run
    will therefore be excluded under the scheme.


    * **correct_peroxide** – Specify whether peroxide/superoxide/ozonide
    corrections are to be applied or not.


    * **check_potcar_hash** (*bool*) – Use potcar hash to verify potcars are correct.



### _class_ pymatgen.entries.compatibility.MITCompatibility(compat_type: Literal['GGA', 'Advanced'] = 'Advanced', correct_peroxide: bool = True, check_potcar_hash: bool = False)
Bases: `CorrectionsList`

This class implements the GGA/GGA+U mixing scheme, which allows mixing of
entries. Note that this should only be used for VASP calculations using the
MIT parameters (see pymatgen.io.vasp.sets MITVaspInputSet). Using
this compatibility scheme on runs with different parameters is not valid.


* **Parameters**


    * **compat_type** – Two options, GGA or Advanced. GGA means all GGA+U
    entries are excluded. Advanced means mixing scheme is
    implemented to make entries compatible with each other,
    but entries which are supposed to be done in GGA+U will have the
    equivalent GGA entries excluded. For example, Fe oxides should
    have a U value under the Advanced scheme. A GGA Fe oxide run
    will therefore be excluded under the scheme.


    * **correct_peroxide** – Specify whether peroxide/superoxide/ozonide
    corrections are to be applied or not.


    * **check_potcar_hash** (*bool*) – Use potcar hash to verify potcars are correct.



### _class_ pymatgen.entries.compatibility.MaterialsProject2020Compatibility(\*args, \*\*kwargs)
Bases: `MaterialsProject2020Compatibility`

This class implements the Materials Project 2020 energy correction scheme, which
incorporates uncertainty quantification and allows for mixing of GGA and GGA+U entries
(see References).

Note that this scheme should only be applied to VASP calculations that use the
Materials Project input set parameters (see pymatgen.io.vasp.sets.MPRelaxSet). Using
this compatibility scheme on calculations with different parameters is not valid.

Note: While the correction scheme is largely composition-based, the energy corrections
applied to ComputedEntry and ComputedStructureEntry can differ for O and S-containing
structures if entry.data[‘oxidation_states’] is not populated or explicitly set. This
occurs because pymatgen will use atomic distances to classify O and S anions as
superoxide/peroxide/oxide and sulfide/polysulfide, resp. when oxidation states are not
provided. If you want the most accurate corrections possible, supply pre-defined
oxidation states to entry.data or pass ComputedStructureEntry.


* **Parameters**


    * **compat_type** – Two options, GGA or Advanced. GGA means all GGA+U
    entries are excluded. Advanced means the GGA/GGA+U mixing scheme
    of Jain et al. (see References) is implemented. In this case,
    entries which are supposed to be calculated in GGA+U (i.e.,
    transition metal oxides and fluorides) will have the corresponding
    GGA entries excluded. For example, Fe oxides should
    have a U value under the Advanced scheme. An Fe oxide run in GGA
    will therefore be excluded.

    To use the “Advanced” type, Entry.parameters must contain a “hubbards”
    key which is a dict of all non-zero Hubbard U values used in the
    calculation. For example, if you ran a Fe2O3 calculation with
    Materials Project parameters, this would look like
    entry.parameters[“hubbards”] = {“Fe”: 5.3}. If the “hubbards” key
    is missing, a GGA run is assumed. Entries obtained from the
    MaterialsProject database will automatically have these fields
    populated. Default: “Advanced”



    * **correct_peroxide** – Specify whether peroxide/superoxide/ozonide
    corrections are to be applied or not. If false, all oxygen-containing
    compounds are assigned the ‘oxide’ correction. Default: True


    * **check_potcar** (*bool*) – Check that the POTCARs used in the calculation are consistent
    with the Materials Project parameters. False bypasses this check altogether. Default: True
    Can also be disabled globally by running pmg config –add PMG_POTCAR_CHECKS false.


    * **check_potcar_hash** (*bool*) – Use potcar hash to verify POTCAR settings are
    consistent with MPRelaxSet. If False, only the POTCAR symbols will
    be used. Default: False


    * **config_file** (*Path*) – Path to the selected compatibility.yaml config file.
    If None, defaults to MP2020Compatibility.yaml distributed with
    pymatgen.


### References

Wang, A., Kingsbury, R., McDermott, M., Horton, M., Jain. A., Ong, S.P.,

    Dwaraknath, S., Persson, K. A framework for quantifying uncertainty
    in DFT energy corrections. Scientific Reports 11: 15496, 2021.
    [https://doi.org/10.1038/s41598-021-94550-5](https://doi.org/10.1038/s41598-021-94550-5)

Jain, A. et al. Formation enthalpies by mixing GGA and GGA + U calculations.

    Phys. Rev. B - Condens. Matter Mater. Phys. 84, 1-10 (2011).


### _class_ pymatgen.entries.compatibility.MaterialsProjectAqueousCompatibility(\*args, \*\*kwargs)
Bases: `MaterialsProjectAqueousCompatibility`

This class implements the Aqueous energy referencing scheme for constructing
Pourbaix diagrams from DFT energies, as described in Persson et al.

This scheme applies various energy adjustments to convert DFT energies into
Gibbs free energies of formation at 298 K and to guarantee that the experimental
formation free energy of H2O is reproduced. Briefly, the steps are:

>
> 1. Beginning with the DFT energy of O2, adjust the energy of H2 so that
> the experimental reaction energy of -2.458 eV/H2O is reproduced.


> 2. Add entropy to the DFT energy of any compounds that are liquid or
> gaseous at room temperature


> 3. Adjust the DFT energies of solid hydrate compounds (compounds that
> contain water, e.g. FeO.nH2O) such that the energies of the embedded
> H2O molecules are equal to the experimental free energy

The above energy adjustments are computed dynamically based on the input
Entries.

### References

K.A. Persson, B. Waldwick, P. Lazic, G. Ceder, Prediction of solid-aqueous
equilibria: Scheme to combine first-principles calculations of solids with
experimental aqueous states, Phys. Rev. B - Condens. Matter Mater. Phys.
85 (2012) 1-12. doi:10.1103/PhysRevB.85.235438.

Initialize the MaterialsProjectAqueousCompatibility class.

Note that this class requires as inputs the ground-state DFT energies of O2 and H2O, plus the value of any
energy adjustments applied to an H2O molecule. If these parameters are not provided in __init__, they can
be automatically populated by including ComputedEntry for the ground state of O2 and H2O in a list of
entries passed to process_entries. process_entries will fail if one or the other is not provided.


* **Parameters**


    * **solid_compat** – Compatibility scheme used to pre-process solid DFT energies prior to applying aqueous
    energy adjustments. May be passed as a class (e.g. MaterialsProject2020Compatibility) or an instance
    (e.g., MaterialsProject2020Compatibility()). If None, solid DFT energies are used as-is.
    Default: MaterialsProject2020Compatibility


    * **o2_energy** – The ground-state DFT energy of oxygen gas, including any adjustments or corrections, in eV/atom.
    If not set, this value will be determined from any O2 entries passed to process_entries.
    Default: None


    * **h2o_energy** – The ground-state DFT energy of water, including any adjustments or corrections, in eV/atom.
    If not set, this value will be determined from any H2O entries passed to process_entries.
    Default: None


    * **h2o_adjustments** – Total energy adjustments applied to one water molecule, in eV/atom.
    If not set, this value will be determined from any H2O entries passed to process_entries.
    Default: None



### _class_ pymatgen.entries.compatibility.MaterialsProjectCompatibility(compat_type: str = 'Advanced', correct_peroxide: bool = True, check_potcar_hash: bool = False)
Bases: `CorrectionsList`

This class implements the GGA/GGA+U mixing scheme, which allows mixing of
entries. Note that this should only be used for VASP calculations using the
MaterialsProject parameters (see pymatgen.io.vasp.sets.MPVaspInputSet).
Using this compatibility scheme on runs with different parameters is not
valid.


* **Parameters**


    * **compat_type** – Two options, GGA or Advanced. GGA means all GGA+U
    entries are excluded. Advanced means mixing scheme is
    implemented to make entries compatible with each other,
    but entries which are supposed to be done in GGA+U will have the
    equivalent GGA entries excluded. For example, Fe oxides should
    have a U value under the Advanced scheme. A GGA Fe oxide run
    will therefore be excluded under the scheme.


    * **correct_peroxide** – Specify whether peroxide/superoxide/ozonide
    corrections are to be applied or not.


    * **check_potcar_hash** (*bool*) – Use potcar hash to verify potcars are correct.


    * **silence_deprecation** (*bool*) – Silence deprecation warning. Defaults to False.



### _class_ pymatgen.entries.compatibility.PotcarCorrection(input_set: type[[pymatgen.io.vasp.sets.VaspInputSet](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.VaspInputSet)], check_potcar: bool = True, check_hash: bool = False)
Bases: `Correction`

Checks that POTCARs are valid within a pre-defined input set. This
ensures that calculations performed using different InputSets are not
compared against each other.

Entry.parameters must contain a “potcar_symbols” key that is a list of
all POTCARs used in the run. Again, using the example of an Fe2O3 run
using Materials Project parameters, this would look like
entry.parameters[“potcar_symbols”] = [‘PAW_PBE Fe_pv 06Sep2000’,
‘PAW_PBE O 08Apr2002’].


* **Parameters**


    * **input_set** ([*InputSet*](pymatgen.io.md#pymatgen.io.core.InputSet)) – object used to generate the runs (used to check
    for correct potcar symbols).


    * **check_potcar** (*bool*) – If False, bypass the POTCAR check altogether. Defaults to True.
    Can also be disabled globally by running pmg config –add PMG_POTCAR_CHECKS false.


    * **check_hash** (*bool*) – If True, uses the potcar hash to check for valid
    potcars. If false, uses the potcar symbol (less reliable). Defaults to False.



* **Raises**

    **ValueError** – if check_potcar=True and entry does not contain “potcar_symbols” key.



#### get_correction(entry: ComputedEntry | ComputedStructureEntry)

* **Parameters**

    **entry** (*AnyComputedEntry*) – ComputedEntry or ComputedStructureEntry.



* **Raises**


    * **ValueError** – If entry does not contain “potcar_symbols” key.


    * **CompatibilityError** – If entry has wrong potcar hash/symbols.



* **Returns**

    0.0 +/- 0.0 (from uncertainties package)



* **Return type**

    ufloat



### _class_ pymatgen.entries.compatibility.UCorrection(\*args, \*\*kwargs)
Bases: `UCorrection`

This class implements the GGA/GGA+U mixing scheme, which allows mixing of
entries. Entry.parameters must contain a “hubbards” key which is a dict
of all non-zero Hubbard U values used in the calculation. For example,
if you ran a Fe2O3 calculation with Materials Project parameters,
this would look like entry.parameters[“hubbards”] = {“Fe”: 5.3}
If the “hubbards” key is missing, a GGA run is assumed.

It should be noted that ComputedEntries assimilated using the
pymatgen.apps.borg package and obtained via the MaterialsProject REST
interface using the pymatgen.matproj.rest package will automatically have
these fields populated.


* **Parameters**


    * **config_file** – Path to the selected compatibility.yaml config file.


    * **input_set** – InputSet object (to check for the +U settings)


    * **compat_type** – Two options, GGA or Advanced. GGA means all GGA+U
    entries are excluded. Advanced means mixing scheme is
    implemented to make entries compatible with each other,
    but entries which are supposed to be done in GGA+U will have the
    equivalent GGA entries excluded. For example, Fe oxides should
    have a U value under the Advanced scheme. A GGA Fe oxide run
    will therefore be excluded under the scheme.


    * **error_file** – Path to the selected compatibilityErrors.yaml config file.


## pymatgen.entries.computed_entries module

This module implements equivalents of the basic ComputedEntry objects, which
is the basic entity that can be used to perform many analyses. ComputedEntries
contain calculated information, typically from VASP or other electronic
structure codes. For example, ComputedEntries can be used as inputs for phase
diagram analysis.


### _class_ pymatgen.entries.computed_entries.CompositionEnergyAdjustment(adj_per_atom, n_atoms, uncertainty_per_atom=nan, name='', cls=None, description='Composition-based energy adjustment')
Bases: `EnergyAdjustment`

An energy adjustment applied to a ComputedEntry based on the atomic composition.
Used in various DFT energy correction schemes.


* **Parameters**


    * **adj_per_atom** – float, energy adjustment to apply per atom, in eV/atom


    * **n_atoms** – float or int, number of atoms.


    * **uncertainty_per_atom** – float, uncertainty in energy adjustment to apply per atom, in eV/atom.
    (Default: np.nan)


    * **name** – str, human-readable name of the energy adjustment.
    (Default: “”)


    * **cls** – dict, Serialized Compatibility class used to generate the energy
    adjustment. (Default: None)


    * **description** – str, human-readable explanation of the energy adjustment.



#### _property_ explain()
Return an explanation of how the energy adjustment is calculated.


#### normalize(factor)
Normalize energy adjustment (in place), dividing value/uncertainty by a
factor.
:param factor: factor to divide by.


#### _property_ uncertainty()
Return the value of the energy adjustment in eV.


#### _property_ value()
Return the value of the energy adjustment in eV.


### _class_ pymatgen.entries.computed_entries.ComputedEntry(composition: [Composition](pymatgen.core.md#pymatgen.core.composition.Composition) | str | dict[str, float], energy: float, correction: float = 0.0, energy_adjustments: list | None = None, parameters: dict | None = None, data: dict | None = None, entry_id: object | None = None)
Bases: `Entry`

Lightweight Entry object for computed data. Contains facilities
for applying corrections to the energy attribute and for storing
calculation parameters.

Initializes a ComputedEntry.


* **Parameters**


    * **composition** ([*Composition*](pymatgen.core.md#pymatgen.core.composition.Composition)) – Composition of the entry. For
    flexibility, this can take the form of all the typical input
    taken by a Composition, including a {symbol: amt} dict,
    a string formula, and others.


    * **energy** (*float*) – Energy of the entry. Usually the final calculated
    energy from VASP or other electronic structure codes.


    * **correction** (*float*) – Manually set an energy correction, will ignore
    energy_adjustments if specified.


    * **energy_adjustments** – An optional list of EnergyAdjustment to
    be applied to the energy. This is used to modify the energy for
    certain analyses. Defaults to None.


    * **parameters** – An optional dict of parameters associated with
    the entry. Defaults to None.


    * **data** – An optional dict of any additional data associated
    with the entry. Defaults to None.


    * **entry_id** – An optional id to uniquely identify the entry.



#### as_dict()

* **Returns**

    MSONable dict.



#### copy()
Returns a copy of the ComputedEntry.


#### _property_ correction(_: floa_ )
Returns:
float: the total energy correction / adjustment applied to the entry,

> in eV.


#### _property_ correction_per_atom(_: floa_ )
Returns:
float: the total energy correction / adjustment applied to the entry,

> normalized by atoms (units of eV/atom).


#### _property_ correction_uncertainty(_: floa_ )
Returns:
float: the uncertainty of the energy adjustments applied to the entry, in eV.


#### _property_ correction_uncertainty_per_atom(_: floa_ )
Returns:
float: the uncertainty of the energy adjustments applied to the entry,

> normalized by atoms (units of eV/atom).


#### _property_ energy(_: floa_ )
the *corrected* energy of the entry.


* **Type**

    return



#### _classmethod_ from_dict(d)

* **Parameters**

    **d** – Dict representation.



* **Returns**

    ComputedEntry



#### normalize(mode: Literal['formula_unit', 'atom'] = 'formula_unit')
Normalize the entry’s composition and energy.


* **Parameters**

    **mode** (*"formula_unit"** | **"atom"*) – “formula_unit” (the default) normalizes to composition.reduced_formula.
    “atom” normalizes such that the composition amounts sum to 1.



#### _property_ uncorrected_energy(_: floa_ )
Returns:
float: the *uncorrected* energy of the entry.


#### _property_ uncorrected_energy_per_atom(_: floa_ )
Returns:
float: the *uncorrected* energy of the entry, normalized by atoms

> (units of eV/atom).


### _class_ pymatgen.entries.computed_entries.ComputedStructureEntry(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), energy: float, correction: float = 0.0, composition: [Composition](pymatgen.core.md#pymatgen.core.composition.Composition) | str | dict[str, float] | None = None, energy_adjustments: list | None = None, parameters: dict | None = None, data: dict | None = None, entry_id: object | None = None)
Bases: `ComputedEntry`

A heavier version of ComputedEntry which contains a structure as well. The
structure is needed for some analyses.

Initializes a ComputedStructureEntry.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – The actual structure of an entry.


    * **energy** (*float*) – Energy of the entry. Usually the final calculated
    energy from VASP or other electronic structure codes.


    * **correction** (*float**, **optional*) – A correction to the energy. This is mutually exclusive with
    energy_adjustments, i.e. pass either or neither but not both. Defaults to 0.


    * **composition** ([*Composition*](pymatgen.core.md#pymatgen.core.composition.Composition)) – Composition of the entry. For
    flexibility, this can take the form of all the typical input
    taken by a Composition, including a {symbol: amt} dict,
    a string formula, and others.


    * **energy_adjustments** – An optional list of EnergyAdjustment to
    be applied to the energy. This is used to modify the energy for
    certain analyses. Defaults to None.


    * **parameters** – An optional dict of parameters associated with
    the entry. Defaults to None.


    * **data** – An optional dict of any additional data associated
    with the entry. Defaults to None.


    * **entry_id** – An optional id to uniquely identify the entry.



#### as_dict()

* **Returns**

    MSONable dict.



#### copy()
Returns a copy of the ComputedStructureEntry.


#### _classmethod_ from_dict(d)

* **Parameters**

    **d** – Dict representation.



* **Returns**

    ComputedStructureEntry



#### normalize(mode: Literal['formula_unit', 'atom'] = 'formula_unit')
Normalize the entry’s composition and energy. The structure remains unchanged.


* **Parameters**

    **mode** (*"formula_unit"** | **"atom"*) – “formula_unit” (the default) normalizes to composition.reduced_formula.
    “atom” normalizes such that the composition amounts sum to 1.



#### _property_ structure(_: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure_ )
the structure of the entry.


* **Type**

    return



### _class_ pymatgen.entries.computed_entries.ConstantEnergyAdjustment(value, uncertainty=nan, name='Constant energy adjustment', cls=None, description='Constant energy adjustment')
Bases: `EnergyAdjustment`

A constant energy adjustment applied to a ComputedEntry. Useful in energy referencing
schemes such as the Aqueous energy referencing scheme.


* **Parameters**


    * **value** – float, value of the energy adjustment in eV


    * **uncertainty** – float, uncertainty of the energy adjustment in eV. (Default: np.nan)


    * **name** – str, human-readable name of the energy adjustment.
    (Default: Constant energy adjustment)


    * **cls** – dict, Serialized Compatibility class used to generate the energy
    adjustment. (Default: None)


    * **description** – str, human-readable explanation of the energy adjustment.



#### _property_ explain()
Return an explanation of how the energy adjustment is calculated.


#### normalize(factor)
Normalize energy adjustment (in place), dividing value/uncertainty by a
factor.
:param factor: factor to divide by.


### _class_ pymatgen.entries.computed_entries.EnergyAdjustment(value, uncertainty=nan, name='Manual adjustment', cls=None, description='')
Bases: `MSONable`

Lightweight class to contain information about an energy adjustment or
energy correction.


* **Parameters**


    * **value** – float, value of the energy adjustment in eV


    * **uncertainty** – float, uncertainty of the energy adjustment in eV. Default: np.nan


    * **name** – str, human-readable name of the energy adjustment.
    (Default: Manual adjustment)


    * **cls** – dict, Serialized Compatibility class used to generate the energy adjustment. (Default: None)


    * **description** – str, human-readable explanation of the energy adjustment.



#### _abstract property_ explain()
Return an explanation of how the energy adjustment is calculated.


#### _abstract_ normalize(factor)
Scale the value of the current energy adjustment by factor in-place.

This method is utilized in ComputedEntry.normalize() to scale the energies to a formula unit basis
(e.g. E_Fe6O9 = 3 x E_Fe2O3).


#### _property_ uncertainty()
Return the uncertainty in the value of the energy adjustment in eV.


#### _property_ value()
Return the value of the energy correction in eV.


### _class_ pymatgen.entries.computed_entries.GibbsComputedStructureEntry(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), formation_enthalpy_per_atom: float, temp: float = 300, gibbs_model: Literal['SISSO'] = 'SISSO', composition: [Composition](pymatgen.core.md#pymatgen.core.composition.Composition) | None = None, correction: float = 0.0, energy_adjustments: list | None = None, parameters: dict | None = None, data: dict | None = None, entry_id: object | None = None)
Bases: `ComputedStructureEntry`

An extension to ComputedStructureEntry which includes the estimated Gibbs
free energy of formation via a machine-learned model.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – The pymatgen Structure object of an entry.


    * **formation_enthalpy_per_atom** (*float*) – Formation enthalpy of the entry;


    * **be** (*must*) – calculated using phase diagram construction (eV)


    * **temp** (*float*) – Temperature in Kelvin. If temperature is not selected from
    one of [300, 400, 500, … 2000 K], then free energies will
    be interpolated. Defaults to 300 K.


    * **gibbs_model** (*'SISSO'*) – Model for Gibbs Free energy. “SISSO”, the descriptor
    created by Bartel et al. (2018) – see reference in documentation, is
    currently the only supported option.


    * **composition** ([*Composition*](pymatgen.core.md#pymatgen.core.composition.Composition)) – The composition of the entry. Defaults to None.


    * **correction** (*float*) – A correction to be applied to the energy. Defaults to 0.


    * **energy_adjustments** (*list*) – A list of energy adjustments to be applied to
    the energy. Defaults to None.


    * **parameters** (*dict*) – An optional dict of parameters associated with
    the entry. Defaults to None.


    * **data** (*dict*) – An optional dict of any additional data associated
    with the entry. Defaults to None.


    * **entry_id** – An optional id to uniquely identify the entry.



#### as_dict()

* **Returns**

    MSONable dict.



#### _classmethod_ from_dict(d)

* **Parameters**

    **d** – Dict representation.



* **Returns**

    GibbsComputedStructureEntry



#### _classmethod_ from_entries(entries, temp=300, gibbs_model='SISSO')
Constructor method for initializing GibbsComputedStructureEntry objects from
T = 0 K ComputedStructureEntry objects, as acquired from a thermochemical
database e.g. The Materials Project.


* **Parameters**


    * **entries** (*[**ComputedStructureEntry**]*) – List of ComputedStructureEntry objects,
    as downloaded from The Materials Project API.


    * **temp** (*int*) – Temperature [K] for estimating Gibbs free energy of formation.


    * **gibbs_model** (*str*) – Gibbs model to use; currently the only option is “SISSO”.



* **Returns**

    list of new entries which replace the orig.

        entries with inclusion of Gibbs free energy of formation at the
        specified temperature.




* **Return type**

    [GibbsComputedStructureEntry]



#### _classmethod_ from_pd(pd, temp=300, gibbs_model='SISSO')
Constructor method for initializing a list of GibbsComputedStructureEntry
objects from an existing T = 0 K phase diagram composed of
ComputedStructureEntry objects, as acquired from a thermochemical database;
(e.g.. The Materials Project).


* **Parameters**


    * **pd** ([*PhaseDiagram*](pymatgen.analysis.md#pymatgen.analysis.phase_diagram.PhaseDiagram)) – T = 0 K phase diagram as created in pymatgen. Must
    contain ComputedStructureEntry objects.


    * **temp** (*int*) – Temperature [K] for estimating Gibbs free energy of formation.


    * **gibbs_model** (*str*) – Gibbs model to use; currently the only option is “SISSO”.



* **Returns**

    list of new entries which replace the orig.

        entries with inclusion of Gibbs free energy of formation at the
        specified temperature.




* **Return type**

    [GibbsComputedStructureEntry]



#### gf_sisso()
Gibbs Free Energy of formation as calculated by SISSO descriptor from Bartel
et al. (2018). Units: eV (not normalized).

WARNING: This descriptor only applies to solids. The implementation here
attempts to detect and use downloaded NIST-JANAF data for common
experimental gases (e.g. CO2) where possible. Note that experimental data is
only for Gibbs Free Energy of formation, so expt. entries will register as
having a formation enthalpy of 0.

Reference: Bartel, C. J., Millican, S. L., Deml, A. M., Rumptz, J. R.,
Tumas, W., Weimer, A. W., … Holder, A. M. (2018). Physical descriptor for
the Gibbs energy of inorganic crystalline solids and
temperature-dependent materials chemistry. Nature Communications, 9(1),
4168. [https://doi.org/10.1038/s41467-018-06682-4](https://doi.org/10.1038/s41467-018-06682-4)


* **Returns**

    the difference between formation enthalpy (T=0 K, Materials
    Project) and the predicted Gibbs free energy of formation  (eV)



* **Return type**

    float



### _class_ pymatgen.entries.computed_entries.ManualEnergyAdjustment(value)
Bases: `ConstantEnergyAdjustment`

A manual energy adjustment applied to a ComputedEntry.


* **Parameters**

    **value** – float, value of the energy adjustment in eV.



### _class_ pymatgen.entries.computed_entries.TemperatureEnergyAdjustment(adj_per_deg, temp, n_atoms, uncertainty_per_deg=nan, name='', cls=None, description='Temperature-based energy adjustment')
Bases: `EnergyAdjustment`

An energy adjustment applied to a ComputedEntry based on the temperature.
Used, for example, to add entropy to DFT energies.


* **Parameters**


    * **adj_per_deg** – float, energy adjustment to apply per degree K, in eV/atom


    * **temp** – float, temperature in Kelvin


    * **n_atoms** – float or int, number of atoms


    * **uncertainty_per_deg** – float, uncertainty in energy adjustment to apply per degree K,
    in eV/atom. (Default: np.nan)


    * **name** – str, human-readable name of the energy adjustment.
    (Default: “”)


    * **cls** – dict, Serialized Compatibility class used to generate the energy
    adjustment. (Default: None)


    * **description** – str, human-readable explanation of the energy adjustment.



#### _property_ explain()
Return an explanation of how the energy adjustment is calculated.


#### normalize(factor)
Normalize energy adjustment (in place), dividing value/uncertainty by a
factor.
:param factor: factor to divide by.


#### _property_ uncertainty()
Return the value of the energy adjustment in eV.


#### _property_ value()
Return the value of the energy correction in eV.

## pymatgen.entries.correction_calculator module

This module calculates corrections for the species listed below, fitted to the experimental and computed
entries given to the CorrectionCalculator constructor.


### _class_ pymatgen.entries.correction_calculator.CorrectionCalculator(species: list[str] | None = None, max_error: float = 0.1, allow_unstable: float | bool = 0.1, exclude_polyanions: list[str] | None = None)
Bases: `object`

A CorrectionCalculator contains experimental and computed entries which it uses to compute corrections.

It graphs residual errors after applying the computed corrections and creates the MPCompatibility.yaml
file the Correction classes use.


#### species()
list of species that corrections are being calculated for


#### exp_compounds()
list of dictionaries which each contain a compound’s formula and experimental data


#### calc_compounds()
dictionary of ComputedEntry objects


#### corrections()
list of corrections in same order as species list


#### corrections_std_error()
list of the variances of the corrections in same order as species list


#### corrections_dict()
dictionary of format {‘species’: (value, uncertainty)} for easier correction lookup

Initializes a CorrectionCalculator.


* **Parameters**


    * **species** – list of species to calculate corrections for


    * **max_error** – maximum tolerable relative uncertainty in experimental energy.
    Compounds with relative uncertainty greater than this value will be excluded from the fit


    * **allow_unstable** – whether unstable entries are to be included in the fit. If True, all compounds will
    be included regardless of their energy above hull. If False or a float, compounds with
    energy above hull greater than the given value (defaults to 0.1 eV/atom) will be
    excluded


    * **exclude_polyanions** – a list of polyanions that contain additional sources of error that may negatively
    influence the quality of the fitted corrections. Compounds with these polyanions
    will be excluded from the fit



#### compute_corrections(exp_entries: list, calc_entries: dict)
Computes the corrections and fills in correction, corrections_std_error, and corrections_dict.


* **Parameters**


    * **exp_entries** – list of dictionary objects with the following keys/values:
    {“formula”: chemical formula, “exp energy”: formation energy in eV/formula unit,
    “uncertainty”: uncertainty in formation energy}


    * **calc_entries** – dictionary of computed entries, of the form {chemical formula: ComputedEntry}



* **Raises**

    **ValueError** – calc_compounds is missing an entry



#### compute_from_files(exp_gz: str, comp_gz: str)

* **Parameters**


    * **exp_gz** – name of .json.gz file that contains experimental data
    data in .json.gz file should be a list of dictionary objects with the following keys/values:
    {“formula”: chemical formula, “exp energy”: formation energy in eV/formula unit,
    “uncertainty”: uncertainty in formation energy}


    * **comp_gz** – name of .json.gz file that contains computed entries
    data in .json.gz file should be a dictionary of {chemical formula: ComputedEntry}.



#### graph_residual_error()
Graphs the residual errors for all compounds after applying computed corrections.


#### graph_residual_error_per_species(specie: str)
Graphs the residual errors for each compound that contains specie after applying computed corrections.


* **Parameters**

    **specie** – the specie/group that residual errors are being plotted for



* **Raises**

    **ValueError** – the specie is not a valid specie that this class fits corrections for



#### make_yaml(name: str = 'MP2020', dir: str | None = None)
Creates the _name_Compatibility.yaml that stores corrections as well as _name_CompatibilityUncertainties.yaml
for correction uncertainties.


* **Parameters**


    * **name** – str, alternate name for the created .yaml file.
    Default: “MP2020”


    * **dir** – str, directory in which to save the file. Pass None (default) to
    save the file in the current working directory.


## pymatgen.entries.entry_tools module

This module implements functions to perform various useful operations on
entries, such as grouping entries by structure.


### _class_ pymatgen.entries.entry_tools.EntrySet(entries: Iterable[[PDEntry](pymatgen.analysis.md#pymatgen.analysis.phase_diagram.PDEntry) | ComputedEntry | ComputedStructureEntry])
Bases: `MutableSet`, `MSONable`

A convenient container for manipulating entries. Allows for generating
subsets, dumping into files, etc.


* **Parameters**

    **entries** – All the entries.



#### add(element)
Add an entry.


* **Parameters**

    **element** – Entry



#### as_dict()
Returns MSONable dict.


#### _property_ chemsys(_: se_ )
Returns:
set representing the chemical system, e.g., {“Li”, “Fe”, “P”, “O”}.


#### discard(element)
Discard an entry.


* **Parameters**

    **element** – Entry



#### _classmethod_ from_csv(filename: str)
Imports PDEntries from a csv.


* **Parameters**

    **filename** – Filename to import from.



* **Returns**

    List of Elements, List of PDEntries



#### get_subset_in_chemsys(chemsys: list[str])
Returns an EntrySet containing only the set of entries belonging to
a particular chemical system (in this definition, it includes all sub
systems). For example, if the entries are from the
Li-Fe-P-O system, and chemsys=[“Li”, “O”], only the Li, O,
and Li-O entries are returned.


* **Parameters**

    **chemsys** – Chemical system specified as list of elements. E.g.,
    [“Li”, “O”]



* **Returns**

    EntrySet



#### _property_ ground_states(_: se_ )
A set containing only the entries that are ground states, i.e., the lowest energy
per atom entry at each composition.


#### is_ground_state(entry)
Boolean indicating whether a given Entry is a ground state.


#### remove_non_ground_states()
Removes all non-ground state entries, i.e., only keep the lowest energy
per atom entry at each composition.


#### to_csv(filename: str, latexify_names: bool = False)
Exports PDEntries to a csv.


* **Parameters**


    * **filename** – Filename to write to.


    * **entries** – PDEntries to export.


    * **latexify_names** – Format entry names to be LaTex compatible,
    e.g., Li_{2}O



### pymatgen.entries.entry_tools.group_entries_by_composition(entries, sort_by_e_per_atom=True)
Given a sequence of Entry-like objects, group them by composition and

    optionally sort by energy above hull.


* **Parameters**


    * **entries** (*List*) – Sequence of Entry-like objects.


    * **sort_by_e_per_atom** (*bool*) – Whether to sort the grouped entries by
    energy per atom (lowest energy first). Default True.



* **Returns**

    Sequence of sequence of entries by composition. e.g,
    [[ entry1, entry2], [entry3, entry4, entry5]]



### pymatgen.entries.entry_tools.group_entries_by_structure(entries, species_to_remove=None, ltol=0.2, stol=0.4, angle_tol=5, primitive_cell=True, scale=True, comparator=None, ncpus=None)
Given a sequence of ComputedStructureEntries, use structure fitter to group
them by structural similarity.


* **Parameters**


    * **entries** – Sequence of ComputedStructureEntries.


    * **species_to_remove** – Sometimes you want to compare a host framework
    (e.g., in Li-ion battery analysis). This allows you to specify
    species to remove before structural comparison.


    * **ltol** (*float*) – Fractional length tolerance. Default is 0.2.


    * **stol** (*float*) – Site tolerance in Angstrom. Default is 0.4 Angstrom.


    * **angle_tol** (*float*) – Angle tolerance in degrees. Default is 5 degrees.


    * **primitive_cell** (*bool*) – If true: input structures will be reduced to
    primitive cells prior to matching. Defaults to True.


    * **scale** – Input structures are scaled to equivalent volume if true;
    For exact matching, set to False.


    * **comparator** – A comparator object implementing an equals method that
    declares equivalency of sites. Default is SpeciesComparator,
    which implies rigid species mapping.


    * **ncpus** – Number of cpus to use. Use of multiple cpus can greatly improve
    fitting speed. Default of None means serial processing.



* **Returns**

    Sequence of sequence of entries by structural similarity. e.g,
    [[ entry1, entry2], [entry3, entry4, entry5]]


## pymatgen.entries.exp_entries module

This module defines Entry classes for containing experimental data.


### _class_ pymatgen.entries.exp_entries.ExpEntry(composition, thermodata, temperature=298)
Bases: [`PDEntry`](pymatgen.analysis.md#pymatgen.analysis.phase_diagram.PDEntry), `MSONable`

An lightweight ExpEntry object containing experimental data for a
composition for many purposes. Extends a PDEntry so that it can be used for
phase diagram generation and reaction calculation.

Current version works only with solid phases and at 298K. Further
extensions for temperature dependence are planned.


* **Parameters**


    * **composition** – Composition of the entry. For flexibility, this can take
    the form of all the typical input taken by a Composition, including
    a {symbol: amt} dict, a string formula, and others.


    * **thermodata** – A sequence of ThermoData associated with the entry.


    * **temperature** – A temperature for the entry in Kelvin. Defaults to 298K.



#### as_dict()

* **Returns**

    MSONable dict



#### _classmethod_ from_dict(d)

* **Parameters**

    **d** – Dict representation.



* **Returns**

    ExpEntry


## pymatgen.entries.mixing_scheme module

This module implements Compatibility corrections for mixing runs of different
functionals.


### _class_ pymatgen.entries.mixing_scheme.MaterialsProjectDFTMixingScheme(structure_matcher: StructureMatcher | None = None, run_type_1: str = 'GGA(+U)', run_type_2: str = 'R2SCAN', compat_1: Compatibility | None = <pymatgen.entries.compatibility.cached_class.<locals>._decorated object>, compat_2: Compatibility | None = None, fuzzy_matching: bool = True, check_potcar: bool = True)
Bases: `Compatibility`

This class implements the Materials Project mixing scheme, which allows mixing of
energies from different DFT functionals. Note that this should only be used for
VASP calculations using the MaterialsProject parameters (e.g. MPRelaxSet or
MPScanRelaxSet). Using this compatibility scheme on runs with different parameters
may lead to unexpected results.

This is the scheme used by the Materials Project to generate Phase Diagrams containing
a mixture of GGA(+U) and R2SCAN calculations. However in principle it can be used to
mix energies from any two functionals.

Instantiate the mixing scheme. The init method creates a generator class that
contains relevant settings (e.g., StructureMatcher instance, Compatibility settings
for each functional) for processing groups of entries.


* **Parameters**


    * **structure_matcher** ([*StructureMatcher*](pymatgen.analysis.md#pymatgen.analysis.structure_matcher.StructureMatcher)) – StructureMatcher object used to determine
    whether calculations from different functionals describe the same material.


    * **run_type_1** – The first DFT run_type. Typically this is the majority or run type or
    the “base case” onto which the other calculations are referenced. Valid choices
    are any run_type recognized by Vasprun.run_type, such as “LDA”, “GGA”, “GGA+U”,
    “PBEsol”, “SCAN”, or “R2SCAN”. The class will ignore any entries that have a
    run_type different than run_type_1 or run_type_2.

    The list of run_type_1 entries provided to process_entries MUST form a complete
    Phase Diagram in order for the mixing scheme to work. If this condition is not
    satisfied, processing the entries will fail.

    Note that the special string “GGA(+U)” (default) will treat both GGA and GGA+U
    calculations as a single type. This option exists because GGA/GGA+U mixing is
    already handled by MaterialsProject2020Compatibility.



    * **run_type_2** – The second DFT run_type. Typically this is the run_type that is ‘preferred’
    but has fewer calculations. If run_type_1 and run_type_2 calculations exist for all
    materials, run_type_2 energies will be used (hence the ‘preferred’ status). The class
    will ignore any entries that have a run_type different than run_type_1 or run_type_2.


    * **compat_1** – Compatibility class used to pre-process entries of run_type_1.
    Defaults to MaterialsProjectCompatibility2020.


    * **compat_2** – Compatibility class used to pre-process entries of run_type_2.
    Defaults to None.


    * **fuzzy_matching** – Whether to use less strict structure matching logic for
    diatomic elements O2, N2, F2, H2, and Cl2 as well as I and Br. Outputs of DFT
    relaxations using
    different functionals frequently fail to structure match for these elements
    even though they come from the same original material. Fuzzy structure matching
    considers the materials equivalent if the formula, number of sites, and
    space group are all identical. If there are multiple materials of run_type_2
    that satisfy these criteria, the one with lowest energy is considered to
    match.


    * **check_potcar** – Whether to ensure the POTCARs used for the run_type_1 and run_type_2 calculations
    are the same. This is useful for ensuring that the mixing scheme is not used on calculations
    that used different POTCARs, which can lead to unphysical results. Defaults to True.
    Has no effect if neither compat_1 nor compat_2 have a check_potcar attribute.
    Can also be disabled globally by running pmg config –add PMG_POTCAR_CHECKS false.



#### _static_ display_entries(entries)
Generate a pretty printout of key properties of a list of ComputedEntry.


#### get_adjustments(entry, mixing_state_data: pd.DataFrame | None = None)
Returns the corrections applied to a particular entry. Note that get_adjustments is not
intended to be called directly in the R2SCAN mixing scheme. Call process_entries instead,
and it will pass the required arguments to get_adjustments.


* **Parameters**


    * **entry** – A ComputedEntry object. The entry must be a member of the list of entries
    used to create mixing_state_data.


    * **mixing_state_data** – A DataFrame containing information about which Entries
    correspond to the same materials, which are stable on the phase diagrams of
    the respective run_types, etc. Can be generated from a list of entries using
    MaterialsProjectDFTMixingScheme.get_mixing_state_data. This argument is included to
    facilitate use of the mixing scheme in high-throughput databases where an alternative
    to get_mixing_state_data is desirable for performance reasons. In general, it should
    always be left at the default value (None) to avoid inconsistencies between the mixing
    state data and the properties of the ComputedStructureEntry.



* **Returns**

    Energy adjustments to be applied to entry.



* **Return type**

    [EnergyAdjustment]



* **Raises**

    **CompatibilityError if the DFT mixing scheme cannot be applied to the entry.** –



#### get_mixing_state_data(entries: list[pymatgen.entries.computed_entries.ComputedStructureEntry])
Generate internal state data to be passed to get_adjustments.


* **Parameters**

    **entries** – The list of ComputedStructureEntry to process. It is assumed that the entries have
    already been filtered using _filter_and_sort_entries() to remove any irrelevant run types,
    apply compat_1 and compat_2, and confirm that all have unique entry_id.



* **Returns**

    A pandas DataFrame that contains information associating structures from

        different functionals with specific materials and establishing how many run_type_1
        ground states have been computed with run_type_2. The DataFrame contains one row
        for each distinct material (Structure), with the following columns:

        > formula: str the reduced_formula
        > spacegroup: int the spacegroup
        > num_sites: int the number of sites in the Structure
        > entry_id_1: the entry_id of the run_type_1 entry
        > entry_id_2: the entry_id of the run_type_2 entry
        > run_type_1: Optional[str] the run_type_1 value
        > run_type_2: Optional[str] the run_type_2 value
        > energy_1: float or nan the ground state energy in run_type_1 in eV/atom
        > energy_2: float or nan the ground state energy in run_type_2 in eV/atom
        > is_stable_1: bool whether this material is stable on the run_type_1 PhaseDiagram
        > hull_energy_1: float or nan the energy of the run_type_1 hull at this composition in eV/atom
        > hull_energy_2: float or nan the energy of the run_type_1 hull at this composition in eV/atom

    None: Returns None if the supplied ComputedStructureEntry are insufficient for applying

        the mixing scheme.




* **Return type**

    DataFrame



#### process_entries(entries: AnyComputedEntry | list[AnyComputedEntry], clean: bool = True, verbose: bool = True, inplace: bool = True, mixing_state_data=None)
Process a sequence of entries with the DFT mixing scheme. Note
that this method will change the data of the original entries.


* **Parameters**


    * **entries** – ComputedEntry or [ComputedEntry]. Pass all entries as a single list, even if they are
    computed with different functionals or require different preprocessing. This list will
    automatically be filtered based on run_type_1 and run_type_2, and processed according to
    compat_1 and compat_2.

    Note that under typical use, when mixing_state_data=None, the entries MUST be
    ComputedStructureEntry. They will be matched using structure_matcher.



    * **clean** (*bool*) – Whether to remove any previously-applied energy adjustments.
    If True, all EnergyAdjustment are removed prior to processing the Entry.
    Default is True.


    * **verbose** (*bool*) – Whether to print verbose error messages about the mixing scheme. Default is True.


    * **inplace** (*bool*) – Whether to adjust input entries in place. Default is True.


    * **mixing_state_data** – A DataFrame containing information about which Entries
    correspond to the same materials, which are stable on the phase diagrams of
    the respective run_types, etc. If None (default), it will be generated from the
    list of entries using MaterialsProjectDFTMixingScheme.get_mixing_state_data.
    This argument is included to facilitate use of the mixing scheme in high-throughput
    databases where an alternative to get_mixing_state_data is desirable for performance
    reasons. In general, it should always be left at the default value (None) to avoid
    inconsistencies between the mixing state data and the properties of the
    ComputedStructureEntry in entries.



* **Returns**

    Adjusted entries. Entries in the original list incompatible with

        chosen correction scheme are excluded from the returned list.




* **Return type**

    list[AnyComputedEntry]