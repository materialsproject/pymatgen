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


### _class_ pymatgen.entries.Entry(composition: [Composition](pymatgen.core.composition.md#pymatgen.core.composition.Composition) | str | dict[str, float], energy: float)
Bases: `MSONable`

A lightweight object containing the energy associated with
a specific chemical composition. This base class is not
intended to be instantiated directly. Note that classes
which inherit from Entry must define a .energy property.

Initializes an Entry.


* **Parameters**


    * **composition** ([*Composition*](pymatgen.core.composition.md#pymatgen.core.composition.Composition)) – Composition of the entry. For
    flexibility, this can take the form of all the typical input
    taken by a Composition, including a {symbol: amt} dict,
    a string formula, and others.


    * **energy** (*float*) – Energy of the entry.



#### as_dict()
MSONable dict.


#### _property_ composition(_: [Composition](pymatgen.core.composition.md#pymatgen.core.composition.Composition_ )
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




* [pymatgen.entries.compatibility module](pymatgen.entries.compatibility.md)


    * [`AnionCorrection`](pymatgen.entries.compatibility.md#pymatgen.entries.compatibility.AnionCorrection)


    * [`AqueousCorrection`](pymatgen.entries.compatibility.md#pymatgen.entries.compatibility.AqueousCorrection)


    * [`Compatibility`](pymatgen.entries.compatibility.md#pymatgen.entries.compatibility.Compatibility)


        * [`Compatibility.explain()`](pymatgen.entries.compatibility.md#pymatgen.entries.compatibility.Compatibility.explain)


        * [`Compatibility.get_adjustments()`](pymatgen.entries.compatibility.md#pymatgen.entries.compatibility.Compatibility.get_adjustments)


        * [`Compatibility.process_entries()`](pymatgen.entries.compatibility.md#pymatgen.entries.compatibility.Compatibility.process_entries)


        * [`Compatibility.process_entry()`](pymatgen.entries.compatibility.md#pymatgen.entries.compatibility.Compatibility.process_entry)


    * [`CompatibilityError`](pymatgen.entries.compatibility.md#pymatgen.entries.compatibility.CompatibilityError)


    * [`Correction`](pymatgen.entries.compatibility.md#pymatgen.entries.compatibility.Correction)


        * [`Correction.correct_entry()`](pymatgen.entries.compatibility.md#pymatgen.entries.compatibility.Correction.correct_entry)


        * [`Correction.get_correction()`](pymatgen.entries.compatibility.md#pymatgen.entries.compatibility.Correction.get_correction)


    * [`CorrectionsList`](pymatgen.entries.compatibility.md#pymatgen.entries.compatibility.CorrectionsList)


        * [`CorrectionsList.explain()`](pymatgen.entries.compatibility.md#pymatgen.entries.compatibility.CorrectionsList.explain)


        * [`CorrectionsList.get_adjustments()`](pymatgen.entries.compatibility.md#pymatgen.entries.compatibility.CorrectionsList.get_adjustments)


        * [`CorrectionsList.get_corrections_dict()`](pymatgen.entries.compatibility.md#pymatgen.entries.compatibility.CorrectionsList.get_corrections_dict)


        * [`CorrectionsList.get_explanation_dict()`](pymatgen.entries.compatibility.md#pymatgen.entries.compatibility.CorrectionsList.get_explanation_dict)


    * [`GasCorrection`](pymatgen.entries.compatibility.md#pymatgen.entries.compatibility.GasCorrection)


    * [`MITAqueousCompatibility`](pymatgen.entries.compatibility.md#pymatgen.entries.compatibility.MITAqueousCompatibility)


    * [`MITCompatibility`](pymatgen.entries.compatibility.md#pymatgen.entries.compatibility.MITCompatibility)


    * [`MaterialsProject2020Compatibility`](pymatgen.entries.compatibility.md#pymatgen.entries.compatibility.MaterialsProject2020Compatibility)


    * [`MaterialsProjectAqueousCompatibility`](pymatgen.entries.compatibility.md#pymatgen.entries.compatibility.MaterialsProjectAqueousCompatibility)


    * [`MaterialsProjectCompatibility`](pymatgen.entries.compatibility.md#pymatgen.entries.compatibility.MaterialsProjectCompatibility)


    * [`PotcarCorrection`](pymatgen.entries.compatibility.md#pymatgen.entries.compatibility.PotcarCorrection)


        * [`PotcarCorrection.get_correction()`](pymatgen.entries.compatibility.md#pymatgen.entries.compatibility.PotcarCorrection.get_correction)


    * [`UCorrection`](pymatgen.entries.compatibility.md#pymatgen.entries.compatibility.UCorrection)


* [pymatgen.entries.computed_entries module](pymatgen.entries.computed_entries.md)


    * [`CompositionEnergyAdjustment`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.CompositionEnergyAdjustment)


        * [`CompositionEnergyAdjustment.explain`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.CompositionEnergyAdjustment.explain)


        * [`CompositionEnergyAdjustment.normalize()`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.CompositionEnergyAdjustment.normalize)


        * [`CompositionEnergyAdjustment.uncertainty`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.CompositionEnergyAdjustment.uncertainty)


        * [`CompositionEnergyAdjustment.value`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.CompositionEnergyAdjustment.value)


    * [`ComputedEntry`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry)


        * [`ComputedEntry.as_dict()`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry.as_dict)


        * [`ComputedEntry.copy()`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry.copy)


        * [`ComputedEntry.correction`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry.correction)


        * [`ComputedEntry.correction_per_atom`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry.correction_per_atom)


        * [`ComputedEntry.correction_uncertainty`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry.correction_uncertainty)


        * [`ComputedEntry.correction_uncertainty_per_atom`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry.correction_uncertainty_per_atom)


        * [`ComputedEntry.energy`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry.energy)


        * [`ComputedEntry.from_dict()`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry.from_dict)


        * [`ComputedEntry.normalize()`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry.normalize)


        * [`ComputedEntry.uncorrected_energy`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry.uncorrected_energy)


        * [`ComputedEntry.uncorrected_energy_per_atom`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry.uncorrected_energy_per_atom)


    * [`ComputedStructureEntry`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry)


        * [`ComputedStructureEntry.as_dict()`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry.as_dict)


        * [`ComputedStructureEntry.copy()`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry.copy)


        * [`ComputedStructureEntry.from_dict()`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry.from_dict)


        * [`ComputedStructureEntry.normalize()`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry.normalize)


        * [`ComputedStructureEntry.structure`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry.structure)


    * [`ConstantEnergyAdjustment`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ConstantEnergyAdjustment)


        * [`ConstantEnergyAdjustment.explain`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ConstantEnergyAdjustment.explain)


        * [`ConstantEnergyAdjustment.normalize()`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ConstantEnergyAdjustment.normalize)


    * [`EnergyAdjustment`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.EnergyAdjustment)


        * [`EnergyAdjustment.explain`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.EnergyAdjustment.explain)


        * [`EnergyAdjustment.normalize()`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.EnergyAdjustment.normalize)


        * [`EnergyAdjustment.uncertainty`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.EnergyAdjustment.uncertainty)


        * [`EnergyAdjustment.value`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.EnergyAdjustment.value)


    * [`GibbsComputedStructureEntry`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.GibbsComputedStructureEntry)


        * [`GibbsComputedStructureEntry.as_dict()`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.GibbsComputedStructureEntry.as_dict)


        * [`GibbsComputedStructureEntry.from_dict()`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.GibbsComputedStructureEntry.from_dict)


        * [`GibbsComputedStructureEntry.from_entries()`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.GibbsComputedStructureEntry.from_entries)


        * [`GibbsComputedStructureEntry.from_pd()`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.GibbsComputedStructureEntry.from_pd)


        * [`GibbsComputedStructureEntry.gf_sisso()`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.GibbsComputedStructureEntry.gf_sisso)


    * [`ManualEnergyAdjustment`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ManualEnergyAdjustment)


    * [`TemperatureEnergyAdjustment`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.TemperatureEnergyAdjustment)


        * [`TemperatureEnergyAdjustment.explain`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.TemperatureEnergyAdjustment.explain)


        * [`TemperatureEnergyAdjustment.normalize()`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.TemperatureEnergyAdjustment.normalize)


        * [`TemperatureEnergyAdjustment.uncertainty`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.TemperatureEnergyAdjustment.uncertainty)


        * [`TemperatureEnergyAdjustment.value`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.TemperatureEnergyAdjustment.value)


* [pymatgen.entries.correction_calculator module](pymatgen.entries.correction_calculator.md)


    * [`CorrectionCalculator`](pymatgen.entries.correction_calculator.md#pymatgen.entries.correction_calculator.CorrectionCalculator)


        * [`CorrectionCalculator.species`](pymatgen.entries.correction_calculator.md#pymatgen.entries.correction_calculator.CorrectionCalculator.species)


        * [`CorrectionCalculator.exp_compounds`](pymatgen.entries.correction_calculator.md#pymatgen.entries.correction_calculator.CorrectionCalculator.exp_compounds)


        * [`CorrectionCalculator.calc_compounds`](pymatgen.entries.correction_calculator.md#pymatgen.entries.correction_calculator.CorrectionCalculator.calc_compounds)


        * [`CorrectionCalculator.corrections`](pymatgen.entries.correction_calculator.md#pymatgen.entries.correction_calculator.CorrectionCalculator.corrections)


        * [`CorrectionCalculator.corrections_std_error`](pymatgen.entries.correction_calculator.md#pymatgen.entries.correction_calculator.CorrectionCalculator.corrections_std_error)


        * [`CorrectionCalculator.corrections_dict`](pymatgen.entries.correction_calculator.md#pymatgen.entries.correction_calculator.CorrectionCalculator.corrections_dict)


        * [`CorrectionCalculator.compute_corrections()`](pymatgen.entries.correction_calculator.md#pymatgen.entries.correction_calculator.CorrectionCalculator.compute_corrections)


        * [`CorrectionCalculator.compute_from_files()`](pymatgen.entries.correction_calculator.md#pymatgen.entries.correction_calculator.CorrectionCalculator.compute_from_files)


        * [`CorrectionCalculator.graph_residual_error()`](pymatgen.entries.correction_calculator.md#pymatgen.entries.correction_calculator.CorrectionCalculator.graph_residual_error)


        * [`CorrectionCalculator.graph_residual_error_per_species()`](pymatgen.entries.correction_calculator.md#pymatgen.entries.correction_calculator.CorrectionCalculator.graph_residual_error_per_species)


        * [`CorrectionCalculator.make_yaml()`](pymatgen.entries.correction_calculator.md#pymatgen.entries.correction_calculator.CorrectionCalculator.make_yaml)


* [pymatgen.entries.entry_tools module](pymatgen.entries.entry_tools.md)


    * [`EntrySet`](pymatgen.entries.entry_tools.md#pymatgen.entries.entry_tools.EntrySet)


        * [`EntrySet.add()`](pymatgen.entries.entry_tools.md#pymatgen.entries.entry_tools.EntrySet.add)


        * [`EntrySet.as_dict()`](pymatgen.entries.entry_tools.md#pymatgen.entries.entry_tools.EntrySet.as_dict)


        * [`EntrySet.chemsys`](pymatgen.entries.entry_tools.md#pymatgen.entries.entry_tools.EntrySet.chemsys)


        * [`EntrySet.discard()`](pymatgen.entries.entry_tools.md#pymatgen.entries.entry_tools.EntrySet.discard)


        * [`EntrySet.from_csv()`](pymatgen.entries.entry_tools.md#pymatgen.entries.entry_tools.EntrySet.from_csv)


        * [`EntrySet.get_subset_in_chemsys()`](pymatgen.entries.entry_tools.md#pymatgen.entries.entry_tools.EntrySet.get_subset_in_chemsys)


        * [`EntrySet.ground_states`](pymatgen.entries.entry_tools.md#pymatgen.entries.entry_tools.EntrySet.ground_states)


        * [`EntrySet.is_ground_state()`](pymatgen.entries.entry_tools.md#pymatgen.entries.entry_tools.EntrySet.is_ground_state)


        * [`EntrySet.remove_non_ground_states()`](pymatgen.entries.entry_tools.md#pymatgen.entries.entry_tools.EntrySet.remove_non_ground_states)


        * [`EntrySet.to_csv()`](pymatgen.entries.entry_tools.md#pymatgen.entries.entry_tools.EntrySet.to_csv)


    * [`group_entries_by_composition()`](pymatgen.entries.entry_tools.md#pymatgen.entries.entry_tools.group_entries_by_composition)


    * [`group_entries_by_structure()`](pymatgen.entries.entry_tools.md#pymatgen.entries.entry_tools.group_entries_by_structure)


* [pymatgen.entries.exp_entries module](pymatgen.entries.exp_entries.md)


    * [`ExpEntry`](pymatgen.entries.exp_entries.md#pymatgen.entries.exp_entries.ExpEntry)


        * [`ExpEntry.as_dict()`](pymatgen.entries.exp_entries.md#pymatgen.entries.exp_entries.ExpEntry.as_dict)


        * [`ExpEntry.from_dict()`](pymatgen.entries.exp_entries.md#pymatgen.entries.exp_entries.ExpEntry.from_dict)


* [pymatgen.entries.mixing_scheme module](pymatgen.entries.mixing_scheme.md)


    * [`MaterialsProjectDFTMixingScheme`](pymatgen.entries.mixing_scheme.md#pymatgen.entries.mixing_scheme.MaterialsProjectDFTMixingScheme)


        * [`MaterialsProjectDFTMixingScheme.display_entries()`](pymatgen.entries.mixing_scheme.md#pymatgen.entries.mixing_scheme.MaterialsProjectDFTMixingScheme.display_entries)


        * [`MaterialsProjectDFTMixingScheme.get_adjustments()`](pymatgen.entries.mixing_scheme.md#pymatgen.entries.mixing_scheme.MaterialsProjectDFTMixingScheme.get_adjustments)


        * [`MaterialsProjectDFTMixingScheme.get_mixing_state_data()`](pymatgen.entries.mixing_scheme.md#pymatgen.entries.mixing_scheme.MaterialsProjectDFTMixingScheme.get_mixing_state_data)


        * [`MaterialsProjectDFTMixingScheme.process_entries()`](pymatgen.entries.mixing_scheme.md#pymatgen.entries.mixing_scheme.MaterialsProjectDFTMixingScheme.process_entries)