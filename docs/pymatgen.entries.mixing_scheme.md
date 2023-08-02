---
layout: default
title: pymatgen.entries.mixing_scheme.md
nav_exclude: true
---

# pymatgen.entries.mixing_scheme module

This module implements Compatibility corrections for mixing runs of different
functionals.


### _class_ pymatgen.entries.mixing_scheme.MaterialsProjectDFTMixingScheme(structure_matcher: StructureMatcher | None = None, run_type_1: str = 'GGA(+U)', run_type_2: str = 'R2SCAN', compat_1: Compatibility | None = <pymatgen.entries.compatibility.cached_class.<locals>._decorated object>, compat_2: Compatibility | None = None, fuzzy_matching: bool = True, check_potcar: bool = True)
Bases: [`Compatibility`](pymatgen.entries.compatibility.md#pymatgen.entries.compatibility.Compatibility)

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


    * **structure_matcher** ([*StructureMatcher*](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.StructureMatcher)) – StructureMatcher object used to determine
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

    [[EnergyAdjustment](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.EnergyAdjustment)]



* **Raises**

    **CompatibilityError if the DFT mixing scheme cannot be applied to the entry.** –



#### get_mixing_state_data(entries: list[[pymatgen.entries.computed_entries.ComputedStructureEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry)])
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