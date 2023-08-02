---
layout: default
title: pymatgen.entries.compatibility.md
nav_exclude: true
---

# pymatgen.entries.compatibility module

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



#### _abstract_ get_adjustments(entry: [ComputedEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry) | [ComputedStructureEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry))
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

    list[[EnergyAdjustment](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.EnergyAdjustment)]



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



#### process_entry(entry: [ComputedEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry), \*\*kwargs)
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



#### _abstract_ get_correction(entry: [ComputedEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry) | [ComputedStructureEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry))
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



#### get_adjustments(entry: [ComputedEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry) | [ComputedStructureEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry))
Get the list of energy adjustments to be applied to an entry.


#### get_corrections_dict(entry: [ComputedEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry) | [ComputedStructureEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry))
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



### _class_ pymatgen.entries.compatibility.PotcarCorrection(input_set: type[[pymatgen.io.vasp.sets.VaspInputSet](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.VaspInputSet)], check_potcar: bool = True, check_hash: bool = False)
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


    * **input_set** ([*InputSet*](pymatgen.io.core.md#pymatgen.io.core.InputSet)) – object used to generate the runs (used to check
    for correct potcar symbols).


    * **check_potcar** (*bool*) – If False, bypass the POTCAR check altogether. Defaults to True.
    Can also be disabled globally by running pmg config –add PMG_POTCAR_CHECKS false.


    * **check_hash** (*bool*) – If True, uses the potcar hash to check for valid
    potcars. If false, uses the potcar symbol (less reliable). Defaults to False.



* **Raises**

    **ValueError** – if check_potcar=True and entry does not contain “potcar_symbols” key.



#### get_correction(entry: [ComputedEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry) | [ComputedStructureEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry))

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