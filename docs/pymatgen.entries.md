---
layout: default
title: pymatgen.entries.md
nav_exclude: true
---

1. TOC
{:toc}

# pymatgen.entries package

Entries are containers for calculated information, which is used in
many analyses. This module contains entry related tools and implements
the base Entry class, which is the basic entity that can be used to
store calculated information. Other Entry classes such as ComputedEntry
and PDEntry inherit from this class.


### _class_ Entry(composition: [Composition](pymatgen.core.md#pymatgen.core.composition.Composition) | str | dict[str, float], energy: float)
Bases: `MSONable`

A lightweight object containing the energy associated with
a specific chemical composition. This base class is not
intended to be instantiated directly. Note that classes
which inherit from Entry must define a .energy property.

Initializes an Entry.


* **Parameters**


    * **composition** ([*Composition*](pymatgen.core.md#pymatgen.core.composition.Composition)) – Composition of the entry. For
    flexibility, this can take the form of all the typical input taken by a
    Composition, including a {symbol: amt} dict, a string formula, and others.


    * **energy** (*float*) – Energy of the entry.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _normalization_factor(mode: Literal['formula_unit', 'atom'] = 'formula_unit')

#### as_dict()
MSONable dict.


#### _property_ composition(_: [Composition](pymatgen.core.md#pymatgen.core.composition.Composition_ )
The composition of the entry.


#### _property_ elements(_: list[[Element](pymatgen.core.md#pymatgen.core.periodic_table.Element) | [Species](pymatgen.core.md#pymatgen.core.periodic_table.Species) | [DummySpecies](pymatgen.core.md#pymatgen.core.periodic_table.DummySpecies)_ )
The set of elements in the entry.


#### _abstract property_ energy(_: floa_ )
The energy of the entry.


#### _property_ energy_per_atom(_: floa_ )
The energy per atom of the entry.


#### _property_ is_element(_: boo_ )
Whether composition of entry is an element.


#### normalize(mode: Literal['formula_unit', 'atom'] = 'formula_unit')
Normalize the entry’s composition and energy.


* **Parameters**

    **mode** (*"formula_unit"** | **"atom"*) – “formula_unit” (the default) normalizes to composition.reduced_formula.
    “atom” normalizes such that the composition amounts sum to 1.



## pymatgen.entries.compatibility module

This module implements Compatibility corrections for mixing runs of different
functionals.


### _class_ AnionCorrection(\*args, \*\*kwargs)
Bases: `AnionCorrection`

Correct anion energies to obtain the right formation energies. Note that
this depends on calculations being run within the same input set.

Used by legacy MaterialsProjectCompatibility and MITCompatibility.


* **Parameters**


    * **config_file** – Path to the selected compatibility.yaml config file.


    * **correct_peroxide** – Specify whether peroxide/superoxide/ozonide
    corrections are to be applied or not.



#### _abc_impl(_ = <_abc._abc_data object_ )

### _class_ AqueousCorrection(\*args, \*\*kwargs)
Bases: `AqueousCorrection`

This class implements aqueous phase compound corrections for elements
and H2O.

Used only by MITAqueousCompatibility.


* **Parameters**


    * **config_file** – Path to the selected compatibility.yaml config file.


    * **error_file** – Path to the selected compatibilityErrors.yaml config file.



#### _abc_impl(_ = <_abc._abc_data object_ )

### _class_ Compatibility()
Bases: `MSONable`

Abstract Compatibility class, not intended for direct use.
Compatibility classes are used to correct the energies of an entry or a set
of entries. All Compatibility classes must implement get_adjustments() method.


#### _abc_impl(_ = <_abc._abc_data object_ )

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



### _exception_ CompatibilityError()
Bases: `Exception`

Exception class for Compatibility. Raised by attempting correction
on incompatible calculation.


### _class_ Correction()
Bases: `object`

A Correction class is a pre-defined scheme for correction a computed
entry based on the type and chemistry of the structure and the
calculation parameters. All Correction classes must implement a
correct_entry method.


#### _abc_impl(_ = <_abc._abc_data object_ )

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



### _class_ CorrectionsList(corrections: Sequence[Correction])
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



#### _abc_impl(_ = <_abc._abc_data object_ )

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



### _class_ GasCorrection(\*args, \*\*kwargs)
Bases: `GasCorrection`

Correct gas energies to obtain the right formation energies. Note that
this depends on calculations being run within the same input set.
Used by legacy MaterialsProjectCompatibility and MITCompatibility.


* **Parameters**

    **config_file** – Path to the selected compatibility.yaml config file.



#### _abc_impl(_ = <_abc._abc_data object_ )

### _class_ MITAqueousCompatibility(compat_type: Literal['GGA', 'Advanced'] = 'Advanced', correct_peroxide: bool = True, check_potcar_hash: bool = False)
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



#### _abc_impl(_ = <_abc._abc_data object_ )

### _class_ MITCompatibility(compat_type: Literal['GGA', 'Advanced'] = 'Advanced', correct_peroxide: bool = True, check_potcar_hash: bool = False)
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



#### _abc_impl(_ = <_abc._abc_data object_ )

### _class_ MaterialsProject2020Compatibility(\*args, \*\*kwargs)
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


#### _abc_impl(_ = <_abc._abc_data object_ )

### _class_ MaterialsProjectAqueousCompatibility(\*args, \*\*kwargs)
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



#### _abc_impl(_ = <_abc._abc_data object_ )

### _class_ MaterialsProjectCompatibility(compat_type: str = 'Advanced', correct_peroxide: bool = True, check_potcar_hash: bool = False)
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



#### _abc_impl(_ = <_abc._abc_data object_ )

### _class_ PotcarCorrection(input_set: type[[VaspInputSet](pymatgen.io.vasp.md#pymatgen.io.vasp.sets.VaspInputSet)], check_potcar: bool = True, check_hash: bool = False)
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



#### _abc_impl(_ = <_abc._abc_data object_ )

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



### _class_ UCorrection(\*args, \*\*kwargs)
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



#### _abc_impl(_ = <_abc._abc_data object_ )
## pymatgen.entries.computed_entries module

This module implements equivalents of the basic ComputedEntry objects, which
is the basic entity that can be used to perform many analyses. ComputedEntries
contain calculated information, typically from VASP or other electronic
structure codes. For example, ComputedEntries can be used as inputs for phase
diagram analysis.


### _class_ CompositionEnergyAdjustment(adj_per_atom, n_atoms, uncertainty_per_atom=nan, name='', cls=None, description='Composition-based energy adjustment')
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


### _class_ ComputedEntry(composition: [Composition](pymatgen.core.md#pymatgen.core.composition.Composition) | str | dict[str, float], energy: float, correction: float = 0.0, energy_adjustments: list | None = None, parameters: dict | None = None, data: dict | None = None, entry_id: object | None = None)
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



#### _abc_impl(_ = <_abc._abc_data object_ )

#### as_dict()
MSONable dict.


#### copy()
Returns a copy of the ComputedEntry.


#### _property_ correction(_: floa_ )
Returns:
float: the total energy correction / adjustment applied to the entry in eV.


#### _property_ correction_per_atom(_: floa_ )
Returns:
float: the total energy correction / adjustment applied to the entry in eV/atom.


#### _property_ correction_uncertainty(_: floa_ )
Returns:
float: the uncertainty of the energy adjustments applied to the entry in eV.


#### _property_ correction_uncertainty_per_atom(_: floa_ )
Returns:
float: the uncertainty of the energy adjustments applied to the entry in eV/atom.


#### _property_ energy(_: floa_ )
The *corrected* energy of the entry.


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
float: the *uncorrected* energy of the entry, normalized by atoms in eV/atom.


### _class_ ComputedStructureEntry(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), energy: float, correction: float = 0.0, composition: [Composition](pymatgen.core.md#pymatgen.core.composition.Composition) | str | dict[str, float] | None = None, energy_adjustments: list | None = None, parameters: dict | None = None, data: dict | None = None, entry_id: object | None = None)
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



#### _abc_impl(_ = <_abc._abc_data object_ )

#### as_dict()
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
The structure of the entry.


### _class_ ConstantEnergyAdjustment(value, uncertainty=nan, name='Constant energy adjustment', cls=None, description='Constant energy adjustment')
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


### _class_ EnergyAdjustment(value, uncertainty=nan, name='Manual adjustment', cls=None, description='')
Bases: `MSONable`

Lightweight class to contain information about an energy adjustment or
energy correction.


* **Parameters**


    * **value** (*float*) – value of the energy adjustment in eV


    * **uncertainty** (*float*) – uncertainty of the energy adjustment in eV. Default: np.nan


    * **name** (*str*) – human-readable name of the energy adjustment.
    (Default: Manual adjustment)


    * **cls** (*dict*) – Serialized Compatibility class used to generate the energy adjustment. Defaults to {}.


    * **description** (*str*) – human-readable explanation of the energy adjustment.



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


### _class_ GibbsComputedStructureEntry(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), formation_enthalpy_per_atom: float, temp: float = 300, gibbs_model: Literal['SISSO'] = 'SISSO', composition: [Composition](pymatgen.core.md#pymatgen.core.composition.Composition) | None = None, correction: float = 0.0, energy_adjustments: list | None = None, parameters: dict | None = None, data: dict | None = None, entry_id: object | None = None)
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



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _static_ _g_delta_sisso(vol_per_atom, reduced_mass, temp)
G^delta as predicted by SISSO-learned descriptor from Eq. (4) in
Bartel et al. (2018).


* **Parameters**


    * **vol_per_atom** (*float*) – volume per atom [Å^3/atom]


    * **reduced_mass** (*float*) – as calculated with pair-wise sum formula
    [amu]


    * **temp** (*float*) – Temperature [K]



* **Returns**

    G^delta [eV/atom]



* **Return type**

    float



#### _static_ _reduced_mass(structure)
Reduced mass as calculated via Eq. 6 in Bartel et al. (2018).


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – The pymatgen Structure object of the entry.



* **Returns**

    reduced mass (amu)



* **Return type**

    float



#### _sum_g_i()
Sum of the stoichiometrically weighted chemical potentials of the elements
at specified temperature, as acquired from “g_els.json”.


* **Returns**

    sum of weighted chemical potentials [eV]



* **Return type**

    float



#### as_dict()
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



### _class_ ManualEnergyAdjustment(value)
Bases: `ConstantEnergyAdjustment`

A manual energy adjustment applied to a ComputedEntry.


* **Parameters**

    **value** – float, value of the energy adjustment in eV.



### _class_ TemperatureEnergyAdjustment(adj_per_deg, temp, n_atoms, uncertainty_per_deg=nan, name='', cls=None, description='Temperature-based energy adjustment')
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


### _class_ CorrectionCalculator(species: list[str] | None = None, max_error: float = 0.1, allow_unstable: float | bool = 0.1, exclude_polyanions: list[str] | None = None)
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


### _class_ EntrySet(entries: Iterable[[PDEntry](pymatgen.analysis.md#pymatgen.analysis.phase_diagram.PDEntry) | ComputedEntry | ComputedStructureEntry])
Bases: `MutableSet`, `MSONable`

A convenient container for manipulating entries. Allows for generating
subsets, dumping into files, etc.


* **Parameters**

    **entries** – All the entries.



#### _abc_impl(_ = <_abc._abc_data object_ )

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



### _get_host(structure, species_to_remove)

### _perform_grouping(args)

### group_entries_by_composition(entries, sort_by_e_per_atom=True)
Given a sequence of Entry-like objects, group them by composition and

    optionally sort by energy above hull.


* **Parameters**


    * **entries** (*list*) – Sequence of Entry-like objects.


    * **sort_by_e_per_atom** (*bool*) – Whether to sort the grouped entries by
    energy per atom (lowest energy first). Default True.



* **Returns**

    Sequence of sequence of entries by composition. e.g,
    [[ entry1, entry2], [entry3, entry4, entry5]]



### group_entries_by_structure(entries, species_to_remove=None, ltol=0.2, stol=0.4, angle_tol=5, primitive_cell=True, scale=True, comparator=None, ncpus=None)
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


### _class_ ExpEntry(composition, thermodata, temperature=298)
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



#### _abc_impl(_ = <_abc._abc_data object_ )

#### as_dict()
MSONable dict.


#### _classmethod_ from_dict(d)

* **Parameters**

    **d** – Dict representation.



* **Returns**

    ExpEntry


## pymatgen.entries.mixing_scheme module

This module implements Compatibility corrections for mixing runs of different
functionals.


### _class_ MaterialsProjectDFTMixingScheme(structure_matcher: StructureMatcher | None = None, run_type_1: str = 'GGA(+U)', run_type_2: str = 'R2SCAN', compat_1: Compatibility | None = <pymatgen.entries.compatibility.cached_class.<locals>._decorated object>, compat_2: Compatibility | None = None, fuzzy_matching: bool = True, check_potcar: bool = True)
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



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _filter_and_sort_entries(entries, verbose=True)
Given a single list of entries, separate them by run_type and return two lists, one containing
only entries of each run_type.


#### _populate_df_row(struct_group, comp, sg, n, pd_type_1, pd_type_2, all_entries)
Helper function to populate a row of the mixing state DataFrame, given
a list of matched structures.


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



#### get_mixing_state_data(entries: list[ComputedStructureEntry])
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