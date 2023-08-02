---
layout: default
title: pymatgen.entries.computed_entries.md
nav_exclude: true
---

# pymatgen.entries.computed_entries module

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


### _class_ pymatgen.entries.computed_entries.ComputedEntry(composition: [Composition](pymatgen.core.composition.md#pymatgen.core.composition.Composition) | str | dict[str, float], energy: float, correction: float = 0.0, energy_adjustments: list | None = None, parameters: dict | None = None, data: dict | None = None, entry_id: object | None = None)
Bases: [`Entry`](pymatgen.entries.md#pymatgen.entries.Entry)

Lightweight Entry object for computed data. Contains facilities
for applying corrections to the energy attribute and for storing
calculation parameters.

Initializes a ComputedEntry.


* **Parameters**


    * **composition** ([*Composition*](pymatgen.core.composition.md#pymatgen.core.composition.Composition)) – Composition of the entry. For
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


### _class_ pymatgen.entries.computed_entries.ComputedStructureEntry(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), energy: float, correction: float = 0.0, composition: [Composition](pymatgen.core.composition.md#pymatgen.core.composition.Composition) | str | dict[str, float] | None = None, energy_adjustments: list | None = None, parameters: dict | None = None, data: dict | None = None, entry_id: object | None = None)
Bases: `ComputedEntry`

A heavier version of ComputedEntry which contains a structure as well. The
structure is needed for some analyses.

Initializes a ComputedStructureEntry.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – The actual structure of an entry.


    * **energy** (*float*) – Energy of the entry. Usually the final calculated
    energy from VASP or other electronic structure codes.


    * **correction** (*float**, **optional*) – A correction to the energy. This is mutually exclusive with
    energy_adjustments, i.e. pass either or neither but not both. Defaults to 0.


    * **composition** ([*Composition*](pymatgen.core.composition.md#pymatgen.core.composition.Composition)) – Composition of the entry. For
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



#### _property_ structure(_: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure_ )
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


### _class_ pymatgen.entries.computed_entries.GibbsComputedStructureEntry(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), formation_enthalpy_per_atom: float, temp: float = 300, gibbs_model: Literal['SISSO'] = 'SISSO', composition: [Composition](pymatgen.core.composition.md#pymatgen.core.composition.Composition) | None = None, correction: float = 0.0, energy_adjustments: list | None = None, parameters: dict | None = None, data: dict | None = None, entry_id: object | None = None)
Bases: `ComputedStructureEntry`

An extension to ComputedStructureEntry which includes the estimated Gibbs
free energy of formation via a machine-learned model.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – The pymatgen Structure object of an entry.


    * **formation_enthalpy_per_atom** (*float*) – Formation enthalpy of the entry;


    * **be** (*must*) – calculated using phase diagram construction (eV)


    * **temp** (*float*) – Temperature in Kelvin. If temperature is not selected from
    one of [300, 400, 500, … 2000 K], then free energies will
    be interpolated. Defaults to 300 K.


    * **gibbs_model** (*'SISSO'*) – Model for Gibbs Free energy. “SISSO”, the descriptor
    created by Bartel et al. (2018) – see reference in documentation, is
    currently the only supported option.


    * **composition** ([*Composition*](pymatgen.core.composition.md#pymatgen.core.composition.Composition)) – The composition of the entry. Defaults to None.


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


    * **pd** ([*PhaseDiagram*](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram)) – T = 0 K phase diagram as created in pymatgen. Must
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