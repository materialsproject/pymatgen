---
layout: default
title: pymatgen.apps.battery.md
nav_exclude: true
---

1. TOC
{:toc}

# pymatgen.apps.battery package

This package contains all battery-related application classes, including
representations of InsertionElectrodes and ConversionElectrodes.


## pymatgen.apps.battery.analyzer module

Analysis classes for batteries.


### _class_ BatteryAnalyzer(struc_oxid, working_ion='Li', oxi_override=None)
Bases: `object`

A suite of methods for starting with an oxidized structure and determining its potential as a battery.

Pass in a structure for analysis.


* **Parameters**


    * **struc_oxid** – a Structure object; oxidation states *must* be assigned for this structure; disordered
    structures should be OK


    * **working_ion** – a String symbol or Element for the working ion.


    * **oxi_override** – a dict of String element symbol, Integer oxidation state pairs.
    by default, H, C, N, O, F, S, Cl, Se, Br, Te, I are considered anions.



#### _get_int_removals_helper(spec_amts_oxi, redox_el, redox_els, num_a)
This is a helper method for get_removals_int_oxid!


* **Parameters**


    * **spec_amts_oxi** – a dict of species to their amounts in the structure


    * **redox_el** – the element to oxidize or reduce


    * **redox_els** – the full list of elements that might be oxidized or reduced


    * **num_a** – a running set of numbers of A ion at integer oxidation steps



* **Returns**

    a set of numbers A; steps for oxidizing oxid_el first, then the other oxid_els in this list



#### _get_max_cap_ah(remove, insert)
Give max capacity in mAh for inserting and removing a charged ion
This method does not normalize the capacity and intended as a helper method.


#### get_max_capgrav(remove=True, insert=True)
Give max capacity in mAh/g for inserting and removing a charged ion
Note that the weight is normalized to the most ion-packed state,
thus removal of 1 Li from LiFePO4 gives the same capacity as insertion of 1 Li into FePO4.


* **Parameters**


    * **remove** – (bool) whether to allow ion removal


    * **insert** – (bool) whether to allow ion insertion



* **Returns**

    max grav capacity in mAh/g



#### get_max_capvol(remove=True, insert=True, volume=None)
Give max capacity in mAh/cc for inserting and removing a charged ion into base structure.


* **Parameters**


    * **remove** – (bool) whether to allow ion removal


    * **insert** – (bool) whether to allow ion insertion


    * **volume** – (float) volume to use for normalization (default=volume of initial structure)



* **Returns**

    max vol capacity in mAh/cc



#### get_removals_int_oxid()
Returns a set of ion removal steps, e.g. set([1 2 4]) etc. in order to
produce integer oxidation states of the redox metals.
If multiple redox metals are present, all combinations of reduction/oxidation are tested.
Note that having more than 3 redox metals will likely slow down the algorithm.

### Examples

LiFePO4 will return [1]
Li4Fe3Mn1(PO4)4 will return [1, 2, 3, 4])
Li6V4(PO4)6 will return [4, 6])  *note that this example is not normalized*


* **Returns**

    array of integer ion removals. If you double the unit cell, your answers will be twice as large!



#### _property_ max_ion_insertion()
Maximum number of ion A that can be inserted while maintaining charge-balance.
No consideration is given to whether there (geometrically speaking) are ion sites to actually accommodate the
extra ions.


* **Returns**

    integer amount of ion. Depends on cell size (this is an ‘extrinsic’ function!)



#### _property_ max_ion_removal()
Maximum number of ion A that can be removed while maintaining charge-balance.


* **Returns**

    integer amount of ion. Depends on cell size (this is an ‘extrinsic’ function!)



### is_redox_active_intercalation(element)
True if element is redox active and interesting for intercalation materials.


* **Parameters**

    **element** – Element object


## pymatgen.apps.battery.battery_abc module

This module defines the abstract base classes for battery-related classes.
Regardless of the kind of electrode, conversion or insertion, there are many
common definitions and properties, e.g., average voltage, capacity, etc. which
can be defined in a general way. The Abc for battery classes implements some of
these common definitions to allow sharing of common logic between them.


### _class_ AbstractElectrode(voltage_pairs: tuple[AbstractVoltagePair, ...], working_ion_entry: [ComputedEntry](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedEntry), framework_formula: str)
Bases: `Sequence`, `MSONable`

An Abstract Base Class representing an Electrode. It is essentially a
sequence of VoltagePairs. Generally, subclasses only need to implement
three abstract properties: voltage_pairs, working_ion and
working_ion_entry.

The general concept is that all other battery properties such as capacity,
etc. are derived from voltage pairs.

One of the major challenges with representing battery materials is keeping
track of the normalization between different entries. For example, one
entry might be TiO2 with one unit cell whereas another is LiTi2O4 with two
unit cells. When computing battery properties, it is needed to always use
a universal reference state otherwise you have normalization errors (e.g.,
the energy of LiTi2O4 must be divided by two to be compared with TiO2).

For properties such as volume, mass, or mAh transferred within the voltage
pair, a universal convention is necessary. AbstractElectrode can query for
extrinsic properties of several different AbstractVoltagePairs belonging to
a single charge/discharge path and be confident that the normalization is
being carried out properly throughout, even if more AbstractVoltagePairs
are added later.

The universal normalization is defined by the reduced structural framework
of the entries, which is common along the entire charge/discharge path. For
example, LiTi2O4 has a reduced structural framework of TiO2. Another
example is Li9V6P16O58 which would have a reduced structural framework of
V3P8O29. Note that reduced structural frameworks need not be
charge-balanced or physical, e.g. V3P8O29 is not charge-balanced, they are
just a tool for normalization.

Example: for a LiTi2O4 -> TiO2 AbstractVoltagePair, extrinsic quantities
like mAh or cell volumes are given per TiO2 formula unit.

Developers implementing a new battery (other than the two general ones
already implemented) need to implement a VoltagePair and an Electrode.


#### voltage_pairs()
Objects that represent each voltage step


* **Type**

    tuple[AbstractVoltagePair, …]



#### working_ion()
Representation of the working ion that only contains element type


#### working_ion_entry()
Representation of the working_ion that contains the energy


* **Type**

    [ComputedEntry](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedEntry)



#### framework_formula()
The compositions of one formula unit of the host material


* **Type**

    str



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _select_in_voltage_range(min_voltage=None, max_voltage=None)
Selects VoltagePairs within a certain voltage range.


* **Parameters**


    * **min_voltage** (*float*) – The minimum allowable voltage for a given
    step.


    * **max_voltage** (*float*) – The maximum allowable voltage allowable for a
    given step.



* **Returns**

    A list of VoltagePair objects



#### _property_ framework()
The composition object representing the framework.


#### framework_formula(_: st_ )

#### get_average_voltage(min_voltage=None, max_voltage=None)
Average voltage for path satisfying between a min and max voltage.


* **Parameters**


    * **min_voltage** (*float*) – The minimum allowable voltage for a given
    step.


    * **max_voltage** (*float*) – The maximum allowable voltage allowable for a
    given step.



* **Returns**

    Average voltage in V across the insertion path (a subset of the
    path can be chosen by the optional arguments)



#### get_capacity_grav(min_voltage=None, max_voltage=None, use_overall_normalization=True)
Get the gravimetric capacity of the electrode.


* **Parameters**


    * **min_voltage** (*float*) – The minimum allowable voltage for a given
    step.


    * **max_voltage** (*float*) – The maximum allowable voltage allowable for a
    given step.


    * **use_overall_normalization** (*booL*) – If False, normalize by the
    discharged state of only the voltage pairs matching the voltage
    criteria. if True, use default normalization of the full
    electrode path.



* **Returns**

    Gravimetric capacity in mAh/g across the insertion path (a subset
    of the path can be chosen by the optional arguments).



#### get_capacity_vol(min_voltage=None, max_voltage=None, use_overall_normalization=True)
Get the volumetric capacity of the electrode.


* **Parameters**


    * **min_voltage** (*float*) – The minimum allowable voltage for a given
    step.


    * **max_voltage** (*float*) – The maximum allowable voltage allowable for a
    given step.


    * **use_overall_normalization** (*booL*) – If False, normalize by the
    discharged state of only the voltage pairs matching the voltage
    criteria. if True, use default normalization of the full
    electrode path.



* **Returns**

    Volumetric capacity in mAh/cc across the insertion path (a subset
    of the path can be chosen by the optional arguments)



#### get_energy_density(min_voltage=None, max_voltage=None, use_overall_normalization=True)

* **Parameters**


    * **min_voltage** (*float*) – The minimum allowable voltage for a given
    step.


    * **max_voltage** (*float*) – The maximum allowable voltage allowable for a
    given step.


    * **use_overall_normalization** (*booL*) – If False, normalize by the
    discharged state of only the voltage pairs matching the voltage
    criteria. if True, use default normalization of the full
    electrode path.



* **Returns**

    Energy density in Wh/L across the insertion path (a subset of the
    path can be chosen by the optional arguments).



#### get_specific_energy(min_voltage=None, max_voltage=None, use_overall_normalization=True)
Returns the specific energy of the battery in mAh/g.


* **Parameters**


    * **min_voltage** (*float*) – The minimum allowable voltage for a given
    step.


    * **max_voltage** (*float*) – The maximum allowable voltage allowable for a
    given step.


    * **use_overall_normalization** (*booL*) – If False, normalize by the
    discharged state of only the voltage pairs matching the voltage
    criteria. if True, use default normalization of the full
    electrode path.



* **Returns**

    Specific energy in Wh/kg across the insertion path (a subset of
    the path can be chosen by the optional arguments)



#### get_sub_electrodes(adjacent_only=True)
If this electrode contains multiple voltage steps, then it is possible
to use only a subset of the voltage steps to define other electrodes.
Must be implemented for each electrode object.


* **Parameters**

    **adjacent_only** – Only return electrodes from compounds that are
    adjacent on the convex hull, i.e. no electrodes returned
    will have multiple voltage steps if this is set true



* **Returns**

    A list of Electrode objects



#### get_summary_dict(print_subelectrodes=True)
Generate a summary dict.


* **Parameters**

    **print_subelectrodes** – Also print data on all the possible
    subelectrodes.



* **Returns**

    A summary of this electrode’s properties in dict format.



#### _property_ max_delta_volume()
Maximum volume change along insertion.


#### _property_ max_voltage()
Highest voltage along insertion.


#### _property_ max_voltage_step()
Maximum absolute difference in adjacent voltage steps.


#### _property_ min_voltage()
Lowest voltage along insertion.


#### _property_ normalization_mass()
Mass used for normalization. This is the mass of the discharged
electrode of the last voltage pair.


* **Type**

    Returns



#### _property_ normalization_volume()
Mass used for normalization. This is the vol of the discharged
electrode of the last voltage pair.


* **Type**

    Returns



#### _property_ num_steps()
The number of distinct voltage steps in from fully charge to discharge
based on the stable intermediate states.


#### voltage_pairs(_: tuple[AbstractVoltagePair, ..._ )

#### _property_ working_ion()
Working ion as pymatgen Element object.


#### working_ion_entry(_: [ComputedEntry](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedEntry_ )

#### _property_ x_charge(_: floa_ )
The number of working ions per formula unit of host in the charged state.


#### _property_ x_discharge(_: floa_ )
The number of working ions per formula unit of host in the discharged state.


### _class_ AbstractVoltagePair(voltage: float, mAh: float, mass_charge: float, mass_discharge: float, vol_charge: float, vol_discharge: float, frac_charge: float, frac_discharge: float, working_ion_entry: [ComputedEntry](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedEntry), framework_formula: str)
Bases: `MSONable`

An Abstract Base Class for a Voltage Pair.


#### voltage()
Voltage of voltage pair.


* **Type**

    float



#### mAh()
Energy in mAh.


* **Type**

    float



#### mass_charge()
Mass of charged pair.


* **Type**

    float



#### mass_discharge()
Mass of discharged pair.


* **Type**

    float



#### vol_charge()
Vol of charged pair.


* **Type**

    float



#### vol_discharge()
Vol of discharged pair.


* **Type**

    float



#### frac_charge()
Frac of working ion in charged pair.


* **Type**

    float



#### frac_discharge()
Frac of working ion in discharged pair.


* **Type**

    float



#### working_ion_entry()
Working ion as an entry.


* **Type**

    [ComputedEntry](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedEntry)



#### framework_formula()
The compositions of one formula unit of the host material


* **Type**

    str



#### frac_charge(_: floa_ )

#### frac_discharge(_: floa_ )

#### _property_ framework(_: [Composition](pymatgen.core.md#pymatgen.core.composition.Composition_ )
The composition object representing the framework.


#### framework_formula(_: st_ )

#### mAh(_: floa_ )

#### mass_charge(_: floa_ )

#### mass_discharge(_: floa_ )

#### vol_charge(_: floa_ )

#### vol_discharge(_: floa_ )

#### voltage(_: floa_ )

#### _property_ working_ion(_: [Element](pymatgen.core.md#pymatgen.core.periodic_table.Element_ )
Working ion as pymatgen Element object.


#### working_ion_entry(_: [ComputedEntry](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedEntry_ )

#### _property_ x_charge(_: floa_ )
The number of working ions per formula unit of host in the charged state.


#### _property_ x_discharge(_: floa_ )
The number of working ions per formula unit of host in the discharged state.

## pymatgen.apps.battery.conversion_battery module

This module contains the classes to build a ConversionElectrode.


### _class_ ConversionElectrode(voltage_pairs: tuple[AbstractVoltagePair, ...], working_ion_entry: [ComputedEntry](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedEntry), framework_formula: str, initial_comp_formula: str)
Bases: `AbstractElectrode`

Class representing a ConversionElectrode, since it is dataclass
this object can be constructed for the attributes.
However, it is usually easier to construct a ConversionElectrode using one of the classmethod
constructors provided.

Attribute:

    voltage_pairs: The voltage pairs making up the Conversion Electrode.
    working_ion_entry: A single ComputedEntry or PDEntry

    > representing the element that carries charge across the
    > battery, e.g. Li.

    initial_comp_formula: Starting composition for ConversionElectrode represented

        as a string/formula.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### _classmethod_ from_composition_and_entries(comp, entries_in_chemsys, working_ion_symbol='Li', allow_unstable=False)
Convenience constructor to make a ConversionElectrode from a
composition and all entries in a chemical system.


* **Parameters**


    * **comp** – Starting composition for ConversionElectrode, e.g.,
    Composition(“FeF3”)


    * **entries_in_chemsys** – Sequence containing all entries in a
    chemical system. E.g., all Li-Fe-F containing entries.


    * **working_ion_symbol** – Element symbol of working ion. Defaults to Li.


    * **allow_unstable** – If True, allow any composition to be used as the
    starting point of a conversion voltage curve, this is useful
    for comparing with insertion electrodes



#### _classmethod_ from_composition_and_pd(comp, pd, working_ion_symbol='Li', allow_unstable=False)
Convenience constructor to make a ConversionElectrode from a
composition and a phase diagram.


* **Parameters**


    * **comp** – Starting composition for ConversionElectrode, e.g.,
    Composition(“FeF3”)


    * **pd** – A PhaseDiagram of the relevant system (e.g., Li-Fe-F)


    * **working_ion_symbol** – Element symbol of working ion. Defaults to Li.


    * **allow_unstable** – Allow compositions that are unstable



#### get_sub_electrodes(adjacent_only=True)
If this electrode contains multiple voltage steps, then it is possible
to use only a subset of the voltage steps to define other electrodes.
For example, an LiTiO2 electrode might contain three subelectrodes:
[LiTiO2 –> TiO2, LiTiO2 –> Li0.5TiO2, Li0.5TiO2 –> TiO2]
This method can be used to return all the subelectrodes with some
options.


* **Parameters**

    **adjacent_only** – Only return electrodes from compounds that are
    adjacent on the convex hull, i.e. no electrodes returned
    will have multiple voltage steps if this is set true



* **Returns**

    A list of ConversionElectrode objects



#### get_summary_dict(print_subelectrodes=True)
Generate a summary dict.
Populates the summary dict with the basic information from the parent method then populates more information.
Since the parent method calls self.get_summary_dict(print_subelectrodes=True) for the subelectrodes.
The current method will be called from within super().get_summary_dict.


* **Parameters**

    **print_subelectrodes** – Also print data on all the possible
    subelectrodes.



* **Returns**

    A summary of this electrode’s properties in dict format.



#### _property_ initial_comp(_: [Composition](pymatgen.core.md#pymatgen.core.composition.Composition_ )
The pymatgen Composition representation of the initial composition.


#### initial_comp_formula(_: st_ )

#### is_super_electrode(conversion_electrode)
Checks if a particular conversion electrode is a sub electrode of the
current electrode. Starting from a more lithiated state may result in
a subelectrode that is essentially on the same path. For example, a
ConversionElectrode formed by starting from an FePO4 composition would
be a super_electrode of a ConversionElectrode formed from an LiFePO4
composition.


### _class_ ConversionVoltagePair(voltage: float, mAh: float, mass_charge: float, mass_discharge: float, vol_charge: float, vol_discharge: float, frac_charge: float, frac_discharge: float, working_ion_entry: [ComputedEntry](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedEntry), framework_formula: str, rxn: [BalancedReaction](pymatgen.analysis.md#pymatgen.analysis.reaction_calculator.BalancedReaction), entries_charge: Iterable[[ComputedEntry](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedEntry)], entries_discharge: Iterable[[ComputedEntry](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedEntry)])
Bases: `AbstractVoltagePair`

A VoltagePair representing a Conversion Reaction with a defined voltage.
Typically not initialized directly but rather used by ConversionElectrode.


#### rxn()
BalancedReaction for the step


* **Type**

    [BalancedReaction](pymatgen.analysis.md#pymatgen.analysis.reaction_calculator.BalancedReaction)



#### voltage()
Voltage for the step


* **Type**

    float



#### mAh()
Capacity of the step


* **Type**

    float



#### vol_charge()
Volume of charged state


* **Type**

    float



#### vol_discharge()
Volume of discharged state


* **Type**

    float



#### mass_charge()
Mass of charged state


* **Type**

    float



#### mass_discharge()
Mass of discharged state


* **Type**

    float



#### frac_charge()
Fraction of working ion in the charged state


* **Type**

    float



#### frac_discharge()
Fraction of working ion in the discharged state


* **Type**

    float



#### entries_charge()
Entries representing decompositions products
in the charged state. Enumerates the decompositions products at the tieline,
so the number of entries will be one fewer than the dimensions of the phase
diagram


* **Type**

    [[ComputedEntry](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedEntry)]



#### entries_discharge()
Entries representing decompositions products
in the discharged state. Enumerates the decompositions products at the tieline,
so the number of entries will be one fewer than the dimensions of the phase
diagram


* **Type**

    [[ComputedEntry](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedEntry)]



#### working_ion_entry()
Entry of the working ion.


* **Type**

    [ComputedEntry](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedEntry)



#### entries_charge(_: Iterable[[ComputedEntry](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedEntry)_ )

#### entries_discharge(_: Iterable[[ComputedEntry](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedEntry)_ )

#### _classmethod_ from_steps(step1, step2, normalization_els, framework_formula)
Creates a ConversionVoltagePair from two steps in the element profile
from a PD analysis.


* **Parameters**


    * **step1** – Starting step


    * **step2** – Ending step


    * **normalization_els** – Elements to normalize the reaction by. To
    ensure correct capacities.


    * **framework_formula** – Formula of the framework.



#### rxn(_: [BalancedReaction](pymatgen.analysis.md#pymatgen.analysis.reaction_calculator.BalancedReaction_ )
## pymatgen.apps.battery.insertion_battery module

This module is used for analysis of materials with potential application as
intercalation batteries.


### _class_ InsertionElectrode(voltage_pairs: tuple[AbstractVoltagePair, ...], working_ion_entry: [ComputedEntry](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedEntry), framework_formula: str, stable_entries: Iterable[[ComputedEntry](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedEntry)], unstable_entries: Iterable[[ComputedEntry](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedEntry)])
Bases: `AbstractElectrode`

A set of topotactically related compounds, with different amounts of a
single element, e.g. TiO2 and LiTiO2, that can be used to define an
insertion battery electrode.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### as_dict_legacy()
Returns: MSONable dict.


#### _classmethod_ from_dict_legacy(d)

* **Parameters**

    **d** (*dict*) – Dict representation.



* **Returns**

    InsertionElectrode



#### _classmethod_ from_entries(entries: Iterable[[ComputedEntry](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedEntry) | [ComputedStructureEntry](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry)], working_ion_entry: [ComputedEntry](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedEntry) | [ComputedStructureEntry](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry) | [PDEntry](pymatgen.analysis.md#pymatgen.analysis.phase_diagram.PDEntry), strip_structures: bool = False)
Create a new InsertionElectrode.


* **Parameters**


    * **entries** – A list of ComputedEntries, ComputedStructureEntries, or
    subclasses representing the different topotactic states
    of the battery, e.g. TiO2 and LiTiO2.


    * **working_ion_entry** – A single ComputedEntry or PDEntry
    representing the element that carries charge across the
    battery, e.g. Li.


    * **strip_structures** – Since the electrode document only uses volume we can make the
    electrode object significantly leaner by dropping the structure data.
    If this parameter is set to True, the ComputedStructureEntry will be
    replaced with a ComputedEntry and the volume will be stored in
    ComputedEntry.data[‘volume’]. If entries provided are ComputedEntries,
    must set strip_structures=False.



#### _property_ fully_charged_entry()
The most charged entry along the topotactic path.


#### _property_ fully_discharged_entry()
The most discharged entry along the topotactic path.


#### get_all_entries(charge_to_discharge=True)
Return all entries input for the electrode.


* **Parameters**

    **charge_to_discharge** – order from most charge to most discharged state? Defaults to
    True.



* **Returns**

    A list of all entries in the electrode (both stable and unstable),
    ordered by amount of the working ion.



#### get_max_instability(min_voltage=None, max_voltage=None)
The maximum instability along a path for a specific voltage range.


* **Parameters**


    * **min_voltage** – The minimum allowable voltage.


    * **max_voltage** – The maximum allowable voltage.



* **Returns**

    Maximum decomposition energy of all compounds along the insertion
    path (a subset of the path can be chosen by the optional arguments)



#### get_max_muO2(min_voltage=None, max_voltage=None)
Maximum critical oxygen chemical potential along path.


* **Parameters**


    * **min_voltage** – The minimum allowable voltage.


    * **max_voltage** – The maximum allowable voltage.



* **Returns**

    Maximum critical oxygen chemical of all compounds along the
    insertion path (a subset of the path can be chosen by the optional
    arguments).



#### get_min_instability(min_voltage=None, max_voltage=None)
The minimum instability along a path for a specific voltage range.


* **Parameters**


    * **min_voltage** – The minimum allowable voltage.


    * **max_voltage** – The maximum allowable voltage.



* **Returns**

    Minimum decomposition energy of all compounds along the insertion
    path (a subset of the path can be chosen by the optional arguments)



#### get_min_muO2(min_voltage=None, max_voltage=None)
Minimum critical oxygen chemical potential along path.


* **Parameters**


    * **min_voltage** – The minimum allowable voltage for a given step


    * **max_voltage** – The maximum allowable voltage allowable for a given
    step



* **Returns**

    Minimum critical oxygen chemical of all compounds along the
    insertion path (a subset of the path can be chosen by the optional
    arguments).



#### get_stable_entries(charge_to_discharge=True)
Get the stable entries.


* **Parameters**

    **charge_to_discharge** – order from most charge to most discharged
    state? Default to True.



* **Returns**

    A list of stable entries in the electrode, ordered by amount of the
    working ion.



#### get_sub_electrodes(adjacent_only=True, include_myself=True)
If this electrode contains multiple voltage steps, then it is possible
to use only a subset of the voltage steps to define other electrodes.
For example, an LiTiO2 electrode might contain three subelectrodes:
[LiTiO2 –> TiO2, LiTiO2 –> Li0.5TiO2, Li0.5TiO2 –> TiO2]
This method can be used to return all the subelectrodes with some
options.


* **Parameters**


    * **adjacent_only** – Only return electrodes from compounds that are
    adjacent on the convex hull, i.e. no electrodes returned
    will have multiple voltage steps if this is set True.


    * **include_myself** – Include this identical electrode in the list of
    results.



* **Returns**

    A list of InsertionElectrode objects



#### get_summary_dict(print_subelectrodes=True)
Generate a summary dict.
Populates the summary dict with the basic information from the parent method then populates more information.
Since the parent method calls self.get_summary_dict(print_subelectrodes=True) for the subelectrodes.
The current method will be called from within super().get_summary_dict.


* **Parameters**

    **print_subelectrodes** – Also print data on all the possible
    subelectrodes.



* **Returns**

    A summary of this electrode’s properties in dict format.



#### get_unstable_entries(charge_to_discharge=True)
Returns the unstable entries for the electrode.


* **Parameters**

    **charge_to_discharge** – Order from most charge to most discharged
    state? Defaults to True.



* **Returns**

    A list of unstable entries in the electrode, ordered by amount of
    the working ion.



#### stable_entries(_: Iterable[[ComputedEntry](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedEntry)_ )

#### unstable_entries(_: Iterable[[ComputedEntry](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedEntry)_ )

### _class_ InsertionVoltagePair(voltage: float, mAh: float, mass_charge: float, mass_discharge: float, vol_charge: float, vol_discharge: float, frac_charge: float, frac_discharge: float, working_ion_entry: [ComputedEntry](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedEntry), framework_formula: str, entry_charge: [ComputedEntry](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedEntry), entry_discharge: [ComputedEntry](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedEntry))
Bases: `AbstractVoltagePair`

Defines an Insertion Voltage Pair.


#### entry_charge(_: [ComputedEntry](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedEntry_ )

#### entry_discharge(_: [ComputedEntry](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedEntry_ )

#### _classmethod_ from_entries(entry1, entry2, working_ion_entry)

* **Parameters**


    * **entry1** – Entry corresponding to one of the entries in the voltage step.


    * **entry2** – Entry corresponding to the other entry in the voltage step.


    * **working_ion_entry** – A single ComputedEntry or PDEntry representing
    the element that carries charge across the battery, e.g. Li.


## pymatgen.apps.battery.plotter module

This module provides plotting capabilities for battery related applications.


### _class_ VoltageProfilePlotter(xaxis='capacity', hide_negative=False)
Bases: `object`

A plotter to make voltage profile plots for batteries.


* **Parameters**


    * **xaxis** – The quantity to use as the xaxis. Can be either


    * **capacity_grav** (*-*) – the graviometric capcity


    * **capacity_vol** (*-*) – the volumetric capacity


    * **x_form** (*-*) – the number of working ions per formula unit of the host


    * **frac_x** (*-*) – the atomic fraction of the working ion


    * **hide_negative** – If True only plot the voltage steps above zero.



#### _choose_best_x_label(formula, wion_symbol)

#### add_electrode(electrode, label=None)
Add an electrode to the plot.


* **Parameters**


    * **electrode** – An electrode. All electrodes satisfying the
    AbstractElectrode interface should work.


    * **label** – A label for the electrode. If None, defaults to a counting
    system, i.e. ‘Electrode 1’, ‘Electrode 2’, …



#### get_plot(width=8, height=8, term_zero=True, ax: Axes | None = None)
Returns a plot object.


* **Parameters**


    * **width** – Width of the plot. Defaults to 8 in.


    * **height** – Height of the plot. Defaults to 6 in.


    * **term_zero** – If True append zero voltage point at the end


    * **ax** (*plt.Axes*) – matplotlib axes object. Defaults to None.



* **Returns**

    matplotlib axes object.



* **Return type**

    plt.Axes



#### get_plot_data(electrode, term_zero=True)

* **Parameters**


    * **electrode** – Electrode object


    * **term_zero** – If True append zero voltage point at the end.



* **Returns**

    Plot data in x, y.



#### get_plotly_figure(width=800, height=600, font_dict=None, term_zero=True, \*\*kwargs)
Return plotly Figure object.


* **Parameters**


    * **width** – Width of the plot. Defaults to 800 px.


    * **height** – Height of the plot. Defaults to 600 px.


    * **font_dict** – define the font. Defaults to {“family”: “Arial”, “size”: 24, “color”: “#000000”}


    * **term_zero** – If True append zero voltage point at the end


    * **\*\*kwargs** – passed to plotly.graph_objects.Layout



#### save(filename, image_format='eps', width=8, height=6)
Save the plot to an image file.


* **Parameters**


    * **filename** – Filename to save to.


    * **image_format** – Format to save to. Defaults to eps.


    * **width** – Width of the plot. Defaults to 8 in.


    * **height** – Height of the plot. Defaults to 6 in.



#### show(width=8, height=6)
Show the voltage profile plot.


* **Parameters**


    * **width** – Width of the plot. Defaults to 8 in.


    * **height** – Height of the plot. Defaults to 6 in.