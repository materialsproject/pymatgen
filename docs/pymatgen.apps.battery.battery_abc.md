---
layout: default
title: pymatgen.apps.battery.battery_abc.md
nav_exclude: true
---

# pymatgen.apps.battery.battery_abc module

This module defines the abstract base classes for battery-related classes.
Regardless of the kind of electrode, conversion or insertion, there are many
common definitions and properties, e.g., average voltage, capacity, etc. which
can be defined in a general way. The Abc for battery classes implements some of
these common definitions to allow sharing of common logic between them.


### _class_ pymatgen.apps.battery.battery_abc.AbstractElectrode(voltage_pairs: tuple[AbstractVoltagePair, ...], working_ion_entry: [ComputedEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry), framework_formula: str)
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

    [ComputedEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry)



#### framework_formula()
The compositions of one formula unit of the host material


* **Type**

    str



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


#### working_ion_entry(_: [ComputedEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry_ )

#### _property_ x_charge(_: floa_ )
The number of working ions per formula unit of host in the charged state.


#### _property_ x_discharge(_: floa_ )
The number of working ions per formula unit of host in the discharged state.


### _class_ pymatgen.apps.battery.battery_abc.AbstractVoltagePair(voltage: float, mAh: float, mass_charge: float, mass_discharge: float, vol_charge: float, vol_discharge: float, frac_charge: float, frac_discharge: float, working_ion_entry: [ComputedEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry), framework_formula: str)
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

    [ComputedEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry)



#### framework_formula()
The compositions of one formula unit of the host material


* **Type**

    str



#### frac_charge(_: floa_ )

#### frac_discharge(_: floa_ )

#### _property_ framework(_: [Composition](pymatgen.core.composition.md#pymatgen.core.composition.Composition_ )
The composition object representing the framework.


#### framework_formula(_: st_ )

#### mAh(_: floa_ )

#### mass_charge(_: floa_ )

#### mass_discharge(_: floa_ )

#### vol_charge(_: floa_ )

#### vol_discharge(_: floa_ )

#### voltage(_: floa_ )

#### _property_ working_ion(_: [Element](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element_ )
Working ion as pymatgen Element object.


#### working_ion_entry(_: [ComputedEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry_ )

#### _property_ x_charge(_: floa_ )
The number of working ions per formula unit of host in the charged state.


#### _property_ x_discharge(_: floa_ )
The number of working ions per formula unit of host in the discharged state.