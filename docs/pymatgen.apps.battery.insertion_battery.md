---
layout: default
title: pymatgen.apps.battery.insertion_battery.md
nav_exclude: true
---

# pymatgen.apps.battery.insertion_battery module

This module is used for analysis of materials with potential application as
intercalation batteries.


### _class_ pymatgen.apps.battery.insertion_battery.InsertionElectrode(voltage_pairs: tuple[[pymatgen.apps.battery.battery_abc.AbstractVoltagePair](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractVoltagePair), ...], working_ion_entry: [ComputedEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry), framework_formula: str, stable_entries: Iterable[[ComputedEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry)], unstable_entries: Iterable[[ComputedEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry)])
Bases: [`AbstractElectrode`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractElectrode)

A set of topotactically related compounds, with different amounts of a
single element, e.g. TiO2 and LiTiO2, that can be used to define an
insertion battery electrode.


#### as_dict_legacy()
Returns: MSONable dict.


#### _classmethod_ from_dict_legacy(d)

* **Parameters**

    **d** (*dict*) – Dict representation.



* **Returns**

    InsertionElectrode



#### _classmethod_ from_entries(entries: Iterable[[ComputedEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry) | [ComputedStructureEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry)], working_ion_entry: [ComputedEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry) | [ComputedStructureEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry) | [PDEntry](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PDEntry), strip_structures: bool = False)
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



#### stable_entries(_: Iterable[[ComputedEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry)_ )

#### unstable_entries(_: Iterable[[ComputedEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry)_ )

### _class_ pymatgen.apps.battery.insertion_battery.InsertionVoltagePair(voltage: float, mAh: float, mass_charge: float, mass_discharge: float, vol_charge: float, vol_discharge: float, frac_charge: float, frac_discharge: float, working_ion_entry: [ComputedEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry), framework_formula: str, entry_charge: [ComputedEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry), entry_discharge: [ComputedEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry))
Bases: [`AbstractVoltagePair`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractVoltagePair)

Defines an Insertion Voltage Pair.


#### entry_charge(_: [ComputedEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry_ )

#### entry_discharge(_: [ComputedEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry_ )

#### _classmethod_ from_entries(entry1, entry2, working_ion_entry)

* **Parameters**


    * **entry1** – Entry corresponding to one of the entries in the voltage step.


    * **entry2** – Entry corresponding to the other entry in the voltage step.


    * **working_ion_entry** – A single ComputedEntry or PDEntry representing
    the element that carries charge across the battery, e.g. Li.