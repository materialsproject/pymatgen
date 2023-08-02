---
layout: default
title: pymatgen.apps.battery.conversion_battery.md
nav_exclude: true
---

# pymatgen.apps.battery.conversion_battery module

This module contains the classes to build a ConversionElectrode.


### _class_ pymatgen.apps.battery.conversion_battery.ConversionElectrode(voltage_pairs: tuple[[AbstractVoltagePair](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractVoltagePair), ...], working_ion_entry: [ComputedEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry), framework_formula: str, initial_comp_formula: str)
Bases: [`AbstractElectrode`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractElectrode)

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



#### _property_ initial_comp(_: [Composition](pymatgen.core.composition.md#pymatgen.core.composition.Composition_ )
The pymatgen Composition representation of the initial composition.


#### initial_comp_formula(_: st_ )

#### is_super_electrode(conversion_electrode)
Checks if a particular conversion electrode is a sub electrode of the
current electrode. Starting from a more lithiated state may result in
a subelectrode that is essentially on the same path. For example, a
ConversionElectrode formed by starting from an FePO4 composition would
be a super_electrode of a ConversionElectrode formed from an LiFePO4
composition.


### _class_ pymatgen.apps.battery.conversion_battery.ConversionVoltagePair(voltage: float, mAh: float, mass_charge: float, mass_discharge: float, vol_charge: float, vol_discharge: float, frac_charge: float, frac_discharge: float, working_ion_entry: [ComputedEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry), framework_formula: str, rxn: [BalancedReaction](pymatgen.analysis.reaction_calculator.md#pymatgen.analysis.reaction_calculator.BalancedReaction), entries_charge: Iterable[[ComputedEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry)], entries_discharge: Iterable[[ComputedEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry)])
Bases: [`AbstractVoltagePair`](pymatgen.apps.battery.battery_abc.md#pymatgen.apps.battery.battery_abc.AbstractVoltagePair)

A VoltagePair representing a Conversion Reaction with a defined voltage.
Typically not initialized directly but rather used by ConversionElectrode.


#### rxn()
BalancedReaction for the step


* **Type**

    [BalancedReaction](pymatgen.analysis.reaction_calculator.md#pymatgen.analysis.reaction_calculator.BalancedReaction)



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

    [[ComputedEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry)]



#### entries_discharge()
Entries representing decompositions products
in the discharged state. Enumerates the decompositions products at the tieline,
so the number of entries will be one fewer than the dimensions of the phase
diagram


* **Type**

    [[ComputedEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry)]



#### working_ion_entry()
Entry of the working ion.


* **Type**

    [ComputedEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry)



#### entries_charge(_: Iterable[[ComputedEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry)_ )

#### entries_discharge(_: Iterable[[ComputedEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry)_ )

#### _classmethod_ from_steps(step1, step2, normalization_els, framework_formula)
Creates a ConversionVoltagePair from two steps in the element profile
from a PD analysis.


* **Parameters**


    * **step1** – Starting step


    * **step2** – Ending step


    * **normalization_els** – Elements to normalize the reaction by. To
    ensure correct capacities.


    * **framework_formula** – Formula of the framework.



#### rxn(_: [BalancedReaction](pymatgen.analysis.reaction_calculator.md#pymatgen.analysis.reaction_calculator.BalancedReaction_ )