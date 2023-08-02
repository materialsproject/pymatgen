---
layout: default
title: pymatgen.apps.battery.analyzer.md
nav_exclude: true
---

# pymatgen.apps.battery.analyzer module

Analysis classes for batteries.


### _class_ pymatgen.apps.battery.analyzer.BatteryAnalyzer(struc_oxid, working_ion='Li', oxi_override=None)
Bases: `object`

A suite of methods for starting with an oxidized structure and determining its potential as a battery.

Pass in a structure for analysis.


* **Parameters**


    * **struc_oxid** – a Structure object; oxidation states *must* be assigned for this structure; disordered
    structures should be OK


    * **working_ion** – a String symbol or Element for the working ion.


    * **oxi_override** – a dict of String element symbol, Integer oxidation state pairs.
    by default, H, C, N, O, F, S, Cl, Se, Br, Te, I are considered anions.



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



### pymatgen.apps.battery.analyzer.is_redox_active_intercalation(element)
True if element is redox active and interesting for intercalation materials.


* **Parameters**

    **element** – Element object