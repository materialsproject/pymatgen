---
layout: default
title: pymatgen.core.bonds.md
nav_exclude: true
---

# pymatgen.core.bonds module

This class implements definitions for various kinds of bonds. Typically used in
Molecule analysis.


### _class_ pymatgen.core.bonds.CovalentBond(site1: [Site](pymatgen.core.sites.md#pymatgen.core.sites.Site), site2: [Site](pymatgen.core.sites.md#pymatgen.core.sites.Site))
Bases: `object`

Defines a covalent bond between two sites.

Initializes a covalent bond between two sites.


* **Parameters**


    * **site1** ([*Site*](pymatgen.core.sites.md#pymatgen.core.sites.Site)) – First site.


    * **site2** ([*Site*](pymatgen.core.sites.md#pymatgen.core.sites.Site)) – Second site.



#### get_bond_order(tol: float = 0.2, default_bl: float | None = None)
The bond order according the distance between the two sites
:param tol: Relative tolerance to test.

> (1 + tol) \* the longest bond distance is considered
> to be the threshold length for a bond to exist.
> (1 - tol) \* the shortest bond distance is considered
> to be the shortest possible bond length
> Defaults to 0.2.


* **Parameters**

    **default_bl** – If a particular type of bond does not exist,
    use this bond length as a default value
    (bond order = 1). If None, a ValueError will be thrown.



* **Returns**

    Float value of bond order. For example, for C-C bond in
    benzene, return 1.7.



#### _static_ is_bonded(site1, site2, tol: float = 0.2, bond_order: float | None = None, default_bl: float | None = None)
Test if two sites are bonded, up to a certain limit.


* **Parameters**


    * **site1** ([*Site*](pymatgen.core.sites.md#pymatgen.core.sites.Site)) – First site


    * **site2** ([*Site*](pymatgen.core.sites.md#pymatgen.core.sites.Site)) – Second site


    * **tol** (*float*) – Relative tolerance to test. Basically, the code
    checks if the distance between the sites is less than (1 +
    tol) \* typical bond distances. Defaults to 0.2, i.e.,
    20% longer.


    * **bond_order** – Bond order to test. If None, the code simply checks
    against all possible bond data. Defaults to None.


    * **default_bl** – If a particular type of bond does not exist, use this
    bond length. If None, a ValueError will be thrown.



* **Returns**

    Boolean indicating whether two sites are bonded.



#### _property_ length(_: floa_ )
Length of the bond.


### pymatgen.core.bonds.get_bond_length(sp1: SpeciesLike, sp2: SpeciesLike, bond_order: float = 1)
Get the bond length between two species.


* **Parameters**


    * **sp1** ([*Species*](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species)) – First specie.


    * **sp2** ([*Species*](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species)) – Second specie.


    * **bond_order** – For species with different possible bond orders,
    this allows one to obtain the bond length for a particular bond
    order. For example, to get the C=C bond length instead of the
    C-C bond length, this should be set to 2. Defaults to 1.



* **Returns**

    Bond length in Angstrom. If no data is available, the sum of the atomic
    radius is used.



### pymatgen.core.bonds.get_bond_order(sp1, sp2, dist: float, tol: float = 0.2, default_bl: float | None = None)
Calculate the bond order given the distance of 2 species.


* **Parameters**


    * **sp1** ([*Species*](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species)) – First specie.


    * **sp2** ([*Species*](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species)) – Second specie.


    * **dist** – Their distance in angstrom


    * **tol** (*float*) – Relative tolerance to test. Basically, the code
    checks if the distance between the sites is larger than
    (1 + tol) \* the longest bond distance or smaller than
    (1 - tol) \* the shortest bond distance to determine if
    they are bonded or the distance is too short.
    Defaults to 0.2.


    * **default_bl** – If a particular type of bond does not exist, use this
    bond length (bond order = 1). If None, a ValueError will be thrown.



* **Returns**

    Float value of bond order. For example, for C-C bond in benzene,
    return 1.7.



### pymatgen.core.bonds.obtain_all_bond_lengths(sp1, sp2, default_bl: float | None = None)
Obtain bond lengths for all bond orders from bond length database.


* **Parameters**


    * **sp1** ([*Species*](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species)) – First specie.


    * **sp2** ([*Species*](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species)) – Second specie.


    * **default_bl** – If a particular type of bond does not exist, use this
    bond length as a default value (bond order = 1).
    If None, a ValueError will be thrown.



* **Returns**

    A dict mapping bond order to bond length in angstrom