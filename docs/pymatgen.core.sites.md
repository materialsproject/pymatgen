---
layout: default
title: pymatgen.core.sites.md
nav_exclude: true
---

# pymatgen.core.sites module

This module defines classes representing non-periodic and periodic sites.


### _class_ pymatgen.core.sites.PeriodicSite(species: SpeciesLike | CompositionLike, coords: ArrayLike, lattice: [Lattice](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice), to_unit_cell: bool = False, coords_are_cartesian: bool = False, properties: dict | None = None, label: str | None = None, skip_checks: bool = False)
Bases: `Site`, `MSONable`

Extension of generic Site object to periodic systems.
PeriodicSite includes a lattice system.

Create a periodic site.


* **Parameters**


    * **species** – Species on the site. Can be:
    i.  A Composition-type object (preferred)
    ii. An  element / species specified either as a string

    > symbols, e.g. “Li”, “Fe2+”, “P” or atomic numbers,
    > e.g., 3, 56, or actual Element or Species objects.

    iii.Dict of elements/species and occupancies, e.g.,

        {“Fe” : 0.5, “Mn”:0.5}. This allows the setup of
        disordered structures.



    * **coords** – Cartesian coordinates of site.


    * **lattice** – Lattice associated with the site.


    * **to_unit_cell** – Translates fractional coordinate to the
    basic unit cell, i.e. all fractional coordinates satisfy 0
    <= a < 1. Defaults to False.


    * **coords_are_cartesian** – Set to True if you are providing
    Cartesian coordinates. Defaults to False.


    * **properties** – Properties associated with the site as a dict, e.g.
    {“magmom”: 5}. Defaults to None.


    * **label** – Label for the site. Defaults to None.


    * **skip_checks** – Whether to ignore all the usual checks and just
    create the site. Use this if the PeriodicSite is created in a
    controlled manner and speed is desired.



#### _property_ a(_: floa_ )
Fractional a coordinate.


#### as_dict(verbosity: int = 0)
JSON-serializable dict representation of PeriodicSite.


* **Parameters**

    **verbosity** (*int*) – Verbosity level. Default of 0 only includes the matrix
    representation. Set to 1 for more details such as Cartesian coordinates, etc.



#### _property_ b(_: floa_ )
Fractional b coordinate.


#### _property_ c(_: floa_ )
Fractional c coordinate.


#### _property_ coords(_: ndarra_ )
Cartesian coordinates.


#### distance(other: PeriodicSite, jimage: ArrayLike | None = None)
Get distance between two sites assuming periodic boundary conditions.


* **Parameters**


    * **other** (*PeriodicSite*) – Other site to get distance from.


    * **jimage** (*3x1 array*) – Specific periodic image in terms of lattice
    translations, e.g., [1,0,0] implies to take periodic image
    that is one a-lattice vector away. If jimage is None,
    the image that is nearest to the site is found.



* **Returns**

    Distance between the two sites



* **Return type**

    distance (float)



#### distance_and_image(other: PeriodicSite, jimage: ArrayLike | None = None)
Gets distance and instance between two sites assuming periodic boundary
conditions. If the index jimage of two sites atom j is not specified it
selects the j image nearest to the i atom and returns the distance and
jimage indices in terms of lattice vector translations. If the index
jimage of atom j is specified it returns the distance between the ith
atom and the specified jimage atom, the given jimage is also returned.


* **Parameters**


    * **other** (*PeriodicSite*) – Other site to get distance from.


    * **jimage** (*3x1 array*) – Specific periodic image in terms of lattice
    translations, e.g., [1,0,0] implies to take periodic image
    that is one a-lattice vector away. If jimage is None,
    the image that is nearest to the site is found.



* **Returns**

    distance and periodic lattice translations
    of the other site for which the distance applies.



* **Return type**

    (distance, jimage)



#### distance_and_image_from_frac_coords(fcoords: ArrayLike, jimage: ArrayLike | None = None)
Gets distance between site and a fractional coordinate assuming
periodic boundary conditions. If the index jimage of two sites atom j
is not specified it selects the j image nearest to the i atom and
returns the distance and jimage indices in terms of lattice vector
translations. If the index jimage of atom j is specified it returns the
distance between the i atom and the specified jimage atom, the given
jimage is also returned.


* **Parameters**


    * **fcoords** (*3x1 array*) – fcoords to get distance from.


    * **jimage** (*3x1 array*) – Specific periodic image in terms of
    lattice translations, e.g., [1,0,0] implies to take periodic
    image that is one a-lattice vector away. If jimage is None,
    the image that is nearest to the site is found.



* **Returns**

    distance and periodic lattice translations
    of the other site for which the distance applies.



* **Return type**

    (distance, jimage)



#### _property_ frac_coords(_: ndarra_ )
Fractional coordinates.


#### _classmethod_ from_dict(dct, lattice=None)
Create PeriodicSite from dict representation.


* **Parameters**


    * **dct** (*dict*) – dict representation of PeriodicSite


    * **lattice** – Optional lattice to override lattice specified in d.
    Useful for ensuring all sites in a structure share the same
    lattice.



* **Returns**

    PeriodicSite



#### is_periodic_image(other: PeriodicSite, tolerance: float = 1e-08, check_lattice: bool = True)
Returns True if sites are periodic images of each other.


* **Parameters**


    * **other** (*PeriodicSite*) – Other site


    * **tolerance** (*float*) – Tolerance to compare fractional coordinates


    * **check_lattice** (*bool*) – Whether to check if the two sites have the
    same lattice.



* **Returns**

    True if sites are periodic images of each other.



* **Return type**

    bool



#### _property_ lattice(_: [Lattice](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice_ )
Lattice associated with PeriodicSite.


#### to_unit_cell(in_place=False)
Move frac coords to within the unit cell.


#### _property_ x(_: floa_ )
Cartesian x coordinate.


#### _property_ y(_: floa_ )
Cartesian y coordinate.


#### _property_ z(_: floa_ )
Cartesian z coordinate.


### _class_ pymatgen.core.sites.Site(species: SpeciesLike | CompositionLike, coords: ArrayLike, properties: dict | None = None, label: str | None = None, skip_checks: bool = False)
Bases: `Hashable`, `MSONable`

A generalized *non-periodic* site. This is essentially a composition
at a point in space, with some optional properties associated with it. A
Composition is used to represent the atoms and occupancy, which allows for
disordered site representation. Coords are given in standard Cartesian
coordinates.

Creates a non-periodic Site.


* **Parameters**


    * **species** – Species on the site. Can be:
    i.  A Composition-type object (preferred)
    ii. An  element / species specified either as a string

    > symbols, e.g. “Li”, “Fe2+”, “P” or atomic numbers,
    > e.g., 3, 56, or actual Element or Species objects.

    iii.Dict of elements/species and occupancies, e.g.,

        {“Fe” : 0.5, “Mn”:0.5}. This allows the setup of
        disordered structures.



    * **coords** – Cartesian coordinates of site.


    * **properties** – Properties associated with the site as a dict, e.g.
    {“magmom”: 5}. Defaults to None.


    * **label** – Label for the site. Defaults to None.


    * **skip_checks** – Whether to ignore all the usual checks and just
    create the site. Use this if the Site is created in a controlled
    manner and speed is desired.



#### as_dict()
JSON-serializable dict representation for Site.


#### distance(other)
Get distance between two sites.


* **Parameters**

    **other** – Other site.



* **Returns**

    distance



* **Return type**

    float



#### distance_from_point(pt)
Returns distance between the site and a point in space.


* **Parameters**

    **pt** – Cartesian coordinates of point.



* **Returns**

    distance



* **Return type**

    float



#### _classmethod_ from_dict(dct: dict)
Create Site from dict representation.


#### _property_ is_ordered(_: boo_ )
True if site is an ordered site, i.e., with a single species with
occupancy 1.


#### _property_ label(_: st_ )
Site label.


#### position_atol(_ = 1e-0_ )

#### _property_ specie(_: [Element](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element) | [Species](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species) | [DummySpecies](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.DummySpecies_ )
The Species/Element at the site. Only works for ordered sites. Otherwise
an AttributeError is raised. Use this property sparingly.  Robust
design should make use of the property species instead. Note that the
singular of species is also species. So the choice of this variable
name is governed by programmatic concerns as opposed to grammar.


* **Raises**

    **AttributeError if Site is not ordered.** –



#### _property_ species(_: [Composition](pymatgen.core.composition.md#pymatgen.core.composition.Composition_ )
The species on the site as a composition, e.g., Fe0.5Mn0.5.


* **Type**

    return



#### _property_ species_string(_: st_ )
String representation of species on the site.


#### _property_ x(_: floa_ )
Cartesian x coordinate.


#### _property_ y(_: floa_ )
Cartesian y coordinate.


#### _property_ z(_: floa_ )
Cartesian z coordinate.