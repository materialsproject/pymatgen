---
layout: default
title: pymatgen.core.interface.md
nav_exclude: true
---

# pymatgen.core.interface module

This module provides classes to store, generate, and manipulate material interfaces.


### _class_ pymatgen.core.interface.Interface(lattice, species, coords, site_properties, validate_proximity=False, to_unit_cell=False, coords_are_cartesian=False, in_plane_offset: tuple[float, float] = (0, 0), gap: float = 0, vacuum_over_film: float = 0, interface_properties: dict | None = None)
Bases: [`Structure`](pymatgen.core.structure.md#pymatgen.core.structure.Structure)

This class stores data for defining an interface between two structures.
It is a subclass of pymatgen.core.structure.Structure.

Makes an interface structure, a structure object with additional information
and methods pertaining to interfaces.


* **Parameters**


    * **lattice** (*Lattice/3x3 array*) – The lattice, either as a
    [`pymatgen.core.lattice.Lattice`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice) or
    simply as any 2D array. Each row should correspond to a lattice
    vector. E.g., [[10,0,0], [20,10,0], [0,0,30]] specifies a
    lattice with lattice vectors [10,0,0], [20,10,0] and [0,0,30].


    * **species** (*[*[*Species*](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species)*]*) – Sequence of species on each site. Can take in
    flexible input, including:


        1. A sequence of element / species specified either as string
    symbols, e.g. [“Li”, “Fe2+”, “P”, …] or atomic numbers,
    e.g., (3, 56, …) or actual Element or Species objects.


        2. List of dict of elements/species and occupancies, e.g.,
    [{“Fe” : 0.5, “Mn”:0.5}, …]. This allows the setup of
    disordered structures.



    * **coords** (*Nx3 array*) – list of fractional/cartesian coordinates of
    each species.


    * **validate_proximity** (*bool*) – Whether to check if there are sites
    that are less than 0.01 Ang apart. Defaults to False.


    * **to_unit_cell** (*bool*) – Whether to translate sites into the unit cell. Defaults to False.


    * **coords_are_cartesian** (*bool*) – Set to True if you are providing
    coordinates in Cartesian coordinates. Defaults to False.


    * **site_properties** (*dict*) – Properties associated with the sites as a
    dict of sequences, e.g., {“magmom”:[5,5,5,5]}. The sequences
    have to be the same length as the atomic species and
    fractional_coords. Defaults to None for no properties.


    * **in_plane_offset** – fractional shift in plane for the film with respect
    to the substrate


    * **gap** – gap between substrate and film in Angstroms; zero corresponds to
    the original distance between substrate and film sites


    * **vacuum_over_film** – vacuum space above the film in Angstroms. Defaults to 0.


    * **interface_properties** – properties associated with the interface. Defaults to None.



#### as_dict()

* **Returns**

    MSONable dict



#### copy()

* **Returns**

    A copy of the Interface.



* **Return type**

    Interface



#### _property_ film(_: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure_ )
A pymatgen Structure for just the film.


#### _property_ film_indices(_: list[int_ )
Site indices of the film sites.


#### _property_ film_layers(_: in_ )
Number of layers of the minimum element in the film composition.


#### _property_ film_sites(_: list[[pymatgen.core.sites.Site](pymatgen.core.sites.md#pymatgen.core.sites.Site)_ )
Return the film sites of the interface.


#### _property_ film_termination(_: st_ )
Label for the film termination chemistry.


#### _classmethod_ from_dict(d)

* **Parameters**

    **d** – dict



* **Returns**

    Creates slab from dict.



#### _classmethod_ from_slabs(substrate_slab: [Slab](pymatgen.core.surface.md#pymatgen.core.surface.Slab), film_slab: [Slab](pymatgen.core.surface.md#pymatgen.core.surface.Slab), in_plane_offset: tuple[float, float] = (0, 0), gap: float = 1.6, vacuum_over_film: float = 0, interface_properties: dict | None = None, center_slab: bool = True)
Makes an interface structure by merging a substrate and film slabs
The film a- and b-vectors will be forced to be the substrate slab’s
a- and b-vectors.

For now, it’s suggested to use a factory method that will ensure the
appropriate interface structure is already met.


* **Parameters**


    * **substrate_slab** ([*Slab*](pymatgen.core.surface.md#pymatgen.core.surface.Slab)) – slab for the substrate


    * **film_slab** ([*Slab*](pymatgen.core.surface.md#pymatgen.core.surface.Slab)) – slab for the film


    * **in_plane_offset** (*tuple*) – fractional shift in plane for the film with respect to the substrate.
    For example, (0.5, 0.5) will shift the film by half the substrate’s a- and b-vectors.
    Defaults to (0, 0).


    * **gap** (*float*) – gap between substrate and film in Angstroms. Defaults to 1.6.


    * **vacuum_over_film** (*float*) – vacuum space above the film in Angstroms. Defaults to 0.


    * **interface_properties** (*dict*) – misc properties to assign to the interface. Defaults to None.


    * **center_slab** (*bool*) – center the slab. Defaults to True.



#### _property_ gap(_: floa_ )
The gap in Cartesian units between the film and the substrate.


#### get_shifts_based_on_adsorbate_sites(tolerance: float = 0.1)
Computes possible in-plane shifts based on an adsorbate site  algorithm.


* **Parameters**

    **tolerance** – tolerance for “uniqueness” for shifts in Cartesian unit
    This is usually Angstroms.



#### get_sorted_structure(key=None, reverse=False)
Get a sorted structure for the interface. The parameters have the same
meaning as in list.sort. By default, sites are sorted by the
electronegativity of the species.


* **Parameters**


    * **key** – Specifies a function of one argument that is used to extract
    a comparison key from each list element: key=str.lower. The
    default value is None (compare the elements directly).


    * **reverse** (*bool*) – If set to True, then the list elements are sorted
    as if each comparison were reversed.



#### _property_ in_plane_offset(_: ndarra_ )
The shift between the film and substrate in fractional
coordinates.


#### _property_ substrate(_: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure_ )
A pymatgen Structure for just the substrate.


#### _property_ substrate_indices(_: list[int_ )
Site indices for the substrate atoms.


#### _property_ substrate_layers(_: in_ )
Number of layers of the minimum element in the substrate composition.


#### _property_ substrate_sites(_: list[[pymatgen.core.sites.Site](pymatgen.core.sites.md#pymatgen.core.sites.Site)_ )
The site objects in the substrate.


#### _property_ substrate_termination(_: st_ )
Label for the substrate termination chemistry.


#### _property_ vacuum_over_film(_: floa_ )
The vacuum space over the film in Cartesian units.


### pymatgen.core.interface.count_layers(struct: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), el=None)
Counts the number of ‘layers’ along the c-axis.


### pymatgen.core.interface.label_termination(slab: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Labels the slab surface termination.