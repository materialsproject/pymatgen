---
layout: default
title: pymatgen.analysis.magnetism.analyzer.md
nav_exclude: true
---

# pymatgen.analysis.magnetism.analyzer module

This module provides some useful functions for dealing with magnetic Structures
(e.g. Structures with associated magmom tags).


### _class_ pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), overwrite_magmom_mode: OverwriteMagmomMode | str = 'none', round_magmoms: bool = False, detect_valences: bool = False, make_primitive: bool = True, default_magmoms: dict | None = None, set_net_positive: bool = True, threshold: float = 0, threshold_nonmag: float = 0.1)
Bases: `object`

A class which provides a few helpful methods to analyze
collinear magnetic structures.

If magnetic moments are not defined, moments will be
taken either from default_magmoms.yaml (similar to the
default magmoms in MPRelaxSet, with a few extra definitions)
or from a specie:magmom dict provided by the default_magmoms
kwarg.

Input magmoms can be replaced using the ‘overwrite_magmom_mode’
kwarg. This can be:
\* “none” to do nothing,
\* “respect_sign” which will overwrite existing magmoms with

> those from default_magmoms but will keep sites with positive magmoms
> positive, negative magmoms negative and zero magmoms zero,


* “respect_zeros”, which will give a ferromagnetic structure
(all positive magmoms from default_magmoms) but still keep sites with
zero magmoms as zero,


* “replace_all” which will try to guess initial magmoms for
all sites in the structure irrespective of input structure
(this is most suitable for an initial DFT calculation),


* “replace_all_if_undefined” is the same as “replace_all” but only if
no magmoms are defined in input structure, otherwise it will respect
existing magmoms.


* “normalize” will normalize magmoms to unity, but will respect sign
(used for comparing orderings), magmoms < theshhold will be set to zero


* **Parameters**


    * **structure** – input Structure object


    * **overwrite_magmom_mode** – “respect_sign”, “respect_zeros”, “replace_all”,
    “replace_all_if_undefined”, “normalize” (default “none”)


    * **round_magmoms** – will round input magmoms to
    specified number of decimal places if integer is supplied, if set
    to a float will try and group magmoms together using a kernel density
    estimator of provided width, and extracting peaks of the estimator
    detect_valences: if True, will attempt to assign valences
    to input structure


    * **make_primitive** – if True, will transform to primitive
    magnetic cell


    * **default_magmoms** – (optional) dict specifying default magmoms


    * **set_net_positive** – if True, will change sign of magnetic
    moments such that the net magnetization is positive. Argument will be
    ignored if mode “respect_sign” is used.


    * **threshold** – number (in Bohr magnetons) below which magmoms
    will be rounded to zero


    * **threshold_nonmag** – number (in Bohr magneton)
    below which nonmagnetic ions (with no magmom specified
    in default_magmoms) will be rounded to zero



#### get_exchange_group_info(symprec: float = 0.01, angle_tolerance: float = 5)
Returns the information on the symmetry of the Hamiltonian
describing the exchange energy of the system, taking into
account relative direction of magnetic moments but not their
absolute direction.

This is not strictly accurate (e.g. some/many atoms will
have zero magnetic moments), but defining symmetry this
way is a useful way of keeping track of distinct magnetic
orderings within pymatgen.


* **Parameters**


    * **symprec** – same as SpacegroupAnalyzer (Default value = 1e-2)


    * **angle_tolerance** – same as SpacegroupAnalyzer (Default value = 5)



* **Returns**

    spacegroup_symbol, international_number



#### get_ferromagnetic_structure(make_primitive: bool = True)
Returns a Structure with all magnetic moments positive
or zero.


* **Parameters**

    **make_primitive** – Whether to make structure primitive after
    making all magnetic moments positive (Default value = True)



* **Returns**

    Structure



#### get_nonmagnetic_structure(make_primitive: bool = True)
Returns a Structure without magnetic moments defined.


* **Parameters**

    **make_primitive** – Whether to make structure primitive after
    removing magnetic information (Default value = True)



* **Returns**

    Structure



#### get_structure_with_only_magnetic_atoms(make_primitive: bool = True)
Returns a Structure with only magnetic atoms present.


* **Parameters**

    **make_primitive** – Whether to make structure primitive after
    removing non-magnetic atoms (Default value = True)


Returns: Structure


#### get_structure_with_spin()
Returns a Structure with species decorated with spin values instead
of using magmom site properties.


#### _property_ is_magnetic(_: boo_ )
Convenience property, returns True if any non-zero magmoms present.


#### _property_ magmoms(_: ndarra_ )
Convenience property, returns magmoms as a numpy array.


#### _property_ magnetic_species_and_magmoms(_: dict[str, Any_ )
Returns a dict of magnetic species and the magnitude of
their associated magmoms. Will return a list if there are
multiple magmoms per species.

Returns: dict of magnetic species and magmoms


#### matches_ordering(other: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Compares the magnetic orderings of one structure with another.


* **Parameters**

    **other** – Structure to compare


Returns: True or False


#### _property_ number_of_magnetic_sites(_: in_ )
Number of magnetic sites present in structure.


#### number_of_unique_magnetic_sites(symprec: float = 0.001, angle_tolerance: float = 5)

* **Parameters**


    * **symprec** – same as in SpacegroupAnalyzer (Default value = 1e-3)


    * **angle_tolerance** – same as in SpacegroupAnalyzer (Default value = 5).


Returns: Number of symmetrically-distinct magnetic sites present
in structure.


#### _property_ ordering(_: Orderin_ )
Applies heuristics to return a magnetic ordering for a collinear
magnetic structure. Result is not guaranteed for correctness.

Returns: Ordering Enum (‘FiM’ is used as the abbreviation for
ferrimagnetic)


#### _property_ types_of_magnetic_specie(_: tuple[[Element](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element) | [Species](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species) | [DummySpecies](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.DummySpecies), ..._ )
Specie->Species rename. Used to maintain backwards compatibility.


#### _property_ types_of_magnetic_species(_: tuple[[Element](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element) | [Species](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species) | [DummySpecies](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.DummySpecies), ..._ )
Equivalent to Structure.types_of_specie but only returns
magnetic species.

Returns: types of Species as a list


### _class_ pymatgen.analysis.magnetism.analyzer.MagneticDeformation(type, deformation)
Bases: `tuple`

Create new instance of MagneticDeformation(type, deformation)


#### deformation()
Alias for field number 1


#### type()
Alias for field number 0


### _class_ pymatgen.analysis.magnetism.analyzer.MagneticStructureEnumerator(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), default_magmoms: dict[str, float] | None = None, strategies: list[str] | tuple[str, ...] = ('ferromagnetic', 'antiferromagnetic'), automatic: bool = True, truncate_by_symmetry: bool = True, transformation_kwargs: dict | None = None)
Bases: `object`

Combines MagneticStructureAnalyzer and MagOrderingTransformation to
automatically generate a set of transformations for a given structure
and produce a list of plausible magnetic orderings.

This class will try generated different collinear
magnetic orderings for a given input structure.

If the input structure has magnetic moments defined, it
is possible to use these as a hint as to which elements are
magnetic, otherwise magnetic elements will be guessed
(this can be changed using default_magmoms kwarg).


* **Parameters**


    * **structure** – input structure


    * **default_magmoms** – (optional, defaults provided) dict of
    magnetic elements to their initial magnetic moments in µB, generally
    these are chosen to be high-spin since they can relax to a low-spin
    configuration during a DFT electronic configuration


    * **strategies** – different ordering strategies to use, choose from:
    ferromagnetic, antiferromagnetic, antiferromagnetic_by_motif,
    ferrimagnetic_by_motif and ferrimagnetic_by_species (here, “motif”,
    means to use a different ordering parameter for symmetry inequivalent
    sites)


    * **automatic** – if True, will automatically choose sensible strategies


    * **truncate_by_symmetry** – if True, will remove very unsymmetrical
    orderings that are likely physically implausible


    * **transformation_kwargs** – keyword arguments to pass to
    MagOrderingTransformation, to change automatic cell size limits, etc.



#### available_strategies(_ = ('ferromagnetic', 'antiferromagnetic', 'ferrimagnetic_by_motif', 'ferrimagnetic_by_species', 'antiferromagnetic_by_motif', 'nonmagnetic'_ )

### _class_ pymatgen.analysis.magnetism.analyzer.Ordering(value)
Bases: `Enum`

Enumeration defining possible magnetic orderings.


#### AFM(_ = 'AFM_ )

#### FM(_ = 'FM_ )

#### FiM(_ = 'FiM_ )

#### NM(_ = 'NM_ )

#### Unknown(_ = 'Unknown_ )

### _class_ pymatgen.analysis.magnetism.analyzer.OverwriteMagmomMode(value)
Bases: `Enum`

Enumeration defining different modes for analyzer.


#### none(_ = 'none_ )

#### normalize(_ = 'normalize_ )

#### replace_all(_ = 'replace_all_ )

#### respect_sign(_ = 'respect_sign_ )

#### respect_zero(_ = 'respect_zeros_ )

### pymatgen.analysis.magnetism.analyzer.magnetic_deformation(structure_A: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), structure_B: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Calculates ‘magnetic deformation proxy’,
a measure of deformation (norm of finite strain)
between ‘non-magnetic’ (non-spin-polarized) and
ferromagnetic structures.

Adapted from Bocarsly et al. 2017,
doi: 10.1021/acs.chemmater.6b04729


* **Parameters**


    * **structure_A** – Structure


    * **structure_B** – Structure


Returns: Magnetic deformation