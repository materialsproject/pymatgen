---
layout: default
title: pymatgen.electronic_structure.core.md
nav_exclude: true
---

# pymatgen.electronic_structure.core module

This module provides core classes needed by all define electronic structure,
such as the Spin, Orbital, etc.


### _class_ pymatgen.electronic_structure.core.Magmom(moment, saxis=(0, 0, 1))
Bases: `MSONable`

New class in active development. Use with caution, feedback is
appreciated.

Class to handle magnetic moments. Defines the magnetic moment of a
site or species relative to a spin quantization axis. Designed for
use in electronic structure calculations.


* For the general case, Magmom can be specified by a vector,
e.g. m = Magmom([1.0, 1.0, 2.0]), and subscripts will work as
expected, e.g. m[0] gives 1.0


* For collinear calculations, Magmom can assumed to be scalar-like,
e.g. m = Magmom(5.0) will work as expected, e.g. float(m) gives 5.0

Both of these cases should be safe and shouldn’t give any surprises,
but more advanced functionality is available if required.

There also exist useful static methods for lists of magmoms:


* Magmom.are_collinear(magmoms) - if true, a collinear electronic
structure calculation can be safely initialized, with float(Magmom)
giving the expected scalar magnetic moment value


* Magmom.get_consistent_set_and_saxis(magmoms) - for non-collinear
electronic structure calculations, a global, consistent spin axis
has to be used. This method returns a list of Magmoms which all
share a common spin axis, along with the global spin axis.

All methods that take lists of magmoms will accept magmoms either as
Magmom objects or as scalars/lists and will automatically convert to
a Magmom representation internally.

The following methods are also particularly useful in the context of
VASP calculations:


* Magmom.get_xyz_magmom_with_001_saxis()


* Magmom.get_00t_magmom_with_xyz_saxis()

See VASP documentation for more information:

[https://cms.mpi.univie.ac.at/wiki/index.php/SAXIS](https://cms.mpi.univie.ac.at/wiki/index.php/SAXIS)


* **Parameters**


    * **moment** – magnetic moment, supplied as float or list/np.ndarray


    * **saxis** – spin axis, supplied as list/np.ndarray, parameter will
    be converted to unit vector (default is [0, 0, 1])



* **Returns**

    Magmom object



#### _static_ are_collinear(magmoms)
Method checks to see if a set of magnetic moments are collinear
with each other.
:param magmoms: list of magmoms (Magmoms, scalars or vectors)
:return: bool.


#### _classmethod_ from_global_moment_and_saxis(global_moment, saxis)
Convenience method to initialize Magmom from a given global
magnetic moment, i.e. magnetic moment with saxis=(0,0,1), and
provided saxis.

Method is useful if you do not know the components of your
magnetic moment in frame of your desired saxis.


* **Parameters**


    * **global_moment** –


    * **saxis** – desired saxis



* **Returns**




#### _classmethod_ from_moment_relative_to_crystal_axes(moment, lattice)
Obtaining a Magmom object from a magnetic moment provided
relative to crystal axes.

Used for obtaining moments from magCIF file.
:param moment: list of floats specifying vector magmom
:param lattice: Lattice
:return: Magmom


#### get_00t_magmom_with_xyz_saxis()
For internal implementation reasons, in non-collinear calculations VASP prefers the following.

> MAGMOM = 0 0 total_magnetic_moment
> SAXIS = x y z

to an equivalent:

> MAGMOM = x y z
> SAXIS = 0 0 1

This method returns a Magmom object with magnetic moment [0, 0, t],
where t is the total magnetic moment, and saxis rotated as required.

A consistent direction of saxis is applied such that t might be positive
or negative depending on the direction of the initial moment. This is useful
in the case of collinear structures, rather than constraining assuming
t is always positive.


* **Returns**

    Magmom



#### _static_ get_consistent_set_and_saxis(magmoms, saxis=None)
Method to ensure a list of magmoms use the same spin axis.
Returns a tuple of a list of Magmoms and their global spin axis.


* **Parameters**


    * **magmoms** – list of magmoms (Magmoms, scalars or vectors)


    * **saxis** – can provide a specific global spin axis



* **Returns**

    (list of Magmoms, global spin axis) tuple



#### get_moment(saxis=(0, 0, 1))
Get magnetic moment relative to a given spin quantization axis.
If no axis is provided, moment will be given relative to the
Magmom’s internal spin quantization axis, i.e. equivalent to
Magmom.moment.


* **Parameters**

    **saxis** – (list/numpy array) spin quantization axis



* **Returns**

    np.ndarray of length 3



#### get_moment_relative_to_crystal_axes(lattice)
If scalar magmoms, moments will be given arbitrarily along z.
Used for writing moments to magCIF file.


* **Parameters**

    **lattice** – Lattice



* **Returns**

    vector as list of floats



#### _static_ get_suggested_saxis(magmoms)
This method returns a suggested spin axis for a set of magmoms,
taking the largest magnetic moment as the reference. For calculations
with collinear spins, this would give a sensible saxis for a ncl
calculation.


* **Parameters**

    **magmoms** – list of magmoms (Magmoms, scalars or vectors)



* **Returns**

    np.ndarray of length 3



#### get_xyz_magmom_with_001_saxis()
Returns a Magmom in the default setting of saxis = [0, 0, 1] and
the magnetic moment rotated as required.


* **Returns**

    Magmom



#### _property_ global_moment()
Get the magnetic moment defined in an arbitrary global reference frame.


* **Returns**

    np.ndarray of length 3



#### _static_ have_consistent_saxis(magmoms)
This method checks that all Magmom objects in a list have a
consistent spin quantization axis. To write MAGMOM tags to a
VASP INCAR, a global SAXIS value for all magmoms has to be used.
If saxis are inconsistent, can create consistent set with:
Magmom.get_consistent_set(magmoms).


* **Parameters**

    **magmoms** – list of magmoms (Magmoms, scalars or vectors)



* **Returns**

    bool



#### _property_ projection()
Projects moment along spin quantisation axis. Useful for obtaining
collinear approximation for slightly non-collinear magmoms.


* **Returns**

    float



### _class_ pymatgen.electronic_structure.core.Orbital(value)
Bases: `Enum`

Enum type for specific orbitals. The indices are basically the order in
which the orbitals are reported in VASP and has no special meaning.


#### dx2(_ = _ )

#### dxy(_ = _ )

#### dxz(_ = _ )

#### dyz(_ = _ )

#### dz2(_ = _ )

#### f0(_ = 1_ )

#### f1(_ = 1_ )

#### f2(_ = 1_ )

#### f3(_ = 1_ )

#### f_1(_ = 1_ )

#### f_2(_ = 1_ )

#### f_3(_ = _ )

#### _property_ orbital_type()
Returns OrbitalType of an orbital.


#### px(_ = _ )

#### py(_ = _ )

#### pz(_ = _ )

#### s(_ = _ )

### _class_ pymatgen.electronic_structure.core.OrbitalType(value)
Bases: `Enum`

Enum type for orbital type. Indices are basically the azimuthal quantum
number, l.


#### d(_ = _ )

#### f(_ = _ )

#### p(_ = _ )

#### s(_ = _ )

### _class_ pymatgen.electronic_structure.core.Spin(value)
Bases: `Enum`

Enum type for Spin. Only up and down.
Usage: Spin.up, Spin.down.


#### down(_ = -_ )

#### up(_ = _ )