---
layout: default
title: pymatgen.symmetry.kpath.md
nav_exclude: true
---

# pymatgen.symmetry.kpath module

Provides classes for generating high-symmetry k-paths using different conventions.


### _class_ pymatgen.symmetry.kpath.KPathBase(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), symprec: float = 0.01, angle_tolerance=5, atol=1e-05, \*args, \*\*kwargs)
Bases: `object`

This is the base class for classes used to generate high-symmetry
paths in reciprocal space (k-paths) for band structure calculations.

Args:
structure (Structure): Structure object.
symprec (float): Tolerance for symmetry finding.
angle_tolerance (float): Angle tolerance for symmetry finding.
atol (float): Absolute tolerance used to compare structures

> and determine symmetric equivalence of points and lines in the BZ.



```
*
```

args: Other arguments supported by subclasses.


```
**
```

kwargs: Other keyword arguments supported by subclasses.


#### get_kpoints(line_density=20, coords_are_cartesian=True)
Returns:
kpoints along the path in Cartesian coordinates
together with the critical-point labels.


#### _property_ kpath()
Returns:
The symmetry line path in reciprocal space.


#### _property_ lattice()
Returns:
The real space lattice.


#### _property_ rec_lattice()
Returns:
The reciprocal space lattice.


#### _property_ structure()
Returns:
The input structure.


### _class_ pymatgen.symmetry.kpath.KPathLatimerMunro(structure, has_magmoms=False, magmom_axis=None, symprec=0.01, angle_tolerance=5, atol=1e-05)
Bases: `KPathBase`

This class looks for a path along high-symmetry lines in the
Brillouin zone. It is based on the method outlined in:
npj Comput Mater 6, 112 (2020). 10.1038/s41524-020-00383-7
The user should ensure that the unit cell of the input structure
is as reduced as possible, i.e. that there is no linear
combination of lattice vectors which can produce a vector of
lesser magnitude than the given set (this is required to
obtain the correct Brillouin zone within the current
implementation). This is checked during initialization and a
warning is issued if the condition is not fulfilled.
In the case of magnetic structures, care must also be taken to
provide the magnetic primitive cell (i.e. that which reproduces
the entire crystal, including the correct magnetic ordering,
upon application of lattice translations). There is no algorithm to

> check for this, so if the input structure is

incorrect, the class will output the incorrect k-path without
any warning being issued.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Structure object


    * **has_magmoms** (*bool*) – Whether the input structure contains
    magnetic moments as site properties with the key ‘magmom.’
    Values may be in the form of 3-component vectors given in
    the basis of the input lattice vectors, or as scalars, in
    which case the spin axis will default to a_3, the third
    real-space lattice vector (this triggers a warning).


    * **magmom_axis** (*list** or **numpy array*) – 3-component vector specifying
    direction along which magnetic moments given as scalars
    should point. If all magnetic moments are provided as
    vectors then this argument is not used.


    * **symprec** (*float*) – Tolerance for symmetry finding


    * **angle_tolerance** (*float*) – Angle tolerance for symmetry finding.


    * **atol** (*float*) – Absolute tolerance used to determine symmetric
    equivalence of points and lines in the BZ.



#### _static_ LabelPoints(index)
Axes used in generating labels for Latimer-Munro convention.


#### _static_ LabelSymbol(index)
Letters used in generating labels for the Latimer-Munro convention.


#### _property_ mag_type()
Returns:
The type of magnetic space group as a string.
Current implementation does not distinguish
between types 3 and 4, so return value is ‘3/4’.
If has_magmoms is False, returns ‘0’.


### _class_ pymatgen.symmetry.kpath.KPathSeek(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), symprec: float = 0.01, angle_tolerance=5, atol=1e-05, system_is_tri=True)
Bases: `KPathBase`

This class looks for a path along high-symmetry lines in the Brillouin zone. It is based on
Hinuma, Y., Pizzi, G., Kumagai, Y., Oba, F., & Tanaka, I. (2017). Band structure diagram paths
based on crystallography. Computational Materials Science, 128, 140-184.
[https://doi.org/10.1016/j.commatsci.2016.10.015](https://doi.org/10.1016/j.commatsci.2016.10.015). It should be used with primitive structures that
comply with the definition given in the paper. The symmetry is determined by spglib using the
SpacegroupAnalyzer class. k-points are generated using the get_kpoints() method for the
reciprocal cell basis defined in the paper.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Structure object


    * **symprec** (*float*) – Tolerance for symmetry finding


    * **angle_tolerance** (*float*) – Angle tolerance for symmetry finding.


    * **atol** (*float*) – Absolute tolerance used to determine edge cases
    for settings of structures.


    * **system_is_tri** (*bool*) – Indicates if the system is time-reversal
    invariant.



### _class_ pymatgen.symmetry.kpath.KPathSetyawanCurtarolo(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), symprec: float = 0.01, angle_tolerance=5, atol=1e-05)
Bases: `KPathBase`

This class looks for a path along high-symmetry lines in
the Brillouin zone.
It is based on Setyawan, W., & Curtarolo, S. (2010).
High-throughput electronic band structure calculations:
Challenges and tools. Computational Materials Science,
49(2), 299-312. doi:10.1016/j.commatsci.2010.05.010
It should be used with primitive structures that
comply with the definition given in the paper.
The symmetry is determined by spglib using the
SpacegroupAnalyzer class. The analyzer can be used to
produce the correct primitive structure with the method
get_primitive_standard_structure(international_monoclinic=False).
A warning will signal possible compatibility problems
with the given structure. k-points generated using the get_kpoints() method
are returned for the reciprocal cell basis defined in the paper.

Args:
structure (Structure): Structure object.
symprec (float): Tolerance for symmetry finding.
angle_tolerance (float): Angle tolerance for symmetry finding.
atol (float): Absolute tolerance used to compare the input

> structure with the one expected as primitive standard.
> A warning will be issued if the cells don’t match.


#### bcc()
BCC Path.


#### bctet1(c, a)
BCT1 Path.


#### bctet2(c, a)
BCT2 Path.


#### _property_ conventional()
Returns:
The conventional cell structure.


#### cubic()
CUB Path.


#### fcc()
FCC Path.


#### hex()
HEX Path.


#### mcl(b, c, beta)
MCL Path.


#### mclc1(a, b, c, alpha)
MCLC1 Path.


#### mclc2(a, b, c, alpha)
MCLC2 Path.


#### mclc3(a, b, c, alpha)
MCLC3 Path.


#### mclc4(a, b, c, alpha)
MCLC4 Path.


#### mclc5(a, b, c, alpha)
MCLC5 Path.


#### orc()
ORC Path.


#### orcc(a, b, c)
ORCC Path.


#### orcf1(a, b, c)
ORFC1 Path.


#### orcf2(a, b, c)
ORFC2 Path.


#### orcf3(a, b, c)
ORFC3 Path.


#### orci(a, b, c)
ORCI Path.


#### _property_ prim()
Returns:
The primitive cell structure.


#### _property_ prim_rec()
Returns:
The primitive reciprocal cell structure.


#### rhl1(alpha)
RHL1 Path.


#### rhl2(alpha)
RHL2 Path.


#### tet()
TET Path.


#### tria()
TRI1a Path.


#### trib()
TRI1b Path.