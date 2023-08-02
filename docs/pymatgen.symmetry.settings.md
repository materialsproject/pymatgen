---
layout: default
title: pymatgen.symmetry.settings.md
nav_exclude: true
---

# pymatgen.symmetry.settings module

This module provides classes for non-standard space-group settings.


### _class_ pymatgen.symmetry.settings.JonesFaithfulTransformation(P, p)
Bases: `object`

Transformation for space-groups defined in a non-standard setting.

Transform between settings using matrix P and origin shift vector p,
using same notation as reference.

Should initialize using from_transformation_string in Jones
faithful notation, given by a string specifying both a
transformation matrix and an origin shift, with parts delimited
by a semi-colon. Best shown by example:


* a,b,c;0,0,0 is the identity (no change)


* -b+c,a+c,-a+b+c;0,0,0 is R3:r to R3:h (rhombohedral to
hexagonal setting)


* a,b,c;-1/4,-1/4,-1/4 is Pnnn:1 to Pnnn:2 (change in origin
choice)


* b,c,a;-1/2,-1/2,-1/2 is Bbab:1 to Ccca:2 (change settin
and origin)

Can transform points (coords), lattices and symmetry operations.

Used for transforming magnetic space groups since these are
commonly used in multiple settings, due to needing to transform
between magnetic and non-magnetic settings.

See: International Tables for Crystallography (2016). Vol. A,
Chapter 1.5, pp. 75-106.


#### _property_ P(_: list[list[float]_ )
transformation matrix


* **Type**

    return



#### _classmethod_ from_origin_shift(origin_shift='0,0,0')
Construct SpaceGroupTransformation from its origin shift string.


* **Parameters**

    **origin_shift** (*str**, **optional*) – Defaults to “0,0,0”.



* **Returns**

    JonesFaithfulTransformation



#### _classmethod_ from_transformation_string(transformation_string='a,b,c;0,0,0')
Construct SpaceGroupTransformation from its transformation string.


* **Parameters**

    **transformation_string** (*str**, **optional*) – Defaults to “a,b,c;0,0,0”.



* **Returns**

    JonesFaithfulTransformation



#### _property_ inverse(_: JonesFaithfulTransformatio_ )
JonesFaithfulTransformation


* **Type**

    return



#### _property_ p(_: list[float_ )
translation vector


* **Type**

    return



#### _static_ parse_transformation_string(transformation_string: str = 'a,b,c;0,0,0')

* **Parameters**

    **transformation_string** (*str**, **optional*) – Defaults to “a,b,c;0,0,0”.



* **Raises**

    **ValueError** – When transformation string fails to parse.



* **Returns**

    transformation matrix & vector



* **Return type**

    tuple[list[list[float]] | np.ndarray, list[float]]



#### transform_coords(coords: list[list[float]] | np.ndarray)
Takes a list of coordinates and transforms them.


#### transform_lattice(lattice: [Lattice](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice))
Transforms a lattice.


#### transform_symmop(symmop: [SymmOp](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp) | [MagSymmOp](pymatgen.core.operations.md#pymatgen.core.operations.MagSymmOp))
Takes a symmetry operation and transforms it.


#### _property_ transformation_string(_: st_ )
transformation string


* **Type**

    return