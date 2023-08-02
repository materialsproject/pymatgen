---
layout: default
title: pymatgen.symmetry.groups.md
nav_exclude: true
---

# pymatgen.symmetry.groups module

Defines SymmetryGroup parent class and PointGroup and SpaceGroup classes.
Shyue Ping Ong thanks Marc De Graef for his generous sharing of his
SpaceGroup data as published in his textbook “Structure of Materials”.


### _class_ pymatgen.symmetry.groups.PointGroup(\*args, \*\*kwargs)
Bases: `PointGroup`

Class representing a Point Group, with generators and symmetry operations.


#### symbol()
Full International or Hermann-Mauguin Symbol.


#### generators()
List of generator matrices. Note that 3x3 matrices are used for Point
Groups.


#### symmetry_ops()
Full set of symmetry operations as matrices.

Initializes a Point Group from its international symbol.


* **Parameters**

    **int_symbol** (*str*) – International or Hermann-Mauguin Symbol.



### _class_ pymatgen.symmetry.groups.SpaceGroup(\*args, \*\*kwargs)
Bases: `SpaceGroup`

Class representing a SpaceGroup.


#### symbol()
Full International or Hermann-Mauguin Symbol.


#### int_number()
International number


#### generators()
List of generator matrices. Note that 4x4 matrices are used for Space
Groups.


#### order()
Order of Space Group

Initializes a Space Group from its full or abbreviated international
symbol. Only standard settings are supported.


* **Parameters**

    **int_symbol** (*str*) – Full International (e.g., “P2/m2/m2/m”) or
    Hermann-Mauguin Symbol (“Pmmm”) or abbreviated symbol. The
    notation is a LaTeX-like string, with screw axes being
    represented by an underscore. For example, “P6_3/mmc”.
    Alternative settings can be accessed by adding a “:identifier”.
    For example, the hexagonal setting  for rhombohedral cells can be
    accessed by adding a “:H”, e.g., “R-3m:H”. To find out all
    possible settings for a spacegroup, use the get_settings()
    classmethod. Alternative origin choices can be indicated by a
    translation vector, e.g., ‘Fm-3m(a-1/4,b-1/4,c-1/4)’.



### _class_ pymatgen.symmetry.groups.SymmetryGroup()
Bases: `Sequence`, [`Stringify`](pymatgen.util.string.md#pymatgen.util.string.Stringify)

Abstract class representing a symmetry group.


#### is_subgroup(supergroup: SymmetryGroup)
True if this group is a subgroup of the supplied group.


* **Parameters**

    **supergroup** (*SymmetryGroup*) – Supergroup to test.



* **Returns**

    True if this group is a subgroup of the supplied group.



#### is_supergroup(subgroup: SymmetryGroup)
True if this group is a supergroup of the supplied group.


* **Parameters**

    **subgroup** (*SymmetryGroup*) – Subgroup to test.



* **Returns**

    True if this group is a supergroup of the supplied group.



#### _abstract property_ symmetry_ops(_: set[[SymmOp](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp)_ )
Returns:
List of symmetry operations associated with the group.


#### to_latex_string()

* **Returns**

    A latex formatted group symbol with proper subscripts and overlines.



### pymatgen.symmetry.groups.in_array_list(array_list: list[np.ndarray] | np.ndarray, arr: np.ndarray, tol: float = 1e-05)
Extremely efficient nd-array comparison using numpy’s broadcasting. This
function checks if a particular array a, is present in a list of arrays.
It works for arrays of any size, e.g., even matrix searches.


* **Parameters**


    * **array_list** (*[**array**]*) – A list of arrays to compare to.


    * **arr** (*array*) – The test array for comparison.


    * **tol** (*float*) – The tolerance. Defaults to 1e-5. If 0, an exact match is done.



* **Returns**

    (bool)



### pymatgen.symmetry.groups.sg_symbol_from_int_number(int_number: int, hexagonal: bool = True)
Obtains a SpaceGroup name from its international number.


* **Parameters**


    * **int_number** (*int*) – International number.


    * **hexagonal** (*bool*) – For rhombohedral groups, whether to return the
    hexagonal setting (default) or rhombohedral setting.



* **Returns**

    Spacegroup symbol



* **Return type**

    str