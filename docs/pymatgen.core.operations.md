---
layout: default
title: pymatgen.core.operations.md
nav_exclude: true
---

# pymatgen.core.operations module

This module provides classes that operate on points or vectors in 3D space.


### _class_ pymatgen.core.operations.MagSymmOp(affine_transformation_matrix: ArrayLike, time_reversal: int, tol: float = 0.01)
Bases: `SymmOp`

Thin wrapper around SymmOp to extend it to support magnetic symmetry by including a time
reversal operator. Magnetic symmetry is similar to conventional crystal symmetry, except
symmetry is reduced by the addition of a time reversal operator which acts on an atom’s magnetic
moment.

Initializes the MagSymmOp from a 4x4 affine transformation matrix and time reversal
operator. In general, this constructor should not be used unless you are transferring
rotations. Use the static constructors instead to generate a SymmOp from proper rotations
and translation.


* **Parameters**


    * **affine_transformation_matrix** (*4x4 array*) – Representing an
    affine transformation.


    * **time_reversal** (*int*) – 1 or -1


    * **tol** (*float*) – Tolerance for determining if matrices are equal.



#### as_dict()

* **Returns**

    MSONABle dict



#### as_xyzt_string()
Returns a string of the form ‘x, y, z, +1’, ‘-x, -y, z, -1’,
‘-y+1/2, x+1/2, z+1/2, +1’, etc. Only works for integer rotation matrices.


#### _classmethod_ from_dict(d: dict)

* **Parameters**

    **d** – dict



* **Returns**

    MagneticSymmOp from dict representation.



#### _static_ from_rotation_and_translation_and_time_reversal(rotation_matrix: ArrayLike = ((1, 0, 0), (0, 1, 0), (0, 0, 1)), translation_vec: ArrayLike = (0, 0, 0), time_reversal: int = 1, tol: float = 0.1)
Creates a symmetry operation from a rotation matrix, translation
vector and time reversal operator.


* **Parameters**


    * **rotation_matrix** (*3x3 array*) – Rotation matrix.


    * **translation_vec** (*3x1 array*) – Translation vector.


    * **time_reversal** (*int*) – Time reversal operator, +1 or -1.


    * **tol** (*float*) – Tolerance to determine if rotation matrix is valid.



* **Returns**

    MagSymmOp object



#### _classmethod_ from_symmop(symmop: SymmOp, time_reversal)
Initialize a MagSymmOp from a SymmOp and time reversal operator.


* **Parameters**


    * **symmop** (*SymmOp*) – SymmOp


    * **time_reversal** (*int*) – Time reversal operator, +1 or -1.



* **Returns**

    MagSymmOp object



#### _static_ from_xyzt_string(xyzt_string: str)

* **Parameters**

    **xyzt_string** (*str*) – of the form ‘x, y, z, +1’, ‘-x, -y, z, -1’,
    ‘-2y+1/2, 3x+1/2, z-y+1/2, +1’, etc.



* **Returns**

    MagSymmOp object



#### operate_magmom(magmom)
Apply time reversal operator on the magnetic moment. Note that
magnetic moments transform as axial vectors, not polar vectors.

See ‘Symmetry and magnetic structures’, Rodríguez-Carvajal and
Bourée for a good discussion. DOI: 10.1051/epjconf/20122200010


* **Parameters**


    * **magmom** – Magnetic moment as electronic_structure.core.Magmom


    * **array-like** (*class** or **as list** or **np*) –



* **Returns**

    Magnetic moment after operator applied as Magmom class



### _class_ pymatgen.core.operations.SymmOp(affine_transformation_matrix: ArrayLike, tol: float = 0.01)
Bases: `MSONable`

A symmetry operation in Cartesian space. Consists of a rotation plus a
translation. Implementation is as an affine transformation matrix of rank 4
for efficiency. Read: [http://en.wikipedia.org/wiki/Affine_transformation](http://en.wikipedia.org/wiki/Affine_transformation).


#### affine_matrix()
A 4x4 numpy.array representing the symmetry operation.

Initializes the SymmOp from a 4x4 affine transformation matrix.
In general, this constructor should not be used unless you are
transferring rotations. Use the static constructors instead to
generate a SymmOp from proper rotations and translation.


* **Parameters**


    * **affine_transformation_matrix** (*4x4 array*) – Representing an
    affine transformation.


    * **tol** (*float*) – Tolerance for determining if matrices are equal.



* **Raises**

    **ValueError** – if matrix is not 4x4.



#### apply_rotation_only(vector: ArrayLike)
Vectors should only be operated by the rotation matrix and not the
translation vector.


* **Parameters**

    **vector** (*3x1 array*) – A vector.



#### are_symmetrically_related(point_a: ArrayLike, point_b: ArrayLike, tol: float = 0.001)
Checks if two points are symmetrically related.


* **Parameters**


    * **point_a** (*3x1 array*) – First point.


    * **point_b** (*3x1 array*) – Second point.


    * **tol** (*float*) – Absolute tolerance for checking distance.



* **Returns**

    True if self.operate(point_a) == point_b or vice versa.



#### are_symmetrically_related_vectors(from_a: ArrayLike, to_a: ArrayLike, r_a: ArrayLike, from_b: ArrayLike, to_b: ArrayLike, r_b: ArrayLike, tol: float = 0.001)
Checks if two vectors, or rather two vectors that connect two points
each are symmetrically related. r_a and r_b give the change of unit
cells. Two vectors are also considered symmetrically equivalent if starting
and end point are exchanged.


* **Parameters**


    * **from_a** (*3x1 array*) – Starting point of the first vector.


    * **to_a** (*3x1 array*) – Ending point of the first vector.


    * **from_b** (*3x1 array*) – Starting point of the second vector.


    * **to_b** (*3x1 array*) – Ending point of the second vector.


    * **r_a** (*3x1 array*) – Change of unit cell of the first vector.


    * **r_b** (*3x1 array*) – Change of unit cell of the second vector.


    * **tol** (*float*) – Absolute tolerance for checking distance.



* **Returns**

    (are_related, is_reversed)



#### as_dict()

* **Returns**

    MSONable dict.



#### as_xyz_string()
Returns a string of the form ‘x, y, z’, ‘-x, -y, z’,
‘-y+1/2, x+1/2, z+1/2’, etc. Only works for integer rotation matrices.


#### _static_ from_axis_angle_and_translation(axis: ArrayLike, angle: float, angle_in_radians: bool = False, translation_vec: ArrayLike = (0, 0, 0))
Generates a SymmOp for a rotation about a given axis plus translation.


* **Parameters**


    * **axis** – The axis of rotation in Cartesian space. For example,
    [1, 0, 0]indicates rotation about x-axis.


    * **angle** (*float*) – Angle of rotation.


    * **angle_in_radians** (*bool*) – Set to True if angles are given in
    radians. Or else, units of degrees are assumed.


    * **translation_vec** – A translation vector. Defaults to zero.



* **Returns**

    SymmOp for a rotation about given axis and translation.



#### _classmethod_ from_dict(d)

* **Parameters**

    **d** – dict



* **Returns**

    SymmOp from dict representation.



#### _static_ from_origin_axis_angle(origin: ArrayLike, axis: ArrayLike, angle: float, angle_in_radians: bool = False)
Generates a SymmOp for a rotation about a given axis through an
origin.


* **Parameters**


    * **origin** (*3x1 array*) – The origin which the axis passes through.


    * **axis** (*3x1 array*) – The axis of rotation in Cartesian space. For
    example, [1, 0, 0]indicates rotation about x-axis.


    * **angle** (*float*) – Angle of rotation.


    * **angle_in_radians** (*bool*) – Set to True if angles are given in
    radians. Or else, units of degrees are assumed.



* **Returns**

    SymmOp.



#### _static_ from_rotation_and_translation(rotation_matrix: ArrayLike = ((1, 0, 0), (0, 1, 0), (0, 0, 1)), translation_vec: ArrayLike = (0, 0, 0), tol: float = 0.1)
Creates a symmetry operation from a rotation matrix and a translation
vector.


* **Parameters**


    * **rotation_matrix** (*3x3 array*) – Rotation matrix.


    * **translation_vec** (*3x1 array*) – Translation vector.


    * **tol** (*float*) – Tolerance to determine if rotation matrix is valid.



* **Returns**

    SymmOp object



#### _static_ from_xyz_string(xyz_string: str)

* **Parameters**

    **xyz_string** – string of the form ‘x, y, z’, ‘-x, -y, z’,
    ‘-2y+1/2, 3x+1/2, z-y+1/2’, etc.



* **Returns**

    SymmOp



#### _property_ inverse(_: SymmO_ )
Returns inverse of transformation.


#### _static_ inversion(origin: ArrayLike = (0, 0, 0))
Inversion symmetry operation about axis.


* **Parameters**

    **origin** (*3x1 array*) – Origin of the inversion operation. Defaults
    to [0, 0, 0].



* **Returns**

    SymmOp representing an inversion operation about the origin.



#### operate(point: ArrayLike)
Apply the operation on a point.


* **Parameters**

    **point** – Cartesian coordinate.



* **Returns**

    Coordinates of point after operation.



#### operate_multi(points: ArrayLike)
Apply the operation on a list of points.


* **Parameters**

    **points** – List of Cartesian coordinates



* **Returns**

    Numpy array of coordinates after operation



#### _static_ reflection(normal: ArrayLike, origin: ArrayLike = (0, 0, 0))
Returns reflection symmetry operation.


* **Parameters**


    * **normal** (*3x1 array*) – Vector of the normal to the plane of
    reflection.


    * **origin** (*3x1 array*) – A point in which the mirror plane passes
    through.



* **Returns**

    SymmOp for the reflection about the plane



#### _property_ rotation_matrix(_: ndarra_ )
A 3x3 numpy.array representing the rotation matrix.


#### _static_ rotoreflection(axis: ArrayLike, angle: float, origin: ArrayLike = (0, 0, 0))
Returns a roto-reflection symmetry operation.


* **Parameters**


    * **axis** (*3x1 array*) – Axis of rotation / mirror normal


    * **angle** (*float*) – Angle in degrees


    * **origin** (*3x1 array*) – Point left invariant by roto-reflection.
    Defaults to (0, 0, 0).



* **Returns**

    Roto-reflection operation



#### transform_tensor(tensor: ndarray)
Applies rotation portion to a tensor. Note that tensor has to be in
full form, not the Voigt form.


* **Parameters**

    **tensor** (*numpy array*) – a rank n tensor



* **Returns**

    Transformed tensor.



#### _property_ translation_vector(_: ndarra_ )
A rank 1 numpy.array of dim 3 representing the translation vector.