---
layout: default
title: pymatgen.core.tensors.md
nav_exclude: true
---

# pymatgen.core.tensors module

This module provides a base class for tensor-like objects and methods for
basic tensor manipulation. It also provides a class, SquareTensor,
that provides basic methods for creating and manipulating rank 2 tensors.


### _class_ pymatgen.core.tensors.SquareTensor(input_array, vscale=None)
Bases: `Tensor`

Base class for doing useful general operations on second rank tensors
(stress, strain etc.).

Create a SquareTensor object. Note that the constructor uses __new__ rather than
__init__ according to the standard method of subclassing numpy ndarrays. Error
is thrown when the class is initialized with non-square matrix.


* **Parameters**


    * **input_array** (*3x3 array-like*) – the 3x3 array-like
    representing the content of the tensor


    * **vscale** (*6x1 array-like*) – 6x1 array-like scaling the
    voigt-notation vector with the tensor entries



#### _property_ det()
Shorthand for the determinant of the SquareTensor.


#### get_scaled(scale_factor)
Scales the tensor by a certain multiplicative scale factor.


* **Parameters**

    **scale_factor** (*float*) – scalar multiplier to be applied to the
    SquareTensor object



#### _property_ inv()
Shorthand for matrix inverse on SquareTensor.


#### is_rotation(tol: float = 0.001, include_improper=True)
Test to see if tensor is a valid rotation matrix, performs a
test to check whether the inverse is equal to the transpose
and if the determinant is equal to one within the specified
tolerance.


* **Parameters**


    * **tol** (*float*) – tolerance to both tests of whether the
    the determinant is one and the inverse is equal
    to the transpose


    * **include_improper** (*bool*) – whether to include improper
    rotations in the determination of validity



#### polar_decomposition(side='right')
Calculates matrices for polar decomposition.


#### _property_ principal_invariants()
Returns a list of principal invariants for the tensor,
which are the values of the coefficients of the characteristic
polynomial for the matrix.


#### refine_rotation()
Helper method for refining rotation matrix by ensuring
that second and third rows are perpendicular to the first.
Gets new y vector from an orthogonal projection of x onto y
and the new z vector from a cross product of the new x and y.


* **Parameters**

    **rotation** (*tol to test for*) –



* **Returns**

    new rotation matrix



#### _property_ trans()
Shorthand for transpose on SquareTensor.


### _class_ pymatgen.core.tensors.Tensor(input_array, vscale=None, check_rank=None)
Bases: `ndarray`, `MSONable`

Base class for doing useful general operations on Nth order tensors,
without restrictions on the type (stress, elastic, strain, piezo, etc.).

Create a Tensor object. Note that the constructor uses __new__
rather than __init__ according to the standard method of
subclassing numpy ndarrays.


* **Parameters**


    * **input_array** – (array-like with shape 3^N): array-like representing
    a tensor quantity in standard (i. e. non-voigt) notation


    * **vscale** – (N x M array-like): a matrix corresponding
    to the coefficients of the voigt-notation tensor


    * **check_rank** – (int): If not None, checks that input_array’s rank == check_rank.
    Defaults to None.



#### as_dict(voigt: bool = False)
Serializes the tensor object.


* **Parameters**

    **voigt** (*bool*) – flag for whether to store entries in
    voigt-notation. Defaults to false, as information
    may be lost in conversion.


Returns (dict):

    serialized format tensor object


#### average_over_unit_sphere(quad=None)
Method for averaging the tensor projection over the unit
with option for custom quadrature.


* **Parameters**

    **quad** (*dict*) – quadrature for integration, should be
    dictionary with “points” and “weights” keys defaults
    to quadpy.sphere.Lebedev(19) as read from file



* **Returns**

    Average of tensor projected into vectors on the unit sphere



#### convert_to_ieee(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), initial_fit=True, refine_rotation=True)
Given a structure associated with a tensor, attempts a
calculation of the tensor in IEEE format according to
the 1987 IEEE standards.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – a structure associated with the
    tensor to be converted to the IEEE standard


    * **initial_fit** (*bool*) – flag to indicate whether initial
    tensor is fit to the symmetry of the structure.
    Defaults to true. Note that if false, inconsistent
    results may be obtained due to symmetrically
    equivalent, but distinct transformations
    being used in different versions of spglib.


    * **refine_rotation** (*bool*) – whether to refine the rotation
    produced by the ieee transform generator, default True



#### einsum_sequence(other_arrays, einsum_string=None)
Calculates the result of an einstein summation expression.


#### fit_to_structure(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), symprec: float = 0.1)
Returns a tensor that is invariant with respect to symmetry
operations corresponding to a structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – structure from which to generate
    symmetry operations


    * **symprec** (*float*) – symmetry tolerance for the Spacegroup Analyzer
    used to generate the symmetry operations



#### _classmethod_ from_dict(d)
Instantiate Tensors from dicts (using MSONable API).


* **Returns**

    hydrated tensor object



* **Return type**

    Tensor



#### _classmethod_ from_values_indices(values, indices, populate=False, structure=None, voigt_rank=None, vsym=True, verbose=False)
Creates a tensor from values and indices, with options
for populating the remainder of the tensor.


* **Parameters**


    * **values** (*floats*) – numbers to place at indices


    * **indices** (*array-likes*) – indices to place values at


    * **populate** (*bool*) – whether to populate the tensor


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – structure to base population
    or fit_to_structure on


    * **voigt_rank** (*int*) – full tensor rank to indicate the
    shape of the resulting tensor. This is necessary
    if one provides a set of indices more minimal than
    the shape of the tensor they want, e.g.
    Tensor.from_values_indices((0, 0), 100)


    * **vsym** (*bool*) – whether to voigt symmetrize during the
    optimization procedure


    * **verbose** (*bool*) – whether to populate verbosely



#### _classmethod_ from_voigt(voigt_input)
Constructor based on the voigt notation vector or matrix.


* **Parameters**

    **voigt_input** (*array-like*) – voigt input for a given tensor



#### get_grouped_indices(voigt=False, \*\*kwargs)
Gets index sets for equivalent tensor values.


* **Parameters**


    * **voigt** (*bool*) – whether to get grouped indices
    of voigt or full notation tensor, defaults
    to false


    * **\*\*kwargs** – keyword args for np.isclose. Can take atol
    and rtol for absolute and relative tolerance, e. g.

    ```python
    >>> tensor.group_array_indices(atol=1e-8)
    ```

    or

    ```python
    >>> tensor.group_array_indices(rtol=1e-5)
    ```




* **Returns**

    list of index groups where tensor values are equivalent to
    within tolerances



#### _static_ get_ieee_rotation(structure, refine_rotation=True)
Given a structure associated with a tensor, determines
the rotation matrix for IEEE conversion according to
the 1987 IEEE standards.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – a structure associated with the
    tensor to be converted to the IEEE standard


    * **refine_rotation** (*bool*) – whether to refine the rotation
    using SquareTensor.refine_rotation



#### get_symbol_dict(voigt=True, zero_index=False, \*\*kwargs)
Creates a summary dict for tensor with associated symbol.


* **Parameters**


    * **voigt** (*bool*) – whether to get symbol dict for voigt
    notation tensor, as opposed to full notation,
    defaults to true


    * **zero_index** (*bool*) – whether to set initial index to zero,
    defaults to false, since tensor notations tend to use
    one-indexing, rather than zero indexing like python


    * **\*\*kwargs** – keyword args for np.isclose. Can take atol
    and rtol for absolute and relative tolerance, e. g.

    ```python
    >>> tensor.get_symbol_dict(atol=1e-8)
    ```

    or

    ```python
    >>> tensor.get_symbol_dict(rtol=1e-5)
    ```




* **Returns**

    list of index groups where tensor values are equivalent to
    within tolerances



#### _static_ get_voigt_dict(rank)
Returns a dictionary that maps indices in the tensor to those
in a voigt representation based on input rank.


* **Parameters**

    **rank** (*int*) – Tensor rank to generate the voigt map



#### is_fit_to_structure(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), tol: float = 0.01)
Tests whether a tensor is invariant with respect to the
symmetry operations of a particular structure by testing
whether the residual of the symmetric portion is below a
tolerance.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – structure to be fit to


    * **tol** (*float*) – tolerance for symmetry testing



#### is_symmetric(tol: float = 1e-05)
Tests whether a tensor is symmetric or not based on the residual
with its symmetric part, from self.symmetrized.


* **Parameters**

    **tol** (*float*) – tolerance to test for symmetry



#### is_voigt_symmetric(tol: float = 1e-06)
Tests symmetry of tensor to that necessary for voigt-conversion
by grouping indices into pairs and constructing a sequence of
possible permutations to be used in a tensor transpose.


#### populate(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), prec: float = 1e-05, maxiter: int = 200, verbose: bool = False, precond: bool = True, vsym: bool = True)
Takes a partially populated tensor, and populates the non-zero
entries according to the following procedure, iterated until
the desired convergence (specified via prec) is achieved.


1. Find non-zero entries


2. Symmetrize the tensor with respect to crystal symmetry and
(optionally) voigt symmetry


3. Reset the non-zero entries of the original tensor


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – structure to base population on


    * **prec** (*float*) – precision for determining a non-zero value


    * **maxiter** (*int*) – maximum iterations for populating the tensor


    * **verbose** (*bool*) – whether to populate verbosely


    * **precond** (*bool*) – whether to precondition by cycling through
    all symmops and storing new nonzero values, default True


    * **vsym** (*bool*) – whether to enforce voigt symmetry, defaults
    to True



* **Returns**

    Populated tensor



* **Return type**

    Tensor



#### project(n)
Convenience method for projection of a tensor into a
vector. Returns the tensor dotted into a unit vector
along the input n.


* **Parameters**

    **n** (*3x1 array-like*) – direction to project onto


Returns (float):

    scalar value corresponding to the projection of
    the tensor into the vector


#### rotate(matrix, tol: float = 0.001)
Applies a rotation directly, and tests input matrix to ensure a valid
rotation.


* **Parameters**


    * **matrix** (*3x3 array-like*) – rotation matrix to be applied to tensor


    * **tol** (*float*) – tolerance for testing rotation matrix validity



#### round(decimals=0)
Wrapper around numpy.round to ensure object
of same type is returned.


* **Parameters**

    **decimals** – Number of decimal places to round to (default: 0).
    If decimals is negative, it specifies the number of
    positions to the left of the decimal point.


Returns (Tensor):

    rounded tensor of same type


#### structure_transform(original_structure, new_structure, refine_rotation=True)
Transforms a tensor from one basis for an original structure
into a new basis defined by a new structure.


* **Parameters**


    * **original_structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – structure corresponding
    to the basis of the current tensor


    * **new_structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – structure corresponding to the
    desired basis


    * **refine_rotation** (*bool*) – whether to refine the rotations
    generated in get_ieee_rotation



* **Returns**

    Tensor that has been transformed such that its basis
    corresponds to the new_structure’s basis



#### symbol(_ = 'T_ )

#### _property_ symmetrized()
Returns a generally symmetrized tensor, calculated by taking
the sum of the tensor and its transpose with respect to all
possible permutations of indices.


#### transform(symm_op)
Applies a transformation (via a symmetry operation) to a tensor.


* **Parameters**

    **symm_op** ([*SymmOp*](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp)) – a symmetry operation to apply to the tensor



#### _property_ voigt()
Returns the tensor in Voigt notation.


#### _property_ voigt_symmetrized()
Returns a “voigt”-symmetrized tensor, i. e. a voigt-notation
tensor such that it is invariant wrt permutation of indices.


#### zeroed(tol: float = 0.001)
Returns the matrix with all entries below a certain threshold (i.e. tol) set to zero.


### _class_ pymatgen.core.tensors.TensorCollection(tensor_list, base_class=<class 'pymatgen.core.tensors.Tensor'>)
Bases: `Sequence`, `MSONable`

A sequence of tensors that can be used for fitting data
or for having a tensor expansion.


* **Parameters**


    * **tensor_list** – List of tensors.


    * **base_class** – Class to be used.



#### as_dict(voigt=False)

* **Parameters**

    **voigt** – Whether to use voight form.



* **Returns**

    Dict representation of TensorCollection.



#### convert_to_ieee(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), initial_fit=True, refine_rotation=True)
Convert all tensors to IEEE.


* **Parameters**


    * **structure** – Structure


    * **initial_fit** – Whether to perform an initial fit.


    * **refine_rotation** – Whether to refine the rotation.



* **Returns**

    TensorCollection.



#### fit_to_structure(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), symprec: float = 0.1)
Fits all tensors to a Structure.


* **Parameters**


    * **structure** – Structure


    * **symprec** – symmetry precision.



* **Returns**

    TensorCollection.



#### _classmethod_ from_dict(d)
Creates TensorCollection from dict.


* **Parameters**

    **d** – dict



* **Returns**

    TensorCollection



#### _classmethod_ from_voigt(voigt_input_list, base_class=<class 'pymatgen.core.tensors.Tensor'>)
Creates TensorCollection from voigt form.


* **Parameters**


    * **voigt_input_list** – List of voigt tensors


    * **base_class** – Class for tensor.



* **Returns**

    TensorCollection.



#### is_fit_to_structure(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), tol: float = 0.01)

* **Parameters**


    * **structure** – Structure


    * **tol** – tolerance



* **Returns**

    Whether all tensors are fitted to Structure.



#### is_symmetric(tol: float = 1e-05)

* **Parameters**

    **tol** – tolerance



* **Returns**

    Whether all tensors are symmetric.



#### is_voigt_symmetric(tol: float = 1e-06)

* **Parameters**

    **tol** – tolerance



* **Returns**

    Whether all tensors are voigt symmetric.



#### _property_ ranks()
Ranks for all tensors.


* **Type**

    return



#### rotate(matrix, tol: float = 0.001)
Rotates TensorCollection.


* **Parameters**


    * **matrix** – Rotation matrix.


    * **tol** – tolerance.



* **Returns**

    TensorCollection.



#### round(\*args, \*\*kwargs)
Round all tensors.


* **Parameters**


    * **args** – Passthrough to Tensor.round


    * **kwargs** – Passthrough to Tensor.round



* **Returns**

    TensorCollection.



#### _property_ symmetrized()
TensorCollection where all tensors are symmetrized.


* **Type**

    return



#### transform(symm_op)
Transforms TensorCollection with a symmetry operation.


* **Parameters**

    **symm_op** – SymmetryOperation.



* **Returns**

    TensorCollection.



#### _property_ voigt()
TensorCollection where all tensors are in voight form.


* **Type**

    return



#### _property_ voigt_symmetrized()
TensorCollection where all tensors are voigt symmetrized.


* **Type**

    return



#### zeroed(tol: float = 0.001)

* **Parameters**

    **tol** – Tolerance



* **Returns**

    TensorCollection where small values are set to 0.



### _class_ pymatgen.core.tensors.TensorMapping(tensors: Sequence[Tensor] = (), values: Sequence = (), tol: float = 1e-05)
Bases: `MutableMapping`

Base class for tensor mappings, which function much like
a dictionary, but use numpy routines to determine approximate
equality to keys for getting and setting items.

This is intended primarily for convenience with things like
stress-strain pairs and fitting data manipulation. In general,
it is significantly less robust than a typical hashing
and should be used with care.

Initialize a TensorMapping.


* **Parameters**


    * **tensors** (*Sequence**[**Tensor**]**, **optional*) – Defaults to (,).


    * **values** (*Sequence**, **optional*) – Values to be associated with tensors. Defaults to (,).


    * **tol** (*float**, **optional*) – an absolute tolerance for getting and setting items in the mapping.
    Defaults to 1e-5.



* **Raises**

    **ValueError** – if tensors and values are not the same length



#### items()

* **Returns**

    Items in mapping.



#### values()

* **Returns**

    Values in mapping.



### pymatgen.core.tensors.get_uvec(vec)
Gets a unit vector parallel to input vector.


### pymatgen.core.tensors.symmetry_reduce(tensors, structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), tol: float = 1e-08, \*\*kwargs)
Function that converts a list of tensors corresponding to a structure
and returns a dictionary consisting of unique tensor keys with symmop
values corresponding to transformations that will result in derivative
tensors from the original list.


* **Parameters**


    * **tensors** (*list** of **tensors*) – list of Tensor objects to test for
    symmetrically-equivalent duplicates


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – structure from which to get symmetry


    * **tol** (*float*) – tolerance for tensor equivalence


    * **kwargs** – keyword arguments for the SpacegroupAnalyzer



* **Returns**

    dictionary consisting of unique tensors with symmetry operations
    corresponding to those which will reconstruct the remaining
    tensors as values