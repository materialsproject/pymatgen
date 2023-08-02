---
layout: default
title: pymatgen.analysis.elasticity.strain.md
nav_exclude: true
---

# pymatgen.analysis.elasticity.strain module

This module provides classes and methods used to describe deformations and
strains, including applying those deformations to structure objects and
generating deformed structure sets for further calculations.


### _class_ pymatgen.analysis.elasticity.strain.Deformation(deformation_gradient)
Bases: [`SquareTensor`](pymatgen.core.tensors.md#pymatgen.core.tensors.SquareTensor)

Subclass of SquareTensor that describes the deformation gradient tensor.

Create a Deformation object. Note that the constructor uses __new__
rather than __init__ according to the standard method of subclassing
numpy ndarrays.


* **Parameters**

    **deformation_gradient** (*3x3 array-like*) – the 3x3 array-like
    representing the deformation gradient



#### apply_to_structure(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Apply the deformation gradient to a structure.


* **Parameters**

    **structure** (*Structure object*) – the structure object to
    be modified by the deformation



#### _classmethod_ from_index_amount(matrixpos, amt)
Factory method for constructing a Deformation object
from a matrix position and amount.


* **Parameters**


    * **matrixpos** (*tuple*) – tuple corresponding the matrix position to
    have a perturbation added


    * **amt** (*float*) – amount to add to the identity matrix at position
    matrixpos



#### get_perturbed_indices(tol: float = 1e-08)
Gets indices of perturbed elements of the deformation gradient,
i. e. those that differ from the identity.


#### _property_ green_lagrange_strain()
Calculates the Euler-Lagrange strain from the deformation gradient.


#### is_independent(tol: float = 1e-08)
Checks to determine whether the deformation is independent.


#### symbol(_ = 'd_ )

### _class_ pymatgen.analysis.elasticity.strain.DeformedStructureSet(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), norm_strains=None, shear_strains=None, symmetry=False)
Bases: `Sequence`

class that generates a set of independently deformed structures that
can be used to calculate linear stress-strain response.

Construct the deformed geometries of a structure. Generates m + n deformed structures
according to the supplied parameters.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – structure to undergo deformation


    * **norm_strains** (*list** of **floats*) – strain values to apply
    to each normal mode.


    * **shear_strains** (*list** of **floats*) – strain values to apply
    to each shear mode.


    * **symmetry** (*bool*) – whether or not to use symmetry reduction.



### _class_ pymatgen.analysis.elasticity.strain.Strain(strain_matrix)
Bases: [`SquareTensor`](pymatgen.core.tensors.md#pymatgen.core.tensors.SquareTensor)

Subclass of SquareTensor that describes the Green-Lagrange strain tensor.

Create a Strain object. Note that the constructor uses __new__
rather than __init__ according to the standard method of
subclassing numpy ndarrays. Note also that the default constructor
does not include the deformation gradient.


* **Parameters**

    **strain_matrix** (*ArrayLike*) – 3x3 matrix or length-6 Voigt notation vector
    representing the Green-Lagrange strain



#### _classmethod_ from_deformation(deformation: ArrayLike)
Factory method that returns a Strain object from a deformation
gradient.


* **Parameters**

    **deformation** (*ArrayLike*) – 3x3 array defining the deformation



#### _classmethod_ from_index_amount(idx, amount)
Like Deformation.from_index_amount, except generates
a strain from the zero 3x3 tensor or Voigt vector with
the amount specified in the index location. Ensures
symmetric strain.


* **Parameters**


    * **idx** (*tuple** or **integer*) – index to be perturbed, can be Voigt or full-tensor notation


    * **amount** (*float*) – amount to perturb selected index



#### get_deformation_matrix(shape: Literal['upper', 'lower', 'symmetric'] = 'upper')
Returns the deformation matrix.


* **Parameters**

    **shape** (*'upper'** | **'lower'** | **'symmetric'*) – method for determining deformation
    ‘upper’ produces an upper triangular defo
    ‘lower’ produces a lower triangular defo
    ‘symmetric’ produces a symmetric defo



#### symbol(_ = 'e_ )

#### _property_ von_mises_strain()
Equivalent strain to Von Mises Stress.


### pymatgen.analysis.elasticity.strain.convert_strain_to_deformation(strain, shape: Literal['upper', 'lower', 'symmetric'])
This function converts a strain to a deformation gradient that will
produce that strain. Supports three methods:


* **Parameters**


    * **strain** (*3x3 array-like*) – strain matrix


    * **shape** – (‘upper’ | ‘lower’ | ‘symmetric’): method for determining deformation
    ‘upper’ produces an upper triangular defo
    ‘lower’ produces a lower triangular defo
    ‘symmetric’ produces a symmetric defo