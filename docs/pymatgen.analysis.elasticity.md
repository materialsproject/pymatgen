---
layout: default
title: pymatgen.analysis.elasticity.md
nav_exclude: true
---

1. TOC
{:toc}

# pymatgen.analysis.elasticity package

Package for analyzing elastic tensors and properties.


## pymatgen.analysis.elasticity.elastic module

This module provides a class used to describe the elastic tensor,
including methods used to fit the elastic tensor from linear response
stress-strain data.


### _class_ ComplianceTensor(s_array)
Bases: [`Tensor`](pymatgen.core.md#pymatgen.core.tensors.Tensor)

This class represents the compliance tensor, and exists
primarily to keep the voigt-conversion scheme consistent
since the compliance tensor has a unique vscale.


* **Parameters**

    **(****)** (*s_array*) –



### _class_ ElasticTensor(input_array, tol: float = 0.0001)
Bases: `NthOrderElasticTensor`

This class extends Tensor to describe the 3x3x3x3
second-order elastic tensor, C_{ijkl}, with various
methods for estimating other properties derived from
the second order elastic tensor.

Create an ElasticTensor object. The constructor throws an error if
the shape of the input_matrix argument is not 3x3x3x3, i. e. in true
tensor notation. Issues a warning if the input_matrix argument does
not satisfy standard symmetries. Note that the constructor uses
__new__ rather than __init__ according to the standard method of
subclassing numpy ndarrays.


* **Parameters**


    * **input_array** (*3x3x3x3 array-like*) – the 3x3x3x3 array-like
    representing the elastic tensor


    * **tol** (*float*) – tolerance for initial symmetry test of tensor



#### cahill_thermalcond(\*args, \*\*kwargs)

#### clarke_thermalcond(\*args, \*\*kwargs)

#### _property_ compliance_tensor()
Returns the Voigt-notation compliance tensor,
which is the matrix inverse of the
Voigt-notation elastic tensor.


#### debye_temperature(\*args, \*\*kwargs)

#### directional_elastic_mod(n)
Calculates directional elastic modulus for a specific vector.


#### directional_poisson_ratio(n, m, tol: float = 1e-08)
Calculates the poisson ratio for a specific direction
relative to a second, orthogonal direction.


* **Parameters**


    * **n** (*3-d vector*) – principal direction


    * **m** (*3-d vector*) – secondary direction orthogonal to n


    * **tol** (*float*) – tolerance for testing of orthogonality



#### _classmethod_ from_independent_strains(strains, stresses, eq_stress=None, vasp=False, tol: float = 1e-10)
Constructs the elastic tensor least-squares fit of independent strains


* **Parameters**


    * **strains** (*list** of **Strains*) – list of strain objects to fit


    * **stresses** (*list** of **Stresses*) – list of stress objects to use in fit
    corresponding to the list of strains


    * **eq_stress** (*Stress*) – equilibrium stress to use in fitting


    * **vasp** (*bool*) – flag for whether the stress tensor should be
    converted based on vasp units/convention for stress


    * **tol** (*float*) – tolerance for removing near-zero elements of the
    resulting tensor.



#### _classmethod_ from_pseudoinverse(strains, stresses)
Class method to fit an elastic tensor from stress/strain
data. Method uses Moore-Penrose pseudo-inverse to invert
the s = C\*e equation with elastic tensor, stress, and
strain in voigt notation.


* **Parameters**


    * **stresses** (*Nx3x3 array-like*) – list or array of stresses


    * **strains** (*Nx3x3 array-like*) – list or array of strains



#### _property_ g_reuss(_: floa_ )
Returns the G_r shear modulus.


#### _property_ g_voigt(_: floa_ )
Returns the G_v shear modulus.


#### _property_ g_vrh(_: floa_ )
Returns the G_vrh (Voigt-Reuss-Hill) average shear modulus.


#### get_structure_property_dict(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), include_base_props: bool = True, ignore_errors: bool = False)
Returns a dictionary of properties derived from the elastic tensor
and an associated structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – structure object for which to calculate
    associated properties


    * **include_base_props** (*bool*) – whether to include base properties,
    like k_vrh, etc.


    * **ignore_errors** (*bool*) – if set to true, will set problem properties
    that depend on a physical tensor to None, defaults to False



#### green_kristoffel(u)
Returns the Green-Kristoffel tensor for a second-order tensor.


#### _property_ homogeneous_poisson()
Returns the homogeneous poisson ratio.


#### _property_ k_reuss(_: floa_ )
Returns the K_r bulk modulus.


#### _property_ k_voigt(_: floa_ )
Returns the K_v bulk modulus.


#### _property_ k_vrh(_: floa_ )
Returns the K_vrh (Voigt-Reuss-Hill) average bulk modulus.


#### long_v(\*args, \*\*kwargs)

#### _property_ property_dict()
Returns a dictionary of properties derived from the elastic tensor.


#### snyder_ac(\*args, \*\*kwargs)

#### snyder_opt(\*args, \*\*kwargs)

#### snyder_total(\*args, \*\*kwargs)

#### trans_v(\*args, \*\*kwargs)

#### _property_ universal_anisotropy()
Returns the universal anisotropy value.


#### _property_ y_mod(_: floa_ )
Calculates Young’s modulus (in SI units) using the
Voigt-Reuss-Hill averages of bulk and shear moduli.


### _class_ ElasticTensorExpansion(c_list)
Bases: [`TensorCollection`](pymatgen.core.md#pymatgen.core.tensors.TensorCollection)

This class is a sequence of elastic tensors corresponding
to an elastic tensor expansion, which can be used to
calculate stress and energy density and inherits all
of the list-based properties of TensorCollection
(e. g. symmetrization, voigt conversion, etc.).

Initialization method for ElasticTensorExpansion.


* **Parameters**

    **c_list** (*list** or **tuple*) – sequence of Tensor inputs
    or tensors from which the elastic tensor
    expansion is constructed.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### calculate_stress(strain)
Calculate’s a given elastic tensor’s contribution to the
stress using Einstein summation.


#### energy_density(strain, convert_GPa_to_eV=True)
Calculates the elastic energy density due to a strain.


#### _classmethod_ from_diff_fit(strains, stresses, eq_stress=None, tol: float = 1e-10, order=3)
Generates an elastic tensor expansion via the fitting function
defined below in diff_fit.


#### get_compliance_expansion()
Gets a compliance tensor expansion from the elastic
tensor expansion.


#### get_effective_ecs(strain, order=2)
Returns the effective elastic constants
from the elastic tensor expansion.


* **Parameters**


    * **strain** (*Strain** or **3x3 array-like*) – strain condition
    under which to calculate the effective constants


    * **order** (*int*) – order of the ecs to be returned



#### get_ggt(n, u)
Gets the Generalized Gruneisen tensor for a given
third-order elastic tensor expansion.


* **Parameters**


    * **n** (*3x1 array-like*) – normal mode direction


    * **u** (*3x1 array-like*) – polarization direction



#### get_gruneisen_parameter(temperature=None, structure=None, quad=None)
Gets the single average gruneisen parameter from the TGT.


* **Parameters**


    * **temperature** (*float*) – Temperature in kelvin, if not specified
    will return non-cv-normalized value


    * **structure** (*float*) – Structure to be used in directional heat
    capacity determination, only necessary if temperature
    is specified


    * **quad** (*dict*) – quadrature for integration, should be
    dictionary with “points” and “weights” keys defaults
    to quadpy.sphere.Lebedev(19) as read from file



#### get_heat_capacity(temperature, structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), n, u, cutoff=100.0)
Gets the directional heat capacity for a higher order tensor
expansion as a function of direction and polarization.


* **Parameters**


    * **temperature** (*float*) – Temperature in kelvin


    * **structure** (*float*) – Structure to be used in directional heat
    capacity determination


    * **n** (*3x1 array-like*) – direction for Cv determination


    * **u** (*3x1 array-like*) – polarization direction, note that
    no attempt for verification of eigenvectors is made


    * **cutoff** (*float*) – cutoff for scale of kt / (hbar \* omega)
    if lower than this value, returns 0



#### get_stability_criteria(s, n)
Gets the stability criteria from the symmetric
Wallace tensor from an input vector and stress
value.


* **Parameters**


    * **s** (*float*) – Stress value at which to evaluate
    the stability criteria


    * **n** (*3x1 array-like*) – direction of the applied
    stress



#### get_strain_from_stress(stress)
Gets the strain from a stress state according
to the compliance expansion corresponding to the
tensor expansion.


#### get_symmetric_wallace_tensor(tau)
Gets the symmetrized wallace tensor for determining
yield strength criteria.


* **Parameters**

    **tau** (*3x3 array-like*) – stress at which to evaluate
    the wallace tensor.



#### get_tgt(temperature=None, structure=None, quad=None)
Gets the thermodynamic Gruneisen tensor (TGT) by via an
integration of the GGT weighted by the directional heat
capacity.

See refs:

    R. N. Thurston and K. Brugger, Phys. Rev. 113, A1604 (1964).
    K. Brugger Phys. Rev. 137, A1826 (1965).


* **Parameters**


    * **temperature** (*float*) – Temperature in kelvin, if not specified
    will return non-cv-normalized value


    * **structure** (*float*) – Structure to be used in directional heat
    capacity determination, only necessary if temperature
    is specified


    * **quad** (*dict*) – quadrature for integration, should be
    dictionary with “points” and “weights” keys defaults
    to quadpy.sphere.Lebedev(19) as read from file



#### get_wallace_tensor(tau)
Gets the Wallace Tensor for determining yield strength
criteria.


* **Parameters**

    **tau** (*3x3 array-like*) – stress at which to evaluate
    the wallace tensor



#### get_yield_stress(n)
Gets the yield stress for a given direction.


* **Parameters**

    **n** (*3x1 array-like*) – direction for which to find the
    yield stress



#### omega(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), n, u)
Finds directional frequency contribution to the heat
capacity from direction and polarization.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Structure to be used in directional heat
    capacity determination


    * **n** (*3x1 array-like*) – direction for Cv determination


    * **u** (*3x1 array-like*) – polarization direction, note that
    no attempt for verification of eigenvectors is made



#### _property_ order()
Order of the elastic tensor expansion, i. e. the order of the
highest included set of elastic constants.


#### thermal_expansion_coeff(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), temperature, mode='debye')
Gets thermal expansion coefficient from third-order constants.


* **Parameters**


    * **temperature** (*float*) – Temperature in kelvin, if not specified
    will return non-cv-normalized value


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Structure to be used in directional heat
    capacity determination, only necessary if temperature
    is specified


    * **mode** (*str*) – mode for finding average heat-capacity,
    current supported modes are ‘debye’ and ‘dulong-petit’



### _class_ NthOrderElasticTensor(input_array, check_rank=None, tol: float = 0.0001)
Bases: [`Tensor`](pymatgen.core.md#pymatgen.core.tensors.Tensor)

An object representing an nth-order tensor expansion
of the stress-strain constitutive equations.


* **Parameters**


    * **(****)** (*tol*) –


    * **(****)** –


    * **(****)** –



#### GPa_to_eV_A3(_ = 0.00624150907446076_ )

#### calculate_stress(strain)
Calculate’s a given elastic tensor’s contribution to the
stress using Einstein summation.


* **Parameters**

    **strain** (*3x3 array-like*) – matrix corresponding to strain



#### energy_density(strain, convert_GPa_to_eV=True)
Calculates the elastic energy density due to a strain.


#### _classmethod_ from_diff_fit(strains, stresses, eq_stress=None, order=2, tol: float = 1e-10)
Takes a list of strains and stresses, and returns a list of coefficients for a
polynomial fit of the given order.


* **Parameters**


    * **strains** – a list of strain values


    * **stresses** – the stress values


    * **eq_stress** – The stress at which the material is assumed to be elastic.


    * **order** – The order of the polynomial to fit. Defaults to 2


    * **tol** (*float*) – tolerance for the fit.



* **Returns**

    the fitted elastic tensor



* **Return type**

    NthOrderElasticTensor



#### _property_ order()
Order of the elastic tensor.


#### symbol(_ = 'C_ )

### diff_fit(strains, stresses, eq_stress=None, order=2, tol: float = 1e-10)
nth order elastic constant fitting function based on
central-difference derivatives with respect to distinct
strain states. The algorithm is summarized as follows:


1. Identify distinct strain states as sets of indices
for which nonzero strain values exist, typically
[(0), (1), (2), (3), (4), (5), (0, 1) etc.]


2. For each strain state, find and sort strains and
stresses by strain value.


3. Find first, second .. nth derivatives of each stress
with respect to scalar variable corresponding to
the smallest perturbation in the strain.


4. Use the pseudo-inverse of a matrix-vector expression
corresponding to the parameterized stress-strain
relationship and multiply that matrix by the respective
calculated first or second derivatives from the
previous step.


5. Place the calculated nth-order elastic
constants appropriately.


* **Parameters**


    * **order** (*int*) – order of the elastic tensor set to return


    * **strains** (*nx3x3 array-like*) – Array of 3x3 strains
    to use in fitting of ECs


    * **stresses** (*nx3x3 array-like*) – Array of 3x3 stresses
    to use in fitting ECs. These should be PK2 stresses.


    * **eq_stress** (*3x3 array-like*) – stress corresponding to
    equilibrium strain (i. e. “0” strain state).
    If not specified, function will try to find
    the state in the list of provided stresses
    and strains. If not found, defaults to 0.


    * **tol** (*float*) – value for which strains below
    are ignored in identifying strain states.



* **Returns**

    Set of tensors corresponding to nth order expansion of
    the stress/strain relation



### find_eq_stress(strains, stresses, tol: float = 1e-10)
Finds stress corresponding to zero strain state in stress-strain list.


* **Parameters**


    * **strains** (*Nx3x3 array-like*) – array corresponding to strains


    * **stresses** (*Nx3x3 array-like*) – array corresponding to stresses


    * **tol** (*float*) – tolerance to find zero strain state



### generate_pseudo(strain_states, order=3)
Generates the pseudo-inverse for a given set of strains.


* **Parameters**


    * **strain_states** (*6xN array like*) – a list of voigt-notation
    “strain-states”, i. e. perturbed indices of the strain
    as a function of the smallest strain e. g. (0, 1, 0, 0, 1, 0)


    * **order** (*int*) – order of pseudo-inverse to calculate



* **Returns**

    for each order tensor, these can be multiplied by the central

        difference derivative of the stress with respect to the strain state

    absent_syms: symbols of the tensor absent from the PI expression




* **Return type**

    pseudo_inverses



### get_diff_coeff(hvec, n=1)
Helper function to find difference coefficients of an
derivative on an arbitrary mesh.


* **Parameters**


    * **hvec** (*1D array-like*) – sampling stencil


    * **n** (*int*) – degree of derivative to find



### get_strain_state_dict(strains, stresses, eq_stress=None, tol: float = 1e-10, add_eq=True, sort=True)
Creates a dictionary of voigt-notation stress-strain sets
keyed by “strain state”, i. e. a tuple corresponding to
the non-zero entries in ratios to the lowest nonzero value,
e.g. [0, 0.1, 0, 0.2, 0, 0] -> (0,1,0,2,0,0)
This allows strains to be collected in stencils as to
evaluate parameterized finite difference derivatives.


* **Parameters**


    * **strains** (*Nx3x3 array-like*) – strain matrices


    * **stresses** (*Nx3x3 array-like*) – stress matrices


    * **eq_stress** (*Nx3x3 array-like*) – equilibrium stress


    * **tol** (*float*) – tolerance for sorting strain states


    * **add_eq** (*bool*) – flag for whether to add eq_strain
    to stress-strain sets for each strain state


    * **sort** (*bool*) – flag for whether to sort strain states



* **Returns**

    strain state keys and dictionaries with stress-strain data corresponding to strain state



* **Return type**

    dict



### get_symbol_list(rank, dim=6)
Returns a symbolic representation of the voigt-notation
tensor that places identical symbols for entries related
by index transposition, i. e. C_1121 = C_1211 etc.


* **Parameters**


    * **dim** (*int*) – dimension of matrix/tensor, e. g. 6 for
    voigt notation and 3 for standard


    * **rank** (*int*) – rank of tensor, e. g. 3 for third-order ECs



* **Returns**

    array representing distinct indices
    c_arr (array): array representing tensor with equivalent

    > indices assigned as above




* **Return type**

    c_vec (array)



### raise_if_unphysical(func)
Wrapper for functions or properties that should raise an error
if tensor is unphysical.


### subs(entry, cmap)
Sympy substitution function, primarily for the purposes
of numpy vectorization.


* **Parameters**


    * **entry** (*symbol** or **exp*) – sympy expr to undergo subs


    * **cmap** (*dict*) – map for symbols to values to use in subs



* **Returns**

    Evaluated expression with substitution


## pymatgen.analysis.elasticity.strain module

This module provides classes and methods used to describe deformations and
strains, including applying those deformations to structure objects and
generating deformed structure sets for further calculations.


### _class_ Deformation(deformation_gradient)
Bases: [`SquareTensor`](pymatgen.core.md#pymatgen.core.tensors.SquareTensor)

Subclass of SquareTensor that describes the deformation gradient tensor.

Create a Deformation object. Note that the constructor uses __new__ rather than
__init__ according to the standard method of subclassing numpy ndarrays.


* **Parameters**

    **deformation_gradient** (*3x3 array-like*) – the 3x3 array-like
    representing the deformation gradient



#### apply_to_structure(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
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

### _class_ DeformedStructureSet(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), norm_strains=None, shear_strains=None, symmetry=False)
Bases: `Sequence`

class that generates a set of independently deformed structures that
can be used to calculate linear stress-strain response.

Construct the deformed geometries of a structure. Generates m + n deformed structures
according to the supplied parameters.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – structure to undergo deformation


    * **norm_strains** (*list** of **floats*) – strain values to apply
    to each normal mode.


    * **shear_strains** (*list** of **floats*) – strain values to apply
    to each shear mode.


    * **symmetry** (*bool*) – whether or not to use symmetry reduction.



#### _abc_impl(_ = <_abc._abc_data object_ )

### _class_ Strain(strain_matrix)
Bases: [`SquareTensor`](pymatgen.core.md#pymatgen.core.tensors.SquareTensor)

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


### convert_strain_to_deformation(strain, shape: Literal['upper', 'lower', 'symmetric'])
This function converts a strain to a deformation gradient that will
produce that strain. Supports three methods:


* **Parameters**


    * **strain** (*3x3 array-like*) – strain matrix


    * **shape** – (‘upper’ | ‘lower’ | ‘symmetric’): method for determining deformation
    ‘upper’ produces an upper triangular defo
    ‘lower’ produces a lower triangular defo
    ‘symmetric’ produces a symmetric defo


## pymatgen.analysis.elasticity.stress module

This module provides the Stress class used to create, manipulate, and
calculate relevant properties of the stress tensor.


### _class_ Stress(stress_matrix)
Bases: [`SquareTensor`](pymatgen.core.md#pymatgen.core.tensors.SquareTensor)

This class extends SquareTensor as a representation of the
stress.

Create a Stress object. Note that the constructor uses __new__
rather than __init__ according to the standard method of
subclassing numpy ndarrays.


* **Parameters**

    **stress_matrix** (*3x3 array-like*) – the 3x3 array-like
    representing the stress



#### _property_ dev_principal_invariants()
Returns the principal invariants of the deviatoric stress tensor,
which is calculated by finding the coefficients of the characteristic
polynomial of the stress tensor minus the identity times the mean
stress.


#### _property_ deviator_stress()
Returns the deviatoric component of the stress.


#### _property_ mean_stress()
Returns the mean stress.


#### piola_kirchoff_1(def_grad)
Calculates the first Piola-Kirchoff stress.


* **Parameters**

    **def_grad** (*3x3 array-like*) – deformation gradient tensor



#### piola_kirchoff_2(def_grad)
Calculates the second Piola-Kirchoff stress.


* **Parameters**

    **def_grad** (*3x3 array-like*) – rate of deformation tensor



#### symbol(_ = 's_ )

#### _property_ von_mises()
Returns the von Mises stress.