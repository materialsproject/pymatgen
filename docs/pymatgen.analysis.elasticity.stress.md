---
layout: default
title: pymatgen.analysis.elasticity.stress.md
nav_exclude: true
---

# pymatgen.analysis.elasticity.stress module

This module provides the Stress class used to create, manipulate, and
calculate relevant properties of the stress tensor.


### _class_ pymatgen.analysis.elasticity.stress.Stress(stress_matrix)
Bases: [`SquareTensor`](pymatgen.core.tensors.md#pymatgen.core.tensors.SquareTensor)

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