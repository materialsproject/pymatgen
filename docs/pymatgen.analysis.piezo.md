---
layout: default
title: pymatgen.analysis.piezo.md
nav_exclude: true
---

# pymatgen.analysis.piezo module

This module provides classes for the Piezoelectric tensor.


### _class_ pymatgen.analysis.piezo.PiezoTensor(input_array, tol: float = 0.001)
Bases: [`Tensor`](pymatgen.core.tensors.md#pymatgen.core.tensors.Tensor)

This class describes the 3x6 piezo tensor in Voigt-notation.

Create an PiezoTensor object. The constructor throws an error if
the shape of the input_matrix argument is not 3x3x3, i. e. in true
tensor notation. Note that the constructor uses __new__ rather than
__init__ according to the standard method of subclassing numpy
ndarrays.


* **Parameters**

    **input_matrix** (*3x3x3 array-like*) – the 3x6 array-like
    representing the piezo tensor



#### _classmethod_ from_vasp_voigt(input_vasp_array)

* **Parameters**

    **input_vasp_array** (*nd.array*) – Voigt form of tensor.



* **Returns**

    PiezoTensor