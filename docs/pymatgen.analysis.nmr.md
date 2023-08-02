---
layout: default
title: pymatgen.analysis.nmr.md
nav_exclude: true
---

# pymatgen.analysis.nmr module

A module for NMR analysis.


### _class_ pymatgen.analysis.nmr.ChemicalShielding(cs_matrix, vscale=None)
Bases: [`SquareTensor`](pymatgen.core.tensors.md#pymatgen.core.tensors.SquareTensor)

This class extends the SquareTensor to perform extra analysis unique to
NMR Chemical shielding tensors.

Three notations to describe chemical shielding tensor (RK Harris; Magn. Resonance
Chem. 2008, 46, 582-598; DOI: 10.1002/mrc.2225) are supported.

Authors: Shyam Dwaraknath, Xiaohui Qu

Create a Chemical Shielding tensor.
Note that the constructor uses __new__
rather than __init__ according to the standard method of
subclassing numpy ndarrays.


* **Parameters**


    * **cs_matrix** (*1x3** or **3x3 array-like*) – the 3x3 array-like
    representing the chemical shielding tensor
    or a 1x3 array of the primary sigma values corresponding
    to the principal axis system


    * **vscale** (*6x1 array-like*) – 6x1 array-like scaling the
    Voigt-notation vector with the tensor entries



#### _class_ HaeberlenNotation(sigma_iso, delta_sigma_iso, zeta, eta)
Bases: `tuple`

Create new instance of HaeberlenNotation(sigma_iso, delta_sigma_iso, zeta, eta)


#### delta_sigma_iso()
Alias for field number 1


#### eta()
Alias for field number 3


#### sigma_iso()
Alias for field number 0


#### zeta()
Alias for field number 2


#### _class_ MarylandNotation(sigma_iso, omega, kappa)
Bases: `tuple`

Create new instance of MarylandNotation(sigma_iso, omega, kappa)


#### kappa()
Alias for field number 2


#### omega()
Alias for field number 1


#### sigma_iso()
Alias for field number 0


#### _class_ MehringNotation(sigma_iso, sigma_11, sigma_22, sigma_33)
Bases: `tuple`

Create new instance of MehringNotation(sigma_iso, sigma_11, sigma_22, sigma_33)


#### sigma_11()
Alias for field number 1


#### sigma_22()
Alias for field number 2


#### sigma_33()
Alias for field number 3


#### sigma_iso()
Alias for field number 0


#### _classmethod_ from_maryland_notation(sigma_iso, omega, kappa)
Initialize from Maryland notation.


* **Parameters**


    * **(****)** (*kappa*) –


    * **(****)** –


    * **(****)** –



* **Returns**

    ChemicalShielding



#### _property_ haeberlen_values()
the Chemical shielding tensor in Haeberlen Notation.


* **Type**

    Returns



#### _property_ maryland_values()
the Chemical shielding tensor in Maryland Notation.


* **Type**

    Returns



#### _property_ mehring_values()
the Chemical shielding tensor in Mehring Notation.


* **Type**

    Returns



#### _property_ principal_axis_system()
Returns a chemical shielding tensor aligned to the principle axis system
so that only the 3 diagonal components are non-zero.


### _class_ pymatgen.analysis.nmr.ElectricFieldGradient(efg_matrix, vscale=None)
Bases: [`SquareTensor`](pymatgen.core.tensors.md#pymatgen.core.tensors.SquareTensor)

This class extends the SquareTensor to perform extra analysis unique to
NMR Electric Field Gradient tensors in units of V/Angstrom^2.

Authors: Shyam Dwaraknath, Xiaohui Qu

Create a Chemical Shielding tensor.
Note that the constructor uses __new__
rather than __init__ according to the standard method of
subclassing numpy ndarrays.


* **Parameters**


    * **efg_matrix** (*1x3** or **3x3 array-like*) – the 3x3 array-like
    representing the electric field tensor
    or a 1x3 array of the primary values corresponding
    to the principal axis system


    * **vscale** (*6x1 array-like*) – 6x1 array-like scaling the
    voigt-notation vector with the tensor entries



#### _property_ V_xx()
First diagonal element.


* **Type**

    Returns



#### _property_ V_yy()
Second diagonal element.


* **Type**

    Returns



#### _property_ V_zz()
Third diagonal element.


* **Type**

    Returns



#### _property_ asymmetry()
Asymmetry of the electric field tensor defined as:
(V_yy - V_xx)/V_zz.


#### coupling_constant(specie)
Computes the coupling constant C_q as defined in:

    Wasylishen R E, Ashbrook S E, Wimperis S. NMR of quadrupolar nuclei
    in solid materials[M]. John Wiley & Sons, 2012. (Chapter 3.2).

C_q for a specific atom type for this electric field tensor:

    > C_q=e\*Q\*V_zz/h

    h: Planck’s constant
    Q: nuclear electric quadrupole moment in mb (millibarn
    e: elementary proton charge


* **Parameters**

    **specie** – flexible input to specify the species at this site.
    Can take a isotope or element string, Species object,
    or Site object



* **Returns**

    the coupling constant as a FloatWithUnit in MHz



#### _property_ principal_axis_system()
Returns a electric field gradient tensor aligned to the principle axis system so that only the 3 diagonal
components are non-zero.