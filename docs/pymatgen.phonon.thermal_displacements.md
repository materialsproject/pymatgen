---
layout: default
title: pymatgen.phonon.thermal_displacements.md
nav_exclude: true
---

# pymatgen.phonon.thermal_displacements module

This module provides classes to handle thermal displacement matrices (anisotropic displacement parameters).


### _class_ pymatgen.phonon.thermal_displacements.ThermalDisplacementMatrices(thermal_displacement_matrix_cart, structure, temperature, thermal_displacement_matrix_cif=None)
Bases: `MSONable`

Class to handle thermal displacement matrices
This class stores thermal displacement matrices in Ucart format.

An earlier implementation based on Matlab can be found here:
[https://github.com/JaGeo/MolecularToolbox](https://github.com/JaGeo/MolecularToolbox)
( J. George, A. Wang, V. L. Deringer, R. Wang, R. Dronskowski, U. Englert, CrystEngComm, 2015, 17, 7414-7422.)


* **Parameters**


    * **thermal_displacement_matrix_cart** – 2D numpy array including the thermal_displacement matrix Ucart
    1st dimension atom types, then compressed thermal displacement matrix will follow
    U11, U22, U33, U23, U13, U12 (xx, yy, zz, yz, xz, xy)
    convention similar to “thermal_displacement_matrices.yaml” in phonopy


    * **structure** – A pymatgen Structure object


    * **temperature** – temperature at which thermal displacement matrix was determined


    * **thermal_displacement_matrix_cif** – 2D numpy array including the thermal_displacement matrix Ucif format
    1st dimension atom types, then compressed thermal displacement matrix will follow
    U11, U22, U33, U23, U13, U12 (xx, yy, zz, yz, xz, xy)
    convention similar to “thermal_displacement_matrices.yaml” in phonopy.



#### _property_ B()
Computation as described in R. W. Grosse-Kunstleve, P. D. Adams, J Appl Cryst 2002, 35, 477-480.
Returns: B as a numpy array, first dimension are the atoms in the structure.


#### _property_ U1U2U3()
Computation as described in R. W. Grosse-Kunstleve, P. D. Adams, J Appl Cryst 2002, 35, 477-480.
Returns: numpy array of eigenvalues of Ucart,  first dimension are the atoms in the structure.


#### _property_ Ucif()
Computation as described in R. W. Grosse-Kunstleve, P. D. Adams, J Appl Cryst 2002, 35, 477-480.
Returns: Ucif as a numpy array, first dimension are the atoms in the structure.


#### _property_ Ustar()
Computation as described in R. W. Grosse-Kunstleve, P. D. Adams, J Appl Cryst 2002, 35, 477-480.
Returns: Ustar as a numpy array, first dimension are the atoms in the structure.


#### _property_ beta()
Computation as described in R. W. Grosse-Kunstleve, P. D. Adams, J Appl Cryst 2002, 35, 477-480.
Returns: beta as a numpy array, first dimension are the atoms in the structure.


#### compute_directionality_quality_criterion(other)
Will compute directionality of prolate displacement ellipsoids as described in
[https://doi.org/10.1039/C9CE00794F](https://doi.org/10.1039/C9CE00794F) with the earlier implementation: [https://github.com/damMroz/Angle/](https://github.com/damMroz/Angle/).


* **Parameters**


    * **other** – ThermalDisplacementMatrix


    * **compared** (*please make sure that the order** of **the atoms in both objects that are*) –


    * **Otherwise** (*is the same.*) –


    * **results** (*this analysis will deliver wrong*) –



* **Returns**

    will return a list including dicts for each atom that include “vector0”
    (largest principal axes of self object),

    > ”vector1” (largest principal axes of the other object), “angle” between both axes,

    >     These vectors can then, for example, be drawn into the structure with VESTA.
    >     Vectors are given in Cartesian coordinates




#### _static_ from_Ucif(thermal_displacement_matrix_cif, structure, temperature)
Starting from a numpy array, it will convert Ucif values into Ucart values and initialize the class.


* **Parameters**


    * **thermal_displacement_matrix_cif** – np.array,
    first dimension are the atoms,
    then reduced form of thermal displacement matrix will follow
    Order as above: U11, U22, U33, U23, U13, U12


    * **structure** – Structure object


    * **temperature** – float
    Corresponding temperature



* **Returns**

    ThermalDisplacementMatrices



#### _static_ from_cif_P1(filename: str)
Reads a cif with P1 symmetry including positions and ADPs.
Currently, no check of symmetry is performed as CifParser methods cannot be easily reused
:param filename: Filename of the cif.

Returns: ThermalDisplacementMatrices Object.


#### _static_ from_structure_with_site_properties_Ucif(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), temperature: float | None = None)
Will create this object with the help of a structure with site properties.


* **Parameters**


    * **structure** – Structure object including U11_cif, U22_cif, U33_cif, U23_cif, U13_cif, U12_cif as site


    * **properties** –


    * **temperature** – temperature for Ucif data


Returns: ThermalDisplacementMatrices Object.


#### _static_ get_full_matrix(thermal_displacement)
Transfers the reduced matrix to the full matrix (order of reduced matrix U11, U22, U33, U23, U13, U12).


* **Parameters**

    **thermal_displacement** – 2d numpy array, first dimension are the atoms



* **Returns**

    3d numpy array including thermal displacements, first dimensions are the atoms



#### _static_ get_reduced_matrix(thermal_displacement)
Transfers the full matrix to reduced matrix (order of reduced matrix U11, U22, U33, U23, U13, U12).


* **Parameters**

    **thermal_displacement** – 2d numpy array, first dimension are the atoms



* **Returns**

    3d numpy array including thermal displacements, first dimensions are the atoms



#### _property_ ratio_prolate()
This will compute ratio between largest and smallest eigenvalue of Ucart.


#### to_structure_with_site_properties_Ucif()
Transfers this object into a structure with site properties (Ucif).
This is useful for sorting the atoms in the structure including site properties.
E.g., with code like this:
def sort_order(site):

> return [site.specie.X, site.frac_coords[0], site.frac_coords[1], site.frac_coords[2]]

new_structure0 = Structure.from_sites(sorted(structure0, key=sort_order)).

Returns: Structure object.


#### visualize_directionality_quality_criterion(other, filename: str = 'visualization.vesta', which_structure: int = 0)
Will create a VESTA file for visualization of the directionality criterion.


* **Parameters**


    * **other** – ThermalDisplacementMatrices


    * **filename** – Filename of the VESTA file


    * **which_structure** – 0 means structure of the self object will be used, 1 means structure of the other
    object will be used



#### write_cif(filename)
Writes a cif including thermal displacements.


* **Parameters**

    **filename** – name of the cif file