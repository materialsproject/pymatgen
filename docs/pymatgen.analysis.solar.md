---
layout: default
title: pymatgen.analysis.solar.md
nav_exclude: true
---

1. TOC
{:toc}

# pymatgen.analysis.solar package

Modules for prediciting theoretical solar-cell efficiency.


## pymatgen.analysis.solar.slme module

Calculate spectroscopy limited maximum efficiency (SLME) given dielectric function data.

Forked and adjusted from :
[https://github.com/usnistgov/jarvis](https://github.com/usnistgov/jarvis)

References: 1) [https://doi.org/10.1021/acs.chemmater.9b02166](https://doi.org/10.1021/acs.chemmater.9b02166)  &


    1. [https://doi.org/10.1103/PhysRevLett.108.068701](https://doi.org/10.1103/PhysRevLett.108.068701)


### absorption_coefficient(dielectric)
Calculate the optical absorption coefficient from an input set of
pymatgen vasprun dielectric constant data.


* **Parameters**

    **dielectric** (*list*) – A list containing the dielectric response function
    in the pymatgen vasprun format.
    - element 0: list of energies
    - element 1: real dielectric tensors, in `[xx, yy, zz, xy, xz, yz]` format.
    - element 2: imaginary dielectric tensors, in `[xx, yy, zz, xy, xz, yz]` format.



* **Returns**

    absorption coefficient using eV as frequency units (cm^-1).



* **Return type**

    (np.array)



### get_dir_indir_gap(run='')
Get direct and indirect bandgaps for a vasprun.xml.


### optics(path='')
Helper function to calculate optical absorption coefficient.


### parse_dielectric_data(data)
Convert a set of 2D vasprun formatted dielectric data to
the eigenvalues of each corresponding 3x3 symmetric numpy matrices.


* **Parameters**

    **data** (*list*) – length N list of dielectric data. Each entry should be
    a list of `[xx, yy, zz, xy , xz, yz ]` dielectric tensor elements.



* **Returns**

    a Nx3 numpy array. Each row contains the eigenvalues

        for the corresponding row in data.




* **Return type**

    np.array



### slme(material_energy_for_absorbance_data, material_absorbance_data, material_direct_allowed_gap, material_indirect_gap, thickness=5e-05, temperature=293.15, absorbance_in_inverse_centimeters=False, cut_off_absorbance_below_direct_allowed_gap=True, plot_current_voltage=False)
Calculate the SLME.


* **Parameters**


    * **material_energy_for_absorbance_data** – energy grid for absorbance data


    * **material_absorbance_data** – absorption coefficient in m^-1


    * **material_direct_allowed_gap** – direct bandgap in eV


    * **material_indirect_gap** – indirect bandgap in eV


    * **thickness** – thickness of the material in m


    * **temperature** – working temperature in K


    * **absorbance_in_inverse_centimeters** – whether the absorbance data is in the unit of cm^-1


    * **cut_off_absorbance_below_direct_allowed_gap** – whether to discard all absorption below bandgap


    * **plot_current_voltage** – whether to plot the current-voltage curve



* **Returns**

    The calculated maximum efficiency.



### to_matrix(xx, yy, zz, xy, yz, xz)
Convert a list of matrix components to a symmetric 3x3 matrix.
Inputs should be in the order xx, yy, zz, xy, yz, xz.


* **Parameters**


    * **xx** (*float*) – xx component of the matrix.


    * **yy** (*float*) – yy component of the matrix.


    * **zz** (*float*) – zz component of the matrix.


    * **xy** (*float*) – xy component of the matrix.


    * **yz** (*float*) – yz component of the matrix.


    * **xz** (*float*) – xz component of the matrix.



* **Returns**

    The matrix, as a 3x3 numpy array.



* **Return type**

    (np.array)