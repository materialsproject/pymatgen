---
layout: default
title: pymatgen.io.feff.outputs.md
nav_exclude: true
---

# pymatgen.io.feff.outputs module

This module defines classes for parsing the FEFF output files.

Currently supports the xmu.dat, ldos.dat output files are for non-spin case.


### _class_ pymatgen.io.feff.outputs.Eels(data)
Bases: `MSONable`

Parse’eels.dat’ file.


* **Parameters**

    **(****)** (*data*) – Eels data.



#### as_dict()
Returns dict representations of Xmu object.


#### _property_ atomic_background()
atomic background.


* **Type**

    Returns



#### _property_ energies()
Returns the energies in eV.


#### _property_ fine_structure()
Fine structure of EELS.


* **Type**

    Returns



#### _static_ from_file(eels_dat_file='eels.dat')
Parse eels spectrum.


* **Parameters**

    **eels_dat_file** (*str*) – filename and path for eels.dat



* **Returns**

    Eels object



#### _property_ total_spectrum()
Returns the total eels spectrum.


### _class_ pymatgen.io.feff.outputs.LDos(complete_dos, charge_transfer)
Bases: `MSONable`

Parser for ldos files ldos01, ldos02, …..


* **Parameters**


    * **complete_dos** ([*CompleteDos*](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.CompleteDos)) – complete dos object


    * **charge_transfer** (*dict*) – computed charge transfer between atoms
    dictionary.



#### _static_ charge_transfer_from_file(feff_inp_file, ldos_file)
Get charge transfer from file.


* **Parameters**


    * **feff_inp_file** (*str*) – name of feff.inp file for run


    * **ldos_file** (*str*) – ldos filename for run, assume consecutive order,
    i.e., ldos01.dat, ldos02.dat….



* **Returns**

    dictionary of dictionaries in order of potential sites
    ({“p”: 0.154, “s”: 0.078, “d”: 0.0, “tot”: 0.232}, …)



#### charge_transfer_to_string()
Returns charge transfer as string.


#### _static_ from_file(feff_inp_file='feff.inp', ldos_file='ldos')
Creates LDos object from raw Feff ldos files by
by assuming they are numbered consecutively, i.e. ldos01.dat
ldos02.dat…


* **Parameters**


    * **feff_inp_file** (*str*) – input file of run to obtain structure


    * **ldos_file** (*str*) – output ldos file of run to obtain dos info, etc.



### _class_ pymatgen.io.feff.outputs.Xmu(header, parameters, absorbing_atom, data)
Bases: `MSONable`

Parser for data in ‘xmu.dat’ file.
The file ‘xmu.dat’ contains XANES, EXAFS or NRIXS data depending on the
situation; \\mu, \\mu_0, and \\chi = \\chi \* \\mu_0/ \\mu_0/(edge+50eV) as
functions of absolute energy E, relative energy E - E_f and wave number k.

Default attributes:

    xmu: Photon absorption cross section of absorbing atom in material
    Energies: Energies of data point
    relative_energies: E - E_fermi
    wavenumber: k=\\sqrt(E -E_fermi)
    mu: The total absorption cross-section.
    mu0: The embedded atomic background absorption.
    chi: fine structure.
    Edge: Aborption Edge
    Absorbing atom: Species of absorbing atom
    Material: Formula of material
    Source: Source of structure
    Calculation: Type of Feff calculation performed


* **Parameters**


    * **header** – Header object


    * **parameters** – Tags object


    * **absorbing_atom** (*str/int*) – absorbing atom symbol or index


    * **data** (*numpy.ndarray**, **Nx6*) – cross_sections.



#### as_dict()
Returns dict representations of Xmu object.


#### _property_ calc()
Returns type of Feff calculation, XANES or EXAFS.


#### _property_ chi()
Returns the normalized fine structure.


#### _property_ e_fermi()
Returns the Fermi level in eV.


#### _property_ edge()
Returns excitation edge.


#### _property_ energies()
Returns the absolute energies in eV.


#### _static_ from_file(xmu_dat_file='xmu.dat', feff_inp_file='feff.inp')
Get Xmu from file.


* **Parameters**


    * **xmu_dat_file** (*str*) – filename and path for xmu.dat


    * **feff_inp_file** (*str*) – filename and path of feff.inp input file



* **Returns**

    Xmu object



#### _property_ material_formula()
Returns chemical formula of material from feff.inp file.


#### _property_ mu()
Returns the total absorption cross-section.


#### _property_ mu0()
Returns the embedded atomic background absorption.


#### _property_ relative_energies()
Returns energy with respect to the Fermi level.
E - E_f.


#### _property_ source()
Returns source identification from Header file.


#### _property_ wavenumber()
Returns The wave number in units of \\AA^-1. k=\\sqrt(E - E_f) where E is
the energy and E_f is the Fermi level computed from electron gas theory
at the average interstitial charge density.