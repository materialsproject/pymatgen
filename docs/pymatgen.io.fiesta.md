---
layout: default
title: pymatgen.io.fiesta.md
nav_exclude: true
---

# pymatgen.io.fiesta module

This module implements input and output for Fiesta ([http://perso.neel.cnrs.fr/xavier.blase/fiesta/index.html](http://perso.neel.cnrs.fr/xavier.blase/fiesta/index.html)).

and

-Nwchem2Fiesta class: to create the input files needed for a Fiesta run
-Fiesta_run: run gw_fiesta and bse_fiesta
-Localised Basis set reader


### _class_ pymatgen.io.fiesta.BSEOutput(filename)
Bases: `object`

A bse output file parser. The start…

All energies are in eV.


* **Parameters**

    **filename** – Filename to read.



### _class_ pymatgen.io.fiesta.BasisSetReader(filename)
Bases: `object`

A basis set reader.
Basis set are stored in data as a dict:
:key l_zeta_ng for each nl orbitals which contain list of tuple (alpha, coef) for each of the ng gaussians
in l_zeta orbital.


* **Parameters**

    **filename** – Filename to read.



#### infos_on_basis_set()
Infos on the basis set as in Fiesta log.


#### set_n_nlmo()

* **Returns**

    the number of nlm orbitals for the basis set



### _class_ pymatgen.io.fiesta.FiestaInput(mol, correlation_grid: dict[str, str] | None = None, Exc_DFT_option: dict[str, str] | None = None, COHSEX_options: dict[str, str] | None = None, GW_options: dict[str, str] | None = None, BSE_TDDFT_options: dict[str, str] | None = None)
Bases: `MSONable`

Input File for Fiesta called “cell.in” by default (mandatory in Fiesta for now).


* **Parameters**


    * **mol** – pymatgen mol


    * **correlation_grid** – dict


    * **Exc_DFT_option** – dict


    * **COHSEX_options** – dict


    * **GW_options** – dict


    * **BSE_TDDFT_options** – dict



#### as_dict()

* **Returns**

    MSONable dict



#### dump_BSE_data_in_GW_run(BSE_dump=True)

* **Parameters**

    **BSE_dump** – boolean



* **Returns**

    set the “do_bse” variable to one in cell.in



#### dump_TDDFT_data_in_GW_run(TDDFT_dump=True)

* **Parameters**

    **TDDFT_dump** – boolean



* **Returns**

    set the do_tddft variable to one in cell.in



#### _classmethod_ from_dict(d)

* **Parameters**

    **d** – Dict representation



* **Returns**

    FiestaInput



#### _classmethod_ from_file(filename)
Read an Fiesta input from a file. Currently tested to work with
files generated from this class itself.


* **Parameters**

    **filename** – Filename to parse.



* **Returns**

    FiestaInput object



#### _classmethod_ from_str(string_input)
Read an FiestaInput from a string. Currently tested to work with
files generated from this class itself.


* **Parameters**

    **string_input** – string_input to parse.



* **Returns**

    FiestaInput object



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### _property_ infos_on_system()
Returns infos on initial parameters as in the log file of Fiesta.


#### _static_ make_FULL_BSE_Densities_folder(folder)
Mkdir “FULL_BSE_Densities” folder (needed for bse run) in the desired folder.


#### _property_ molecule()
Returns molecule associated with this FiestaInput.


#### set_BSE_options(n_excitations=10, nit_bse=200)
Set parameters in cell.in for a BSE computation
:param nv_bse: number of valence bands
:param nc_bse: number of conduction bands
:param n_excitations: number of excitations
:param nit_bse: number of iterations.


#### set_GW_options(nv_band=10, nc_band=10, n_iteration=5, n_grid=6, dE_grid=0.5)
Set parameters in cell.in for a GW computation
:param nv__band: number of valence bands to correct with GW
:param nc_band: number of conduction bands to correct with GW
:param n_iteration: number of iteration
:param n_grid and dE_grid:: number of points and spacing in eV for correlation grid.


#### set_auxiliary_basis_set(folder, auxiliary_folder, auxiliary_basis_set_type='aug_cc_pvtz')
copy in the desired folder the needed auxiliary basis set “X2.ion” where X is a specie.
:param auxiliary_folder: folder where the auxiliary basis sets are stored
:param auxiliary_basis_set_type: type of basis set (string to be found in the extension of the file name; must

> be in lower case). ex: C2.ion_aug_cc_pvtz_RI_Weigend find “aug_cc_pvtz”.


#### write_file(filename)
Write FiestaInput to a file
:param filename: Filename.


### _class_ pymatgen.io.fiesta.FiestaOutput(filename)
Bases: `object`

A Fiesta output file parser.

All energies are in eV.


* **Parameters**

    **filename** – Filename to read.



### _class_ pymatgen.io.fiesta.FiestaRun(folder: str | None = None, grid: tuple[int, int, int] = (2, 2, 2), log_file: str = 'log')
Bases: `MSONable`

To run FIESTA inside python:

    if grid is [x,x] then bse runs
    if grid is [x,x,y] the fiesta(gw) runs
    otherwise it breaks.


* **Parameters**


    * **folder** – Folder to look for runs.


    * **grid** –


    * **log_file** – logfile of Fiesta.



#### as_dict()

* **Returns**

    MSONable dict



#### bse_run()
Performs BSE run.


#### _classmethod_ from_dict(d)

* **Parameters**

    **d** – Dict representation



* **Returns**

    FiestaRun



#### run()
Performs FIESTA (gw) run.


### _class_ pymatgen.io.fiesta.Nwchem2Fiesta(folder, filename='nwchem', log_file='log_n2f')
Bases: `MSONable`

To run NWCHEM2FIESTA inside python:

If nwchem.nw is the input, nwchem.out the output, and structure.movecs the
“movecs” file, the syntax to run NWCHEM2FIESTA is: NWCHEM2FIESTA
nwchem.nw  nwchem.nwout  structure.movecs > log_n2f

folder: where are stored the nwchem
filename: name of nwchem files read by NWCHEM2FIESTA (filename.nw, filename.nwout and filename.movecs)
logfile: logfile of NWCHEM2FIESTA.

the run method launches NWCHEM2FIESTA


#### as_dict()

* **Returns**

    MSONable dict



#### _classmethod_ from_dict(d)

* **Parameters**

    **d** – Dict representation.



* **Returns**

    Nwchem2Fiesta



#### run()
Performs actual NWCHEM2FIESTA run.