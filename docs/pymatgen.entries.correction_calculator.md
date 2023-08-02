---
layout: default
title: pymatgen.entries.correction_calculator.md
nav_exclude: true
---

# pymatgen.entries.correction_calculator module

This module calculates corrections for the species listed below, fitted to the experimental and computed
entries given to the CorrectionCalculator constructor.


### _class_ pymatgen.entries.correction_calculator.CorrectionCalculator(species: list[str] | None = None, max_error: float = 0.1, allow_unstable: float | bool = 0.1, exclude_polyanions: list[str] | None = None)
Bases: `object`

A CorrectionCalculator contains experimental and computed entries which it uses to compute corrections.

It graphs residual errors after applying the computed corrections and creates the MPCompatibility.yaml
file the Correction classes use.


#### species()
list of species that corrections are being calculated for


#### exp_compounds()
list of dictionaries which each contain a compound’s formula and experimental data


#### calc_compounds()
dictionary of ComputedEntry objects


#### corrections()
list of corrections in same order as species list


#### corrections_std_error()
list of the variances of the corrections in same order as species list


#### corrections_dict()
dictionary of format {‘species’: (value, uncertainty)} for easier correction lookup

Initializes a CorrectionCalculator.


* **Parameters**


    * **species** – list of species to calculate corrections for


    * **max_error** – maximum tolerable relative uncertainty in experimental energy.
    Compounds with relative uncertainty greater than this value will be excluded from the fit


    * **allow_unstable** – whether unstable entries are to be included in the fit. If True, all compounds will
    be included regardless of their energy above hull. If False or a float, compounds with
    energy above hull greater than the given value (defaults to 0.1 eV/atom) will be
    excluded


    * **exclude_polyanions** – a list of polyanions that contain additional sources of error that may negatively
    influence the quality of the fitted corrections. Compounds with these polyanions
    will be excluded from the fit



#### compute_corrections(exp_entries: list, calc_entries: dict)
Computes the corrections and fills in correction, corrections_std_error, and corrections_dict.


* **Parameters**


    * **exp_entries** – list of dictionary objects with the following keys/values:
    {“formula”: chemical formula, “exp energy”: formation energy in eV/formula unit,
    “uncertainty”: uncertainty in formation energy}


    * **calc_entries** – dictionary of computed entries, of the form {chemical formula: ComputedEntry}



* **Raises**

    **ValueError** – calc_compounds is missing an entry



#### compute_from_files(exp_gz: str, comp_gz: str)

* **Parameters**


    * **exp_gz** – name of .json.gz file that contains experimental data
    data in .json.gz file should be a list of dictionary objects with the following keys/values:
    {“formula”: chemical formula, “exp energy”: formation energy in eV/formula unit,
    “uncertainty”: uncertainty in formation energy}


    * **comp_gz** – name of .json.gz file that contains computed entries
    data in .json.gz file should be a dictionary of {chemical formula: ComputedEntry}.



#### graph_residual_error()
Graphs the residual errors for all compounds after applying computed corrections.


#### graph_residual_error_per_species(specie: str)
Graphs the residual errors for each compound that contains specie after applying computed corrections.


* **Parameters**

    **specie** – the specie/group that residual errors are being plotted for



* **Raises**

    **ValueError** – the specie is not a valid specie that this class fits corrections for



#### make_yaml(name: str = 'MP2020', dir: str | None = None)
Creates the _name_Compatibility.yaml that stores corrections as well as _name_CompatibilityUncertainties.yaml
for correction uncertainties.


* **Parameters**


    * **name** – str, alternate name for the created .yaml file.
    Default: “MP2020”


    * **dir** – str, directory in which to save the file. Pass None (default) to
    save the file in the current working directory.