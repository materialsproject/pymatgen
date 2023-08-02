---
layout: default
title: pymatgen.io.cp2k.utils.md
nav_exclude: true
---

# pymatgen.io.cp2k.utils module

Utility functions for assisting with cp2k IO.


### pymatgen.io.cp2k.utils.chunk(string: str)
Chunk the string from a cp2k basis or potential file.


### pymatgen.io.cp2k.utils.get_truncated_coulomb_cutoff(inp_struct: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Get the truncated Coulomb cutoff for a given structure.


### pymatgen.io.cp2k.utils.get_unique_site_indices(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure) | [Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule))
Get unique site indices for a structure according to site properties. Whatever site-property
has the most unique values is used for indexing.

For example, if you have magnetic CoO with half Co atoms having a positive moment, and the
other half having a negative moment. Then this function will create a dict of sites for
Co_1, Co_2, O. This function also deals with “Species” properties like oxi_state and spin by
pushing them to site properties.

This creates unique sites, based on site properties, but does not have anything to do with
turning those site properties into CP2K input parameters. This will only be done for properties
which can be turned into CP2K input parameters, which are stored in parsable_site_properties.


### pymatgen.io.cp2k.utils.natural_keys(text: str)
Sort text by numbers coming after an underscore with natural number
convention,
Ex: [file_1, file_12, file_2] becomes [file_1, file_2, file_12].


### pymatgen.io.cp2k.utils.postprocessor(data: str)
Helper function to post process the results of the pattern matching functions in Cp2kOutput
and turn them to Python types.


* **Parameters**

    **data** (*str*) – The data to be post processed.



* **Raises**

    **ValueError** – If the data cannot be parsed.



* **Returns**

    The post processed data.



* **Return type**

    str | int | float | bool | None



### pymatgen.io.cp2k.utils.preprocessor(data: str, dir: str = '.')
Cp2k contains internal preprocessor flags that are evaluated before execution. This helper
function recognizes those preprocessor flags and replaces them with an equivalent cp2k input
(this way everything is contained neatly in the cp2k input structure, even if the user preferred
to use the flags.

CP2K preprocessor flags (with arguments) are:

> @INCLUDE FILENAME: Insert the contents of FILENAME into the file at

>     this location.

> @SET VAR VALUE: set a variable, VAR, to have the value, VALUE.
> $VAR or ${VAR}: replace these with the value of the variable, as set

> > by the @SET flag.

> @IF/@ELIF: Not implemented yet.


* **Parameters**


    * **data** (*str*) – cp2k input to preprocess


    * **dir** (*str**, **optional*) – Path for include files. Default is ‘.’ (current directory).



* **Returns**

    Preprocessed string