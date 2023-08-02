---
layout: default
title: pymatgen.io.cif.md
nav_exclude: true
---

# pymatgen.io.cif module

Wrapper classes for Cif input and output from Structures.


### _class_ pymatgen.io.cif.CifBlock(data, loops, header)
Bases: `object`

Object for storing cif data. All data is stored in a single dictionary.
Data inside loops are stored in lists in the data dictionary, and
information on which keys are grouped together are stored in the loops
attribute.


* **Parameters**


    * **data** – dict of data to go into the cif. Values should be convertible to string,
    or lists of these if the key is in a loop


    * **loops** – list of lists of keys, grouped by which loop they should appear in


    * **header** – name of the block (appears after the

    ```
    data_
    ```

     on the first line).



#### _classmethod_ from_str(string)
Reads CifBlock from string.


* **Parameters**

    **string** – String representation.



* **Returns**

    CifBlock



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### maxlen(_ = 7_ )

### _class_ pymatgen.io.cif.CifFile(data, orig_string=None, comment=None)
Bases: `object`

Reads and parses CifBlocks from a .cif file or string.


* **Parameters**


    * **data** (*dict*) – Of CifBlock objects.


    * **orig_string** (*str*) – The original cif string.


    * **comment** (*str*) – Comment string.



#### _classmethod_ from_file(filename)
Reads CifFile from a filename.


* **Parameters**

    **filename** – Filename



* **Returns**

    CifFile



#### _classmethod_ from_str(string)
Reads CifFile from a string.


* **Parameters**

    **string** – String representation.



* **Returns**

    CifFile



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


### _class_ pymatgen.io.cif.CifParser(filename: str | StringIO, occupancy_tolerance: float = 1.0, site_tolerance: float = 0.0001, frac_tolerance: float = 0.0001)
Bases: `object`

Parses a CIF file. Attempts to fix CIFs that are out-of-spec, but will
issue warnings if corrections applied. These are also stored in the
CifParser’s errors attribute.


* **Parameters**


    * **filename** (*str*) – CIF filename, gzipped or bzipped CIF files are fine too.


    * **occupancy_tolerance** (*float*) – If total occupancy of a site is between 1 and occupancy_tolerance, the
    occupancies will be scaled down to 1.


    * **site_tolerance** (*float*) – This tolerance is used to determine if two sites are sitting in the same position,
    in which case they will be combined to a single disordered site. Defaults to 1e-4.


    * **frac_tolerance** (*float*) – This tolerance is used to determine is a coordinate should be rounded to an ideal
    value. E.g., 0.6667 is rounded to 2/3. This is desired if symmetry operations are going to be applied.
    However, for very large CIF files, this may need to be set to 0.



#### as_dict()

* **Returns**

    MSONable dict



#### _static_ from_str(cif_string: str, \*\*kwargs)
Creates a CifParser from a string.


* **Parameters**


    * **cif_string** (*str*) – String representation of a CIF.


    * **\*\*kwargs** – Passthrough of all kwargs supported by CifParser.



* **Returns**

    CifParser



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### get_bibtex_string()
Get BibTeX reference from CIF file.
:param data:
:return: BibTeX string.


#### get_lattice(data, length_strings=('a', 'b', 'c'), angle_strings=('alpha', 'beta', 'gamma'), lattice_type=None)
Generate the lattice from the provided lattice parameters. In
the absence of all six lattice parameters, the crystal system
and necessary parameters are parsed.


#### _static_ get_lattice_no_exception(data, length_strings=('a', 'b', 'c'), angle_strings=('alpha', 'beta', 'gamma'), lattice_type=None)
Take a dictionary of CIF data and returns a pymatgen Lattice object.


* **Parameters**


    * **data** – a dictionary of the CIF file


    * **length_strings** – The strings that are used to identify the length parameters in the CIF file.


    * **angle_strings** – The strings that are used to identify the angles in the CIF file.


    * **lattice_type** – The type of lattice.  This is a string, and can be any of the following:



* **Returns**

    Lattice object



#### get_magsymops(data)
Equivalent to get_symops except for magnetic symmetry groups.
Separate function since additional operation for time reversal symmetry
(which changes magnetic moments on sites) needs to be returned.


#### get_structures(primitive: bool = True, symmetrized: bool = False, on_error: Literal['ignore', 'warn', 'raise'] = 'warn')
Return list of structures in CIF file.


* **Parameters**


    * **primitive** (*bool*) – Set to False to return conventional unit cells.
    Defaults to True. With magnetic CIF files, will return primitive
    magnetic cell which may be larger than nuclear primitive cell.


    * **symmetrized** (*bool*) – If True, return a SymmetrizedStructure which will
    include the equivalent indices and symmetry operations used to
    create the Structure as provided by the CIF (if explicit symmetry
    operations are included in the CIF) or generated from information
    in the CIF (if only space group labels are provided). Note that
    currently Wyckoff labels and space group labels or numbers are
    not included in the generated SymmetrizedStructure, these will be
    notated as “Not Parsed” or -1 respectively.


    * **on_error** (*'ignore'** | **'warn'** | **'raise'*) – What to do in case of KeyError or ValueError
    while parsing CIF file. Defaults to ‘warn’.



* **Returns**

    All structures in CIF file.



* **Return type**

    list[[Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure)]



#### get_symops(data)
In order to generate symmetry equivalent positions, the symmetry
operations are parsed. If the symops are not present, the space
group symbol is parsed, and symops are generated.


#### _property_ has_errors()
Whether there are errors/warnings detected in CIF parsing.


* **Type**

    return



#### _static_ parse_magmoms(data, lattice=None)
Parse atomic magnetic moments from data dictionary.


#### _static_ parse_oxi_states(data)
Parse oxidation states from data dictionary.


### _class_ pymatgen.io.cif.CifWriter(struct, symprec=None, write_magmoms=False, significant_figures=8, angle_tolerance=5.0, refine_struct=True)
Bases: `object`

A wrapper around CifFile to write CIF files from pymatgen structures.


* **Parameters**


    * **struct** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – structure to write


    * **symprec** (*float*) – If not none, finds the symmetry of the structure
    and writes the cif with symmetry information. Passes symprec
    to the SpacegroupAnalyzer. See also refine_struct.


    * **write_magmoms** (*bool*) – If True, will write magCIF file. Incompatible
    with symprec


    * **significant_figures** (*int*) – Specifies precision for formatting of floats.
    Defaults to 8.


    * **angle_tolerance** (*float*) – Angle tolerance for symmetry finding. Passes
    angle_tolerance to the SpacegroupAnalyzer. Used only if symprec
    is not None.


    * **refine_struct** – Used only if symprec is not None. If True, get_refined_structure
    is invoked to convert input structure from primitive to conventional.



#### _property_ ciffile()
CifFile associated with the CifWriter.


* **Type**

    Returns



#### write_file(filename)
Write the cif file.


### pymatgen.io.cif.str2float(text)
Remove uncertainty brackets from strings and return the float.