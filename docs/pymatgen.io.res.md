---
layout: default
title: pymatgen.io.res.md
nav_exclude: true
---

# pymatgen.io.res module

Provides parsing and read/write support for ShelX .res files as produced by the AIRSS code.

Converting from and back to pymatgen objects is expected to be reversible, i.e. you
should get the same Structure or ComputedStructureEntry back. On the other hand, converting
from and back to a string/file is not guaranteed to be reversible, i.e. a diff on the output
would not be empty. The difference should be limited to whitespace, float precision, and the
REM entries.


### _class_ pymatgen.io.res.AirssProvider(res: Res, parse_rems: Literal['gentle', 'strict'] = 'gentle')
Bases: `ResProvider`

Provides access to the res file as does `ResProvider`. This class additionally provides
access to fields in the TITL entry and various other fields found in the REM entries
that AIRSS puts in the file. Values in the TITL entry that AIRSS could not get end up as 0.
If the TITL entry is malformed, empty, or missing then attempting to construct this class
from a res file will raise a ResError.

While AIRSS supports a number of geometry and energy solvers, CASTEP is the default. As such,
fetching the information from the REM entries is only supported if AIRSS was used with CASTEP.
The other properties that get put in the TITL should still be accessible even if CASTEP was
not used.

The `parse_rems` attribute controls whether functions that fail to retrieve information
from the REM entries should return `None`. If this is set to `"strict"`,
then a `ParseError` may be raised, but the return value will not be `None`.
If it is set to `"gentle"`, then `None` will be returned instead of raising an
exception. This setting applies to all methods of this class that are typed to return
an Optional type. Default is `"gentle"`.

The `from_str()` and `from_file()` methods should be used instead of constructing this directly.


#### _property_ appearances(_: in_ )
This is sometimes the number of times a structure was found in an AIRSS search.
Using the cryan tool that comes with AIRSS may be a better approach than relying
on this property.


#### as_dict(verbose: bool = True)
Get dict with title fields, structure and rems of this AirssProvider.


#### _property_ energy(_: floa_ )
Energy of the structure. With CASTEP, this is usually the enthalpy and is in eV.


#### _property_ entry(_: [ComputedStructureEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry_ )
Get this res file as a ComputedStructureEntry.


#### _classmethod_ from_file(filename: str, parse_rems: Literal['gentle', 'strict'] = 'gentle')
Construct a Provider from a file.


#### _classmethod_ from_str(string: str, parse_rems: Literal['gentle', 'strict'] = 'gentle')
Construct a Provider from a string.


#### get_airss_version()
Retrieves the version of AIRSS that was used along with the build date (not compile date).


* **Returns**

    (version string, date)



#### get_castep_version()
Retrieves the version of CASTEP that the res file was computed with from the REM entries.


* **Returns**

    version string



#### get_cut_grid_gmax_fsbc()
Retrieves the cut-off energy, grid scale, Gmax, and finite basis set correction setting
from the REM entries.


* **Returns**

    (cut-off, grid scale, Gmax, fsbc)



#### get_func_rel_disp()
Retrieves the functional, relativity scheme, and dispersion correction from the REM entries.


* **Returns**

    (functional, relativity, dispersion)



#### get_mpgrid_offset_nkpts_spacing()
Retrieves the MP grid, the grid offsets, number of kpoints, and maximum kpoint spacing.


* **Returns**

    (MP grid), (offsets), No. kpts, max spacing)



#### get_pspots()
Retrieves the OTFG pseudopotential string that can be used to generate the
pseudopotentials used in the calculation.


* **Returns**

    dict[specie, potential]



#### get_run_start_info()
Retrieves the run start date and the path it was started in from the REM entries.


* **Returns**

    (date, path)



#### _property_ integrated_absolute_spin_density(_: floa_ )
Corresponds to the last `Integrated |Spin Density|` in the CASTEP file.


#### _property_ integrated_spin_density(_: floa_ )
Corresponds to the last `Integrated Spin Density` in the CASTEP file.


#### _property_ pressure(_: floa_ )
Pressure for the run. This is in GPa if CASTEP was used.


#### _property_ seed(_: st_ )
The seed name, typically also the name of the res file.


#### _property_ spacegroup_label(_: st_ )
The Hermann-Mauguin notation of the spacegroup with ascii characters.
So no. 225 would be Fm-3m, and no. 194 would be P6_3/mmc.


#### _property_ volume(_: floa_ )
Volume of the structure. This is in cubic Angstroms if CASTEP was used.


### _exception_ pymatgen.io.res.ResError()
Bases: `ValueError`

This exception indicates a problem was encountered while trying to retrieve a value or
perform an action that a provider for the res file does not support.


### _class_ pymatgen.io.res.ResIO()
Bases: `object`

Class providing convenience methods for converting a Structure or ComputedStructureEntry
to/from a string or file in the res format as used by AIRSS.

Note: Converting from and back to pymatgen objects is expected to be reversible, i.e. you
should get the same Structure or ComputedStructureEntry back. On the other hand, converting
from and back to a string/file is not guaranteed to be reversible, i.e. a diff on the output
would not be empty. The difference should be limited to whitespace, float precision, and the
REM entries.

If the TITL entry doesn’t exist or is malformed or empty, then you can only get
a Structure. Attempting to get an Entry will raise a ResError.


#### _classmethod_ entry_from_file(filename: str)
Produce a pymatgen ComputedStructureEntry from a res file.


#### _classmethod_ entry_from_str(string: str)
Produce a pymatgen ComputedStructureEntry from contents of a res file.


#### _classmethod_ entry_to_file(entry: [ComputedStructureEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry), filename: str)
Write a pymatgen ComputedStructureEntry to a res file.


#### _classmethod_ entry_to_str(entry: [ComputedStructureEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry))
Produce the contents of a res file from a pymatgen ComputedStructureEntry.


#### _classmethod_ structure_from_file(filename: str)
Produces a pymatgen Structure from a res file.


#### _classmethod_ structure_from_str(string: str)
Produces a pymatgen Structure from contents of a res file.


#### _classmethod_ structure_to_file(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), filename: str)
Write a pymatgen Structure to a res file.


#### _classmethod_ structure_to_str(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Produce the contents of a res file from a pymatgen Structure.


### _exception_ pymatgen.io.res.ResParseError()
Bases: [`ParseError`](pymatgen.io.core.md#pymatgen.io.core.ParseError)

This exception indicates a problem was encountered during parsing due to unexpected formatting.


### _class_ pymatgen.io.res.ResProvider(res: Res)
Bases: `MSONable`

Provides access to elements of the res file in the form of familiar pymatgen objects.

The `from_str()` and `from_file()` methods should be used instead of constructing this directly.


#### _classmethod_ from_file(filename: str)
Construct a Provider from a file.


#### _classmethod_ from_str(string: str)
Construct a Provider from a string.


#### _property_ lattice(_: [Lattice](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice_ )
Construct a Lattice from the res file.


#### _property_ rems(_: list[str_ )
The full list of REM entries contained within the res file.


#### _property_ sites(_: list[[pymatgen.core.sites.PeriodicSite](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite)_ )
Construct a list of PeriodicSites from the res file.


#### _property_ structure(_: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure_ )
Construct a Structure from the res file.


### _class_ pymatgen.io.res.ResWriter(entry: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure) | [ComputedStructureEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry))
Bases: `object`

This class provides a means to write a Structure or ComputedStructureEntry to a res file.

This class can be constructed from either a pymatgen Structure or ComputedStructureEntry object.


#### _property_ string(_: st_ )
The contents of the res file.


#### write(filename: str)
Write the res data to a file.