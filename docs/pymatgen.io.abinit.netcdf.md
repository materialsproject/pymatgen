---
layout: default
title: pymatgen.io.abinit.netcdf.md
nav_exclude: true
---

# pymatgen.io.abinit.netcdf module

Wrapper for netCDF readers.


### _class_ pymatgen.io.abinit.netcdf.ETSF_Reader(path)
Bases: `NetcdfReader`

This object reads data from a file written according to the ETSF-IO specifications.

We assume that the netcdf file contains at least the crystallographic section.

Open the Netcdf file specified by path (read mode).


#### chemical_symbols()
Chemical symbols char [number of atom species][symbol length].


#### read_abinit_hdr()
Read the variables associated to the Abinit header.

Return `AbinitHeader`


#### read_abinit_xcfunc()
Read ixc from an Abinit file. Return `XcFunc` object.


#### read_structure(cls=<class 'pymatgen.core.structure.Structure'>)
Returns the crystalline structure stored in the rootgrp.


#### typeidx_from_symbol(symbol)
Returns the type index from the chemical symbol. Note python convention.


### _class_ pymatgen.io.abinit.netcdf.NO_DEFAULT()
Bases: `object`

Signal that read_value should raise an Error.


### _class_ pymatgen.io.abinit.netcdf.NetcdfReader(path)
Bases: `object`

Wraps and extends netCDF4.Dataset. Read only mode. Supports with statements.

Additional documentation available at:

    [http://netcdf4-python.googlecode.com/svn/trunk/docs/netCDF4-module.html](http://netcdf4-python.googlecode.com/svn/trunk/docs/netCDF4-module.html)

Open the Netcdf file specified by path (read mode).


#### Error()
alias of `NetcdfReaderError`


#### close()
Close the file.


#### print_tree()
Print all the groups in the file.


#### read_dimvalue(dimname, path='/', default=<class 'pymatgen.io.abinit.netcdf.NO_DEFAULT'>)
Returns the value of a dimension.


* **Parameters**


    * **dimname** – Name of the variable


    * **path** – path to the group.


    * **default** – return default if dimname is not present and
    default is not NO_DEFAULT else raise self.Error.



#### read_keys(keys, dict_cls=<class 'monty.collections.AttrDict'>, path='/')
Read a list of variables/dimensions from file. If a key is not present the corresponding
entry in the output dictionary is set to None.


#### read_value(varname, path='/', cmode=None, default=<class 'pymatgen.io.abinit.netcdf.NO_DEFAULT'>)
Returns the values of variable with name varname in the group specified by path.


* **Parameters**


    * **varname** – Name of the variable


    * **path** – path to the group.


    * **cmode** – if cmode==”c”, a complex ndarrays is constructed and returned
    (netcdf does not provide native support from complex datatype).


    * **default** – returns default if varname is not present.
    self.Error is raised if default is set to NO_DEFAULT



* **Returns**

    numpy array if varname represents an array, scalar otherwise.



#### read_variable(varname, path='/')
Returns the variable with name varname in the group specified by path.


#### read_varnames(path='/')
List of variable names stored in the group specified by path.


#### walk_tree(top=None)
Navigate all the groups in the file starting from top.
If top is None, the root group is used.


### pymatgen.io.abinit.netcdf.as_etsfreader(file)
Return an ETSF_Reader. Accepts filename or ETSF_Reader.


### pymatgen.io.abinit.netcdf.as_ncreader(file)
Convert file into a NetcdfReader instance.
Returns reader, closeit where closeit is set to True
if we have to close the file before leaving the procedure.


### pymatgen.io.abinit.netcdf.structure_from_ncdata(ncdata, site_properties=None, cls=<class 'pymatgen.core.structure.Structure'>)
Reads and returns a pymatgen structure from a NetCDF file
containing crystallographic data in the ETSF-IO format.


* **Parameters**


    * **ncdata** – filename or NetcdfReader instance.


    * **site_properties** – Dictionary with site properties.


    * **cls** – The Structure class to instantiate.