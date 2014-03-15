"""Wrapper for netCDF readers."""
from __future__ import division, print_function

import os.path

from pymatgen.core.units import ArrayWithUnit
from pymatgen.core.structure import Structure
from monty.dev import requires


__author__ = "Matteo Giantomassi"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"
__email__ = "gmatteo at gmail.com"
__status__ = "Development"
__date__ = "$Feb 21, 2013M$"

__all__ = [
    "as_ncreader",
    "as_etsfreader",
    "NetcdfReader",
    "ETSF_Reader",
    "structure_from_etsf_file",
]

try:
    import netCDF4
except ImportError:
    netCDF4 = None


def _asreader(file, cls):
    closeit = False
    if not isinstance(file, cls):
        file, closeit = cls(file), True
    return file, closeit


def as_ncreader(file):
    """
    Convert file into a NetcdfReader instance.
    Returns reader, closeit where closeit is set to True
    if we have to close the file before leaving the procedure.
    """
    return _asreader(file, NetcdfReader)


def as_etsfreader(file):
    return _asreader(file, ETSF_Reader)


class NetcdfReaderError(Exception):
    """Base error class for NetcdfReader"""


class NetcdfReader(object):
    """
    Wraps and extends netCDF4.Dataset. Read only mode. Supports with statements.

    Additional documentation available at:
        http://netcdf4-python.googlecode.com/svn/trunk/docs/netCDF4-module.html
    """
    Error = NetcdfReaderError

    @requires(netCDF4 is not None, "netCDF4 must be installed to use this class")
    def __init__(self, path):
        """Open the Netcdf file specified by path (read mode)."""
        self.path = os.path.abspath(path)

        try:
            self.rootgrp = netCDF4.Dataset(self.path, mode="r")
        except Exception as exc:
            raise self.Error("%s: %s" % (self.path, str(exc)))

        self.ngroups = len(list(self.walk_tree()))

        #self.path2group = collections.OrderedDict()
        #for children in self.walk_tree():
        #   for child in children:
        #       #print child.group,  child.path
        #       self.path2group[child.path] = child.group

    def __enter__(self):
        """Activated when used in the with statement."""
        return self

    def __exit__(self, type, value, traceback):
        """
        Activated at the end of the with statement. It automatically closes the file.
        """
        self.rootgrp.close()

    def close(self):
        self.rootgrp.close()

    #@staticmethod
    #def pathjoin(*args):
    #    return "/".join(args)

    def walk_tree(self, top=None):
        """
        Navigate all the groups in the file starting from top.
        If top is None, the root group is used.
        """
        if top is None:
            top = self.rootgrp

        values = top.groups.values()
        yield values
        for value in top.groups.values():
            for children in self.walk_tree(value):
                yield children

    def print_tree(self):
        for children in self.walk_tree():
            for child in children:
                print(child)

    def read_dimvalue(self, dimname, path="/"):
        """Returns the value of a dimension."""
        dim = self._read_dimensions(dimname, path=path)[0]
        return len(dim)

    def read_varnames(self, path="/"):
        """List of variable names stored in the group specified by path."""
        if path == "/":
            return self.rootgrp.variables.keys()
        else:
            group = self.path2group[path]
            return group.variables.keys()

    def read_value(self, varname, path="/", cmode=None):
        """
        Returns the values of variable with name varname in the group specified by path.

        Args:
            varname:
                Name of the variable
            path:
                path to the group.
            cmode:
                if cmode=="c", a complex ndarrays is constructed and returned
                (netcdf does not provide native support from complex datatype).
        """
        try:
            var = self.read_variable(varname, path=path)
        except:
            raise

        if cmode is None:
            # scalar or array
            return var[0] if not var.shape else var[:]
        else:
            assert var.shape[-1] == 2
            if cmode == "c":
                return var[...,0] + 1j*var[...,1]
            else:
                raise ValueError("Wrong value for cmode %s" % cmode)

    def read_variable(self, varname, path="/"):
        """
        Returns the variable with name varname in the group specified by path.
        """
        return self._read_variables(varname, path=path)[0]

    def _read_dimensions(self, *dimnames, **kwargs):
        path = kwargs.get("path", "/")
        try:
            if path == "/":
                return [self.rootgrp.dimensions[dname] for dname in dimnames]
            else:
                group = self.path2group[path]
                return [group.dimensions[dname] for dname in dimnames]

        except KeyError:
            raise self.Error("dimnames %s, kwargs %s" % (dimnames, kwargs))

    def _read_variables(self, *varnames, **kwargs):
        path = kwargs.get("path", "/")
        try:
            if path == "/":
                return [self.rootgrp.variables[vname] for vname in varnames]
            else:
                group = self.path2group[path]
                return [group.variables[vname] for vname in varnames]

        except KeyError:
            raise self.Error("varnames %s, kwargs %s" % (varnames, kwargs))

    def read_values_with_map(self, names, map_names=None, path="/"):
        """
        Read (dimensions, variables) with a mapping.

        Args:
            names:
                list of netCDF keywords to read.
            map_names:
                dictionary used to map names to the netCDF keywords used to access data on file.
            path:
                Used to access groups.

        returns: od, missing
            od is the dictionary. Values are stored in d[name] for name in names.
            missing is a list of 2-d tuple with the keywords that are not found.
        """
        if map_names is None:
            map_names = {}

        od, missing = {}, []
        for k in names:
            try:
                key = map_names[k]
            except KeyError:
                # Read k.
                key = k

            try:
                # Try to read a variable.
                od[k] = self.read_value(key, path=path)
            except self.Error:
                try:
                    # Try to read a dimension.
                    od[k] = self.read_dimvalue(key, path=path)
                except self.Error:
                    # key is missing!
                    missing.append((k, key))

        return od, missing


class ETSF_Reader(NetcdfReader):
    """
    This object reads data from a file written according to the
    ETSF-IO specifications.

    We assume that the netcdf file contains at least the crystallographic section.
    """
    @property
    def chemical_symbols(self):
        """Chemical symbols char [number of atom species][symbol length]."""
        if not hasattr(self, "_chemical_symbols"):
            symbols = self.read_value("chemical_symbols")
            self._chemical_symbols = []
            for s in symbols:
                self._chemical_symbols.append("".join(s))

        return self._chemical_symbols

    def typeidx_from_symbol(self, symbol):
        """Returns the type index from the chemical symbol. Note python convention."""
        return self._chemical_symbols.index(symbol)

    def read_structure(self):
        """
        Returns the crystalline structure.

        Args:
            site_properties:
                Optional dictionary with site properties.
        """
        if self.ngroups != 1:
            raise NotImplementedError("ngroups != 1")

        return structure_from_etsf_file(self)


def structure_from_etsf_file(ncdata, site_properties=None):
    """
    Reads and returns a pymatgen structure from a NetCDF file
    containing crystallographic data in the ETSF-IO format.

    Args:
        ncdata:
            filename or NetcdfReader instance.
        site_properties:
            Dictionary with site properties.
    """
    ncdata, closeit = as_ncreader(ncdata)

    # TODO check whether atomic units are used
    lattice = ArrayWithUnit(ncdata.read_value("primitive_vectors"),
                            "bohr").to("ang")

    red_coords = ncdata.read_value("reduced_atom_positions")
    natom = len(red_coords)

    znucl_type = ncdata.read_value("atomic_numbers")

    # type_atom[0:natom] --> index Between 1 and number of atom species
    type_atom = ncdata.read_value("atom_species")

    # Fortran to C index and float --> int conversion.
    species = natom * [None]
    for atom in range(natom):
        type_idx = type_atom[atom] - 1
        species[atom] = int(znucl_type[type_idx])

    d = {}
    if site_properties is not None:
        for prop in site_properties:
            d[property] = ncdata.read_value(prop)

    structure = Structure(lattice, species, red_coords, site_properties=d)

    # Quick and dirty hack.
    # I need an abipy structure since I need to_abivars and other methods.
    #from pymatgen.io.abinitio.abiobjects import AbiStructure
    #structure.__class__ = AbiStructure
    try:
        from abipy.core.structure import Structure as AbipyStructure
        structure.__class__ = AbipyStructure
    except ImportError:
        pass

    if closeit:
        ncdata.close()

    return structure
