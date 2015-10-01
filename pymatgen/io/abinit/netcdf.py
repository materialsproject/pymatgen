# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""Wrapper for netCDF readers."""
from __future__ import unicode_literals, division, print_function

import os.path

from monty.dev import requires, deprecated
from monty.collections import AttrDict
from monty.functools import lazy_property
from pymatgen.core.units import ArrayWithUnit
from pymatgen.core.structure import Structure

import logging
logger = logging.getLogger(__name__)


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
    "structure_from_ncdata",
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


class NO_DEFAULT(object):
    """Signal that read_value should raise an Error"""


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
            raise self.Error("In file %s: %s" % (self.path, str(exc)))

        self.ngroups = len(list(self.walk_tree()))

        #self.path2group = collections.OrderedDict()
        #for children in self.walk_tree():
        #   for child in children:
        #       #print(child.group,  child.path)
        #       self.path2group[child.path] = child.group

    def __enter__(self):
        """Activated when used in the with statement."""
        return self

    def __exit__(self, type, value, traceback):
        """Activated at the end of the with statement. It automatically closes the file."""
        self.rootgrp.close()

    def close(self):
        try:
            self.rootgrp.close()
        except Exception as exc:
            logger.warning("Exception %s while trying to close %s" % (exc, self.path))

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

    def read_value(self, varname, path="/", cmode=None, default=NO_DEFAULT):
        """
        Returns the values of variable with name varname in the group specified by path.

        Args:
            varname: Name of the variable
            path: path to the group.
            cmode: if cmode=="c", a complex ndarrays is constructed and returned
                (netcdf does not provide native support from complex datatype).
            default: read_value returns default if varname is not present.

        Returns:
            numpy array if varname represents an array, scalar otherwise.
        """
        try:
            var = self.read_variable(varname, path=path)
        except self.Error:
            if default is NO_DEFAULT: raise
            return default

        if cmode is None:
            # scalar or array
            # getValue is not portable!
            try:
                return var.getValue()[0] if not var.shape else var[:]
            except IndexError:
                return var.getValue() if not var.shape else var[:]

        else:
            assert var.shape[-1] == 2
            if cmode == "c":
                return var[...,0] + 1j*var[...,1]
            else:
                raise ValueError("Wrong value for cmode %s" % cmode)

    def read_variable(self, varname, path="/"):
        """Returns the variable with name varname in the group specified by path."""
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
            raise self.Error("In file %s:\ndimnames %s, kwargs %s" % (self.path, dimnames, kwargs))

    def _read_variables(self, *varnames, **kwargs):
        path = kwargs.get("path", "/")
        try:
            if path == "/":
                return [self.rootgrp.variables[vname] for vname in varnames]
            else:
                group = self.path2group[path]
                return [group.variables[vname] for vname in varnames]

        except KeyError:
            raise self.Error("In file %s:\nvarnames %s, kwargs %s" % (self.path, varnames, kwargs))

    def read_keys(self, keys, dict_cls=AttrDict, path="/"):
        """
        Read a list of variables/dimensions from file. If a key is not present the corresponding
        entry in the output dictionary is set to None.
        """
        od = dict_cls()
        for k in keys:
            try:
                # Try to read a variable.
                od[k] = self.read_value(k, path=path)
            except self.Error:
                try:
                    # Try to read a dimension.
                    od[k] = self.read_dimvalue(k, path=path)
                except self.Error:
                    od[k] = None

        return od


class ETSF_Reader(NetcdfReader):
    """
    This object reads data from a file written according to the ETSF-IO specifications.

    We assume that the netcdf file contains at least the crystallographic section.
    """
    @lazy_property
    def chemical_symbols(self):
        """Chemical symbols char [number of atom species][symbol length]."""
        charr = self.read_value("chemical_symbols")
        symbols = []
        for v in charr:
            symbols.append("".join(c for c in v))

        #symbols = ["".join(str(c)) for symb in symbols for c in symb]
        #symbols = [s.decode("ascii") for s in symbols]
        #chemical_symbols = [str("".join(s)) for s in symbols]
        #print(symbols)
        return symbols

    def typeidx_from_symbol(self, symbol):
        """Returns the type index from the chemical symbol. Note python convention."""
        return self.chemical_symbols.index(symbol)

    def read_structure(self, cls=Structure):
        """Returns the crystalline structure."""
        if self.ngroups != 1:
            raise NotImplementedError("In file %s: ngroups != 1" % self.path)

        return structure_from_ncdata(self, cls=cls)


def structure_from_ncdata(ncdata, site_properties=None, cls=Structure):
    """
    Reads and returns a pymatgen structure from a NetCDF file
    containing crystallographic data in the ETSF-IO format.

    Args:
        ncdata: filename or NetcdfReader instance.
        site_properties: Dictionary with site properties.
        cls: The Structure class to instanciate.
    """
    ncdata, closeit = as_ncreader(ncdata)

    # TODO check whether atomic units are used
    lattice = ArrayWithUnit(ncdata.read_value("primitive_vectors"), "bohr").to("ang")

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

    structure = cls(lattice, species, red_coords, site_properties=d)

    # Quick and dirty hack.
    # I need an abipy structure since I need to_abivars and other methods.
    try:
        from abipy.core.structure import Structure as AbipyStructure
        structure.__class__ = AbipyStructure
    except ImportError:
        pass

    if closeit:
        ncdata.close()

    return structure
