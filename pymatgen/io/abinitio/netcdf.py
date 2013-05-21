"""Wrapper for netCDF readers."""
from __future__ import division, print_function

import os.path

from pymatgen.core.physical_constants import Bohr2Ang, Ha2eV
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from pymatgen.electronic_structure.core import Spin
from pymatgen.util.decorators import requires

__author__ = "Matteo Giantomassi"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"
__email__ = "gmatteo at gmail.com"
__status__ = "Development"
__date__ = "$Feb 21, 2013M$"

__all__ = [
    "NetcdfReader",
    "GSR_Reader",
    "structure_from_etsf_file",
]

try:
    import netCDF4
except ImportError:
    netCDF4 = None


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
    #    return "/".join([arg for arg in args])

    def walk_tree(self, top=None):
        """
        Navigate all the groups in the file starting from top.
        If top is None, the root group is used.
        """
        if top is None: top = self.rootgrp

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

    def read_value(self, varname, path="/"):
        """
        Returns the values of variable with name varname in the group specified by path.
        """
        var = self.read_variable(varname, path=path)
        # scalar or array
        return var[0] if not var.shape else var[:]

    def read_variable(self, varname, path="/"):
        """
        Returns the variable with name varname in the group specified by path.
        """
        return self._read_variables(varname, path=path)[0]

    def read_values(self, *varnames, **kwargs):
        """Retunrs a list with the values of the variables."""
        vars = self._read_variables(*varnames, **kwargs)
        values = []
        # scalar or array
        for var in vars:
            v = var[0] if not var.shape else var[:]
            values.append(v)
        return values

    def _read_dimensions(self, *dimnames, **kwargs):
        path = kwargs.get("path", "/")
        if path == "/":
            return [self.rootgrp.dimensions[dname] for dname in dimnames]
        else:
            group = self.path2group[path]
            return [group.dimensions[dname] for dname in dimnames]

    def _read_variables(self, *varnames, **kwargs):
        path = kwargs.get("path", "/")
        if path == "/":
            return [self.rootgrp.variables[vname] for vname in varnames]
        else:
            group = self.path2group[path]
            return [group.variables[vname] for vname in varnames]

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
            except KeyError:
                try:
                    # Try to read a dimension.
                    od[k] = self.read_dimvalue(key, path=path)
                except KeyError:
                    # key is missing!
                    missing.append((k, key))

        return od, missing


    # This quantities are Abinit specific.
    # We assume that the netcdf file contains the crystallographic section.

    @property
    def chemical_symbols(self):
        """Chemical symbols char [number of atom species][symbol length]."""
        if not hasattr(self, "_chemical_symbols"):
            symbols = self.read_value("chemical_symbols")
            self._chemical_symbols = []
            for s in symbols:
                self._chemical_symbols.append("".join(c for c in s))

        return self._chemical_symbols

    def typeidx_from_symbol(self, symbol):
        """Returns the type index from the chemical symbol. Note python convention."""
        return self._chemical_symbols.index(symbol)

    def read_structure(self, site_properties=None):
        """
        Returns the crystalline structure.

        Args:
            site_properties:
                Optional dictionary with site properties.
        """
        if self.ngroups != 1:
            raise NotImplementedError("ngroups != 1")

        return structure_from_etsf_file(self, site_properties=site_properties)

################################################################################


class GSR_Reader(NetcdfReader):
    """
    This object reads the results stored in the _GSR (Ground-State Results)
    file. produced by ABINIT. It provides helper function to access the most
    important quantities.
    """

    def read_band_structure(self):
        raise NotImplementedError("")
        structure = self.read_structure()
        from pprint import pprint

        kpoints = self.read_value("reduced_coordinates_of_kpoints")
        efermi = Ha2eV(self.read_value("fermie"))
        np_eigvals = Ha2eV(self.read_value("eigenvalues"))
        # TODO
        #assert np_eigvals.units == "atomic units"
        nsppol = np_eigvals.shape[0]

        # FIXME: Here I need the labels
        labels_dict = {}
        for (i, kpoint) in enumerate(kpoints):
            labels_dict[str(i)] = kpoint

        eigenvals = {}
        for isp in range(nsppol):
            spin = Spin.up
            if isp == 1: spin = Spin.down
            eigenvals[spin] = np_eigvals[isp,:,:].transpose()
            print(eigenvals[spin].shape)
            #tmp = np_eigvals[isp,:,:].transpose()

        #bands = BandStructure(kpoints, eigenvals, structure.lattice, efermi,
        # labels_dict=None, structure=structure)

        bands = BandStructureSymmLine(kpoints, eigenvals, structure.lattice,
                                      efermi, labels_dict, structure=structure)
        return bands

##########################################################################################


def structure_from_etsf_file(ncdata, site_properties=None):
    """
    Reads and returns a pymatgen structure from a NetCDF file
    containing crystallographic data in the ETSF-IO format.

    Args:
        ncdata:
            filename or NetcdfReader instance.
        site_properties:
    """
    open_and_close = isinstance(ncdata, str)
    if open_and_close:
        ncdata = NetcdfReader(ncdata)

    # TODO check whether atomic units are used
    lattice = Bohr2Ang(ncdata.read_value("primitive_vectors"))

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

    if open_and_close:
        ncdata.close()

    return structure
