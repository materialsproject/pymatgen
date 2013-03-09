"""Wrapper for netCDF readers."""
from __future__ import division, print_function

import numpy as np
import os.path 
import collections 

from pymatgen.core.structure import Structure
from pymatgen.core.physical_constants import Bohr2Ang

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

##########################################################################################

class NetcdfReader(object):
    "Wraps and extend netCDF4.Dataset. Read only mode"

    def __init__(self, filename):
        self.path = os.path.abspath(filename)

        try:
            #raise ImportError
            import netCDF4
            try:
                self.rootgrp = netCDF4.Dataset(self.path, mode="r")
            except Exception as exc:
                raise RuntimeError("%s : %s" % (self.path, str(exc)))

            self.ngroups = len( list(self.walk_tree()) )
            self._have_netcdf4 = True

        except ImportError:
            from scipy.io import netcdf
            try:
                self.rootgrp = netcdf.netcdf_file(self.path, mode="r")
            except Exception as exc:
                raise RuntimeError("%s : %s" % (self.path, str(exc)))

            self.ngroups = 1
            self._have_netcdf4 = False

        #self.path2group = collections.OrderedDict()
        #for children in self.walk_tree():
        #   for child in children:
        #       #print child.group,  child.path
        #       self.path2group[child.path] = child.group

    #@staticmethod
    #def join(*args):
    #    return "/".join([arg for arg in args])

    def __enter__(self):
        "Activated when used in the with statement."
        return self

    def __exit__(self, type, value, traceback):
        "Activated at the end of the with statement. It automatically close the file."
        self.rootgrp.close()

    def close(self):
        self.rootgrp.close()

    def walk_tree(self, top=None):
        """
        Navigate all the groups in the file starting from top.
        If top is None, the root group is used.
        """
        if top is None: top = self.rootgrp
        values = top.groups.values()
        yield values
        for value in top.groups.values():
            for children in walktree(value):
                yield children

    def print_tree(self, top=None):
        for children in self.walk_tree():
            for child in children:
                print(child)

    def get_varnames(self, path="/"):
        if path == "/":
            return self.rootgrp.variables.keys()
        else:
            group = self.path2group[path]
            return group.variables.keys()

    def get_value(self, varname, path="/"):
        var = self.get_variable(varname, path=path)
        # scalar or array
        return var[0] if not var.shape else var[:]

    def get_variable(self, varname, path="/"):
        return self.get_variables(varname, path=path)[0]

    def get_values(self, *varnames, **kwargs):
        vars = self.get_variables(*varnames, **kwargs)
        values = []
        # scalar or array
        for var in vars:
            v = var[0] if not var.shape else var[:]
            values.append(v)
        return values

    def get_variables(self, *varnames, **kwargs):
        path = kwargs.get("path","/")
        if path == "/":
            return [self.rootgrp.variables[vname] for vname in varnames]
        else:
            group = self.path2group[path]
            return [group.variables[vname] for vname in varnames]

##########################################################################################

class GSR_Reader(NetcdfReader):

    def get_structure(self, site_properties=None):
        if self.ngroups != 1:
            raise NotImplementedError("ngroups != 1")

        return structure_from_etsf_file(self, site_properties=site_properties)

    #def isconverged(self):
    #    "True if calculation is converged"

##########################################################################################

def structure_from_etsf_file(ncdata, site_properties=None):
    """
    Return a new instance from a NetCDF file containing crystallographic data in the ETSF-IO format.

    Args:
        ncdata: filename or NetcdfReader instance.
    """

    open_and_close = isinstance(ncdata, str)
    if open_and_close:
        ncdata = NetcdfReader(ncdata)

    # TODO check whether atomic units are used

    lattice = Bohr2Ang(ncdata.get_value("primitive_vectors"))

    red_coords = ncdata.get_value("reduced_atom_positions")
    natom = len(red_coords)

    znucl_type = ncdata.get_value("atomic_numbers")

    typat = ncdata.get_value("atom_species")

    species = np.zeros(natom, np.int)
    for atom in range(natom):
        # Fortran to C index and float --> int conversion.
        type_idx = typat[atom] - 1
        species[atom] = znucl_type[type_idx-1]

    d = {}
    if site_properties is not None:
        for property in site_properties:
                d[property] = ncdata.get_value(property)
    
    new_structure = Structure(lattice, species, red_coords, site_properties=d)

    if open_and_close:
        ncdata.close()

    return new_structure
