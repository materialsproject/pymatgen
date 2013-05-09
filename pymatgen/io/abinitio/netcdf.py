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

###############################################################################


class NetcdfReader(object):
    "Wraps and extends netCDF4.Dataset. Read only mode"

    @requires(netCDF4 is not None,
              "netCDF4 library must be installed to use this class")
    def __init__(self, filename):
        self.path = os.path.abspath(filename)

        try:
            self.rootgrp = netCDF4.Dataset(self.path, mode="r")
        except Exception as exc:
            raise RuntimeError("%s : %s" % (self.path, str(exc)))

        self.ngroups = len( list(self.walk_tree()) )

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
        "Activated at the end of the with statement. It automatically closes the file."
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
            for children in self.walk_tree(value):
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

################################################################################
class GSR_Reader(NetcdfReader):
    "Netcdf reader for the Ground-State Results file."

    def get_structure(self, site_properties=None):
        if self.ngroups != 1:
            raise NotImplementedError("ngroups != 1")

        return structure_from_etsf_file(self, site_properties=site_properties)

    def get_band_structure(self):
        raise NotImplementedError("")
        structure = self.get_structure()
        from pprint import pprint

        kpoints = self.get_value("reduced_coordinates_of_kpoints")
        efermi = Ha2eV(self.get_value("fermie"))
        np_eigvals = Ha2eV(self.get_value("eigenvalues"))
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

        #bands = BandStructure(kpoints, eigenvals, structure.lattice, efermi, labels_dict=None, structure=structure)

        bands = BandStructureSymmLine(kpoints, eigenvals, structure.lattice,
                                      efermi, labels_dict, structure=structure)
        return bands

##########################################################################################

def structure_from_etsf_file(ncdata, site_properties=None):
    """
    Return a new instance from a NetCDF file containing crystallographic data
    in the ETSF-IO format.

    Args:
        ncdata:
            filename or NetcdfReader instance.
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

    species = natom * [None,]

    for atom in range(natom):
        # Fortran to C index and float --> int conversion.
        type_idx = typat[atom] - 1
        species[atom] = int(znucl_type[type_idx-1])

    d = {}
    if site_properties is not None:
        for property in site_properties:
            d[property] = ncdata.get_value(property)

    new = Structure(lattice, species, red_coords, site_properties=d)

    if open_and_close:
        ncdata.close()

    return new
