"""Wrapper for netCDF readers."""
from __future__ import division, print_function

import numpy as np
import os.path 
from collections import OrderedDict

try:
    import scipy.io.netcdf as nc
except ImportError:
    pass

__author__ = "Matteo Giantomassi"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"
__email__ = "gmatteo at gmail.com"
__status__ = "Development"
__date__ = "$Feb 21, 2013M$"

#__all__ = [
#"NetcdfReader",
#"GSR_Reader",
#]

##########################################################################################

def wrap_ncreader(fname):
    return nc.netcdf_file(fname)

##########################################################################################

def ncread_keys(ncfile):
    """
    Return the keywords defined in the netCDF file ncfile.

    :arg ncfile: string with the filename or NC file-object.
    :return: (dim_names, var_name) list with the names of the dimensions and of the variables found in ncfile.
    """

    do_close = False
    if isinstance(ncfile, str):
        ncdata = wrap_ncreader(ncfile)
        do_close = True
    else:
        ncdata = ncfile

    dim_names = ncdata.dimensions.keys()
    var_names = ncdata.variables.keys()

    if do_close: ncdata.close()
    return dim_names, var_names

##########################################################################################

def ncread_key(ncfile, key):
    """
    Return the value of the keyword defined in the netCDF file ncfile.

    :arg ncfile: string with the filename or NC file-object.
    :arg key: NC keyword
    :return: the value of key, raise KeyError if key is not in ncfile.
    """

    do_close = False
    if isinstance(ncfile, str):
        ncdata = wrap_ncreader(ncfile)
        do_close = True
    else:
        ncdata = ncfile

    varobj = None
    if key in ncdata.dimensions:
        varobj = ncdata.dimensions[key]

    elif key in ncdata.variables:
        try:
            varobj   = ncdata.variables[key][:] # Read array.
            varshape = np.shape(varobj)
            varobj   = np.reshape(np.array(varobj), varshape)

        except IndexError:
            #print("IndexError for key: ", key," type: ",type(varobj))
            varobj = ncdata.variables[key].getValue() # Read scalar variable.

    if do_close: ncdata.close()
    if varobj is None: raise KeyError
    return varobj

#########################################################################################

def ncread_varsdims(obj, ncfile, names, map_names=None, overwrite=False):
    """
    Read names from the netCDF file ncfile.

    :arg obj: object with a __dict__ dictionary.
    :arg ncfile: string with the filename or NC file-object.
    :arg names: list of netCDF keywords to read. Values are stored in obj.name for name in names.
    :arg map_names: dictionary used to map name to the netCDF keyword used to access data on file.
    :arg overwrite: if overwrite==False and obj.name already exists, ValueError is raised.
    :return missing: list of tuple with the keywords that are not found in ncfile.
    """
    if map_names is None: map_names = {}

    do_close = False
    if isinstance(ncfile, str):
        ncdata = wrap_ncreader(ncfile)
        do_close = True
    else:
        ncdata = ncfile

    #print("DIMS ",ncdata.dimensions.keys())
    #print("VARS ",ncdata.variables.keys())
    #
    # Add nc dimensions and variables to obj.__dict__
    missing = list()
    for k in names:

        try:
            etsf_k = map_names[k]
        except KeyError:
            etsf_k = k  # Read k.

        if not overwrite and k in obj.__dict__:
            raise ValueError("Cannot overwrite key: "+ str(k))

        if etsf_k in ncdata.dimensions:
            obj.__dict__[k] = ncdata.dimensions[etsf_k]

        elif etsf_k in ncdata.variables:
            try:
                varobj = ncdata.variables[etsf_k][:] # Read array.
                varshape = np.shape(varobj)
                obj.__dict__[k] = np.reshape(np.array(varobj), varshape)

            except IndexError:
                #print("IndexError for key: %s, type: %s" % (etsf_k, type(varobj)))
                varobj = ncdata.variables[etsf_k].getValue() # Read scalar variable.
                obj.__dict__[k] = varobj

        else:
            missing.append((k,etsf_k))

    if do_close: ncdata.close()

    return missing

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

        #self.path2group = OrderedDict()
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
        for children in gsfile.walk_tree():
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

    def get_structure(self):
        if self.ngroups != 1:
            raise NotImplementedError("ngroups != 1")
        return Structure.from_etsf_file(self.rootgrp)

    #def isconverged(self):

##########################################################################################

if __name__ == "__main__":
    filename = "hello_crystal.nc"
    with GSR_Reader(filename) as gsr_file:
    #print gsfile.variables
    #print dir(gsfile.variables["etotal"])
        print(gsr_file.get_varnames())
        gsr_file.print_tree()
        etotal, cforces = gsr_file.get_values("etotal", "cartesian_forces")
        print(etotal, cforces)
        etotal, fermie = gsr_file.get_values("etotal", "fermie")
        print("etotal %s, fermi %s" % (etotal, fermie))
        #print etotal,type(etotal)
        #print gsr_file.etotal
        #print gsr_file.getncattr("etotal")
        for group in gsr_file.walk_tree():
            print("group: " + str(group))
