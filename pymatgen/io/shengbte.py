# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License
import warnings

from monty.dev import requires
from monty.json import MSONable

try:
    import f90nml
except:
    f90nml = None

"""
This module defines IO ShengBTE for reading, updating, and writing the
CONTROL input file
"""

__author__ = "Rees Chang"
__copyright__ = "Copyright 2019, The Materials Project"
__version__ = "0.1"
__email__ = "rc564@cornell.edu"
__date__ = "June 27, 2019"


class Control(dict, MSONable):
    """
    Class for reading, updating, and writing ShengBTE CONTROL files.
    See  https://bitbucket.org/sousaw/shengbte/src/master/ for more
    detailed description and default values of CONTROL arguments.

    Args:
        ngrid (size 3 list): number of grid planes along each axis in
            reciprocal space
        lfactor (float): unit of measurement for lattice vectors (nm)
        scalebroad (float): scale parameter for Gaussian smearing. A value
            of 1.0 is theoretically guaranteed to work, but significant
            speedups can sometimes be achieved by reducing it with
            negligible loss of precision.
        t (int or float): temperature (Kelvin)
        **kwargs: Other ShengBTE parameters. Several parameters are required
            for ShengBTE to run - we have listed these parameters below:
            - nelements (int): number of different elements in the compound
            - natoms (int): number of atoms in the unit cell
            - lattvec (size 3x3 array): real-space lattice vectors, in units
              of lfactor
            - types (size natom list): a vector of natom integers, ranging
              from 1 to nelements, assigning an element to each atom in the
              system
            - elements (size natom list): a vector of element names
            - positions (size natomx3 array): atomic positions in lattice
              coordinates
            - scell (size 3 list): supercell sizes along each crystal axis
              used for the 2nd-order force constant calculation
    """

    required_params = ["nelements", "natoms", "ngrid", "lattvec", "types",
                       "elements", "positions", "scell"]
    allocations_keys = ["nelements", "natoms", "ngrid", "norientations"]
    crystal_keys = ["lfactor", "lattvec", "types", "elements", "positions",
                    "masses", "gfactors", "epsilon", "born", "scell",
                    "orientations"]
    params_keys = ["t", "t_min", "t_max", "t_step", "omega_max", "scalebroad",
                   "rmin", "rmax", "dr", "maxiter", "nticks", "eps"]
    flags_keys = ["nonanalytic", "convergence", "isotopes", "autoisotopes",
                  "nanowires", "onlyharmonic", "espresso"]

    @requires(f90nml,
              "ShengBTE Control object requires f90nml to be installed. "
              "Please get it at https://pypi.org/project/f90nml.")
    def __init__(self, ngrid=None, lfactor=0.1,
                 scalebroad=0.5, t=300, **kwargs):
        if ngrid is None:
            ngrid = [25, 25, 25]

        self["ngrid"] = ngrid
        self["lfactor"] = lfactor
        self["scalebroad"] = scalebroad
        self["t"] = t
        self.update(kwargs)

    @classmethod
    def from_file(cls, filepath):
        """
        Read a CONTROL namelist file and output a 'Control' object

        Args:
            filepath (String): Path of the CONTROL file.

        Returns:
            'Control' object with parameters instantiated.
        """
        nml = f90nml.read(filepath)
        sdict = nml.todict()

        all_dict = {}
        all_dict.update(sdict["allocations"])
        all_dict.update(sdict["crystal"])
        all_dict.update(sdict["parameters"])
        all_dict.update(sdict["flags"])

        return cls.from_dict(all_dict)

    @classmethod
    def from_dict(cls, sdict):
        """
        Write a CONTROL file from a Python dictionary.
        Description and default parameters can be found at
        https://bitbucket.org/sousaw/shengbte/src/master/.
        Note some parameters are mandatory. Optional parameters
        default here to None and will not be written to file.

        Args:
            dict: A Python dictionary of ShengBTE input parameters.
        """
        return cls(**sdict)

    def to_file(self, filename='CONTROL'):
        """
        Writes ShengBTE CONTROL file from 'Control' object
        """

        for param in self.required_params:
            if param not in self.as_dict():
                warnings.warn(
                    "Required parameter '{}' not specified!".format(param))

        alloc_dict = {k: self[k] for k in self.allocations_keys
                      if k in self and self[k] is not None}
        alloc_nml = f90nml.Namelist({"allocations": alloc_dict})
        control_str = str(alloc_nml) + "\n"

        crystal_dict = {k: self[k] for k in self.crystal_keys
                      if k in self and self[k] is not None}
        crystal_nml = f90nml.Namelist({"crystal": crystal_dict})
        control_str += str(crystal_nml) + "\n"

        params_dict = {k: self[k] for k in self.params_keys
                      if k in self and self[k] is not None}
        params_nml = f90nml.Namelist({"parameters": params_dict})
        control_str += str(params_nml) + "\n"

        flags_dict = {k: self[k] for k in self.flags_keys
                      if k in self and self[k] is not None}
        flags_nml = f90nml.Namelist({"flags": flags_dict})
        control_str += str(flags_nml)

        with open(filename, "w") as file:
            file.write(control_str)

    def as_dict(self):
        return dict(self)
