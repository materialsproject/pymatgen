# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License

"""
This module implements reading and writing of ShengBTE CONTROL files.
"""

import warnings
from typing import Any, Dict, List, Optional, Union

import numpy as np
from monty.dev import requires
from monty.json import MSONable

from pymatgen.core.structure import Structure
from pymatgen.io.vasp import Kpoints

try:
    import f90nml
except ImportError:
    f90nml = None

__author__ = "Rees Chang, Alex Ganose"
__copyright__ = "Copyright 2019, The Materials Project"
__version__ = "0.1"
__email__ = "rc564@cornell.edu, aganose@lbl.gov"
__date__ = "June 27, 2019"


class Control(MSONable, dict):
    """
    Class for reading, updating, and writing ShengBTE CONTROL files.
    See  https://bitbucket.org/sousaw/shengbte/src/master/ for more
    detailed description and default values of CONTROL arguments.
    """

    required_params = [
        "nelements",
        "natoms",
        "ngrid",
        "lattvec",
        "types",
        "elements",
        "positions",
        "scell",
    ]
    allocations_keys = ["nelements", "natoms", "ngrid", "norientations"]
    crystal_keys = [
        "lfactor",
        "lattvec",
        "types",
        "elements",
        "positions",
        "masses",
        "gfactors",
        "epsilon",
        "born",
        "scell",
        "orientations",
    ]
    params_keys = [
        "t",
        "t_min",
        "t_max",
        "t_step",
        "omega_max",
        "scalebroad",
        "rmin",
        "rmax",
        "dr",
        "maxiter",
        "nticks",
        "eps",
    ]
    flags_keys = [
        "nonanalytic",
        "convergence",
        "isotopes",
        "autoisotopes",
        "nanowires",
        "onlyharmonic",
        "espresso",
    ]

    def __init__(self, ngrid: Optional[List[int]] = None, temperature: Union[float, Dict[str, float]] = 300, **kwargs):
        """
        Args:
            ngrid: Reciprocal space grid density as a list of 3 ints.
            temperature: The temperature to calculate the lattice thermal
                conductivity for. Can be given as a single float, or a dictionary
                with the keys "min", "max", "step".
            **kwargs: Other ShengBTE parameters. Several parameters are required
                for ShengBTE to run - we have listed these parameters below:

                - nelements (int): number of different elements in the compound
                - natoms (int): number of atoms in the unit cell
                - lattvec (size 3x3 array): real-space lattice vectors, in units
                  of lfactor
                - lfactor (float): unit of measurement for lattice vectors (nm).
                    I.e., set to 0.1 if lattvec given in Angstrom.
                - types (size natom list): a vector of natom integers, ranging
                  from 1 to nelements, assigning an element to each atom in the
                  system
                - elements (size natom list): a vector of element names
                - positions (size natomx3 array): atomic positions in lattice
                  coordinates
                - scell (size 3 list): supercell sizes along each crystal axis
                  used for the 2nd-order force constant calculation
        """
        super().__init__()
        if ngrid is None:
            ngrid = [25, 25, 25]

        self["ngrid"] = ngrid

        if isinstance(temperature, (int, float)):
            self["t"] = temperature

        elif isinstance(temperature, dict):
            self["t_min"] = temperature["min"]
            self["t_max"] = temperature["max"]
            self["t_step"] = temperature["step"]
        else:
            raise ValueError("Unsupported temperature type, must be float or dict")

        self.update(kwargs)

    @classmethod
    @requires(
        f90nml,
        "ShengBTE Control object requires f90nml to be installed. Please get it at https://pypi.org/project/f90nml.",
    )
    def from_file(cls, filepath: str):
        """
        Read a CONTROL namelist file and output a 'Control' object

        Args:
            filepath: Path of the CONTROL file.

        Returns:
            'Control' object with parameters instantiated.
        """
        nml = f90nml.read(filepath)
        sdict = nml.todict()

        all_dict: Dict[str, Any] = {}
        all_dict.update(sdict["allocations"])
        all_dict.update(sdict["crystal"])
        all_dict.update(sdict["parameters"])
        all_dict.update(sdict["flags"])
        all_dict.pop("_start_index")  # remove unnecessary cruft

        return cls.from_dict(all_dict)

    @classmethod
    def from_dict(cls, control_dict: Dict):
        """
        Write a CONTROL file from a Python dictionary. Description and default
        parameters can be found at
        https://bitbucket.org/sousaw/shengbte/src/master/.
        Note some parameters are mandatory. Optional parameters default here to
        None and will not be written to file.

        Args:
            control_dict: A Python dictionary of ShengBTE input parameters.
        """
        return cls(**control_dict)

    @requires(
        f90nml,
        "ShengBTE Control object requires f90nml to be installed. Please get it at https://pypi.org/project/f90nml.",
    )
    def to_file(self, filename: str = "CONTROL"):
        """
        Writes ShengBTE CONTROL file from 'Control' object

        Args:
            filename: A file name.
        """

        for param in self.required_params:
            if param not in self.as_dict():
                warnings.warn(f"Required parameter '{param}' not specified!")

        alloc_dict = _get_subdict(self, self.allocations_keys)
        alloc_nml = f90nml.Namelist({"allocations": alloc_dict})
        control_str = str(alloc_nml) + "\n"

        crystal_dict = _get_subdict(self, self.crystal_keys)
        crystal_nml = f90nml.Namelist({"crystal": crystal_dict})
        control_str += str(crystal_nml) + "\n"

        params_dict = _get_subdict(self, self.params_keys)
        params_nml = f90nml.Namelist({"parameters": params_dict})
        control_str += str(params_nml) + "\n"

        flags_dict = _get_subdict(self, self.flags_keys)
        flags_nml = f90nml.Namelist({"flags": flags_dict})
        control_str += str(flags_nml) + "\n"

        with open(filename, "w") as file:
            file.write(control_str)

    @classmethod
    def from_structure(cls, structure: Structure, reciprocal_density: Optional[int] = 50000, **kwargs):
        """
        Get a ShengBTE control object from a structure.

        Args:
            structure: A structure object.
            reciprocal_density: If not None, the q-point grid ("ngrid") will be
                set using this density.
            kwargs: Additional options to be passed to the Control constructor.
                See the docstring of the __init__ method for more details

        Returns:
            A ShengBTE control object.
        """

        elements = list(map(str, structure.composition.elements))

        unique_nums = np.unique(structure.atomic_numbers)
        types_dict = dict(zip(unique_nums, range(len(unique_nums))))
        types = [types_dict[i] + 1 for i in structure.atomic_numbers]

        control_dict = {
            "nelements": structure.ntypesp,
            "natoms": structure.num_sites,
            "norientations": 0,
            "lfactor": 0.1,
            "lattvec": structure.lattice.matrix.tolist(),
            "elements": elements,
            "types": types,
            "positions": structure.frac_coords.tolist(),
        }

        if reciprocal_density:
            kpoints = Kpoints.automatic_density(structure, reciprocal_density)
            control_dict["ngrid"] = kpoints.kpts[0]

        control_dict.update(**kwargs)

        return Control(**control_dict)

    def get_structure(self) -> Structure:
        """
        Get a pymatgen Structure from a ShengBTE control object.

        The control object must have the "lattvec", "types", "elements", and
        "positions" settings otherwise an error will be thrown.

        Returns:
            The structure.
        """
        required = ["lattvec", "types", "elements", "positions"]
        if not all(r in self for r in required):
            raise ValueError("All of ['lattvec', 'types', 'elements', 'positions'] must be in control object")

        unique_elements = self["elements"]
        n_unique_elements = len(unique_elements)
        element_map = dict(zip(range(1, n_unique_elements + 1), unique_elements))
        species = [element_map[i] for i in self["types"]]

        cell = np.array(self["lattvec"])

        if "lfactor" in self:
            cell *= self["lfactor"] * 10  # to nm then to Angstrom

        return Structure(cell, species, self["positions"])

    def as_dict(self):
        """
        Returns: MSONAble dict
        """
        return dict(self)


def _get_subdict(master_dict, subkeys):
    """Helper method to get a set of keys from a larger dictionary"""
    return {k: master_dict[k] for k in subkeys if k in master_dict and master_dict[k] is not None}
