# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

"""
This module implements classes for generating Lammps input files.
The input files consist of the main input file with the control parameters and
the data file.
"""

import json
import os
from functools import partial
from collections import OrderedDict

from monty.json import MSONable, MontyDecoder

from pymatgen.io.lammps.data import LammpsData, LammpsForceFieldData

__author__ = "Kiran Mathew"
__credits__ = "Navnidhi Rajput"
__email__ = "kmathew@lbl.gov"

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))


class DictLammpsInput(MSONable):
    """
    Implementation of LammpsInputSet that is initialized from a dict
    settings. It is typically used by other LammpsInputSets for
    initialization from json or yaml source files.

    Args:
        name (str): A name for the input set.
        config_dict (dict): The config dictionary to use.
        lammps_data (LammpsData): LammpsData object
        data_filename (str): name of the the lammps data file
        user_lammps_settings (dict): User lammps settings. This allows a user
            to override lammps settings, e.g., setting a different force field
            or bond type.
    """

    def __init__(self, name, config_dict, lammps_data=None,
                 data_file="in.data", user_lammps_settings={}):
        self.name = name
        self.lines = []
        self.config_dict = config_dict
        self.lammps_data = lammps_data
        self.data_file = data_file
        self.config_dict["read_data"] = data_file
        self.user_lammps_settings = user_lammps_settings
        if self.user_lammps_settings:
            self.config_dict.update(self.user_lammps_settings)

    def __str__(self):
        """
        string representation of the lammps input file with the
        control parameters
        """
        lines = ""
        for k1, v1 in self.config_dict.items():
            if isinstance(v1, dict):
                v1 = v1.values()
            if isinstance(v1, list):
                for x in v1:
                    lines = "".join([lines, "{} ".format(k1)])
                    lines = "".join([lines, str(x), os.linesep])
            else:
                lines = "".join([lines, "{} ".format(k1)])
                lines = "".join([lines, " {}{}".format(str(v1), os.linesep)])
        return lines

    def write_input(self, filename, data_file=None):
        """
        get the string representation of the main input file and write it.
        Also writes the data file if the lammps_data attribute is set.

        Args:
            filename (string): name of the input file
        """
        if data_file:
            self.config_dict["read_data"] = data_file
            self.data_file = data_file
        # write the main input file
        with open(filename, 'w') as f:
            f.write(self.__str__())
        # write the data file if present
        if self.lammps_data:
            print("Writing the data to {}".format(self.data_file))
            self.lammps_data.write_data_file(filename=self.data_file)

    @staticmethod
    def from_file(name, filename, data_obj=None, data_file=None,
                  is_forcefield=False):
        """
        Read in the input settings from json file as ordereddict. Also
        reads in the datafile if provided.
        Note: with monty.serialization.loadfn the order of paramters in the
        json file is not preserved

        Args:
            filename (string): name of the file with the lamps control
                paramters
            data_obj (LammpsData): LammpsData object
            data_file (string): name of the data file name
            is_forcefield (bool): whether the data file has forcefield and
                topology info in it.

        Returns:
            DictLammpsInput
        """
        with open(filename) as f:
            config_dict = json.load(f, object_pairs_hook=OrderedDict)
        lammps_data = None
        if data_file:
            if is_forcefield:
                lammps_data = LammpsForceFieldData.from_file(data_file)
            else:
                lammps_data = LammpsData.from_file(data_file)
        if data_obj and isinstance(data_obj, LammpsData):
            lammps_data = data_obj
        return DictLammpsInput(name, config_dict, lammps_data)

    def as_dict(self):
        d = MSONable.as_dict(self)
        if hasattr(self, "kwargs"):
            d.update(**self.kwargs)
        d["config_dict"] = self.config_dict.items()
        return d

    @classmethod
    def from_dict(cls, d):
        decoded = {k: MontyDecoder().process_decoded(v) for k, v in d.items()
                   if k not in ["@module", "@class", "config_dict"]}
        decoded["config_dict"] = OrderedDict(d["config_dict"])
        return cls(**decoded)


# NVT
NVTLammpsInput = partial(DictLammpsInput.from_file, "NVT",
                         os.path.join(MODULE_DIR, "NVT.json"))

# NPT
NPTLammpsInput = partial(DictLammpsInput.from_file, "NPT",
                         os.path.join(MODULE_DIR, "NPT.json"))

# NPT followed by NVT
NPTNVTLammpsInput = partial(DictLammpsInput.from_file, "NPT_NVT",
                            os.path.join(MODULE_DIR, "NPT_NVT.json"))
