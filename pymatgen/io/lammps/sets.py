# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module implements classes for reading and generating Lammps inputset.

For the ease of management we divide LAMMPS input into 2 files:

    1.Data file: All structure related settings such as the atomic positions,
            bonds, angles, dihedrals, corresponding parametrizations etc are
            set in the data file.

    2. Control/input file: This is the main input file that should be fed to the
            lammps binary. The main input file consists of the path to the
            afore-mentioned data file and the job control parameters such as
            the ensemble type(NVT, NPT etc), max number of iterations etc.
"""

import os
import six

from monty.json import MSONable, MontyDecoder

from pymatgen.io.lammps.data import LammpsData
from pymatgen.io.lammps.input import LammpsInput

__author__ = "Kiran Mathew"
__email__ = "kmathew@lbl.gov"


class LammpsInputSet(MSONable):

    def __init__(self, name, lammps_input, lammps_data=None,
                 data_filename="in.data", user_lammps_settings=None):
        """
        Implementation of LammpsInputSet that is initialized from a dict
        settings. It is typically used by other LammpsInputSets for
        initialization from json or yaml source files.

        Args:
            name (str): A name for the input set.
            lammps_input (LammpsInput): The config dictionary to use.
            lammps_data (LammpsData): LammpsData object
            data_filename (str): name of the the lammps data file.
                Note: this will override the value for 'data_file' key in lammps_input
            user_lammps_settings (dict): User lammps settings. This allows a user
                to override lammps settings, e.g., setting a different force field
                or bond type.
        """
        self.name = name
        self.lines = []
        self.lammps_input = lammps_input
        self.lammps_data = lammps_data
        self.data_filename = data_filename
        self.lammps_input.settings["data_file"] = data_filename
        self.user_lammps_settings = user_lammps_settings or {}
        self.lammps_input.settings.update(self.user_lammps_settings)

    def write_input(self, input_filename, data_filename=None):
        """
        Get the string representation of the main input file and write it.
        Also writes the data file if the lammps_data attribute is set.

        Args:
            input_filename (string): name of the input file
            data_filename (string): override the data file name with this
        """
        if data_filename:
            data_filename = os.path.abspath(os.path.join(os.getcwd(), data_filename))
        if data_filename and ("data_file" in self.lammps_input.settings):
            self.lammps_input.settings["data_file"] = data_filename
            self.data_filename = data_filename
        self.lammps_input.write_file(input_filename)
        # write the data file if present
        if self.lammps_data:
            self.lammps_data.write_file(filename=self.data_filename)

    @classmethod
    def from_file(cls, name, input_template, user_settings,
                  lammps_data=None, data_filename="in.data"):
        """
        Returns LammpsInputSet from  input file template and input data.

        Args:
            name (str)
            input_template (string): path to the input template file.
            user_settings (dict): User lammps settings, the keys must
                correspond to the keys in the template.
            lammps_data (string/LammpsData): path to the
                data file or an appropriate object
            data_filename (string): name of the the lammps data file.

        Returns:
            LammpsInputSet
        """
        user_settings["data_file"] = data_filename
        lammps_input = LammpsInput.from_file(input_template, user_settings)
        if isinstance(lammps_data, six.string_types):
            lammps_data = LammpsData.from_file(lammps_data)
        return cls(name, lammps_input, lammps_data=lammps_data,
                   data_filename=data_filename)

    def as_dict(self):
        d = MSONable.as_dict(self)
        if hasattr(self, "kwargs"):
            d.update(**self.kwargs)
        d["lammps_input"] = self.lammps_input.as_dict()
        return d

    @classmethod
    def from_dict(cls, d):
        decoded = {k: MontyDecoder().process_decoded(v) for k, v in d.items()
                   if k not in ["@module", "@class", "lammps_input"]}
        decoded["lammps_input"] = LammpsInput.from_dict(d["lammps_input"])
        return cls(**decoded)
