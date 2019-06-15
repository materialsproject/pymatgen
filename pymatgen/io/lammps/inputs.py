# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
import abc
import re
import os
import glob
import shutil
import warnings
from itertools import combinations
from copy import deepcopy
import numpy as np
from pathlib import Path
from monty.serialization import loadfn
from monty.io import zopen
import json
from monty.dev import deprecated
from string import Template

from monty.json import MSONable
from pymatgen.io.lammps.data import LammpsData
from pymatgen.io.lammps.utils import Pair
from pymatgen.core.periodic_table import Element

"""
This module implements methods for writing LAMMPS input files.

"""

__author__ = "Kiran Mathew, Brandon Wood, Zhi Deng"
__copyright__ = "Copyright 2018, The Materials Virtual Lab"
__version__ = "1.0"
__maintainer__ = "Zhi Deng"
__email__ = "z4deng@eng.ucsd.edu"
__date__ = "Aug 1, 2018"

with open(str(Path(__file__).absolute().parent / "forcefield.json"), "rt") as f:
    _forcefields = json.load(f)


class LammpsRun(MSONable):
    """
    Examples for various simple LAMMPS runs with given simulation box,
    force field and a few more settings. Experience LAMMPS users should
    consider using write_lammps_inputs method with more sophisticated
    templates.

    """

    template_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "templates")

    def __init__(self, script_template, settings, data,
                 script_filename):
        """
        Base constructor.

        Args:
            script_template (str): String template for input script
                with placeholders. The format for placeholders has to
                be '$variable_name', e.g., '$temperature'
            settings (dict): Contains values to be written to the
                placeholders, e.g., {'temperature': 1}.
            data (LammpsData or str): Data file as a LammpsData
                instance or path to an existing data file. Default to
                None, i.e., no data file supplied. Useful only when
                read_data cmd is in the script.
            script_filename (str): Filename for the input script.

        """
        self.script_template = script_template
        self.settings = settings
        self.data = data
        self.script_filename = script_filename

    def write_inputs(self, output_dir, **kwargs):
        """
        Writes all input files (input script, and data if needed).
        Other supporting files are not handled at this moment.

        Args:
            output_dir (str): Directory to output the input files.
            **kwargs: kwargs supported by LammpsData.write_file.

        """
        write_lammps_inputs(output_dir=output_dir,
                            script_template=self.script_template,
                            settings=self.settings, data=self.data,
                            script_filename=self.script_filename, **kwargs)

    @classmethod
    def md(cls, data, force_field, temperature, nsteps,
           other_settings=None):
        """
        Example for a simple MD run based on template md.txt.

        Args:
            data (LammpsData or str): Data file as a LammpsData
                instance or path to an existing data file.
            force_field (str): Combined force field related cmds. For
                example, 'pair_style eam\npair_coeff * * Cu_u3.eam'.
            temperature (float): Simulation temperature.
            nsteps (int): No. of steps to run.
            other_settings (dict): other settings to be filled into
                placeholders.

        """
        template_path = os.path.join(cls.template_dir, "md.txt")
        with open(template_path) as f:
            script_template = f.read()
        settings = other_settings.copy() if other_settings is not None else {}
        settings.update({'force_field': force_field,
                         "temperature": temperature, "nsteps": nsteps})
        script_filename = "in.md"
        return cls(script_template=script_template,
                   settings=settings, data=data, script_filename=script_filename)


class LammpsInput(MSONable, metaclass=abc.ABCMeta):

    def __init__(self, commands, config_dict={}, lammps_data=None):
        self.commands    = commands
        self.config_dict = config_dict
        self.lammps_data = lammps_data

    @classmethod
    def from_file(cls, filename):
        with zopen(filename=filename) as f:
            return LammpsInput(commands=f.read(), config_dict={}, lammps_data=None)


class LammpsInputSet(MSONable):

    """
    Class representing the input information necessary to run Lammps. This is comprised of an input object
    (representing the lammps input file), optional configuration dictionary, and optional lammps data file
    for input.
    """

    def __init__(self, lammps_input, config_dict={}, lammps_data=None, **kwargs):
        """

        Base constructor.

        Args:
            lammps_input: (LammpsInput or str) the lammps input object or input file path to use. If read from an
                        input file, that file can, and often should, be a template file. In the template file,
                        all variables should be specified with a dollar sign, followed by a unique name.

                        Example: write_data $output_data. Where output_data is the variable defined in config_dict

            config_dict: (dict) a dictionary containing the variables to substitute into a template LammpsInput.
                        Default: {}

            lammps_data: (LammpsData or str) the lammps data object to use. Can either be given as an initialized
                        LammpsData object, or you can supply a filepath to read-in the data. This does not need to be
                        given. If you, for example, have your lammps data file in the same directory as your execution
                        and would simply like to use a variable from the config_dict to let lammps read it,
                        this works fine.
        """
        if isinstance(lammps_input, str):
            if os.path.isfile(lammps_input):
                self.lammps_input = LammpsInput.from_file(filename=lammps_input)
        elif isinstance(lammps_input, dict):
            self.lammps_input = LammpsInput.from_dict(lammps_input)
        else:
            self.lammps_input = lammps_input

        self.config_dict = deepcopy(config_dict)

        self.lammps_data = lammps_data
        if isinstance(lammps_data, str):
            if os.path.isfile(lammps_data):
                atom_style = config_dict.get('atom_style', 'full')
                self.lammps_data = LammpsData.from_file(lammps_data, atom_style=atom_style)
        elif isinstance(lammps_data, dict):
            self.lammps_data = LammpsData.from_dict(lammps_data)
        elif lammps_data is None:
            try:
                atom_style = config_dict.get('atom_style', 'full')
                self.lammps_data = LammpsData.from_file(os.path.join(os.getcwd(),'lammps.data'), atom_style=atom_style)
            except:
                print('Could not located any lammps data')
        self.kwargs      = kwargs

    def write_input(self, input_filename="in.lammps", output_dir="", data_filename="lammps.data"):

        """
        Write the inputs necessary to run Lammps. This usually just means the input file, but could also mean
        a data file for lammps to read.

        Args:
            input_filename: The name of the input file to write. Default: "in.lammps"
            output_dir: Output directory. Defailt: "" e.g: the current directory
            data_filename:
        Return:
            None
        """
        script_template = self.lammps_input.commands
        variables       = {} if self.config_dict is None else self.config_dict

        input_script = script_template

        read_data = bool(re.search(r"read_data\s+(.*)\n", script_template))
        if read_data:
            if isinstance(self.lammps_data, LammpsData):
                print("FOUND LAMMPS DATA OBJECT")
                distance = self.kwargs.get('distance', 6)
                velocity = self.kwargs.get('velocity', 8)
                charge = self.kwargs.get('charge', 3)
                self.lammps_data.write_file(os.path.join(output_dir, data_filename),
                                            distance=distance, velocity=velocity, charge=charge)
            elif isinstance(self.lammps_data, str) and os.path.exists(os.path.join(output_dir, self.lammps_data)):
                print("Located LAMMPS Data File...")
                self.lammps_data = LammpsData.from_file(os.path.join(output_dir, self.lammps_data),
                                                        atom_style=self.config_dict.get('atom_style', 'full'))
                #shutil.copyfile(self.lammps_data, os.path.join(output_dir, data_filename))
            else:
                warnings.warn("No data file supplied. Skip writing.")

        for key, value in variables.items():
            if key == 'pair_coeffs':
                sub = self._get_pair_coeffs(value)
            else:
                sub = value
            var = "$" + str(key)
            regex = re.escape(var)
            input_script = re.sub(regex, str(sub), input_script)

        if bool(output_dir) and not os.path.exists(output_dir):
            os.makedirs(output_dir)

        with open(os.path.join(output_dir, input_filename), "w") as f:
            f.write(input_script)

    @classmethod
    def from_file(cls, filename, config_dict, lammps_data):
        input = LammpsInput.from_file(filename=filename)
        if isinstance(lammps_data, str):
            lammps_data = LammpsData.from_file(lammps_data, atom_style=config_dict.get('atom_style', 'full'))
        return LammpsInputSet(input, config_dict=config_dict, lammps_data=lammps_data)

    def _get_pair_coeffs(self, pairs):

        try:
            masses = self.lammps_data.masses.copy()
            unique_indices = np.unique(masses['mass'], return_index=True)[1]
            unique_masses  = [masses['mass'].values[index] for index in sorted(unique_indices)]

            ref_masses = [el.atomic_mass.real for el in Element]
            diff = [np.abs(np.array(ref_masses) - mass) for mass in unique_masses]
            atomic_numbers = [np.argmin(d) + 1 for d in diff]
            symbols = [Element.from_Z(an).symbol for an in atomic_numbers]

            key = []
            numbers = []
            for i, s in enumerate(symbols):
                key.append(s)
                numbers.append(i)
        except:
            key = self.config_dict.get('key')
            numbers = [v-1 for k, v in key.items()]
            key = [k for k,v in key.items()]

        pair_numbers = (list(combinations(numbers, 2)) + [(i,i) for i in numbers])
        pair_numbers.sort(key=lambda x: x[0])
        string = ''
        for p in pair_numbers:
            string += 'pair_coeff '+str(p[0]+1)+' '+str(p[1]+1)+' '
            try:
                pair = Pair(key[p[0]], key[p[1]])
                params = _forcefields[str(pair.get_id())]['params']
            except:
                pair = Pair(key[p[1]], key[p[0]])
                params = _forcefields[str(pair.get_id())]['params']
            for k, v in params.items():
                string += str(v) + ' '
            string += '\n'

        return string
