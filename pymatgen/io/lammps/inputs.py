# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module implements methods for writing LAMMPS input files.
"""

import os
import re
import shutil
import warnings
from collections import OrderedDict
from pathlib import Path
from string import Template
from typing import Dict, Optional, Union

from monty.dev import deprecated
from monty.io import zopen
from monty.json import MSONable

from pymatgen.io.core import InputFile
from pymatgen.io.lammps.data import LammpsData  # , CombinedData
from pymatgen.io.template import TemplateInputGen

__author__ = "Kiran Mathew, Brandon Wood, Zhi Deng"
__copyright__ = "Copyright 2018, The Materials Virtual Lab"
__version__ = "1.0"
__maintainer__ = "Zhi Deng"
__email__ = "z4deng@eng.ucsd.edu"
__date__ = "Aug 1, 2018"


class LammpsInputFile(InputFile):
    """
    Class representing a LAMMPS input settings file, e.g. in.lammps.
    Allows for LAMMPS input generation in line/stage wise manner. A stage
    here is defined as a block of LAMMPS input commands usually performing a
    specific task during the simulation such as energy minimization or
    NPT/NVT runs. But more broadly, a stage can also be block of LAMMPS
    input where simulation box is setup, a set of variables are declared or
    quantities are computed. A comment beginning with '#' is treated as a
    header for a stage and marks the start of a new stage in the LAMMPS input.
    Other comments starting with '##' are treated as conventional comments.
    """

    def __init__(self, input_settings: OrderedDict = None):
        """
        Args:
            input_settings: Ordered Dictionary of LAMMPS input settings.
        """

        self.nstages = 0
        self.ncomments = 0
        self.input_settings = input_settings if input_settings is not None else OrderedDict()
        self.curr_stage_name = None

    def init_stage(self):
        """
        Use this to initiate new stage/black in LAMMPS input file.
        """

        self.nstages += 1

        stage_name = "stage %s" % self.nstages
        self._add_stage_name(stage_name)

    def add_comment(self, comment: str, is_stage_header: bool = False) -> None:
        """
         Method to add a comment or a stage header.

        Args:
            comment: Comment string to be added. The comment will be
                preceded '##' (double hash) in the generated input.
            is_stage_header: Set this to True, if the comment should be
                treated as stage/block header to be used at the beginning of
                each stage/block. Headers are preceded by '#' (single hash)
                in the generated input.
        """

        if is_stage_header:
            self._add_command("header", "# " + comment)
        else:
            self.ncomments += 1
            self._add_command(f"comment_{self.ncomments}", "## " + comment)

    def add_commands(self, commands: Union[str, list, dict], stage_name: str = None):
        """
        Adds LAMMPS command/s and its arguments to LAMMPS input file.

        Args:
            commands: LAMMPS command/s. Can pass single string or list of LAMMPS command
            with its arguments. Also accepts dictionary of LAMMPS commands are
            corresponding arguments as key value pairs.
            stage_name: If a stage name is mentioned, the command is added
                under that stage block else latest stage is assumed.
        """

        if self.curr_stage_name is None:
            self.init_stage()
        self._add_stage_name(stage_name)

        if isinstance(commands, str):
            self._add_command(commands)
        elif isinstance(commands, list):
            for command in commands:
                self._add_command(command=command)
        elif isinstance(commands, dict):
            for command, args in commands.items():
                self._add_command(command=command, args=args)

    def add_stage(self, stage_commands: Union[list, dict], header: str = None, description: str = None):
        """
        Adds an entire stage or block of LAMMPS input commands and argument.

        Args:
            stage_commands:A dictionary containing LAMMPS commands and arguments
                as key value pairs. Example: {'units': 'metal', 'atom_style': 'charge'}.
                Alternatively, a list of strings containing both LAMMPS commands
                and corresponding arguments can be provided. Example: ['units  metal',
                'atom_style charge']
            header: The header for the stage. The stage number will be used
                as the header if not provided.
            description: Add short description for this stage. This will
                add inline with the stage header in the input
        """

        self.init_stage()

        if header is None:
            header = self.curr_stage_name

        if description is not None:
            self.add_comment(f"{header} : {description}", is_stage_header=True)
        else:
            self.add_comment(header, is_stage_header=True)  # type: ignore

        self.add_commands(commands=stage_commands)

    def generate_lammps_input(self):
        """
        Returns: LammpsInputFile object
        """
        return LammpsInputFile(self.input_settings)

    def get_string(self):
        """
        Generates and returns the string representation of LAMMPS input file.

        Returns: String representation of LAMMPS input file.

        """

        lammps_input = "## LAMMPS input generated from LammpsInputConstructor\n\n"
        for stage, command_dict in self.input_settings.items():
            for command, args in command_dict.items():
                if command.startswith("comment") or command.startswith("header"):
                    if command.startswith("header"):
                        lammps_input += "\n"
                    lammps_input += args + "\n"
                else:
                    lammps_input += f"{command:20} {args}\n"
        return lammps_input

    def from_string(self, s: str):  # type: ignore
        """
        Helper method to parse string representation of LammpsInputFile

        Args:
            s: String representation of LammpsInputFile.

        Returns: LammpsInputFile
        """

        for line in self._clean_lines(s.splitlines()[1:]):
            if line[:2] == "##":
                self.add_comment(comment=line[2:], is_stage_header=False)
            if line[0] == "#":
                self.init_stage()
                self.add_comment(comment=line[2:], is_stage_header=True)
            else:
                self.add_commands(line)

        return LammpsInputFile(self.input_settings)

    def from_file(self, path: Union[str, Path]):  # type: ignore
        """
        Creates an InputFile object from a file.

        Args:
            path: Filename to read, including path.

        Returns:
            InputFile
        """
        filename = path if isinstance(path, Path) else Path(path)
        with zopen(filename, "rt") as f:
            return self.from_string(f.read())

    @classmethod
    def from_dict(cls, d: dict):
        """
        Reads LammpsInputFile from Ordered dictionary.
        Args:
            d: LammpsIncarFile ordered dict.

        Returns:
            LammpsInputFile object
        """
        return LammpsInputFile(OrderedDict({k: v for k, v in d.items() if not k.startswith("@")}))

    def as_dict(self):
        """
        Returns ordered dictionary of LammpsInputFile parameters.

        Returns:
            LammpsInputFile ordered dictionary.
        """
        d = self.input_settings.copy()
        d["@module"] = type(self).__module__
        d["@class"] = type(self).__name__
        return d

    def __str__(self):
        return self.get_string()

    @staticmethod
    def _clean_lines(string_list: list):
        """
        Helper method to strips whitespaces, carriage returns and empty
        lines from a list of strings.

        Args:
            string_list (list): List of strings.
        """
        for s in string_list:
            clean_s = s
            clean_s = clean_s.strip()
            if not (clean_s == "" or clean_s == "\n"):
                yield clean_s

    def _add_stage_name(self, stage_name: str = None) -> None:
        """
        Helper method to generate and add stage name internally.

        Args:
            stage_name: Stage name to be added.
        """

        if stage_name is None:
            stage_name = self.curr_stage_name if self.curr_stage_name is not None else "stage %s" % self.nstages

        if not (stage_name in self.input_settings.keys()):
            self.input_settings[stage_name] = OrderedDict()
            self.curr_stage_name = stage_name  # type: ignore

    def _add_command(self, command: str, args: str = None, stage_name: str = None):
        """
        Helper method to add a single LAMMPS command and its arguments to
        LAMMPS input file.

        Args:
            command: LAMMPS command
            args: Arguments for the LAMMPS command
            stage_name: If a stage name is mentioned, the command is added
                under that stage block else latest stage is assumed.
        """
        if self.curr_stage_name is None:
            self.init_stage()
        self._add_stage_name(stage_name)

        if args is None:
            string_split = command.split()
            command = string_split[0]
            args = " ".join(string_split[1:])
        self.input_settings[self.curr_stage_name][command] = args


class CombinedData(InputFile):
    """
    Class representing a LAMMPS data file, e.g. system.data
    """


class LammpsRun(MSONable):
    """
    Examples for various simple LAMMPS runs with given simulation box,
    force field and a few more settings. Experienced LAMMPS users should
    consider using write_lammps_inputs method with more sophisticated
    templates.

    """

    template_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "templates")

    def __init__(self, script_template, settings, data, script_filename):
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
        write_lammps_inputs(
            output_dir=output_dir,
            script_template=self.script_template,
            settings=self.settings,
            data=self.data,
            script_filename=self.script_filename,
            **kwargs,
        )

    @classmethod
    def md(cls, data, force_field, temperature, nsteps, other_settings=None):
        r"""
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
        settings.update({"force_field": force_field, "temperature": temperature, "nsteps": nsteps})
        script_filename = "in.md"
        return cls(
            script_template=script_template,
            settings=settings,
            data=data,
            script_filename=script_filename,
        )


class LammpsTemplateGen(TemplateInputGen):
    """
    Creates an InputSet object for a LAMMPS run based on a template file.
    The input script is constructed by substituting variables into placeholders
    in the template file using python's Template.safe_substitute() function.
    The data file containing coordinates and topology information can be provided
    as a LammpsData instance. Alternatively, you can include a read_data command
    in the template file that points to an existing data file.
    Other supporting files are not handled at the moment.

    To write the input files to a directory, call LammpsTemplateSet.write_input()
    See pymatgen.io.template.py for additional documentation of this method.
    """

    def get_input_set(  # type: ignore
        self,
        script_template: Union[str, Path],
        settings: Optional[Dict] = None,
        script_filename: str = "in.lammps",
        data: Union[LammpsData, CombinedData] = None,
        data_filename: str = "system.data",
    ):
        """
        Args:
            script_template: String template for input script with
                placeholders. The format for placeholders has to be
                '$variable_name', e.g., '$temperature'
            settings: Contains values to be written to the
                placeholders, e.g., {'temperature': 1}. Default to None.
            data: Data file as a LammpsData instance. Default to None, i.e., no
                data file supplied. Note that a matching 'read_data' command
                must be provided in the script template in order for the data
                file to actually be read.
            script_filename: Filename for the input file.
            data_filename: Filename for the data file, if provided.
        """
        input_set = super().get_input_set(template=script_template, variables=settings, filename=script_filename)

        if data:
            input_set.update({data_filename: data})
        return input_set


@deprecated(LammpsTemplateGen, "This method will be retired in the future. Consider using LammpsTemplateSet instead.")
def write_lammps_inputs(
    output_dir,
    script_template,
    settings=None,
    data=None,
    script_filename="in.lammps",
    make_dir_if_not_present=True,
    **kwargs,
):
    """
    Writes input files for a LAMMPS run. Input script is constructed
    from a str template with placeholders to be filled by custom
    settings. Data file is either written from a LammpsData
    instance or copied from an existing file if read_data cmd is
    inspected in the input script. Other supporting files are not
    handled at the moment.

    Args:
        output_dir (str): Directory to output the input files.
        script_template (str): String template for input script with
            placeholders. The format for placeholders has to be
            '$variable_name', e.g., '$temperature'
        settings (dict): Contains values to be written to the
            placeholders, e.g., {'temperature': 1}. Default to None.
        data (LammpsData or str): Data file as a LammpsData instance or
            path to an existing data file. Default to None, i.e., no
            data file supplied. Useful only when read_data cmd is in
            the script.
        script_filename (str): Filename for the input script.
        make_dir_if_not_present (bool): Set to True if you want the
            directory (and the whole path) to be created if it is not
            present.
        **kwargs: kwargs supported by LammpsData.write_file.

    Examples:
        >>> eam_template = '''units           metal
        ... atom_style      atomic
        ...
        ... lattice         fcc 3.615
        ... region          box block 0 20 0 20 0 20
        ... create_box      1 box
        ... create_atoms    1 box
        ...
        ... pair_style      eam
        ... pair_coeff      1 1 Cu_u3.eam
        ...
        ... velocity        all create $temperature 376847 loop geom
        ...
        ... neighbor        1.0 bin
        ... neigh_modify    delay 5 every 1
        ...
        ... fix             1 all nvt temp $temperature $temperature 0.1
        ...
        ... timestep        0.005
        ...
        ... run             $nsteps'''
        >>> write_lammps_inputs('.', eam_template, settings={'temperature': 1600.0, 'nsteps': 100})
        >>> with open('in.lammps') as f:
        ...     script = f.read()
        ...
        >>> print(script)
        units           metal
        atom_style      atomic

        lattice         fcc 3.615
        region          box block 0 20 0 20 0 20
        create_box      1 box
        create_atoms    1 box

        pair_style      eam
        pair_coeff      1 1 Cu_u3.eam

        velocity        all create 1600.0 376847 loop geom

        neighbor        1.0 bin
        neigh_modify    delay 5 every 1

        fix             1 all nvt temp 1600.0 1600.0 0.1

        timestep        0.005

        run             100


    """
    variables = {} if settings is None else settings
    template = Template(script_template)
    input_script = template.safe_substitute(**variables)
    if make_dir_if_not_present and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    with open(os.path.join(output_dir, script_filename), "w") as f:
        f.write(input_script)
    read_data = re.search(r"read_data\s+(.*)\n", input_script)
    if read_data:
        data_filename = read_data.group(1).split()[0]
        if isinstance(data, LammpsData):
            data.write_file(os.path.join(output_dir, data_filename), **kwargs)
        elif isinstance(data, str) and os.path.exists(data):
            shutil.copyfile(data, os.path.join(output_dir, data_filename))
        else:
            warnings.warn(f"No data file supplied. Skip writing {data_filename}.")
