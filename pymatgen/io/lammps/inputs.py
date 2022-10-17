# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module implements methods for writing LAMMPS input files.
"""

from __future__ import annotations

import os
import re
import shutil
import warnings
from pathlib import Path
from string import Template

from monty.dev import deprecated
from monty.json import MSONable

from pymatgen.io.core import InputFile
from pymatgen.io.lammps.data import CombinedData, LammpsData
from pymatgen.io.template import TemplateInputGen

__author__ = "Kiran Mathew, Brandon Wood, Zhi Deng, Manas Likhit, Guillaume Brunin"
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
    quantities are computed.
    """

    def __init__(self, input_settings: dict = None):
        """
        Args:
            input_settings: dictionary of LAMMPS input settings.
        """
        self.nstages = 0
        self.ncomments = 0
        self.input_settings = input_settings if input_settings else {}

    def add_stage(self, command: str | list | dict, stage_name: str = None):
        """
        Adds LAMMPS command(s) and its arguments to LAMMPS input file.

        Args:
            command: LAMMPS command(s) for this stage of the run.
                Can pass single string or list of LAMMPS commands
                with their arguments. Also accepts dictionary of LAMMPS commands and
                corresponding arguments as key, value pairs.
            stage_name: If a stage name is mentioned, the command is added
                under that stage block else a new stage is created and named from numerotation.
        """
        # Name the stage if not given
        # + 1 so that if a new (list of) command is added without stage name,
        # it goes in a new stage.
        if stage_name is None:
            stage_name = f"Stage {self.nstages + 1}"

        # Initialize the stage if not already present
        if stage_name not in self.input_settings.keys():
            self.input_settings[stage_name] = {}
            self.nstages += 1

        # Handle the different input formats to add commands to the stage
        if isinstance(command, str):
            self.add_command(command=command, stage_name=stage_name)

        elif isinstance(command, list):
            for comm in command:
                if comm[0] == "#":
                    self.add_comment(comment=comm[1:], inline=True, stage_name=stage_name, index_comment=True)
                else:
                    self.add_command(command=comm, stage_name=stage_name)

        elif isinstance(command, dict):
            for comm, args in command.items():
                if comm[0] == "#":
                    self.add_comment(comment=comm[1:], inline=True, stage_name=stage_name, index_comment=True)
                else:
                    self.add_command(command=comm, args=args, stage_name=stage_name)

        else:
            raise TypeError("The command should be a string, list of strings or dictionary.")

    def add_command(self, stage_name: str, command: str, args: str = None):
        """
        Helper method to add a single LAMMPS command and its arguments to
        a LAMMPS input file. The stage name should be provided: a default behavior
        is avoided here to avoid mistakes.

        Args:
            stage_name (str): name of the stage to which the command should be added.
            command (str): LAMMPS command, with or without the arguments.
            args (str): Arguments for the LAMMPS command.
        """
        if args is None:
            string_split = command.split()
            command = string_split[0]
            args = " ".join(string_split[1:])

        if command not in self.input_settings[stage_name]:
            self.input_settings[stage_name][command] = args
        else:
            # In LAMMPS, the same command can be used multiple times
            self.input_settings[stage_name][command] += "\n" + command + " " + args

    def add_comment(self, comment: str, inline: bool = False, stage_name: str = None, index_comment: bool = False):
        """
         Method to add a comment in a stage or as a whole stage (which will do nothing).

        Args:
            comment: Comment string to be added. The comment will be
                preceded by '#' in the generated input.
            inline: True if the comment should be inline within a given block of commands.
            stage_name: set the stage_name to which the comment will be written. Required if inline is True.
            index_comment: True if the comment should start with "Comment x" with x its number in the ordering.
                Used only for inline comments.
        """
        self.ncomments += 1
        # "Stage" of comment only
        if not inline:
            if stage_name is None:
                stage_name = f"Comment {self.ncomments}"
            self.add_stage(command="# " + comment, stage_name=stage_name)
            self.nstages += -1

        # Inline comment
        elif inline and stage_name:
            if index_comment:
                command = f"# Comment {self.ncomments}:"
            else:
                command = "#"
            self.add_command(command=command, args=comment, stage_name=stage_name)

        else:
            raise NotImplementedError("If you want to add an inline comment, please specify the stage name.")

    def get_string(self):
        """
        Generates and returns the string representation of LAMMPS input file.
        Stages are separated by empty lines.
        The headers of the stages will be put in comments preceding each stage.
        Other comments will be put inline within stages, where they have been added.

        Returns: String representation of LAMMPS input file.
        """
        lammps_input = "# LAMMPS input generated from LammpsInputFile\n"

        for stage_name, stage_dict in self.input_settings.items():
            # Print first the name of the stage in a comment
            if "Comment" not in stage_name:
                lammps_input += "\n# " + stage_name + "\n"

            # In case of a block of comment, the header is not printed (useless)
            else:
                lammps_input += "\n"

            # Then print the commands and arguments (including inline comments)
            for command, args in stage_dict.items():
                lammps_input += f"{command} {args}\n"

        return lammps_input

    @classmethod
    def from_string(cls, s: str):
        """
        Helper method to parse string representation of LammpsInputFile.
        If you created the input file by hand, there is no guarantee that the representation
        will be perfect as it is difficult to account for all the cosmetic changes you
        could have done on your input script. Always check that you have what you want !

        Args:
            s: String representation of LammpsInputFile.

        Returns: LammpsInputFile
        """
        LIF = cls()

        # Strip string from starting and/or ending white spaces
        s = s.strip()

        # Remove "&" symbols at the end of lines
        while "&" in s:
            sequence = "&"
            index = s.index("&")
            next_symbol = ""
            i = 0
            while next_symbol != "\n":
                sequence += next_symbol
                i += 1
                next_symbol = s[index + i]
            s = s.replace(sequence + "\n", "")

        # Remove unwanted lines from the string
        lines = cls._clean_lines(s.splitlines())
        # Split the string into blocks based on the empty lines of the input file
        blocks = cls._get_blocks(lines)

        for block in blocks:
            # Block of comment(s)
            if all([line[0] == "#" for line in block]):
                LIF.add_comment(comment=block[0][1:].strip(), inline=False)
                stage_name = f"Comment {LIF.ncomments}"
                if len(block) > 1:
                    for line in block[1:]:
                        LIF.add_comment(comment=line[1:].strip(), inline=True, stage_name=stage_name)

            # Header of a stage
            elif block[0][0] == "#":
                # Find the name of the header.
                # If the comment is on multiple lines, the header will be the whole text
                icomm_max = len(block)
                for i, line in enumerate(block):
                    if line[0] != "#" and i <= icomm_max:
                        icomm_max = i

                comments = block[:icomm_max]
                header = ""
                for line in comments:
                    header += line[1:].strip() + " "

                commands = block[icomm_max:]
                LIF.add_stage(command=commands, stage_name=header)

            # Stage with no header
            else:
                LIF.add_stage(command=block)

        return LIF

    @staticmethod
    def _clean_lines(string_list: list):
        """
        Helper method to strips whitespaces, carriage returns and redundant empty
        lines from a list of strings.
        Transforms "& \n" and "&\n" into "" as the & symbol means the line continues.
        Also removes lines with "# LAMMPS input generated from LammpsInputFile"
        to avoid possible repetitions.

        Args:
            string_list (list): List of strings.
        """
        if len(string_list) == 0 or all([s == "" for s in string_list]):
            raise ValueError("The list of strings should contain some non-empty strings.")

        # Remove the first comment line possibly added by previous input creations
        while "# LAMMPS input generated from LammpsInputFile" in string_list:
            string_list.remove("# LAMMPS input generated from LammpsInputFile")

        # Get rid of possible initial empty lines
        # and final empty lines
        imin = len(string_list)
        imax = 0
        for i, s in enumerate(string_list):
            if s != "" and i <= imin:
                imin = i
            if s != "" and i >= imax:
                imax = i
        string_list = string_list[imin : imax + 1]

        # Keep only a single empty lines when there are multiple ones
        lines = [string_list[0]]

        for i, s in enumerate(string_list[1:-1]):
            if s != "":
                lines.append(s)
            if s == "" and string_list[i + 2] != "":
                lines.append(s)

        lines.append(string_list[-1])

        return lines

    @staticmethod
    def _get_blocks(string_list: list):
        """
        Helper method to return a list of blocks of lammps commands,
        separated from "" in a list of all commands.

        Args:
            string_list (list): List of strings.
        """
        blocks: list[list[str]] = [[]]
        i_block = 0

        for s in string_list:
            if s != "":
                blocks[i_block].append(s)
            if s == "":
                blocks.append([])
                i_block += 1

        return blocks


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
        Example for a simple MD run based on template md.template.

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
        template_path = os.path.join(cls.template_dir, "md.template")
        with open(template_path) as f:
            script_template = f.read()
        settings = other_settings.copy() if other_settings else {}
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
        script_template: str | Path,
        settings: dict | None = None,
        script_filename: str = "in.lammps",
        data: LammpsData | CombinedData | None = None,
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
