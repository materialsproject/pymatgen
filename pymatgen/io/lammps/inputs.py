# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module implements methods for reading/manupilating/writing LAMMPS input files.
It does not implement methods for automatically creating inputs based on a structure
and computation type. For this, see the InputSet and InputGenerator in sets.py, or
https://github.com/Matgenix/atomate2-lammps
"""

from __future__ import annotations

import os
import re
import shutil
import warnings
from pathlib import Path
from string import Template

import numpy as np
from monty.dev import deprecated
from monty.io import zopen
from monty.json import MSONable

from pymatgen.io.core import InputFile
from pymatgen.io.lammps.data import CombinedData, LammpsData
from pymatgen.io.template import TemplateInputGen

__author__ = "Kiran Mathew, Brandon Wood, Zhi Deng, Manas Likhit, Guillaume Brunin (Matgenix)"
__copyright__ = "Copyright 2018, The Materials Virtual Lab"
__version__ = "2.0"
__maintainer__ = "Zhi Deng, Guillaume Brunin (Matgenix)"
__email__ = "z4deng@eng.ucsd.edu, info@matgenix.com"
__date__ = "Nov 2022"


class LammpsInputFile(InputFile):
    """
    Class representing a LAMMPS input settings file, e.g. in.lammps.
    Allows for LAMMPS input generation in line/stage wise manner. A stage
    here is defined as a block of LAMMPS input commands usually performing a
    specific task during the simulation such as energy minimization or
    NPT/NVT runs. But more broadly, a stage can also be a block of LAMMPS
    input where the simulation box is set up, a set of variables are declared or
    quantities are computed.

    The LammpsInputFile is defined by the attribute `list_of_commands`,
    i.e. a list of stages (strings) together with the corresponding LAMMPS input settings (strings).
    The structure is the following:
    ```
    list_of_commands = [
        ["Stage 1", [[cmd1, args1], [cmd2, args2]]],
        ["Stage 2", [[cmd3, args3]]]
    ]
    ```
    where cmd's are the LAMMPS command names (e.g., "units", or "pair_coeff")
    and the args are the corresponding arguments.
    "Stage 1" and "Stage 2" are examples of stage names.

    """

    def __init__(self, list_of_commands: list | None = None):
        """
        Args:
            list_of_commands: list of LAMMPS input settings.
        """
        self.list_of_commands = list_of_commands if list_of_commands else []
        self.nstages = self.get_nstages()
        self.ncomments = self.get_ncomments()

    @property
    def stages_names(self) -> list:
        """List of names for all the stages present in list_of_commands."""
        return [stage_name for stage_name, _ in self.list_of_commands] if self.list_of_commands else []

    def get_nstages(self) -> int:
        """Returns the number of stages in the current LammpsInputFile."""
        return len(self.stages_names)

    def get_ncomments(self) -> int:
        """
        Returns the number of comments in the current LammpsInputFile. Includes the blocks of comments as well
        as inline comments (comment lines within blocks of LAMMPS commands).
        """
        ncomments = 0
        for _, stage in self.list_of_commands:
            # Block of comment = 1 comment
            if all(cmd.startswith("#") for cmd, args in stage):
                ncomments += 1
            else:
                # Else, inline comment each count as one
                ncomments += sum(1 for cmd, args in stage if cmd.startswith("#"))

        return ncomments

    def get_args(self, command: str, stage_name: str | None = None) -> list | str:
        """
        Given a command, returns the corresponding arguments (or list of arguments) in the LammpsInputFile.
        A stage name can be given; in this case the search will happen only for this stage.
        If the command is not found, None is returned.

        Args:
            command (str): String with the command to find in the input file (e.g., "units").
            stage_name (str): String giving the stage name where the change should take place.

        Returns:
            Value of the argument corresponding to the command.
            List if the same command is used multiple times.
        """
        args = []
        stages_to_look = [stage_name] if stage_name else self.stages_names

        for stage, cmd_list in self.list_of_commands:
            if stage in stages_to_look:
                for cmd, arg in cmd_list:
                    if command == cmd:
                        args.append(arg)

        if len(args) == 1:
            return args[0]
        else:
            return args

    def contains_command(self, command: str, stage_name: str | None = None) -> bool:
        """
        Returns whether a given command is present in the LammpsInputFile.
        A stage name can be given; in this case the search will happen only for this stage.

        Args:
            command (str): String with the command to find in the input file (e.g., "units").
            stage_name (str): String giving the stage name where the change should take place.

        Returns:
            True if the command is present, False is not.
        """
        return bool(self.get_args(command, stage_name))

    def set_args(self, command: str, argument: str, stage_name: str | None = None, how: str | int | list[int] = "all"):
        """
        Sets the arguments for the given command to the given string.
        If a stage name is specified, it will be replaced or set only for this stage.
        If the command is set multiple times in the file/stage, it will be replaced based on "how":
        either the first occurrence, all of them, or the index of the occurrence.

        Args:
            command (str): String representing the command to change, e.g., "units".
            argument (str): String with the new value for the command, e.g., "atomic".
            stage_name (str): String giving the stage name where the change should take place.
            how (str or int or list): "all" for changing all occurrences of the command within the stage_name
                or the whole input file, "first" for the first occurrence, int i for the i-th time the command
                is present in the stage_name or the whole input file, starting at 0. Can be a list of indexes as well.
        """
        # Get the stages to look in
        stages_to_look = [stage_name] if stage_name else self.stages_names

        # Translates how into range of indices
        if how == "first":
            how = [0]
        elif how == "all":
            getargs = self.get_args(command, stage_name)
            N = 1 if isinstance(getargs, str) else len(getargs)
            how = [i for i in range(N)]
        elif isinstance(how, int):
            how = [how]

        # Look for occurrences in the relevant stages
        i = 0
        for i_stage, (stage, cmd_list) in enumerate(self.list_of_commands):
            if stage in stages_to_look:
                for i_cmd, (cmd, _) in enumerate(cmd_list):
                    if command == cmd:
                        if i in how:
                            self.list_of_commands[i_stage][1][i_cmd][1] = argument
                        i += 1

    def add_stage(
        self,
        command: str | list[str] | dict[str, str | float],
        stage_name: str | None = None,
        after_stage: str | None = None,
    ):
        r"""
        Adds LAMMPS command(s) and its arguments to LAMMPS input file.

        Examples:
            1) In order to add a stage defining the potential to be used, you can use:
            ```
            your_input_file.add_stage(
                command=["pair_coeff 1 1 morse 0.0580 3.987 3.404", "pair_coeff 1 4 morse 0.0408 1.399 3.204"],
                stage_name="Definition of the potential"
            )
            ```

            2) Another stage could consist in an energy minimization. In that case, the commands could look like
            ```
            command = [
                "thermo 1",
                "thermo_style custom step lx ly lz press pxx pyy pzz pe",
                "dump dmp all atom 5 run.dump",
                "min_style cg",
                "fix 1 all box/relax iso 0.0 vmax 0.001",
                "minimize 1.0e-16 1.0e-16 5000 10000",
                "write_data run.data"
            ]
            ```
            or a similar version with a long string containing all (separated by \n), or a dictionary such as
            `{"thermo": 1, ...}`.

        Args:
            command (str or list or dict): LAMMPS command(s) for this stage of the run.
                Can pass a single string or a list of LAMMPS commands
                with their arguments. Also accepts a dictionary of LAMMPS commands and
                corresponding arguments as key, value pairs.
            stage_name (str): If a stage name is mentioned, the command is added
                under that stage block, else a new stage is created and named from numbering.
            after_stage (str): Name of the stage after which this stage should be added.
                If None, the stage is added at the end of the LammpsInputFile.
        """
        # Name the stage if not given
        # + 1 so that if a new (list of) command is added without stage name,
        # it goes in a new stage.
        if stage_name is None:
            stage_name = f"Stage {self.nstages + 1}"

        # Initialize the stage if not already present
        if stage_name not in self.stages_names:
            if not self.list_of_commands:
                self.list_of_commands = [[stage_name, []]]
            else:
                if after_stage is None:
                    index_insert = -1
                elif after_stage in self.stages_names:
                    index_insert = self.stages_names.index(after_stage) + 1
                    if index_insert == len(self.stages_names):
                        index_insert = -1
                else:
                    raise ValueError("The stage after which this one should be added does not exist.")
                if index_insert == -1:
                    self.list_of_commands.append([stage_name, []])
                else:
                    self.list_of_commands.insert(index_insert, [stage_name, []])
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

    def remove_stage(self, stage_name: str):
        """
        Removes a whole stage from the LammpsInputFile.

        Args:
            stage_name (str): name of the stage to remove.
        """
        if stage_name in self.stages_names:
            self.list_of_commands.pop(self.stages_names.index(stage_name))
            self.nstages = self.get_nstages()
            self.ncomments = self.get_ncomments()
        else:
            raise LookupError("The given stage name is not present in this LammpsInputFile.")

    def merge_stages(self, stage_names: list[str]):
        """
        Merges multiple stages of a LammpsInputFile together.
        The merged stage will be at the same index as the first of the stages to be merged.
        The others will appear in the same order as provided in the list. Other non-merged stages will follow.

        Args:
             stage_names (list): list of strings giving the names of the stages to be merged.
        """
        if not all([stage in self.stages_names for stage in stage_names]):
            raise ValueError("At least one of the stages to be merged is not in the LammpsInputFile.")

        indices_stages_to_merge = [self.stages_names.index(stage) for stage in stage_names]
        if not np.all([np.array(indices_stages_to_merge[1:]) >= np.array(indices_stages_to_merge[:-1])]):
            raise ValueError(
                """The provided stages are not in the right order. You should merge stages
            order of appearance in your LammpsInputFile. If you want to reorder stages, use add_stage and
            remove_stage."""
            )

        list_of_commands = self.list_of_commands[: indices_stages_to_merge[0]]
        stages_names = self.stages_names[: indices_stages_to_merge[0]]

        merge_name = "Merge of: " + ", ".join([self.stages_names[i] for i in indices_stages_to_merge])
        merged_commands = []
        for i in indices_stages_to_merge:
            for j in range(len(self.list_of_commands[i][1])):
                merged_commands.append(self.list_of_commands[i][1][j])

        merged_stages = [merge_name, merged_commands]
        stages_names.append(merge_name)
        list_of_commands.append(merged_stages)

        for i_stage, stage in enumerate(self.list_of_commands):
            if i_stage > indices_stages_to_merge[0] and i_stage not in indices_stages_to_merge:
                list_of_commands.append(stage)

        self.list_of_commands = list_of_commands
        self.ncomments = self.get_ncomments()
        self.nstages = self.get_nstages()

    def add_command(self, stage_name: str, command: str, args: str | float | None = None):
        """
        Helper method to add a single LAMMPS command and its arguments to
        a LAMMPS input file. The stage name should be provided: a default behavior
        is avoided here to avoid mistakes. To add a command at the end of the input file,
        use add_stage directly.

        Example:
            In order to add the command ``pair_coeff 1 1 morse 0.0580 3.987 3.404``
            to the stage "Definition of the potential", simply use
            ```
            your_input_file.add_command(
                stage_name="Definition of the potential",
                command="pair_coeff 1 1 morse 0.0580 3.987 3.404"
            )
            ```
            or
            ```
            your_input_file.add_command(
                stage_name="Definition of the potential",
                command="pair_coeff",
                args="1 1 morse 0.0580 3.987 3.404"
            )
            ```

        Args:
            stage_name (str): name of the stage to which the command should be added.
            command (str): LAMMPS command, with or without the arguments.
            args (str): Arguments for the LAMMPS command.
        """
        if args is None:
            string_split = command.split()
            command = string_split[0]
            args = " ".join(string_split[1:])

        # Find where the command should be added (stage index instead of stage name)
        idx = self.stages_names.index(stage_name)
        # Add the command
        if not self.list_of_commands[idx][1]:
            self.list_of_commands[idx][1] = [[command, args]]
        else:
            self.list_of_commands[idx][1].append([command, args])

    def remove_command(self, command: str, stage_name: str | list[str] | None = None, remove_empty_stages: bool = True):
        """
        Removes a given command from a given stage. If no stage is given, removes all occurrences of the command.
        In case removing a command completely empties a stage, the choice whether to keep this stage in the
        LammpsInputFile is given by remove_empty_stages.

        Args:
            command (str): command to be removed.
            stage_name (str or list): names of the stages where the command should be removed.
            remove_empty_stages (bool): whether to remove the stages emptied by removing the command or not.
        """
        if stage_name is None:
            stage_name = self.stages_names
        elif isinstance(stage_name, str):
            stage_name = [stage_name]
        elif not isinstance(stage_name, list):
            raise ValueError("If given, stage_name should be a string or a list of strings.")

        n_removed = 0
        indices_to_remove = []
        new_list_of_stages = []
        for i_stage, (stage, cmd_list) in enumerate(self.list_of_commands):
            if stage in stage_name:
                new_list_of_commands = []
                for i_cmd, (cmd, arg) in enumerate(cmd_list):
                    if cmd == command:
                        n_removed += 1
                    else:
                        new_list_of_commands.append([cmd, arg])
                        indices_to_remove.append([i_stage, i_cmd])
                if new_list_of_commands or not remove_empty_stages:
                    new_list_of_stages.append([stage, new_list_of_commands])
            else:
                new_list_of_stages.append([stage, cmd_list])

        self.list_of_commands = new_list_of_stages
        self.ncomments = self.get_ncomments()
        self.nstages = self.get_nstages()

        if n_removed == 0:
            warnings.warn(f"{command} not found in the LammpsInputFile.")

    def add_comment(
        self, comment: str, inline: bool = False, stage_name: str | None = None, index_comment: bool = False
    ):
        """
         Method to add a comment inside a stage (between actual commands)
         or as a whole stage (which will do nothing when LAMMPS runs).

        Args:
            comment (str): Comment string to be added. The comment will be
                preceded by '#' in the generated input.
            inline (bool): True if the comment should be inline within a given block of commands.
            stage_name (str): set the stage_name to which the comment will be written. Required if inline is True.
            index_comment (bool): True if the comment should start with "Comment x" with x its number in the ordering.
                Used only for inline comments.
        """
        self.ncomments += 1
        # "Stage" of comment only
        if not inline:
            if stage_name is None:
                stage_name = f"Comment {self.ncomments}"
            self.add_stage(command="# " + comment, stage_name=stage_name)

        # Inline comment
        elif inline and stage_name:
            command = "#"
            if index_comment:
                if "Comment" in comment and comment.strip()[9] == ":":
                    args = ":".join(comment.split(":")[1:])
                else:
                    args = comment
            else:
                args = comment
            self.add_command(command=command, args=args, stage_name=stage_name)

        else:
            raise NotImplementedError("If you want to add an inline comment, please specify the stage name.")

    def append(self, lmp_input_file: LammpsInputFile, renumbering: bool = True):
        """
        Appends a LammpsInputFile to another. The list_of_commands are merged,
        and the numbering of stages/comments is either kept the same or updated.

        Args:
            lmp_input_file (LammpsInputFile): LammpsInputFile to append.
            renumbering (bool): whether to renumber stages and comments for the added LammpsInputFile.
        """
        # Renumbering comments and stages of the lmp_input_file.list_of_commands
        new_list_to_add = lmp_input_file.list_of_commands
        if renumbering:
            for i_stage, stage in enumerate(lmp_input_file.list_of_commands):
                if "Comment" == stage[0].split()[0] and stage[0].split()[1].isdigit():
                    i_comment = int(stage[0].split()[1]) + self.ncomments
                    new_list_to_add[i_stage][0] = f"Comment {i_comment}"
                if "Stage" == stage[0].split()[0] and stage[0].split()[1].isdigit():
                    this_stage = int(stage[0].split()[1]) + self.nstages
                    new_list_to_add[i_stage][0] = f"Stage {this_stage}"
        self.list_of_commands += new_list_to_add
        self.ncomments = self.get_ncomments()
        self.nstages = self.get_nstages()

    def get_string(self, ignore_comments: bool = False, keep_stages: bool = True) -> str:
        """
        Generates and returns the string representation of LAMMPS input file.
        Stages are separated by empty lines.
        The headers of the stages will be put in comments preceding each stage.
        Other comments will be put inline within stages, where they have been added.

        Args:
            ignore_comments (bool): True if only the commands should be kept from the input file.
            keep_stages (bool): True if the block structure from the input file should be kept.
                                If False, a single block is assumed.

        Returns: String representation of LAMMPS input file.
        """
        lammps_input = "# LAMMPS input generated from LammpsInputFile\n"
        if not keep_stages:
            lammps_input += "\n"

        for stage_name, cmd_list in self.list_of_commands:
            if keep_stages:
                # Print first the name of the stage in a comment.
                # We print this even if ignore_comments is True.
                if "Comment" not in stage_name and len(self.list_of_commands) > 1:
                    lammps_input += "\n# " + stage_name + "\n"

                # In case of a block of comment, the header is not printed (useless)
                else:
                    lammps_input += "\n"

            # Then print the commands and arguments (including inline comments)
            for command, args in cmd_list:
                if not (ignore_comments and "#" in command):
                    lammps_input += f"{command} {args.strip()}\n"

        return lammps_input

    def write_file(self, filename: str | Path, ignore_comments: bool = False, keep_stages: bool = True) -> None:
        """
        Writes the input file.

        Args:
            filename (str or path): The filename to output to, including path.
            ignore_comments (bool): True if only the commands should be kept from the input file.
            keep_stages (bool): True if the block structure from the input file should be kept.
                                If False, a single block is assumed.
        """
        filename = filename if isinstance(filename, Path) else Path(filename)
        with zopen(filename, "wt") as f:
            f.write(self.get_string(ignore_comments=ignore_comments, keep_stages=keep_stages))

    @classmethod
    def from_string(cls, contents: str, ignore_comments: bool = False, keep_stages: bool = False) -> LammpsInputFile:
        """
        Helper method to parse string representation of LammpsInputFile.
        If you created the input file by hand, there is no guarantee that the representation
        will be perfect as it is difficult to account for all the cosmetic changes you
        could have done on your input script. Always check that you have what you want !
        By default, a single stage containing all the input settings is created.
        If the block structure of your input file should be kept and stored as
        different stages, set keep_stages to True.

        Args:
            contents (str): String representation of LammpsInputFile.
            ignore_comments (bool): True if only the commands should be kept from the input file.
            keep_stages (bool): True if the block structure from the input file should be kept.
                                If False, a single block is assumed.

        Returns:
            LammpsInputFile
        """
        LIF = cls()

        # Strip string from starting and/or ending white spaces
        s = contents.strip()

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
        lines = cls._clean_lines(s.splitlines(), ignore_comments=ignore_comments)
        # Split the string into blocks based on the empty lines of the input file
        blocks = cls._get_blocks(lines, keep_stages=keep_stages)

        for block in blocks:
            keep_block = True
            # Block of comment(s)
            if all([line[0] == "#" for line in block]):
                if ignore_comments:
                    keep_block = False
                else:
                    LIF.add_comment(comment=block[0][1:].strip(), inline=False)
                    stage_name = f"Comment {LIF.ncomments}"
                    if len(block) > 1:
                        for line in block[1:]:
                            LIF.add_comment(comment=line[1:].strip(), inline=True, stage_name=stage_name)

            # Header of a stage
            elif block[0][0] == "#" and keep_block:
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

                header = header.strip()
                commands = block[icomm_max:]
                LIF.add_stage(command=commands, stage_name=None if (ignore_comments or not keep_stages) else header)

            # Stage with no header
            else:
                LIF.add_stage(command=block)

        return LIF

    @classmethod
    def from_file(cls, path: str | Path, ignore_comments: bool = False, keep_stages: bool = False) -> LammpsInputFile:
        """
        Creates an InputFile object from a file.

        Args:
            path (str or path): Filename to read, including path.
            ignore_comments (bool): True if only the commands should be kept from the input file.
            keep_stages (bool): True if the block structure from the input file should be kept.
                                If False, a single block is assumed.

        Returns:
            LammpsInputFile
        """
        filename = path if isinstance(path, Path) else Path(path)
        with zopen(filename, "rt") as f:
            return cls.from_string(f.read(), ignore_comments=ignore_comments, keep_stages=keep_stages)

    def __repr__(self):
        return self.get_string()

    @staticmethod
    def _clean_lines(string_list: list, ignore_comments: bool = False) -> list:
        r"""
        Helper method to strips whitespaces, carriage returns and redundant empty
        lines from a list of strings.
        Transforms "& \n" and "&\n" into "" as the & symbol means the line continues.
        Also removes lines with "# LAMMPS input generated from LammpsInputFile"
        to avoid possible repetitions.

        Args:
            string_list (list): List of strings.
            ignore_comments (bool): True if the strings starting with # should be ignored.

        Returns:
            List of strings
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

        # Get rid of empty comments that are there just for cosmetic reasons
        new_list = []
        for s in string_list:
            if len(s) > 1:
                new_list.append(s)
            elif not (len(s.strip()) == 1 and s[0] == "#"):
                new_list.append(s)
        string_list = new_list

        # Keep only a single empty lines when there are multiple ones
        lines = [string_list[0]]

        for i, s in enumerate(string_list[1:-1]):
            if s != "" and not (s[0] == "#" and ignore_comments):
                lines.append(s)
            elif s == "" and string_list[i + 2] != "":
                lines.append(s)

        lines.append(string_list[-1])

        return lines

    @staticmethod
    def _get_blocks(string_list: list[str], keep_stages: bool = False) -> list[list[str]]:
        """
        Helper method to return a list of blocks of LAMMPS commands,
        separated from "" in a list of all commands.

        Args:
            string_list (list): List of strings.
            keep_stages (bool): True if the block structure from the input file should be kept.
                         If False, a single block is assumed.

        Returns:
            List of list of strings containing the different blocks
        """
        blocks: list[list[str]] = [[]]
        i_block = 0

        for s in string_list:
            if s != "":
                blocks[i_block].append(s)
            if s == "" and keep_stages:
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
