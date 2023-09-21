"""
This module implements methods for reading/manupilating/writing LAMMPS input files.
It does not implement methods for automatically creating inputs based on a structure
and computation type. For this, see the InputSet and InputGenerator in sets.py, or
https://github.com/Matgenix/atomate2-lammps.
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

from pymatgen.core import __version__ as CURRENT_VER
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

    The LammpsInputFile is defined by the attribute `stages`,
    i.e. a list of dicts each with keys `stage_name` and `commands`,
    defining the stage names and the corresponding LAMMPS input settings (list of tuples of two strings each).
    The structure is the following:
    ```
    stages = [
        {"stage_name": "Stage 1", "commands": [(cmd1, args1), (cmd2, args2)]},
        {"stage_name": "Stage 2", "commands": [(cmd3, args3)]}
    ]
    ```
    where cmd's are the LAMMPS command names (e.g., "units", or "pair_coeff")
    and the args are the corresponding arguments.
    "Stage 1" and "Stage 2" are examples of stage names.
    """

    def __init__(self, stages: list | None = None):
        """
        Args:
            stages: list of LAMMPS input settings.
        """
        self.stages = stages or []

        # Enforce that all stage names are unique
        if len(self.stages_names) != len(set(self.stages_names)):
            raise ValueError("Stage names should be unique.")

        # Check the format of each stage:
        for stage in self.stages:
            self._check_stage_format(stage)

    @property
    def stages_names(self) -> list:
        """List of names for all the stages present in `stages`."""
        return [stage["stage_name"] for stage in self.stages]

    @property
    def nstages(self) -> int:
        """Returns the number of stages in the current LammpsInputFile."""
        return len(self.stages)

    @property
    def ncomments(self) -> int:
        """
        Returns the number of comments in the current LammpsInputFile. Includes the blocks of comments as well
        as inline comments (comment lines within blocks of LAMMPS commands).
        """
        ncomments = 0
        for stage in self.stages:
            # Block of comment = 1 comment
            if all(cmd.strip().startswith("#") for (cmd, args) in stage["commands"]):
                ncomments += 1
            else:
                # Else, inline comment each count as one
                ncomments += sum(1 for cmd, args in stage["commands"] if cmd.strip().startswith("#"))

        return ncomments

    def get_args(self, command: str, stage_name: str | None = None) -> list | str:
        """
        Given a command, returns the corresponding arguments (or list of arguments) in the LammpsInputFile.
        A stage name can be given; in this case the search will happen only for this stage.
        If a stage name is not given, the command will be searched for through all of them.
        If the command is not found, an empty list is returned.

        Args:
            command (str): String with the command to find in the input file (e.g., "units").
            stage_name (str): String giving the stage name where the change should take place.

        Returns:
            Value of the argument corresponding to the command.
            List if the same command is used multiple times.
        """
        args = []
        stages_to_look = [stage_name] if stage_name else self.stages_names

        for stage in self.stages:
            if stage["stage_name"] in stages_to_look:
                for cmd, arg in stage["commands"]:
                    if command == cmd:
                        args.append(arg)

        return args if len(args) != 1 else args[0]

    def contains_command(self, command: str, stage_name: str | None = None) -> bool:
        """
        Returns whether a given command is present in the LammpsInputFile.
        A stage name can be given; in this case the search will happen only for this stage.

        Args:
            command (str): String with the command to find in the input file (e.g., "units").
            stage_name (str): String giving the stage name where the change should take place.

        Returns:
            bool: True if the command is present, False if not.
        """
        return bool(self.get_args(command, stage_name))

    def set_args(self, command: str, argument: str, stage_name: str | None = None, how: str | int | list[int] = "all"):
        """
        Sets the arguments for the given command to the given string.
        If the command is not found, nothing is done. Use LammpsInputFile.add_commands instead.
        If a stage name is specified, it will be replaced or set only for this stage.
        If no stage name is given, it will apply the change in all of them that contain the given command.
        If the command is set multiple times in the file/stage, it will be replaced based on "how":
        either the first occurrence, all of them, or the index of the occurrence.

        Args:
            command (str): String representing the command to change, e.g., "units".
            argument (str): String with the new value for the command, e.g., "atomic".
            stage_name (str): String giving the stage name where the change should take place.
            how (str or int or list): "all" for changing all occurrences of the command within the stage_name
                or the whole InputFile, "first" for the first occurrence, int i for the i-th time the command
                is present in the stage_name or the whole InputFile, starting at 0. Can be a list of indexes as well.
        """
        # Get the stages to look in
        stages_to_look = [stage_name] if stage_name else self.stages_names

        # Translates how into range of indices
        if how == "first":
            how = [0]
        elif how == "all":
            getargs = self.get_args(command, stage_name)
            N = 1 if isinstance(getargs, str) else len(getargs)
            how = list(range(N))
        elif isinstance(how, int):
            how = [how]
        elif not (isinstance(how, list) and all(isinstance(h, int) for h in how)):
            raise ValueError("""The argument 'how' should be a 'first', 'all', an integer or a list of integers.""")

        # Look for occurrences in the relevant stages
        i = 0
        for i_stage, stage in enumerate(self.stages):
            if stage["stage_name"] in stages_to_look:
                for i_cmd, (cmd, _) in enumerate(stage["commands"]):
                    if command == cmd:
                        if i in how:
                            self.stages[i_stage]["commands"][i_cmd] = (cmd, argument)
                        i += 1

    def add_stage(
        self,
        stage: dict | None = None,
        commands: str | list[str] | dict[str, str | float] | None = None,
        stage_name: str | None = None,
        after_stage: str | int | None = None,
    ):
        r"""
        Adds a new stage to the LammpsInputFile, either from a whole stage (dict) or
        from a stage_name and commands. Both ways are mutually exclusive.

        Examples:
            1) In order to add a stage defining the force field to be used, you can use:
            ```
            your_input_file.add_stage(
                commands=["pair_coeff 1 1 morse 0.0580 3.987 3.404", "pair_coeff 1 4 morse 0.0408 1.399 3.204"],
                stage_name="Definition of the force field"
            )
            ```
            or
            ```
            your_input_file.add_stage(
                {
                    "stage_name": "Definition of the force field",
                    "commands": [
                        ("pair_coeff", "1 1 morse 0.0580 3.987 3.404"),
                        ("pair_coeff", "1 4 morse 0.0408 1.399 3.204")
                    ],
                }
            )
            ```
            2) Another stage could consist in an energy minimization. In that case, the commands could look like
            ```
            commands = [
                "thermo 1",
                "thermo_style custom step lx ly lz press pxx pyy pzz pe",
                "dump dmp all atom 5 run.dump",
                "min_style cg",
                "fix 1 all box/relax iso 0.0 vmax 0.001",
                "minimize 1.0e-16 1.0e-16 5000 10000",
                "write_data run.data"
            ]
            ```
            or a dictionary such as `{"thermo": 1, ...}`, or a string with a single command (e.g., "units atomic").

        Args:
            stage (dict): if provided, this is the stage that will be added to the LammpsInputFile.stages
            commands (str or list or dict): if provided, these are the LAMMPS command(s)
                that will be included in the stage to add.
                Can pass a list of LAMMPS commands with their arguments.
                Also accepts a dictionary of LAMMPS commands and
                corresponding arguments as key, value pairs.
                A single string can also be passed (single command together with its arguments).
                Not used in case a whole stage is given.
            stage_name (str): If a stage name is mentioned, the commands are added
                under that stage block, else the new stage is named from numbering.
                If given, stage_name cannot be one of the already present stage names.
                Not used in case a whole stage is given.
            after_stage (str): Name of the stage after which this stage should be added.
                If None, the stage is added at the end of the LammpsInputFile.
        """
        # Get the index of the stage to add
        if after_stage is None:
            index_insert = -1
        elif isinstance(after_stage, int):
            index_insert = after_stage + 1
        elif after_stage in self.stages_names:
            index_insert = self.stages_names.index(after_stage) + 1
            if index_insert == len(self.stages_names):
                index_insert = -1
        else:
            raise ValueError("The stage after which this one should be added does not exist.")

        # Adds a stage depending on the given inputs.
        # If a stage is given, we check its format and insert it where it should directly.
        if stage:
            if commands or stage_name:
                warnings.warn(
                    "A stage has been passed together with commands and stage_name. This is incompatible. "
                    "Only the stage will be used."
                )

            # Make sure the given stage has the correct format
            if stage["stage_name"] in self.stages_names:
                raise ValueError("The provided stage name is already present in LammpsInputFile.stages.")
            self._check_stage_format(stage)

            # Insert the stage in the LammpsInputFile.stages
            if index_insert == -1:
                self.stages.append(stage)
            else:
                self.stages.insert(index_insert, stage)
        else:
            # Initialize the stage (even if no commands are given)
            self._initialize_stage(stage_name=stage_name, stage_index=index_insert)
            if commands:
                # Adds the commands to the created stage
                self.add_commands(stage_name=self.stages_names[index_insert], commands=commands)

    def remove_stage(self, stage_name: str):
        """
        Removes a whole stage from the LammpsInputFile.

        Args:
            stage_name (str): name of the stage to remove.
        """
        if stage_name in self.stages_names:
            idx = self.stages_names.index(stage_name)
            self.stages.pop(idx)
        else:
            raise LookupError("The given stage name is not present in this LammpsInputFile.")

    def rename_stage(self, stage_name: str, new_name: str):
        """
        Renames a stage `stage_name` from LammpsInputFile into `new_name`.
        First checks that the stage to rename is present, and that
        the new name is not already a stage name.

        Args:
            stage_name (str): name of the stage to rename.
            new_name (str): new name of the stage.
        """
        if stage_name in self.stages_names:
            if new_name in self.stages_names:
                raise ValueError("The provided stage name is already present in LammpsInputFile.stages.")
            idx = self.stages_names.index(stage_name)
            self.stages[idx]["stage_name"] = new_name
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
        if not all(stage in self.stages_names for stage in stage_names):
            raise ValueError("At least one of the stages to be merged is not in the LammpsInputFile.")

        indices_stages_to_merge = [self.stages_names.index(stage) for stage in stage_names]
        if not np.all([np.array(indices_stages_to_merge[1:]) >= np.array(indices_stages_to_merge[:-1])]):
            raise ValueError(
                """The provided stages are not in the right order. You should merge stages in the order of appearance
                in your LammpsInputFile. If you want to reorder stages, modify LammpsInputFile.stages directly. """
            )

        stages = self.stages[: indices_stages_to_merge[0]]

        merge_name = "Merge of: " + ", ".join([self.stages_names[i] for i in indices_stages_to_merge])
        merged_commands = []
        for i in indices_stages_to_merge:
            for j in range(len(self.stages[i]["commands"])):
                merged_commands.append(self.stages[i]["commands"][j])

        merged_stages = {"stage_name": merge_name, "commands": merged_commands}
        stages.append(merged_stages)

        for i_stage, stage in enumerate(self.stages):
            if i_stage > indices_stages_to_merge[0] and i_stage not in indices_stages_to_merge:
                stages.append(stage)

        self.stages = stages

    def add_commands(self, stage_name: str, commands: str | list[str] | dict):
        """
        Method to add a LAMMPS commands and their arguments to a stage of
        the LammpsInputFile. The stage name should be provided: a default behavior
        is avoided here to avoid mistakes (e.g., the commands are added to the wrong stage).

        Example:
            In order to add the command ``pair_coeff 1 1 morse 0.0580 3.987 3.404``
            to the stage "Definition of the potential", simply use
            ```
            your_input_file.add_commands(
                stage_name="Definition of the potential",
                commands="pair_coeff 1 1 morse 0.0580 3.987 3.404"
            )
            ```
            To add multiple commands, use a dict or a list, e.g.,
            ```
            your_input_file.add_commands(
                stage_name="Definition of the potential",
                commands=["pair_coeff 1 1 morse 0.0580 3.987 3.404", "units atomic"]
            )
            your_input_file.add_commands(
                stage_name="Definition of the potential",
                commands={"pair_coeff": "1 1 morse 0.0580 3.987 3.404", "units": "atomic"}
            )
            ```

        Args:
            stage_name (str): name of the stage to which the command should be added.
            commands (str or list or dict): LAMMPS command, with or without the arguments.
        """
        if stage_name not in self.stages_names:
            raise ValueError("The provided stage name does not correspond to one of the LammpsInputFile.stages.")

        # Handle the different input formats to add commands to the stage
        if isinstance(commands, str):
            self._add_command(command=commands, stage_name=stage_name)

        elif isinstance(commands, list):
            for comm in commands:
                if comm[0] == "#":
                    self._add_comment(comment=comm[1:].strip(), inline=True, stage_name=stage_name, index_comment=True)
                else:
                    self._add_command(command=comm, stage_name=stage_name)

        elif isinstance(commands, dict):
            for comm, args in commands.items():
                if comm[0] == "#":
                    self._add_comment(comment=comm[1:].strip(), inline=True, stage_name=stage_name, index_comment=True)
                else:
                    self._add_command(command=comm, args=args, stage_name=stage_name)

        else:
            raise TypeError("The command should be a string, list of strings or dictionary.")

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
        for i_stage, stage in enumerate(self.stages):
            if stage["stage_name"] in stage_name:
                new_commands = []
                for i_cmd, (cmd, arg) in enumerate(stage["commands"]):
                    if cmd == command:
                        n_removed += 1
                    else:
                        new_commands.append((cmd, arg))
                        indices_to_remove.append([i_stage, i_cmd])
                if new_commands or not remove_empty_stages:
                    new_list_of_stages.append({"stage_name": stage["stage_name"], "commands": new_commands})
            else:
                new_list_of_stages.append(stage)

        self.stages = new_list_of_stages

        if n_removed == 0:
            warnings.warn(f"{command} not found in the LammpsInputFile.")

    def append(self, lmp_input_file: LammpsInputFile):
        """
        Appends a LammpsInputFile to another. The `stages` are merged,
        and the numbering of stages/comments is either kept the same or updated.

        Args:
            lmp_input_file (LammpsInputFile): LammpsInputFile to append.
        """
        # Renumbering comments and stages of the lmp_input_file.stages
        new_list_to_add = lmp_input_file.stages
        for i_stage, stage in enumerate(lmp_input_file.stages):
            stage_name = stage["stage_name"]
            if stage_name.split()[0].strip() == "Comment" and stage_name.split()[1].isdigit():
                i_comment = int(stage_name.split()[1]) + self.ncomments
                new_list_to_add[i_stage]["stage_name"] = f"Comment {i_comment}"
            if stage_name.split()[0] == "Stage" and stage_name.split()[1].isdigit():
                this_stage = int(stage_name.split()[1]) + self.nstages
                new_list_to_add[i_stage]["stage_name"] = f"Stage {this_stage}"

        # Making sure no stage_name of lmp_input_file clash with those from self.
        # If it is the case, we rename them.
        for i_stage, stage in enumerate(lmp_input_file.stages):
            if stage["stage_name"] in self.stages_names:
                stage["stage_name"] = f"Stage {self.nstages + i_stage + 1} (previously {stage['stage_name']})"

        # Append the two list of stages
        self.stages += new_list_to_add

    @np.deprecate(message="Use get_str instead")
    def get_string(self, *args, **kwargs) -> str:
        return self.get_str(*args, **kwargs)

    def get_str(self, ignore_comments: bool = False, keep_stages: bool = True) -> str:
        """
        Generates and ² the string representation of the LammpsInputFile.
        Stages are separated by empty lines.
        The headers of the stages will be put in comments preceding each stage.
        Other comments will be put inline within stages, where they have been added.

        Args:
            ignore_comments (bool): True if only the commands should be kept from the InputFile.
            keep_stages (bool): If True, the string is formatted in a block structure with stage names
                and newlines that differentiate commands in the respective stages of the InputFile.
                If False, stage names are not printed and all commands appear in a single block.

        Returns:
            str: String representation of the LammpsInputFile.
        """
        lammps_input = f"# LAMMPS input generated from LammpsInputFile with pymatgen v{CURRENT_VER}\n"
        if not keep_stages:
            lammps_input += "\n"

        for stage in self.stages:
            if keep_stages:
                # Print first the name of the stage in a comment.
                # We print this even if ignore_comments is True.
                if "Comment" not in stage["stage_name"] and len(self.stages) > 1:
                    lammps_input += "\n# " + stage["stage_name"] + "\n"

                # In case of a block of comment, the header is not printed (useless)
                else:
                    lammps_input += "\n"

            # Then print the commands and arguments (including inline comments)
            for command, args in stage["commands"]:
                if not (ignore_comments and "#" in command):
                    lammps_input += f"{command} {args.strip()}\n"

        return lammps_input

    def write_file(self, filename: str | Path, ignore_comments: bool = False, keep_stages: bool = True) -> None:
        """
        Writes the input file.

        Args:
            filename (str or path): The filename to output to, including path.
            ignore_comments (bool): True if only the commands should be kept from the InputFile.
            keep_stages (bool): True if the block structure from the InputFile should be kept.
                                If False, a single block is assumed.
        """
        filename = filename if isinstance(filename, Path) else Path(filename)
        with zopen(filename, "wt") as f:
            f.write(self.get_str(ignore_comments=ignore_comments, keep_stages=keep_stages))

    @classmethod
    @np.deprecate(message="Use from_str instead")
    def from_string(cls, *args, **kwargs):
        return cls.from_str(*args, **kwargs)

    @classmethod
    def from_str(cls, contents: str, ignore_comments: bool = False, keep_stages: bool = False) -> LammpsInputFile:
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
            if all(line[0] == "#" for line in block):
                if ignore_comments:
                    keep_block = False
                else:
                    LIF._add_comment(comment=block[0][1:].strip(), inline=False)
                    stage_name = f"Comment {LIF.ncomments}"
                    if len(block) > 1:
                        for line in block[1:]:
                            LIF._add_comment(comment=line[1:].strip(), inline=True, stage_name=stage_name)

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
                stage_name = f"Stage {LIF.nstages + 1}" if (ignore_comments or not keep_stages) else header
                commands = block[icomm_max:]
                LIF.add_stage(commands=commands, stage_name=stage_name)

            # Stage with no header
            else:
                stage_name = f"Stage {LIF.nstages + 1}"
                LIF.add_stage(commands=block, stage_name=stage_name)
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
            return cls.from_str(f.read(), ignore_comments=ignore_comments, keep_stages=keep_stages)

    def __repr__(self):
        return self.get_str()

    def _initialize_stage(self, stage_name: str | None = None, stage_index: int | None = None):
        """
        Initialize an empty stage with the given name in the LammpsInputFile.

        Args:
            stage_name (str): If a stage name is mentioned, the command is added
                under that stage block, else the new stage is named from numbering.
                If given, stage_name cannot be one of the already present stage names.
            stage_index (int): Index of the stage where it should be added.
                If None, the stage is added at the end of the LammpsInputFile.
        """
        if stage_name is None:
            stage_name = f"Stage {self.nstages + 1}"
            if stage_name in self.stages_names:
                # This situation can happen when some stages have been removed, then others added
                stage_numbers = [
                    int(stage_name.split()[1])
                    for stage_name in self.stages_names
                    if stage_name.split()[0] == "Stage" and stage_name.split()[1].isdigit()
                ]
                stage_name = f"Stage {np.max(stage_numbers) + 1}"

        if not isinstance(stage_name, str):
            raise TypeError("Stage names should be strings.")

        if stage_name in self.stages_names:
            raise ValueError("The provided stage name is already present in LammpsInputFile.stages.")

        # Initialize the stage
        if stage_index is None or stage_index == -1:
            self.stages.append({"stage_name": stage_name, "commands": []})
        elif stage_index > len(self.stages):
            raise IndexError("The provided index is too large with respect to the current number of stages.")
        else:
            self.stages.insert(stage_index, {"stage_name": stage_name, "commands": []})

    @staticmethod
    def _check_stage_format(stage: dict):
        if list(stage) != ["stage_name", "commands"]:
            raise KeyError(
                "The provided stage does not have the correct keys. It should be 'stage_name' and 'commands'."
            )
        if not isinstance(stage["stage_name"], str):
            raise TypeError("The value of 'stage_name' should be a string.")
        if not isinstance(stage["commands"], list):
            raise TypeError("The provided commands should be a list.")
        if len(stage["commands"]) >= 1 and (
            not all(len(cmdargs) == 2 for cmdargs in stage["commands"])
            or not all(isinstance(cmd, str) and isinstance(arg, str) for (cmd, arg) in stage["commands"])
        ):
            raise ValueError("The provided commands should be a list of 2-strings tuples.")

    def _add_command(self, stage_name: str, command: str, args: str | float | None = None):
        """
        Helper method to add a single LAMMPS command and its arguments to
        the LammpsInputFile. The stage name should be provided: a default behavior
        is avoided here to avoid mistakes.

        Example:
            In order to add the command ``pair_coeff 1 1 morse 0.0580 3.987 3.404``
            to the stage "Definition of the potential", simply use
            ```
            your_input_file._add_command(
                stage_name="Definition of the potential",
                command="pair_coeff 1 1 morse 0.0580 3.987 3.404"
            )
            ```
            or
            ```
            your_input_file._add_command(
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
        if not self.stages[idx]["commands"]:
            self.stages[idx]["commands"] = [(command, args)]
        else:
            self.stages[idx]["commands"].append((command, args))

    def _add_comment(
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
        # "Stage" of comments only
        if not inline:
            if stage_name is None:
                stage_name = f"Comment {self.ncomments + 1}"
                self.stages.append({"stage_name": stage_name, "commands": [("#", comment)]})
            elif stage_name in self.stages_names:
                self._add_command(command="#", args=comment, stage_name=stage_name)
            else:
                self.stages.append({"stage_name": stage_name, "commands": [("#", comment)]})

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
            self._add_command(command=command, args=args, stage_name=stage_name)

        else:
            raise NotImplementedError("If you want to add an inline comment, please specify the stage name.")

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
        if len(string_list) == 0 or all(s == "" for s in string_list):
            raise ValueError("The list of strings should contain some non-empty strings.")

        # Remove the first comment line possibly added by previous input creations
        while "# LAMMPS input generated from LammpsInputFile" in string_list:
            string_list.remove("# LAMMPS input generated from LammpsInputFile")

        # Get rid of possible initial empty lines
        # and final empty lines
        imin = len(string_list)
        imax = 0
        for idx, string in enumerate(string_list):
            if string != "" and idx <= imin:
                imin = idx
            if string != "" and idx >= imax:
                imax = idx
        string_list = string_list[imin : imax + 1]

        # Get rid of empty comments that are there just for cosmetic reasons
        new_list = []
        for string in string_list:
            if len(string) > 1 or not (len(string.strip()) == 1 and string[0] == "#"):
                new_list.append(string)
        string_list = new_list

        # Keep only a single empty lines when there are multiple ones
        lines = [string_list[0]]

        for idx, string in enumerate(string_list[1:-1]):
            if (
                string != ""
                and not (string[0] == "#" and ignore_comments)
                or string == ""
                and string_list[idx + 2] != ""
            ):
                lines.append(string)

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
