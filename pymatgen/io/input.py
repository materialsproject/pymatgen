"""
This module defines the abstract interface for pymatgen Input classes.
"""

import abc
from pathlib import Path
from string import Template
from typing import Union, Optional, Dict
from monty.json import MSONable
from monty.io import zopen
from pymatgen.core import Structure, Molecule


__author__ = "Ryan Kingsbury"
__email__ = "RKingsbury@lbl.gov"
__status__ = "Development"
__date__ = "August 2021"


class Input(MSONable):
    """
    Abstract base class for all Input classes. Input classes serve
    as containers for all calculation input data.

    All Input must implement a write_inputs method. Implementing
    the from_files method and is_valid property are optional.
    """

    @abc.abstractmethod
    def write_inputs(
        self,
        dir: Union[str, Path],
        make_dir_if_not_present: bool = True,
        overwrite_if_present: bool = True,
    ):
        """
        Write Inputs to one or more files

        Args:
            dir: Directory to write input files to
            make_dir_if_not_present: Whether to create the directory if it does not already exist.
            overwrite_if_present: Whether to overwrite an input file if it already exists.
        """
        pass

    @classmethod
    @abc.abstractmethod
    def from_files(cls, dir: Union[str, Path]):
        """
        Construct an Input from a directory of one or more files.

        Args:
            dir: Directory to read input files from
        """
        raise NotImplementedError(f"from_files has not been implemented in {cls}")

    @property
    def is_valid(self) -> bool:
        """
        A place to implement basic checks to verify the validity of an
        input set. Can be as simple or as complex as desired.

        Will raise a NotImplementedError unless overloaded by the inheriting class.
        """
        raise NotImplementedError(f"is_valid has not been implemented in {self.__class__}")


class InputSettings(MSONable):
    """
    InputSettings classes serve as generators for Input objects. They contain
    settings or sets of instructions for how to create Input from a set of
    coordinates or a previous calculation directory.
    """

    @staticmethod
    @abc.abstractmethod
    def get_input(coordinates: Union[Structure, Molecule], prev_dir: Union[str, Path] = None) -> Input:
        """
        Generate an Input object

        Args:
            coordinates: Atomic coordinates. May be passed as a Structure, Molecule, or path to
                any file that can be converted to one via Structure.from_file() or Molecule.from_file().
            prev_dir: Location of a previous calculation to inherit settings from.
        """
        pass


class TemplateInput(Input):
    """
    Concrete implementation of Input that is based on a single template input
    file with variables. The .write_inputs() method simply substitutes values into
    the variables and writes out a new file.

    This class is provided as a low-barrier way to support new codes and to provide
    an intuitive way for users to transition from manual scripts to pymatgen I/O
    classes.
    """

    def __init__(self, template: Union[str, Path], variables: Optional[Dict] = None):
        """
        Args:
            template: the input file template containing variable strings to be
                replaced.
            variables: dict of variables to replace in the template. Keys are the
                text to replaced with the values, e.g. {"TEMPERATURE": 298} will
                replace the text $TEMPERATURE in the template. See Python's
                Template.safe_substitute() method documentation for more details.
        """
        self.template = template
        self.variables = variables if variables else {}

    def write_inputs(
        self,
        dir: Union[str, Path],
        make_dir_if_not_present: bool = True,
        overwrite_if_present: bool = True,
        filename: str = "input.txt",
    ):
        """
        Args:
            dir: Directory to write input files to
            make_dir_if_not_present: Whether to create the directory if it does not already exist.
            overwrite_if_present: Whether to overwrite an input file if it already exists.
            filename: name of the file to be written
        """
        path = dir if isinstance(dir, Path) else Path(dir)
        file = path / filename

        if not path.exists():
            if make_dir_if_not_present:
                path.parent.mkdir(parents=True, exist_ok=True)

        if file.exists() and not overwrite_if_present:
            raise FileExistsError(f"File {str(filename)} already exists!")
        file.touch()

        # load the template
        with zopen(self.template, "r") as f:
            template_str = f.read()

        # replace all variables
        to_write = Template(template_str).safe_substitute(**self.variables)

        # write the file
        with zopen(file, "wt") as f:
            f.write(to_write)
