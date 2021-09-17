"""
This module defines the abstract interface for pymatgen InputSet and InputSetGenerator classes.
"""

import abc
from pathlib import Path
from string import Template
from typing import Union, Optional, Dict
from monty.json import MSONable
from monty.io import zopen


__author__ = "Ryan Kingsbury"
__email__ = "RKingsbury@lbl.gov"
__status__ = "Development"
__date__ = "August 2021"


class InputSet(MSONable):
    """
    Abstract base class for all InputSet classes. InputSet classes serve
    as containers for all calculation input data.

    All InputSet must implement a write_inputs method. Implementing
    the from_files and validate methods is optional.
    """

    @abc.abstractmethod
    def write_inputs(
        self,
        directory: Union[str, Path],
        make_dir: bool = True,
        overwrite: bool = True,
    ):
        """
        Write Inputs to one or more files

        Args:
            directory: Directory to write input files to
            make_dir: Whether to create the directory if it does not already exist.
            overwrite: Whether to overwrite an input file if it already exists.
        """
        pass

    @classmethod
    @abc.abstractmethod
    def from_directory(cls, directory: Union[str, Path]):
        """
        Construct an InputSet from a directory of one or more files.

        Args:
            directory: Directory to read input files from
        """
        raise NotImplementedError(f"from_directory has not been implemented in {cls}")

    def validate(self) -> bool:
        """
        A place to implement basic checks to verify the validity of an
        input set. Can be as simple or as complex as desired.

        Will raise a NotImplementedError unless overloaded by the inheriting class.
        """
        raise NotImplementedError(f".validate() has not been implemented in {self.__class__}")


class InputSetGenerator(MSONable):
    """
    InputSetGenerator classes serve as generators for Input objects. They contain
    settings or sets of instructions for how to create Input from a set of
    coordinates or a previous calculation directory.
    """

    @abc.abstractmethod
    def get_input_set(self) -> InputSet:
        """
        Generate an InputSet object. Typically the first argument to this method
        will be a Structure or other form of atomic coordinates.
        """
        pass


class TemplateInputSet(InputSet):
    """
    Concrete implementation of InputSet that is based on a single template input
    file with variables.

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

        # load the template
        with zopen(self.template, "r") as f:
            template_str = f.read()

        # replace all variables
        self.data = Template(template_str).safe_substitute(**self.variables)

    def write_inputs(
        self,
        directory: Union[str, Path],
        make_dir: bool = True,
        overwrite: bool = True,
        filename: str = "input.txt",
    ):
        """
        Args:
            directory: Directory to write input files to
            make_dir: Whether to create the directory if it does not already exist.
            overwrite: Whether to overwrite an input file if it already exists.
            filename: name of the file to be written
        """
        path = directory if isinstance(directory, Path) else Path(directory)
        file = path / filename

        if not path.exists():
            if make_dir:
                path.mkdir(parents=True, exist_ok=True)

        if file.exists() and not overwrite:
            raise FileExistsError(f"File {str(filename)} already exists!")
        file.touch()

        # write the file
        with zopen(file, "wt") as f:
            f.write(self.data)
    
    @classmethod
    def from_directory(cls, directory: Union[str, Path]):
        """
        Construct an InputSet from a directory of one or more files.

        Args:
            directory: Directory to read input files from
        """
        raise NotImplementedError(f"from_directory has not been implemented in {cls}")
