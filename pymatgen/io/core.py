"""
This module defines the abstract interface for pymatgen InputSet and InputSetGenerator classes.
"""

import abc
import os
from pathlib import Path
from collections.abc import Set
from string import Template
from typing import Union, Optional, Dict
from zipfile import ZipFile
from monty.json import MSONable
from monty.io import zopen


__author__ = "Ryan Kingsbury"
__email__ = "RKingsbury@lbl.gov"
__status__ = "Development"
__date__ = "September 2021"


class InputFile(MSONable):
    """
    Abstract base class to represent a single input file. Note that use
    of this class is optional; it is possible create an InputSet that
    does not rely on underlying Inputfile objects.

    All InputFile classes must implement a __str__ method, which
    is called by write_file.
    """

    @abc.abstractmethod
    def __str__(self):
        """
        String representation of an entire input file.
        """

    def write_file(self, filename: Union[str, Path]):
        """
        Write the input file.

        Args:
            filename: The filename to output to, including path.
        """
        filename = filename if isinstance(filename, Path) else Path(filename)
        with open(filename, "wt") as f:
            f.write(self.__str__())

    @abc.abstractmethod
    @staticmethod
    def from_string(contents):
        """
        Create an InputFile object from a string

        Args:
            contents: The contents of the file as a single string

        Returns:
            InputFile
        """

    @classmethod
    def from_file(cls, filename: Union[str, Path]):
        """
        Creates an InputFile object from a file.

        Args:
            filename: Filename to read, including path.

        Returns:
            InputFile
        """
        filename = filename if isinstance(filename, Path) else Path(filename)
        with zopen(filename, "rt") as f:
            return cls.from_string(f.read())


class InputSet(MSONable, Set):
    """
    Abstract base class for all InputSet classes. InputSet classes serve
    as containers for all calculation input data.

    All InputSet must implement a _get_inputs and from_directory method.
    Implementing the validate method is optional.
    """

    @abc.abstractmethod
    def _get_inputs(self) -> Dict[str, Union[str, InputFile]]:
        """
        Generate a dictionary of one or more input files to be written. Keys
        are filenames, values are InputFile objects or strings representing
        the entire contents of the file.

        This method is called by write_input(), which performs the actual file
        write operations.
        """
        pass

    def write_input(
        self,
        directory: Union[str, Path],
        make_dir: bool = True,
        overwrite: bool = True,
        zip_inputs: bool = False,
        **kwargs,
    ):
        """
        Write Inputs to one or more files

        Args:
            directory: Directory to write input files to
            make_dir: Whether to create the directory if it does not already exist.
            overwrite: Whether to overwrite an input file if it already exists.
            Additional kwargs are passed to generate_inputs
            zip_inputs: If True, inputs will be zipped into a file with the
                same name as the InputSet (e.g., InputSet.zip)
        """
        path = directory if isinstance(directory, Path) else Path(directory)
        # the following line will trigger a mypy error due to a bug in mypy
        # will be fixed soon. See https://github.com/python/mypy/commit/ea7fed1b5e1965f949525e918aa98889fb59aebf
        files = self._get_inputs(**kwargs)  # type: ignore
        for fname, contents in files.items():
            file = path / fname

            if not path.exists():
                if make_dir:
                    path.mkdir(parents=True, exist_ok=True)

            if file.exists() and not overwrite:
                raise FileExistsError(f"File {str(fname)} already exists!")
            file.touch()

            # write the file
            if isinstance(contents, str):
                with zopen(file, "wt") as f:
                    f.write(contents)
            else:
                contents.write_file(file)

        if zip_inputs:
            zipfilename = path / f"{self.__class__.__name__}.zip"
            with ZipFile(zipfilename, "w") as zip:
                for fname, contents in files.items():
                    file = path / fname
                    try:
                        zip.write(file)
                        os.remove(file)
                    except FileNotFoundError:
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

    def __len__(self):
        return len(self._get_inputs().keys())

    def __contains__(self, item):
        return item in self._get_inputs().values()

    def __iter__(self):
        for k, v in self._get_inputs().items():
            yield k, v


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

    def _get_inputs(self, filename: str = "input.txt"):
        """
        Args:
            filename: name of the file to be written
        """
        return {filename: self.data}

    @classmethod
    def from_directory(cls, directory: Union[str, Path]):
        """
        Construct an InputSet from a directory of one or more files.

        Args:
            directory: Directory to read input files from
        """
        raise NotImplementedError(f"from_directory has not been implemented in {cls}")
