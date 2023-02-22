# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module defines the abstract interface for reading and writing calculation
inputs in pymatgen. The interface comprises a 3-tiered hierarchy of classes.

1. An InputFile object represents the contents of a single input file, e.g.
   the INCAR. This class standardizes file read and write operations.
2. An InputSet is a dict-like container that maps filenames (keys) to file
   contents (either strings or InputFile objects). This class provides a standard
   write_input() method.
3. InputGenerator classes implement a get_input_set method that, when provided
   with a structure, return an InputSet object with all parameters set correctly.
   Calculation input files can be written to disk with the write_inputs method.

If you want to implement a new InputGenerator, please take note of the following:

1. You must implement a get_input_set method that returns an InputSet
2. All customization of calculation parameters should be done in the __init__
   method of the InputGenerator. The idea is that the generator contains
   the "recipe", but nothing that is specific to a particular system. get_input_set
   takes system-specific information (such as structure) and applies the recipe.
3. All InputGenerator must save all supplied args and kwargs as instance variables.
   E.g., self.my_arg = my_arg and self.kwargs = kwargs in the __init__. This
   ensures the as_dict and from_dict work correctly.
"""

from __future__ import annotations

import abc
import os
from collections.abc import MutableMapping
from pathlib import Path
from zipfile import ZipFile

from monty.io import zopen
from monty.json import MSONable

__author__ = "Ryan Kingsbury"
__email__ = "RKingsbury@lbl.gov"
__status__ = "Development"
__date__ = "October 2021"


class InputFile(MSONable):
    """
    Abstract base class to represent a single input file. Note that use
    of this class is optional; it is possible create an InputSet that
    does not rely on underlying Inputfile objects.

    All InputFile classes must implement a get_string method, which
    is called by write_file.

    If InputFile classes implement an __init__ method, they must assign all
    arguments to __init__ as attributes.
    """

    @abc.abstractmethod
    def get_string(self) -> str:
        """
        Return a string representation of an entire input file.
        """

    def write_file(self, filename: str | Path) -> None:
        """
        Write the input file.

        Args:
            filename: The filename to output to, including path.
            kwargs: Keyword arguments passed to get_string()
        """
        filename = filename if isinstance(filename, Path) else Path(filename)
        with zopen(filename, "wt") as f:
            f.write(self.get_string())

    @classmethod
    @abc.abstractmethod
    def from_string(cls, contents: str):
        """
        Create an InputFile object from a string

        Args:
            contents: The contents of the file as a single string

        Returns:
            InputFile
        """

    @classmethod
    def from_file(cls, path: str | Path):
        """
        Creates an InputFile object from a file.

        Args:
            path: Filename to read, including path.

        Returns:
            InputFile
        """
        filename = path if isinstance(path, Path) else Path(path)
        with zopen(filename, "rt") as f:
            return cls.from_string(f.read())


class InputSet(MSONable, MutableMapping):
    """
    Abstract base class for all InputSet classes. InputSet are dict-like
    containers for all calculation input data.

    Since InputSet inherits dict, it can be instantiated in the same manner,
    or a custom __init__ can be provided. Either way, `self` should be
    populated with keys that are filenames to be written, and values that are
    InputFile objects or strings representing the entire contents of the file.

    All InputSet must implement from_directory. Implementing the validate method
    is optional.
    """

    def __init__(self, inputs: dict[str | Path, str | InputFile] | None = None, **kwargs):
        """
        Instantiate an InputSet.

        Args:
            inputs: The core mapping of filename: file contents that defines the InputSet data.
                This should be a dict where keys are filenames and values are InputFile objects
                or strings representing the entire contents of the file. If a value is not an
                InputFile object nor a str, but has a __str__ method, this str representation
                of the object will be written to the corresponding file. This mapping will
                become the .inputs attribute of the InputSet.
            **kwargs: Any kwargs passed will be set as class attributes e.g.
                InputSet(inputs={}, foo='bar') will make InputSet.foo == 'bar'.
        """
        self.inputs = inputs or {}
        self._kwargs = kwargs
        self.__dict__.update(**kwargs)

    def __getattr__(self, k):
        # allow accessing keys as attributes
        if k in self._kwargs:
            return self.get(k)
        raise AttributeError(f"'{type(self).__name__}' object has no attribute {k!r}")

    def __len__(self):
        return len(self.inputs)

    def __iter__(self):
        return iter(self.inputs)

    def __getitem__(self, key):
        return self.inputs[key]

    def __setitem__(self, key, value):
        self.inputs[key] = value

    def __delitem__(self, key):
        del self.inputs[key]

    def write_input(
        self,
        directory: str | Path,
        make_dir: bool = True,
        overwrite: bool = True,
        zip_inputs: bool = False,
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

        for fname, contents in self.inputs.items():
            file = path / fname

            if not path.exists() and make_dir:
                path.mkdir(parents=True, exist_ok=True)

            if file.exists() and not overwrite:
                raise FileExistsError(f"File {str(fname)} already exists!")
            file.touch()

            # write the file
            if isinstance(contents, InputFile):
                contents.write_file(file)
            else:
                with zopen(file, "wt") as f:
                    f.write(str(contents))

        if zip_inputs:
            zipfilename = path / f"{type(self).__name__}.zip"
            with ZipFile(zipfilename, "w") as zip:
                for fname in self.inputs:
                    file = path / fname
                    try:
                        zip.write(file)
                        os.remove(file)
                    except FileNotFoundError:
                        pass

    @classmethod
    def from_directory(cls, directory: str | Path):
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


class InputGenerator(MSONable):
    """
    InputGenerator classes serve as generators for Input objects. They contain
    settings or sets of instructions for how to create Input from a set of
    coordinates or a previous calculation directory.
    """

    @abc.abstractmethod
    def get_input_set(self) -> InputSet:
        """
        Generate an InputSet object. Typically the first argument to this method
        will be a Structure or other form of atomic coordinates.
        """
