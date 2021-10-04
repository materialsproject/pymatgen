# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Input sets for LAMMPS
"""

import logging
from pathlib import Path
from typing import Union, Dict

from pymatgen.io.core import InputSet, InputSetGenerator

__author__ = "Ryan Kingsbury, ..."
__copyright__ = "Copyright 2021, The Materials Project"
__version__ = "0.1"

logger = logging.getLogger(__name__)


class LammpsInputSet(InputSet):
    """
    Container class for all LAMMPS inputs. This class is intended to provide
    general functionality that can be customized to many purposes.

    InputSetGenerator classes elsewhere in this module are used to create
    specific instances of LammpsInputSet that are tailored to specific purposes.
    """

    def get_inputs(self) -> Dict[str, str]:
        """
        Generate a dictionary of one or more input files to be written. Keys
        are filenames, values are the contents of each file.

        This method is called by write_input(), which performs the actual file
        write operations.
        """
        pass

    @classmethod
    def from_directory(cls, directory: Union[str, Path]):
        """
        Construct an InputSet from a directory of one or more files.

        Args:
            directory: Directory to read input files from
        """
        raise NotImplementedError(f"from_directory has not been implemented in {cls}")


class LammpsAqueousSet(InputSetGenerator):
    """
    Yields a LammpsInputSet tailored for simulating aqueous electrolytes
    """

    def get_input_set(self) -> InputSet:
        """
        Generate a LammpsInputSet tailored for simulating aqueous electrolytes

        Typically the first argument to this method
        will be a Structure or other form of atomic coordinates.
        """
        pass
