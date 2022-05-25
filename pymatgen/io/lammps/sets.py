# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Input sets for LAMMPS
"""

import logging
from pathlib import Path
from typing import Dict, List, Union

from pymatgen.io.core import InputGenerator, InputSet
from pymatgen.io.lammps.inputs import CombinedData, LammpsInputFile

__author__ = "Ryan Kingsbury, ..."
__copyright__ = "Copyright 2021, The Materials Project"
__version__ = "0.1"

logger = logging.getLogger(__name__)


class LammpsInputSet(InputSet):
    """
    Container class for all LAMMPS inputs. This class is intended to provide
    general functionality that can be customized to many purposes.

    InputGenerator classes elsewhere in this module are used to create
    specific instances of LammpsInputSet that are tailored to specific purposes.
    """

    def __init__(self, inputfile: LammpsInputFile, data: CombinedData):
        """
        Args:
            inputfile: The input file containing settings
            data: the data file containing structure and topology information
        """
        self.inputfile = inputfile
        self.datafile = data

    def get_inputs(self) -> Dict[str, Union[str, LammpsInputFile]]:
        """
        Generate a dictionary of one or more input files to be written. Keys
        are filenames, values are the contents of each file.

        This method is called by write_input(), which performs the actual file
        write operations.
        """
        return {"in.lammps": self.inputfile, "system.data": self.datafile}

    @classmethod
    def from_directory(cls, directory: Union[str, Path]):
        """
        Construct an InputSet from a directory of one or more files.

        Args:
            directory: Directory to read input files from
        """
        raise NotImplementedError(f"from_directory has not been implemented in {cls}")


class LammpsAqueousSet(InputGenerator):
    """
    Yields a LammpsInputSet tailored for simulating aqueous electrolytes
    """

    def get_input_set(self, mols: List, numbers: List[int]) -> LammpsInputSet:  # type: ignore
        """
        Generate a LammpsInputSet tailored for simulating aqueous electrolytes

        Typically the first argument to this method
        will be a Structure or other form of atomic coordinates.

        Args:
            mols: list of Molecule objects
            numbers: list of int describing the number of each type of  molecule.
                Numbers must be passed in the order corresponding to the order
                of Molecule in mol.
        """
