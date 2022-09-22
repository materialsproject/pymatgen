# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Input sets for LAMMPS
"""
import logging
import os
from pathlib import Path
from typing import Dict, List

from pymatgen.core import Structure
from pymatgen.io.core import InputGenerator, InputSet
from pymatgen.io.lammps.data import CombinedData, LammpsData
from pymatgen.io.lammps.inputs import LammpsInputFile
from pymatgen.io.template import TemplateInputGen

__author__ = "Ryan Kingsbury, ..."
__copyright__ = "Copyright 2021, The Materials Project"
__version__ = "0.1"

logger = logging.getLogger(__name__)
template_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "templates")


class LammpsInputSet(InputSet):
    """
    Container class for all LAMMPS inputs. This class is intended to provide
    general functionality that can be customized to many purposes.

    InputGenerator classes elsewhere in this module are used to create
    specific instances of LammpsInputSet that are tailored to specific purposes.
    """

    def __init__(
        self,
        inputfile: LammpsInputFile,
        data: LammpsData | CombinedData,  # pylint: disable=E1131
        calc_type: str = None,
        template_file: str = None,
    ):
        """
        Args:
            inputfile: The input file containing settings
            data: the data file containing structure and topology information
        """
        self.inputfile = inputfile
        self.data = data
        self.calc_type = calc_type
        self.template_file = template_file

        super().__init__(inputs={"in.lammps": self.inputfile, "system.data": self.data})

    def get_inputs(self) -> Dict[str, str | LammpsInputFile]:  # pylint: disable=E1131
        """
        Generate a dictionary of one or more input files to be written. Keys
        are filenames, values are the contents of each file.

        This method is called by write_input(), which performs the actual file
        write operations.
        """
        return {"in.lammps": self.inputfile, "system.data": self.data}

    @classmethod
    def from_directory(cls, directory: str | Path):  # pylint: disable=E1131
        """
        Construct an InputSet from a directory of one or more files.

        Args:
            directory: Directory to read input files from
        """
        raise NotImplementedError(f"from_directory has not been implemented in {cls}")


class LammpsMinimization(InputGenerator):
    """
    Yields a LammpsInputSet tailored for minimizing the energy of a system by iteratively
    adjusting atom coordinates.
    """

    template = os.path.join(template_dir, "minimization.template")
    calc_type = "minimization"

    def __init__(
        self, structure: Structure | LammpsData | CombinedData | None, settings: dict | None = None
    ):  # pylint: disable=E1131
        """ """
        if isinstance(structure, Structure):
            self.data = LammpsData.from_structure(structure)
        else:
            self.data = structure

        self.settings = settings

        # Get the input file using Templates
        template_set = TemplateInputGen().get_input_set(
            template=self.template, variables=self.settings, filename="in.lammps"
        )
        self.input_file = LammpsInputFile.from_string(template_set.inputs["in.lammps"])

    def get_input_set(self) -> LammpsInputSet:
        """
        Generate a LammpsInputSet tailored for minimizing the energy of a system
        """

        # Get the LammpsInputSet from the input file and data
        input_set = LammpsInputSet(
            inputfile=self.input_file, data=self.data, calc_type=self.calc_type, template_file=self.template
        )

        return input_set


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
