# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Input sets for LAMMPS
"""

from __future__ import annotations

import logging
import os
from pathlib import Path
from string import Template

from monty.io import zopen

from pymatgen.core import Structure
from pymatgen.io.core import InputGenerator, InputSet
from pymatgen.io.lammps.data import CombinedData, LammpsData
from pymatgen.io.lammps.inputs import LammpsInputFile

__author__ = "Ryan Kingsbury, Guillaume Brunin (Matgenix)"
__copyright__ = "Copyright 2021, The Materials Project"
__version__ = "0.2"

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
        inputfile: LammpsInputFile | str,  # pylint: disable=E1131
        data: LammpsData | CombinedData,  # pylint: disable=E1131
        calc_type: str = None,
        template_file: str = None,
    ):
        """
        Args:
            inputfile: The input file containing settings.
                       It can be a LammpsInputFile object or a string representation.
            data: The data file containing structure and topology information.
                  It can be a LammpsData or a CombinedData object.
            calc_type: String used to shortly describe the type of computations performed by LAMMPS.
            template_file: Path (string) to the template file used to create the input file for LAMMPS.
        """
        self.inputfile = inputfile
        self.data = data
        self.calc_type = calc_type
        self.template_file = template_file

        super().__init__(inputs={"in.lammps": self.inputfile, "system.data": self.data})

    def get_inputs(self) -> dict[str, str | LammpsInputFile]:  # pylint: disable=E1131
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
        Construct a LammpsInputSet from a directory of two or more files.
        TODO: accept directories with only the input file, that should include the structure as well.

        Args:
            directory: Directory to read input files from. It should contain at least two files:
                       in.lammps for the LAMMPS input file, and system.data with the system information.
        """
        input_file = LammpsInputFile.from_file(os.path.join(directory, "in.lammps"))
        atom_style = input_file.get_args("atom_style")
        data_file = LammpsData.from_file(os.path.join(directory, "system.data"), atom_style=atom_style)
        return LammpsInputSet(inputfile=input_file, data=data_file, calc_type="read_from_dir")

    def validate(self) -> bool:
        """
        A place to implement basic checks to verify the validity of an
        input set. Can be as simple or as complex as desired.

        Will raise a NotImplementedError unless overloaded by the inheriting class.
        """
        raise NotImplementedError(f".validate() has not been implemented in {self.__class__}")


class BaseLammpsGenerator(InputGenerator):
    """
    Base class to generate LAMMPS input sets.
    Uses template files for the input. The variables that can be changed
    in the input template file are those starting with a $ sign, e.g., $nsteps.
    This generic class is specialized for each template in subclasses, e.g. LammpsMinimization.
    You can create a template for your own task following those present in pymatgen/io/lammps/templates.
    The parameters are then replaced based on the values found
    in the settings dictionary that you provide, e.g., {"nsteps": 1000}.

    Parameters:
        template: Path (string) to the template file used to create the input file for LAMMPS.
        calc_type: String used to shortly describe the type of computations performed by LAMMPS.

    Args:
        settings: Dictionary containing the values of the parameters to replace in the template.
    """

    template: str = os.path.join(template_dir, "md.template")
    calc_type: str = "lammps"

    def __init__(self, settings: dict | None = None):  # pylint: disable=E1131
        self.settings = settings or {}

        # load the template
        with zopen(self.template, "r") as f:
            template_str = f.read()

        # replace all variables
        self.input_str = Template(template_str).safe_substitute(**self.settings)
        self.input_file = LammpsInputFile.from_string(self.input_str)

    def get_input_set(  # type: ignore
        self, structure: Structure | LammpsData | CombinedData | None  # pylint: disable=E1131
    ) -> LammpsInputSet:
        """
        Generate a LammpsInputSet from the structure/data file, tailored to the template file.
        """
        if isinstance(structure, Structure):
            data = LammpsData.from_structure(structure)
        else:
            data = structure

        # Get the LammpsInputSet from the input file and data
        input_set = LammpsInputSet(
            inputfile=self.input_file, data=data, calc_type=self.calc_type, template_file=self.template
        )

        return input_set


class LammpsMinimization(BaseLammpsGenerator):
    """
    Yields a LammpsInputSet tailored for minimizing the energy of a system by iteratively
    adjusting atom coordinates.
    """

    template = os.path.join(template_dir, "minimization.template")
    calc_type = "minimization"

    def __init__(
        self,
        units: str = "metal",
        atom_style: str = "full",
        dimension: int = 3,
        boundary: str = "p p p",
        read_data: str = "system.data",
    ):
        self.units = units
        self.atom_style = atom_style
        self.dimension = dimension
        self.boundary = boundary
        self.read_data = read_data
        self.settings = {
            "units": units,
            "atom_style": atom_style,
            "dimension": dimension,
            "boundary": boundary,
            "read_data": read_data,
        }
        super().__init__(settings=self.settings)


class LammpsAqueousSet(InputGenerator):
    """
    Yields a LammpsInputSet tailored for simulating aqueous electrolytes
    """

    def get_input_set(self, mols: list, numbers: list[int]) -> LammpsInputSet:  # type: ignore
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
