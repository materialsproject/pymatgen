"""
Input sets for LAMMPS.
This InputSet and InputGenerator implementation is based on templates and is not intended to be very flexible.
For instance, pymatgen will not detect whether a given variable should be adapted based on others
(e.g., the number of steps from the temperature), it will not check for convergence nor will it actually run LAMMPS.
For additional flexibility and automation, use the atomate2-lammps implementation
(https://github.com/Matgenix/atomate2-lammps).
"""

from __future__ import annotations

import os
from typing import TYPE_CHECKING

from pymatgen.io.core import InputSet
from pymatgen.io.lammps.data import CombinedData, LammpsData
from pymatgen.io.lammps.inputs import LammpsInputFile

if TYPE_CHECKING:
    from typing_extensions import Self

    from pymatgen.util.typing import PathLike

__author__ = "Ryan Kingsbury, Guillaume Brunin (Matgenix)"
__copyright__ = "Copyright 2021, The Materials Project"
__version__ = "0.2"

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))


class LammpsInputSet(InputSet):
    r"""
    Container class for all LAMMPS inputs. This class is intended to provide
    general functionality that can be customized to many purposes.
    InputGenerator classes elsewhere in this module are used to create
    specific instances of LammpsInputSet that are tailored to specific purposes.

    /!\ This InputSet and InputGenerator implementation is based on templates and is not intended to be very flexible.
    For instance, pymatgen will not detect whether a given variable should be adapted based on others
    (e.g., the number of steps from the temperature), it will not check for convergence nor will it actually run LAMMPS.
    For additional flexibility and automation, use the atomate2-lammps implementation
    (https://github.com/Matgenix/atomate2-lammps).
    """

    def __init__(
        self,
        inputfile: LammpsInputFile | str,
        data: LammpsData | CombinedData,
        calc_type: str = "",
        template_file: PathLike = "",
        keep_stages: bool = False,
    ) -> None:
        """
        Args:
            inputfile: The input file containing settings. It can be a LammpsInputFile object
                or a string representation.
            data: The data file containing structure and topology information.
                It can be a LammpsData or a CombinedData object.
            calc_type: Human-readable string used to briefly describe the type of computations performed by LAMMPS.
            template_file: Path (string) to the template file used to create the input file for LAMMPS.
            keep_stages: Whether to keep the stage structure of the LammpsInputFile or not.
        """
        if isinstance(inputfile, LammpsInputFile):
            self.inputfile = inputfile
        else:
            self.inputfile = LammpsInputFile.from_str(inputfile, keep_stages=keep_stages)
        self.data = data
        self.calc_type = calc_type
        self.template_file = template_file
        self.keep_stages = keep_stages

        super().__init__(inputs={"in.lammps": self.inputfile, "system.data": self.data})

    @classmethod
    def from_directory(cls, directory: PathLike, keep_stages: bool = False) -> Self:
        """Construct a LammpsInputSet from a directory of two or more files.

        TODO: accept directories with only the input file, that should include the structure as well.

        Args:
            directory: Directory to read input files from. It should contain at least two files:
                       in.lammps for the LAMMPS input file, and system.data with the system information.
            keep_stages: Whether to keep the stage structure of the LammpsInputFile or not.
        """
        input_file = LammpsInputFile.from_file(f"{directory}/in.lammps", keep_stages=keep_stages)
        atom_style = input_file.get_args("atom_style")
        if isinstance(atom_style, list):
            raise ValueError("Variable atom_style is specified multiple times in the input file.")  # noqa: TRY004
        data_file = LammpsData.from_file(f"{directory}/system.data", atom_style=atom_style)
        return cls(inputfile=input_file, data=data_file, calc_type="read_from_dir")

    def validate(self) -> bool:
        """
        A place to implement basic checks to verify the validity of an
        input set. Can be as simple or as complex as desired.
        Will raise a NotImplementedError unless overloaded by the inheriting class.
        """
        raise NotImplementedError(f".validate() has not been implemented in {type(self).__name__}")
