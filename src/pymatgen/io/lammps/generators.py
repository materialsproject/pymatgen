"""
Input set generators for LAMMPS.
This InputSet and InputGenerator implementation is based on templates and is not intended to be very flexible.
For instance, pymatgen will not detect whether a given variable should be adapted based on others
(e.g., the number of steps from the temperature), it will not check for convergence nor will it actually run LAMMPS.
For additional flexibility and automation, use the atomate2-lammps implementation
(https://github.com/Matgenix/atomate2-lammps).
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from string import Template

from monty.io import zopen

from pymatgen.core import Structure
from pymatgen.io.core import InputGenerator
from pymatgen.io.lammps.data import CombinedData, LammpsData
from pymatgen.io.lammps.inputs import LammpsInputFile
from pymatgen.io.lammps.sets import LammpsInputSet

__author__ = "Ryan Kingsbury, Guillaume Brunin (Matgenix)"
__copyright__ = "Copyright 2021, The Materials Project"
__version__ = "0.2"

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
TEMPLATE_DIR = f"{MODULE_DIR}/templates"


@dataclass
class BaseLammpsGenerator(InputGenerator):
    r"""
    Base class to generate LAMMPS input sets.
    Uses template files for the input. The variables that can be changed
    in the input template file are those starting with a $ sign, e.g. $nsteps.
    This generic class is specialized for each template in subclasses, e.g. LammpsMinimization.
    You can create a template for your own task following those present in pymatgen/io/lammps/templates.
    The parameters are then replaced based on the values found
    in the settings dictionary that you provide, e.g. `{"nsteps": 1000}`.

    Attributes:
        template: Path (string) to the template file used to create the InputFile for LAMMPS.
        calc_type: Human-readable string used to briefly describe the type of computations performed by LAMMPS.
        settings: Dictionary containing the values of the parameters to replace in the template.
        keep_stages: If True, the string is formatted in a block structure with stage names
        and newlines that differentiate commands in the respective stages of the InputFile.
        If False, stage names are not printed and all commands appear in a single block.

    /!\ This InputSet and InputGenerator implementation is based on templates and is not intended to be very flexible.
    For instance, pymatgen will not detect whether a given variable should be adapted based on others
    (e.g., the number of steps from the temperature), it will not check for convergence nor will it actually run LAMMPS.
    For additional flexibility and automation, use the atomate2-lammps implementation
    (https://github.com/Matgenix/atomate2-lammps).
    """

    inputfile: LammpsInputFile | None = field(default=None)
    template: str = field(default_factory=str)
    data: LammpsData | CombinedData | None = field(default=None)
    settings: dict = field(default_factory=dict)
    calc_type: str = field(default="lammps")
    keep_stages: bool = field(default=True)

    def get_input_set(self, structure: Structure | LammpsData | CombinedData) -> LammpsInputSet:
        """Generate a LammpsInputSet from the structure/data, tailored to the template file."""
        data: LammpsData = LammpsData.from_structure(structure) if isinstance(structure, Structure) else structure

        # Load the template
        with zopen(self.template, mode="rt", encoding="utf-8") as file:
            template_str = file.read()

        # Replace all variables
        input_str = Template(template_str).safe_substitute(**self.settings)
        # Get LammpsInputFile
        input_file = LammpsInputFile.from_str(input_str, keep_stages=self.keep_stages)

        # Get the LammpsInputSet from the InputFile and data
        return LammpsInputSet(
            inputfile=input_file,
            data=data,
            calc_type=self.calc_type,
            template_file=self.template,
        )


class LammpsMinimization(BaseLammpsGenerator):
    """
    Generator that yields a LammpsInputSet tailored for minimizing the energy of a system by iteratively
    adjusting atom coordinates.
    Example usage:
    ```
    structure = Structure.from_file("mp-149.cif")
    lmp_minimization = LammpsMinimization(units="atomic").get_input_set(structure)
    ```.

    Do not forget to specify the force field, otherwise LAMMPS will not be able to run!

    This InputSet and InputGenerator implementation is based on templates and is not intended to be very flexible.
    For instance, pymatgen will not detect whether a given variable should be adapted based on others
    (e.g., the number of steps from the temperature), it will not check for convergence nor will it actually run LAMMPS.
    For additional flexibility and automation, use the atomate2-lammps implementation
    (https://github.com/Matgenix/atomate2-lammps).
    """

    def __init__(
        self,
        template: str | None = None,
        units: str = "metal",
        atom_style: str = "full",
        dimension: int = 3,
        boundary: str = "p p p",
        read_data: str = "system.data",
        force_field: str = "Unspecified force field!",
        keep_stages: bool = False,
    ) -> None:
        r"""
        Args:
            template: Path (string) to the template file used to create the InputFile for LAMMPS.
            units: units to be used for the LAMMPS calculation (see LAMMPS docs).
            atom_style: atom_style to be used for the LAMMPS calculation (see LAMMPS docs).
            dimension: dimension to be used for the LAMMPS calculation (see LAMMPS docs).
            boundary: boundary to be used for the LAMMPS calculation (see LAMMPS docs).
            read_data: read_data to be used for the LAMMPS calculation (see LAMMPS docs).
            force_field: force field to be used for the LAMMPS calculation (see LAMMPS docs).
                Note that you should provide all the required information as a single string.
                In case of multiple lines expected in the input file,
                separate them with '\n' in force_field.
            keep_stages: If True, the string is formatted in a block structure with stage names
                and newlines that differentiate commands in the respective stages of the InputFile.
                If False, stage names are not printed and all commands appear in a single block.
        """
        if template is None:
            template = f"{TEMPLATE_DIR}/minimization.template"
        settings = {
            "units": units,
            "atom_style": atom_style,
            "dimension": dimension,
            "boundary": boundary,
            "read_data": read_data,
            "force_field": force_field,
        }

        super().__init__(
            template=template,
            settings=settings,
            calc_type="minimization",
            keep_stages=keep_stages,
        )

    @property
    def units(self) -> str:
        """The argument of the command 'units' passed to the generator."""
        return self.settings["units"]

    @property
    def atom_style(self) -> str:
        """The argument of the command 'atom_style' passed to the generator."""
        return self.settings["atom_style"]

    @property
    def dimension(self) -> int:
        """The argument of the command 'dimension' passed to the generator."""
        return self.settings["dimension"]

    @property
    def boundary(self) -> str:
        """The argument of the command 'boundary' passed to the generator."""
        return self.settings["boundary"]

    @property
    def read_data(self) -> str:
        """The argument of the command 'read_data' passed to the generator."""
        return self.settings["read_data"]

    @property
    def force_field(self) -> str:
        """The details of the force field commands passed to the generator."""
        return self.settings["force_field"]
