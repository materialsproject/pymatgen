"""
This module defines a simple concrete implementation of the InputGenerator class that can be
used to facilitate writing large numbers of input files based on a template.
"""

from __future__ import annotations

from string import Template
from typing import TYPE_CHECKING

from monty.io import zopen

from pymatgen.io.core import InputGenerator, InputSet

if TYPE_CHECKING:
    from pymatgen.util.typing import PathLike

__author__ = "Ryan Kingsbury"
__email__ = "RKingsbury@lbl.gov"
__status__ = "Development"
__date__ = "October 2021"


class TemplateInputGen(InputGenerator):
    """
    Concrete implementation of InputGenerator that is based on a single template input
    file with variables.

    This class is provided as a low-barrier way to support new codes and to provide
    an intuitive way for users to transition from manual scripts to pymatgen I/O
    classes.
    """

    def get_input_set(
        self,
        template: PathLike,
        variables: dict | None = None,
        filename: PathLike = "input.txt",
    ) -> InputSet:
        """
        Args:
            template: the input file template containing variable strings to be
                replaced.
            variables: dict of variables to replace in the template. Keys are the
                text to replaced with the values, e.g. {"TEMPERATURE": 298} will
                replace the text $TEMPERATURE in the template. See Python's
                Template.safe_substitute() method documentation for more details.
            filename: the file to be written.
        """
        self.template = template
        self.variables = variables or {}
        self.filename = str(filename)

        # Load the template
        with zopen(self.template, mode="rt", encoding="utf-8") as file:
            template_str = file.read()

        # Replace all variables
        self.data = Template(template_str).safe_substitute(**self.variables)  # type:ignore[arg-type]
        return InputSet({self.filename: self.data})
