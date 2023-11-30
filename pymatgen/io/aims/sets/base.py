"""Module defining base FHI-aims input set and generator."""
from __future__ import annotations

import contextlib
import copy
import json
import logging
import os
import shutil
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from monty.json import MontyEncoder

from pymatgen.io.aims.inputs import AimsControlIn, AimsGeometryIn
from pymatgen.io.core import InputFile, InputSet

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence

    from pymatgen.core import Molecule, Structure

TMPDIR_NAME: str = "tmpdir"
OUTPUT_FILE_NAME: str = "aims.out"
CONTROL_FILE_NAME: str = "control.in"
PARAMS_JSON_FILE_NAME: str = "parameters.json"
GEOMETRY_FILE_NAME: str = "geometry.in"


@contextlib.contextmanager
def cwd(path: str, mkdir: bool = False, rmdir: bool = False):
    """Change cwd intermediately.

    Example
    -------
    >>> with cwd(some_path):
    >>>     do so some stuff in some_path
    >>> do so some other stuff in old cwd

    Args:
        path (str) or Path: Path to change working directory to
        mkdir (bool): If True make path if it does not exist
        rmdir (bool): If True remove the working directory upon exiting
    """
    CWD = os.getcwd()

    if os.path.exists(path) is False and mkdir:
        os.makedirs(path)

    os.chdir(path)
    yield

    os.chdir(CWD)
    if rmdir:
        shutil.rmtree(path)


DEFAULT_AIMS_PROPERTIES = [
    "energy",
    "free_energy",
    "forces",
    "stress",
    "stresses",
    "dipole",
    "magmom",
]

logger = logging.getLogger(__name__)


@dataclass
class AimsInputFile(InputFile):
    """The input file for an FHI-aims calculation.

    Args:
        _content_str (str): The contents of the input file as a string
    """

    _content_str: str = ""

    def get_string(self) -> str:
        """Get the contents of the input file.

        Returns:
            str: The contents of the input file
        """
        return self._content_str

    def get_str(self) -> str:
        """Get the contents of the input file.

        Returns:
            str: The contents of the input file
        """
        return self._content_str

    @classmethod
    def from_string(cls, contents: str):
        """Create an input file from the contents string.

        Args:
            contents (str): The contents of the input file
        """
        return cls(contents)

    @classmethod
    def from_str(cls, contents: str):
        """Create an input file from the contents string.

        Args:
            contents (str): The contents of the input file
        """
        return cls(contents)


class AimsInputSet(InputSet):
    """A class to represent a set of Aims inputs."""

    def __init__(
        self,
        parameters: dict[str, Any],
        structure: Structure | Molecule,
        properties: Sequence[str] = ("energy", "free_energy"),
    ):
        """Construct the AimsInputSet.

        Args:
            parameters (Dict[str, Any]): The ASE parameters object for the calculation
            structure (Structure or Molecule): The Structure/Molecule objects to
                create the inputs for
            properties (Sequence[str]): The properties to calculate for the calculation
        """
        self._parameters = parameters
        self._structure = structure
        self._properties = properties

        aims_control_in, aims_geometry_in = self.get_input_files()
        super().__init__(
            inputs={
                CONTROL_FILE_NAME: aims_control_in,
                GEOMETRY_FILE_NAME: aims_geometry_in,
                PARAMS_JSON_FILE_NAME: json.dumps(self._parameters, indent=2, cls=MontyEncoder),
            }
        )

    def get_input_files(self) -> tuple[str, str]:
        """Get the input file contents for the calculation.

        Returns:
            tuple[str, str]: The contents of the control.in and geometry.in file
        """
        property_flags = {
            "forces": "compute_forces",
            "stress": "compute_analytical_stress",
            "stresses": "compute_heat_flux",
        }
        updated_params = dict(**self._parameters)
        for prop in self._properties:
            aims_name = property_flags.get(prop, None)
            if aims_name is not None:
                updated_params[aims_name] = True

        aims_geometry_in = AimsGeometryIn.from_structure(self._structure)
        aims_control_in = AimsControlIn(updated_params)

        with cwd(TMPDIR_NAME, mkdir=True, rmdir=True):
            aims_control_in.write_file(self._structure)
            aims_control_in_content = AimsInputFile.from_file("control.in")

            aims_geometry_in.write_file()
            aims_geometry_in_content = AimsInputFile.from_file("geometry.in")

        return aims_control_in_content, aims_geometry_in_content

    @property
    def control_in(self) -> str | InputFile | slice:
        """Get the control.in file contents."""
        return self[CONTROL_FILE_NAME]

    @property
    def geometry_in(self) -> str | InputFile | slice:
        """Get the geometry.in file contents."""
        return self[GEOMETRY_FILE_NAME]

    @property
    def parameters_json(self) -> str | InputFile | slice:
        """Get the JSON representation of the parameters dict."""
        return self[PARAMS_JSON_FILE_NAME]

    def set_parameters(self, *args, **kwargs) -> dict[str, Any]:
        """Set the parameters object for the AimsTemplate.

        This sets the parameters object that is passed to an AimsTempalte and
        resets the control.in file

        One can pass a dictionary mapping the aims variables to their values or
        the aims variables as keyword arguments. A combination of the two
        options is also allowed.

        Returns:
            dict[str, Any]: dictionary with the variables that have been added.
        """
        self._parameters.clear()
        for arg in args:
            self._parameters.update(arg)

        self._parameters.update(kwargs)

        aims_control_in, _ = self.get_input_files()

        self.inputs[CONTROL_FILE_NAME] = aims_control_in
        self.inputs[PARAMS_JSON_FILE_NAME] = json.dumps(self._parameters, indent=2, cls=MontyEncoder)

        inputs = {str(key): val for key, val in self.inputs.items()}
        self.__dict__.update(inputs)

        return self._parameters

    def remove_parameters(self, keys: Iterable[str] | str, strict: bool = True) -> dict[str, Any]:
        """Remove the aims parameters listed in keys.

        This removes the aims variables from the parameters object.

        Args:
            keys (Iterable[str] or str): string or list of strings with the names of
                the aims parameters to be removed.
            strict (bool): whether to raise a KeyError if one of the aims parameters
                to be removed is not present.

        Returns:
            dict[str, Any]: Dictionary with the variables that have been removed.
        """
        if isinstance(keys, str):
            keys = [keys]
        for key in keys:
            if strict and key not in self._parameters:
                raise ValueError(f"The key ({key}) is not in self._parameters")

            if key not in self._parameters:
                continue

            del self._parameters[key]

        return self.set_parameters(**self._parameters)

    def set_structure(self, structure: Structure | Molecule):
        """Set the structure object for this input set.

        Args:
            structure (Structure or Molecule): The new Structure or Molecule
                for the calculation
        """
        self._structure = structure

        aims_control_in, aims_geometry_in = self.get_input_files()
        self.inputs[GEOMETRY_FILE_NAME] = aims_geometry_in
        self.inputs[CONTROL_FILE_NAME] = aims_control_in
        inputs = {str(key): val for key, val in self.inputs.items()}
        self.__dict__.update(inputs)

    def deepcopy(self):
        """Deep copy of the input set."""
        return copy.deepcopy(self)
