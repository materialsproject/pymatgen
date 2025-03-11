"""Classes for reading/manipulating/writing FHI-aims input files.

Works for aims cube objects, geometry.in and control.in
"""

from __future__ import annotations

import os
import re
import textwrap
import time
from copy import deepcopy
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING
from warnings import warn

import numpy as np
from monty.io import zopen
from monty.json import MontyDecoder, MSONable
from monty.os.path import zpath
from pyfhiaims.control.control import AimsControl
from pyfhiaims.geometry.geometry import AimsGeometry
from pyfhiaims.species_defaults.species import SpeciesDefaults as PyFHIAimsSD

from pymatgen.core import SETTINGS, Element, Molecule, Structure

if TYPE_CHECKING:
    from typing import Any, Self


__author__ = "Thomas A. R. Purcell"
__version__ = "1.0"
__email__ = "purcellt@arizona.edu"
__date__ = "November 2023"


@dataclass
class AimsGeometryIn(MSONable):
    """Representation of an aims geometry.in file.

    Attributes:
        _content (str): The content of the input file
        _structure (Structure | Molecule): The structure or molecule
            representation of the file
    """

    _content: str
    _structure: Structure | Molecule

    @classmethod
    def from_str(cls, contents: str) -> Self:
        """Create an input from the content of an input file.

        Args:
            contents (str): The content of the file

        Returns:
            The AimsGeometryIn file for the string contents
        """
        content_lines = [
            line.strip() for line in contents.split("\n") if len(line.strip()) > 0 and line.strip()[0] != "#"
        ]
        geometry = AimsGeometry.from_strings(content_lines)

        return cls(_content="\n".join(content_lines), _structure=geometry.structure)

    @classmethod
    def from_file(cls, filepath: str | Path) -> Self:
        """Create an AimsGeometryIn from an input file.

        Args:
            filepath (str | Path): The path to the input file (either plain text of gzipped)

        Returns:
            AimsGeometryIn: The input object represented in the file
        """
        with zopen(filepath, mode="rt", encoding="utf-8") as in_file:
            content = in_file.read()
        return cls.from_str(content)

    @classmethod
    def from_structure(cls, structure: Structure | Molecule) -> Self:
        """Construct an input file from an input structure.

        Args:
            structure (Structure | Molecule): The structure for the file

        Returns:
            AimsGeometryIn: The input object for the structure
        """
        geometry = AimsGeometry.from_structure(structure)

        content = textwrap.dedent(
            f"""\
            #{"=" * 79}
            # FHI-aims geometry file: geometry.in
            # File generated from pymatgen
            # {time.asctime()}
            #{"=" * 79}
            """
        )
        content += geometry.to_string()
        return cls(_content=content, _structure=structure)

    @property
    def structure(self) -> Structure | Molecule:
        """Access structure for the file."""
        return self._structure

    @property
    def content(self) -> str:
        """Access the contents of the file."""
        return self._content

    def write_file(self, directory: str | Path | None = None, overwrite: bool = False) -> None:
        """Write the geometry.in file.

        Args:
            directory (str | Path | None): The directory to write the geometry.in file
            overwrite (bool): If True allow to overwrite existing files
        """
        directory = directory or Path.cwd()
        file_name = Path(directory) / "geometry.in"

        if not overwrite and file_name.exists():
            raise ValueError(f"geometry.in file exists in {directory}")

        with open(f"{directory}/geometry.in", mode="w", encoding="utf-8") as file:
            file.write(self.get_header(file_name.as_posix()))
            file.write(self.content)
            file.write("\n")

    def as_dict(self) -> dict[str, Any]:
        """Get a dictionary representation of the geometry.in file."""
        dct = {}
        dct["@module"] = type(self).__module__
        dct["@class"] = type(self).__name__
        dct["content"] = self.content
        dct["structure"] = self.structure
        return dct

    @classmethod
    def from_dict(cls, dct: dict[str, Any]) -> Self:
        """Initialize from dictionary.

        Args:
            dct (dict[str, Any]): The MontyEncoded dictionary of the AimsGeometryIn object

        Returns:
            The input object represented by the dict
        """
        decoded = {key: MontyDecoder().process_decoded(val) for key, val in dct.items() if not key.startswith("@")}

        return cls(_content=decoded["content"], _structure=decoded["structure"])


@dataclass
class AimsControlIn(MSONable):
    """An FHI-aims control.in file.

    Attributes:
        _parameters (dict[str, Any]): The parameters dictionary containing all input
            flags (key) and values for the control.in file
    """

    _parameters: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Initialize the output list of _parameters."""
        self._parameters.setdefault("output", [])

    def __getitem__(self, key: str) -> Any:
        """Get an input parameter.

        Args:
            key (str): The parameter to get

        Returns:
            The setting for that parameter

        Raises:
            KeyError: If the key is not in self._parameters
        """
        if key not in self._parameters:
            raise KeyError(f"{key} not set in AimsControl")
        return self._parameters[key]

    def __setitem__(self, key: str, value: Any) -> None:
        """Set an attribute of the class.

        Args:
            key (str): The parameter to get
            value (Any): The value for that parameter
        """
        if key == "output":
            if value in self._parameters[key]:
                return

            if isinstance(value, str):
                value = [value]
            self._parameters[key] += value
        else:
            self._parameters[key] = value

    def __delitem__(self, key: str) -> Any:
        """Delete a parameter from the input object.

        Args:
        key (str): The key in the parameter to remove

        Returns:
            Either the value of the deleted parameter or None if key is
            not in self._parameters
        """
        return self._parameters.pop(key, None)

    @property
    def parameters(self) -> dict[str, Any]:
        """The dictionary of input parameters for control.in."""
        return self._parameters

    @parameters.setter
    def parameters(self, parameters: dict[str, Any]) -> None:
        """Reset a control.in inputs from a parameters dictionary.

        Args:
            parameters (dict[str, Any]): The new set of parameters to use
        """
        self._parameters = parameters
        self._parameters.setdefault("output", [])

    def get_content(
        self,
        structure: Structure | Molecule,
        verbose_header: bool = False,
        directory: str | Path | None = None,
    ) -> str:
        """Get the content of the file.

        Args:
            structure (Structure | Molecule): The structure to write the input
                file for
            verbose_header (bool): If True print the input option dictionary
            directory: str | Path | None = The directory for the calculation,

        Returns:
            str: The content of the file for a given structure
        """
        parameters = deepcopy(self._parameters)

        if directory is None:
            directory = ""

        if parameters["xc"] == "LDA":
            parameters["xc"] = "pw-lda"

        geometry = AimsGeometry.from_structure(structure)
        magmom = np.array([atom.initial_moment for atom in geometry.atoms])
        if (
            parameters.get("spin", "") == "collinear"
            and np.allclose(magmom, 0.0)
            and ("default_initial_moment" not in parameters)
        ):
            warn(
                "Removing spin from parameters since no spin information is in the structure",
                RuntimeWarning,
                stacklevel=1,
            )
            parameters.pop("spin")

        outputs = parameters.pop("output", [])
        control_in = AimsControl(parameters=parameters, outputs=outputs)

        species_defaults_map = self._parameters.get("species_dir", SETTINGS.get("AIMS_SPECIES_DIR", ""))
        if not species_defaults_map:
            raise KeyError("Species' defaults not specified in the parameters")

        species_dir = None
        if isinstance(species_defaults_map, str):
            species_defaults_map = Path(species_defaults_map)
            if species_defaults_map.is_absolute():
                species_dir = species_defaults_map.parent
                basis_set = species_defaults_map.name
            else:
                basis_set = str(species_defaults_map)
        else:
            basis_set = species_defaults_map

        if species_dir is None:
            species_dir = SETTINGS.get("AIMS_SPECIES_DIR", "")

        elements = {site.species_string: site.specie.name for site in structure[::-1]}
        elements_map = {}
        for label in elements:
            clean_label = re.sub(r",\s*spin\s*=\s*[+-]?([0-9]*[.])?[0-9]+", "", label)
            elements_map[clean_label] = elements.get(clean_label, clean_label)

        for label, el in elements_map.items():
            if isinstance(basis_set, dict):
                basis_name = basis_set.get(label, None)
                if basis_name is None:
                    raise ValueError(f"Basis set not found for specie {label} (represented by element {el})")
            else:
                basis_name = basis_set
            if el != "Emptium":
                if not hasattr(Element, el):
                    raise ValueError(f"{el} is not a valid element name.")
                el_obj = Element(el)
                species_file_name = f"{el_obj.Z:02}_{el}_default"
            else:
                species_file_name = "00_Emptium_default"

            paths_to_try = [
                (Path(species_dir) / basis_name / species_file_name).expanduser().as_posix(),
                (Path(species_dir) / "defaults_2020" / basis_name / species_file_name).expanduser().as_posix(),
            ]
            for path in paths_to_try:
                path = zpath(path)
                if os.path.isfile(path):
                    with zopen(path, mode="rt") as file:
                        sdf_content = [line.strip() for line in file.readlines()]
                    geometry.set_species(label, PyFHIAimsSD.from_strings(sdf_content))

        content = "\n".join(
            [
                f"#{'=' * 79}",
                f"# FHI-aims geometry file: {directory or '.'}/geometry.in",
                "# File generated from pymatgen",
                f"# {time.asctime()}",
                f"#{'=' * 79}\n",
            ]
        )
        content += control_in.get_content(geometry=geometry, verbose_header=verbose_header)
        return content

    def write_file(
        self,
        structure: Structure | Molecule,
        directory: str | Path | None = None,
        verbose_header: bool = False,
        overwrite: bool = False,
    ) -> None:
        """Write the control.in file.

        Args:
            structure (Structure | Molecule): The structure to write the input
                file for
            directory (str or Path): The directory to write the control.in file.
                If None use cwd
            verbose_header (bool): If True print the input option dictionary
            overwrite (bool): If True allow to overwrite existing files

        Raises:
            ValueError: If a file must be overwritten and overwrite is False
            ValueError: If k-grid is not provided for the periodic structures
        """
        directory = directory or Path.cwd()

        if (Path(directory) / "control.in").exists() and not overwrite:
            raise ValueError(f"control.in file already in {directory}")

        if isinstance(structure, Structure) and (
            "k_grid" not in self._parameters and "k_grid_density" not in self._parameters
        ):
            raise ValueError("k-grid must be defined for periodic systems")

        content = self.get_content(structure, verbose_header)

        with open(f"{directory}/control.in", mode="w") as file:
            file.write(content)

    def as_dict(self) -> dict[str, Any]:
        """Get a dictionary representation of the geometry.in file."""
        dct: dict[str, Any] = {}
        dct["@module"] = type(self).__module__
        dct["@class"] = type(self).__name__
        dct["parameters"] = self.parameters
        return dct

    @classmethod
    def from_dict(cls, dct: dict[str, Any]) -> Self:
        """Initialize from dictionary.

        Args:
            dct (dict[str, Any]): The MontyEncoded dictionary

        Returns:
            The AimsControl for dct

        """
        decoded = {key: MontyDecoder().process_decoded(val) for key, val in dct.items() if not key.startswith("@")}

        return cls(_parameters=decoded["parameters"])
