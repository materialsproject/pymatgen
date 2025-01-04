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
from pyfhiaims.control.control import AimsControlIn as PyFHIAimsControl
from pyfhiaims.geometry.atom import FHIAimsAtom
from pyfhiaims.geometry.geometry import AimsGeometry
from pyfhiaims.species_defaults.species import SpeciesDefaults as PyFHIAimsSD

from pymatgen.core import SETTINGS, Element, Lattice, Molecule, Species, Structure

if TYPE_CHECKING:
    from typing import Any

    from typing_extensions import Self


__author__ = "Thomas A. R. Purcell"
__version__ = "1.0"
__email__ = "purcellt@arizona.edu"
__date__ = "November 2023"


def structure2aimsgeo(structure: Structure | Molecule) -> AimsGeometry:
    """Convert a structure into an AimsGeometry object

    Args:
        structure (Structure | Molecule): structure to convert

    Returns:
        AimsGeometry: The resulting AimsGeometry
    """
    lattice_vectors = getattr(structure, "lattice", None)
    if lattice_vectors is not None:
        lattice_vectors = lattice_vectors.matrix

    atoms = []
    for site in structure:
        element = site.species_string.split(",spin=")[0]
        charge = site.properties.get("charge", None)
        spin = site.properties.get("magmom", None)
        coord = site.coords
        v = site.properties.get("velocity", None)

        if isinstance(site.specie, Species) and site.specie.spin is not None:
            if spin is not None and spin != site.specie.spin:
                raise ValueError("species.spin and magnetic moments don't agree. Please only define one")
            spin = site.specie.spin
        elif spin is None:
            spin = 0.0

        if isinstance(site.specie, Species) and site.specie.oxi_state is not None:
            if charge is not None and charge != site.specie.oxi_state:
                raise ValueError("species.oxi_state and charge don't agree. Please only define one")
            charge = site.specie.oxi_state
        elif charge is None:
            charge = 0.0

        atoms.append(
            FHIAimsAtom(
                symbol=element,
                position=coord,
                velocity=v,
                initial_charge=charge,
                initial_moment=spin,
            )
        )

    return AimsGeometry(atoms=atoms, lattice_vectors=lattice_vectors)


def aimsgeo2structure(geometry: AimsGeometry) -> Structure | Molecule:
    """Convert an AimsGeometry object into a Structure of Molecule

    Args:
        structure (Structure | Molecule): structure to convert

    Returns:
        Structure | Molecule: The resulting Strucucture/Molecule
    """
    charge = np.array([atom.initial_charge for atom in geometry.atoms])
    magmom = np.array([atom.initial_moment for atom in geometry.atoms])
    velocity: list[None | list[float]] = [atom.velocity for atom in geometry.atoms]
    species = [atom.symbol for atom in geometry.atoms]
    coords = [atom.position for atom in geometry.atoms]

    for vv, vel in enumerate(velocity):
        if vel is not None and np.sum(np.abs(vel)) < 1e-10:
            velocity[vv] = None

    lattice = Lattice(geometry.lattice_vectors) if geometry.lattice_vectors is not None else None

    site_props = {"charge": charge, "magmom": magmom}
    if any(vel is not None for vel in velocity):
        site_props["velocity"] = velocity

    if lattice is None:
        structure = Molecule(species, coords, np.sum(charge), site_properties=site_props)
    else:
        structure = Structure(
            lattice,
            species,
            coords,
            np.sum(charge),
            coords_are_cartesian=True,
            site_properties=site_props,
        )

    return structure


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
        structure = aimsgeo2structure(geometry)

        return cls(_content="\n".join(content_lines), _structure=structure)

    @classmethod
    def from_file(cls, filepath: str | Path) -> Self:
        """Create an AimsGeometryIn from an input file.

        Args:
            filepath (str | Path): The path to the input file (either plain text of gzipped)

        Returns:
            AimsGeometryIn: The input object represented in the file
        """
        with zopen(filepath, mode="rt") as in_file:
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
        geometry = structure2aimsgeo(structure)

        content = textwrap.dedent(
            f"""\
            #{'=' * 79}
            # FHI-aims geometry file: geometry.in
            # File generated from pymatgen
            # {time.asctime()}
            #{'=' * 79}
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

        with open(f"{directory}/geometry.in", mode="w") as file:
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
            raise KeyError(f"{key} not set in AimsControlIn")
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

        geometry = structure2aimsgeo(structure)
        magmom = np.array([atom.initial_moment for atom in geometry.atoms])
        if (
            parameters.get("spin", "") == "collinear"
            and np.all(magmom == 0.0)
            and ("default_initial_moment" not in parameters)
        ):
            warn(
                "Removing spin from parameters since no spin information is in the structure",
                RuntimeWarning,
                stacklevel=1,
            )
            parameters.pop("spin")

        outputs = parameters.pop("output", [])
        control_in = PyFHIAimsControl(parameters=parameters, outputs=outputs)

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

    # def get_species_block(
    #     self,
    #     structure: Structure | Molecule,
    #     basis_set: str | dict[str, str],
    #     species_dir: str | Path | None = None,
    # ) -> str:
    #     """Get the basis set information for a structure

    #     Args:
    #         structure (Molecule or Structure): The structure to get the basis set information for
    #         basis_set (str | dict[str, str]):
    #             a name of a basis set (`light`, `tight`...) or a mapping from site labels to basis set names.
    #             The name of a basis set can either correspond to the subfolder in `defaults_2020` folder
    #             or be a full path from the `FHI-aims/species_defaults` directory.
    #         species_dir (str | Path | None): The base species directory

    #     Returns:
    #         The block to add to the control.in file for the species

    #     Raises:
    #         ValueError: If a file for the species is not found
    #     """
    #     species_defaults = SpeciesDefaults.from_structure(structure, basis_set, species_dir)
    #     return str(species_defaults)

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
            The AimsControlIn for dct
        """
        decoded = {key: MontyDecoder().process_decoded(val) for key, val in dct.items() if not key.startswith("@")}

        return cls(_parameters=decoded["parameters"])


# @dataclass
# class AimsSpeciesFile:
#     """An FHI-aims single species' defaults file.

#     Attributes:
#         data (str): A string of the complete species defaults file
#         label (str): A string representing the name of species
#     """

#     data: str = ""
#     label: str | None = None

#     def __post_init__(self) -> None:
#         """Set default label"""
#         if self.label is None:
#             for line in self.data.splitlines():
#                 if "species" in line:
#                     self.label = line.split()[1]

#     @classmethod
#     def from_file(cls, filename: str, label: str | None = None) -> Self:
#         """Initialize from file.

#         Args:
#             filename (str): The filename of the species' defaults file
#             label (str): A string representing the name of species

#         Returns:
#             AimsSpeciesFile
#         """
#         with zopen(filename, mode="rt") as file:
#             return cls(data=file.read(), label=label)

#     @classmethod
#     def from_element_and_basis_name(
#         cls,
#         element: str,
#         basis: str,
#         *,
#         species_dir: str | Path | None = None,
#         label: str | None = None,
#     ) -> Self:
#         """Initialize from element and basis names.

#         Args:
#             element (str): the element name (not to confuse with the species)
#             basis (str): the directory in which the species' defaults file is located relative to the
#                 root `species_defaults` (or `species_defaults/defaults_2020`) directory.`.
#             label (str): A string representing the name of species. If not specified,
#                 then equal to element

#         Returns:
#             AimsSpeciesFile
#         """
#         # check if element is in the Periodic Table (+ Emptium)
#         if element != "Emptium":
#             if not hasattr(Element, element):
#                 raise ValueError(f"{element} is not a valid element name.")
#             el_obj = Element(element)
#             species_file_name = f"{el_obj.Z:02}_{element}_default"
#         else:
#             species_file_name = "00_Emptium_default"

#         aims_species_dir = species_dir or SETTINGS.get("AIMS_SPECIES_DIR")

#         if aims_species_dir is None:
#             raise ValueError(
#                 "No AIMS_SPECIES_DIR variable found in the config file. "
#                 "Please set the variable in ~/.config/.pmgrc.yaml to the root of `species_defaults` "
#                 "folder in FHIaims/ directory."
#             )
#         paths_to_try = [
#             (Path(aims_species_dir) / basis / species_file_name).expanduser().as_posix(),
#             (Path(aims_species_dir) / "defaults_2020" / basis / species_file_name).expanduser().as_posix(),
#         ]
#         for path in paths_to_try:
#             path = zpath(path)
#             if os.path.isfile(path):
#                 return cls.from_file(path, label)

#         raise RuntimeError(
#             f"Can't find the species' defaults file for {element} in {basis} basis set. Paths tried: {paths_to_try}"
#         )

#     def __str__(self) -> str:
#         """String representation of the species' defaults file"""
#         return re.sub(
#             r"^ *species +\w+",
#             f"  species        {self.label}",
#             self.data,
#             flags=re.MULTILINE,
#         )

#     @property
#     def element(self) -> str:
#         if match := re.search(r"^ *species +(\w+)", self.data, flags=re.MULTILINE):
#             return match[1]
#         raise ValueError("Can't find element in species' defaults file")

#     def as_dict(self) -> dict[str, Any]:
#         """Dictionary representation of the species' defaults file."""
#         return {
#             "label": self.label,
#             "data": self.data,
#             "@module": type(self).__module__,
#             "@class": type(self).__name__,
#         }

#     @classmethod
#     def from_dict(cls, dct: dict[str, Any]) -> Self:
#         """Deserialization of the AimsSpeciesFile object"""
#         return cls(**dct)


# class SpeciesDefaults(list, MSONable):
#     """A list containing a set of species' defaults objects with
#     methods to read and write them to files
#     """

#     def __init__(
#         self,
#         labels: Sequence[str],
#         basis_set: str | dict[str, str],
#         *,
#         species_dir: str | Path | None = None,
#         elements: dict[str, str] | None = None,
#     ) -> None:
#         """
#         Args:
#             labels (list[str]): a list of labels, for which to build species' defaults
#             basis_set (str | dict[str, str]):
#                 a name of a basis set (`light`, `tight`...) or a mapping from site labels to basis set names.
#                 The name of a basis set can either correspond to the subfolder in `defaults_2020` folder
#                 or be a full path from the `FHI-aims/species_defaults` directory.
#             species_dir (str | Path | None): The base species directory
#             elements (dict[str, str] | None):
#                 a mapping from site labels to elements. If some label is not in this mapping,
#                 it coincides with an element.
#         """
#         super().__init__()
#         self.labels = labels
#         self.basis_set = basis_set
#         self.species_dir = species_dir

#         if elements is None:
#             elements = {}

#         self.elements = {}
#         for label in self.labels:
#             label = re.sub(r",\s*spin\s*=\s*[+-]?([0-9]*[.])?[0-9]+", "", label)
#             self.elements[label] = elements.get(label, label)
#         self._set_species()

#     def _set_species(self) -> None:
#         """Initialize species defaults from the instance data"""
#         del self[:]

#         for label, el in self.elements.items():
#             if isinstance(self.basis_set, dict):
#                 basis_set = self.basis_set.get(label, None)
#                 if basis_set is None:
#                     raise ValueError(f"Basis set not found for specie {label} (represented by element {el})")
#             else:
#                 basis_set = self.basis_set
#             self.append(
#                 AimsSpeciesFile.from_element_and_basis_name(el, basis_set, species_dir=self.species_dir, label=label)
#             )

#     def __str__(self):
#         """String representation of the species' defaults"""
#         return "".join([str(x) for x in self])

#     @classmethod
#     def from_structure(
#         cls,
#         struct: Structure | Molecule,
#         basis_set: str | dict[str, str],
#         species_dir: str | Path | None = None,
#     ):
#         """Initialize species defaults from a structure."""
#         labels = []
#         elements = {}
#         for site in struct:
#             el = site.specie
#             if site.species_string not in labels:
#                 labels.append(site.species_string)
#                 elements[site.species_string] = el.name
#         return SpeciesDefaults(sorted(labels), basis_set, species_dir=species_dir, elements=elements)

#     def to_dict(self):
#         """Dictionary representation of the species' defaults"""
#         return {
#             "labels": self.labels,
#             "elements": self.elements,
#             "basis_set": self.basis_set,
#             "@module": type(self).__module__,
#             "@class": type(self).__name__,
#         }

#     @classmethod
#     def from_dict(cls, dct: dict[str, Any]) -> SpeciesDefaults:
#         """Deserialization of the SpeciesDefaults object"""
#         return SpeciesDefaults(dct["labels"], dct["basis_set"], elements=dct["elements"])
