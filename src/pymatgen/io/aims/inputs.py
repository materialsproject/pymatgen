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
from monty.dev import deprecated
from monty.io import zopen
from monty.json import MontyDecoder, MSONable
from monty.os.path import zpath

from pymatgen.core import SETTINGS, Element, Lattice, Molecule, Species, Structure

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any

    from typing_extensions import Self

    from pymatgen.core.structure import IStructure


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

        species, coords, is_frac, lattice_vectors = [], [], [], []
        charges_dct, moments_dct, velocities_dct = {}, {}, {}

        for line in content_lines:
            inp = line.split()
            if (inp[0] == "atom") or (inp[0] == "atom_frac"):
                coords.append([float(ii) for ii in line.split()[1:4]])
                species.append(inp[4])
                is_frac.append(inp[0] == "atom_frac")
            if inp[0] == "lattice_vector":
                lattice_vectors.append([float(ii) for ii in line.split()[1:4]])
            if inp[0] == "initial_moment":
                moments_dct[len(coords) - 1] = float(inp[1])
            if inp[0] == "initial_charge":
                charges_dct[len(coords) - 1] = float(inp[1])
            if inp[0] == "velocity":
                velocities_dct[len(coords) - 1] = [float(x) for x in inp[1:]]

        charge = np.zeros(len(coords))
        for key, val in charges_dct.items():
            charge[key] = val

        magmom = np.zeros(len(coords))
        for key, val in moments_dct.items():
            magmom[key] = val

        velocity: list[None | list[float]] = [None for _ in coords]
        for key, v in velocities_dct.items():
            velocity[key] = v

        if len(lattice_vectors) == 3:
            lattice = Lattice(lattice_vectors)
            for cc in range(len(coords)):
                if is_frac[cc]:
                    coords[cc] = lattice.get_cartesian_coords(np.array(coords[cc]).reshape(1, 3)).flatten()
        elif len(lattice_vectors) == 0:
            lattice = None
            if any(is_frac):
                raise ValueError("Fractional coordinates given in file with no lattice vectors.")
        else:
            raise ValueError("Incorrect number of lattice vectors passed.")

        site_props = {"magmom": magmom, "charge": charge}
        if velocities_dct:
            site_props["velocity"] = velocity  # type:ignore[assignment]

        if lattice is None:
            structure = Molecule(species, coords, np.sum(charge), site_properties=site_props)  # type:ignore[arg-type]
        else:
            structure = Structure(  # type:ignore[assignment]
                lattice,
                species,
                coords,
                np.sum(charge),  # type:ignore[arg-type]
                coords_are_cartesian=True,
                site_properties=site_props,
            )

        return cls(_content="\n".join(content_lines), _structure=structure)

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
        return cls.from_str(content)  # type:ignore[arg-type]

    @classmethod
    def from_structure(cls, structure: Structure | IStructure | Molecule) -> Self:
        """Construct an input file from an input structure.

        Args:
            structure (Structure | Molecule): The structure for the file

        Returns:
            AimsGeometryIn: The input object for the structure
        """
        content_lines: list[str] = []

        if isinstance(structure, Structure):
            for lv in structure.lattice.matrix:
                content_lines.append(f"lattice_vector {lv[0]: .12e} {lv[1]: .12e} {lv[2]: .12e}")

        for site in structure:
            element = site.species_string.split(",spin=")[0]
            charge = site.properties.get("charge", 0)
            spin = site.properties.get("magmom", None)
            coord = site.coords
            v = site.properties.get("velocity", [0.0, 0.0, 0.0])

            if isinstance(site.specie, Species) and site.specie.spin is not None:
                if spin is not None and spin != site.specie.spin:
                    raise ValueError("species.spin and magnetic moments don't agree. Please only define one")
                spin = site.specie.spin

            if isinstance(site.specie, Species) and site.specie.oxi_state is not None:
                if charge is not None and charge != site.specie.oxi_state:
                    raise ValueError("species.oxi_state and charge don't agree. Please only define one")
                charge = site.specie.oxi_state

            content_lines.append(f"atom {coord[0]: .12e} {coord[1]: .12e} {coord[2]: .12e} {element}")
            if charge != 0:
                content_lines.append(f"     initial_charge {charge:.12e}")
            if (spin is not None) and (spin != 0):
                content_lines.append(f"     initial_moment {spin:.12e}")
            if (v is not None) and any(v_i != 0.0 for v_i in v):
                content_lines.append(f"     velocity   {'  '.join([f'{v_i:.12e}' for v_i in v])}")

        return cls(_content="\n".join(content_lines), _structure=structure)  # type:ignore[arg-type]

    @property
    def structure(self) -> Structure | Molecule:
        """Access structure for the file."""
        return self._structure

    @property
    def content(self) -> str:
        """Access the contents of the file."""
        return self._content

    def get_header(self, filename: str) -> str:
        """A header of geometry.in file

        Args:
            filename (str): A name of the file for the header
        """
        return textwrap.dedent(
            f"""\
        #{"=" * 72}
        # FHI-aims geometry file: {filename}
        # File generated from pymatgen
        # {time.asctime()}
        #{"=" * 72}
        """
        )

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
        dct["structure"] = self.structure  # type:ignore[assignment]
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


ALLOWED_AIMS_CUBE_TYPES = (
    "delta_density",
    "spin_density",
    "stm",
    "total_density",
    "total_density_integrable",
    "long_range_potential",
    "hartree_potential",
    "xc_potential",
    "delta_v",
    "ion_dens",
    "dielec_func",
    "elf",
)

ALLOWED_AIMS_CUBE_TYPES_STATE = (
    "first_order_density",
    "eigenstate",
    "eigenstate_imag",
    "eigenstate_density",
)

ALLOWED_AIMS_CUBE_FORMATS = (
    "cube",
    "gOpenMol",
    "xsf",
)


@dataclass
class AimsCube(MSONable):
    """The FHI-aims cubes.

    Attributes:
        type (str): The value to be outputted as a cube file
        origin (Sequence[float] or tuple[float, float, float]): The origin of the cube
        edges (Sequence[Sequence[float]]): Specifies the edges of a cube: dx, dy, dz
            dx (float): The length of the step in the x direction
            dy (float): The length of the step in the y direction
            dx (float): The length of the step in the x direction
        points (Sequence[int] or tuple[int, int, int]): The number of points
            along each edge
        spin_state (int): The spin-channel to use either 1 or 2
        kpoint (int): The k-point to use (the index of the list printed from
            `output k_point_list`)
        filename (str): The filename to use
        format (str): The format to output the cube file in: cube, gOpenMol, or xsf
        elf_type (int): The type of electron localization function to use (
            see FHI-aims manual)
    """

    type: str = field(default_factory=str)
    origin: Sequence[float] | tuple[float, float, float] = field(default_factory=lambda: [0.0, 0.0, 0.0])
    edges: Sequence[Sequence[float]] = field(default_factory=lambda: 0.1 * np.eye(3))  # type:ignore[return-value, arg-type]
    points: Sequence[int] | tuple[int, int, int] = field(default_factory=lambda: [0, 0, 0])
    format: str = "cube"
    spin_state: int | None = None
    kpoint: int | None = None
    filename: str | None = None
    elf_type: int | None = None

    def __eq__(self, other: object) -> bool:
        """Check if two cubes are equal to each other."""
        if not isinstance(other, AimsCube):
            return NotImplemented

        if self.type != other.type:
            return False

        if not np.allclose(self.origin, other.origin):
            return False

        if not np.allclose(self.edges, other.edges):
            return False

        if not np.allclose(self.points, other.points):
            return False

        if self.format != other.format:
            return False

        if self.spin_state != other.spin_state:
            return False

        if self.kpoint != other.kpoint:
            return False

        if self.filename != other.filename:
            return False

        return self.elf_type == other.elf_type

    def __post_init__(self) -> None:
        """Check the inputted variables to make sure they are correct.

        Raises:
            ValueError: If any of the inputs is invalid
        """
        split_type = self.type.split()
        cube_type = split_type[0]
        if split_type[0] in ALLOWED_AIMS_CUBE_TYPES:
            if len(split_type) > 1:
                raise ValueError(f"{cube_type=} can not have a state associated with it")
        elif split_type[0] in ALLOWED_AIMS_CUBE_TYPES_STATE:
            if len(split_type) != 2:
                raise ValueError(f"{cube_type=} must have a state associated with it")
        else:
            raise ValueError("Cube type undefined")

        if self.format not in ALLOWED_AIMS_CUBE_FORMATS:
            raise ValueError(f"{self.format} is invalid. Cube files must have a format of {ALLOWED_AIMS_CUBE_FORMATS}")

        valid_spins = (1, 2, None)
        if self.spin_state not in valid_spins:
            raise ValueError(f"Spin state must be one of {valid_spins}")

        if len(self.origin) != 3:
            raise ValueError("The cube origin must have 3 components")

        if len(self.points) != 3:
            raise ValueError("The number of points per edge must have 3 components")

        if len(self.edges) != 3:
            raise ValueError("Only three cube edges can be passed")

        for edge in self.edges:
            if len(edge) != 3:
                raise ValueError("Each cube edge must have 3 components")

        if self.elf_type is not None and self.type != "elf":
            raise ValueError("elf_type is only used when the cube type is elf. Otherwise it must be None")

    @property
    def control_block(self) -> str:
        """The block of text for the control.in file of the Cube."""
        cb = f"output cube {self.type}\n"
        cb += f"    cube origin {self.origin[0]: .12e} {self.origin[1]: .12e} {self.origin[2]: .12e}\n"
        for idx in range(3):
            cb += f"    cube edge {self.points[idx]} "
            cb += f"{self.edges[idx][0]: .12e} "
            cb += f"{self.edges[idx][1]: .12e} "
            cb += f"{self.edges[idx][2]: .12e}\n"
        cb += f"    cube format {self.format}\n"
        if self.spin_state is not None:
            cb += f"    cube spinstate {self.spin_state}\n"
        if self.kpoint is not None:
            cb += f"    cube kpoint {self.kpoint}\n"
        if self.filename is not None:
            cb += f"    cube filename {self.filename}\n"
        if self.elf_type is not None:
            cb += f"    cube elf_type {self.elf_type}\n"

        return cb

    def as_dict(self) -> dict[str, Any]:
        """Get a dictionary representation of the geometry.in file."""
        dct: dict[str, Any] = {}
        dct["@module"] = type(self).__module__
        dct["@class"] = type(self).__name__
        dct["type"] = self.type
        dct["origin"] = self.origin
        dct["edges"] = self.edges
        dct["points"] = self.points
        dct["format"] = self.format
        dct["spin_state"] = self.spin_state
        dct["kpoint"] = self.kpoint
        dct["filename"] = self.filename
        dct["elf_type"] = self.elf_type
        return dct

    @classmethod
    def from_dict(cls, dct: dict[str, Any]) -> Self:
        """Initialize from dictionary.

        Args:
            dct (dict[str, Any]): The MontyEncoded dictionary

        Returns:
            AimsCube
        """
        attrs = (
            "type",
            "origin",
            "edges",
            "points",
            "format",
            "spin_state",
            "kpoint",
            "filename",
            "elf_type",
        )
        decoded = {key: MontyDecoder().process_decoded(dct[key]) for key in attrs}

        return cls(**decoded)


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

    def get_aims_control_parameter_str(self, key: str, value: Any, fmt: str) -> str:
        """Get the string needed to add a parameter to the control.in file.

        Args:
            key (str): The name of the input flag
            value (Any): The value to be set for the flag
            fmt (str): The format string to apply to the value

        Returns:
            str: The line to add to the control.in file
        """
        if value is None:
            return ""
        return f"{key:35s}{fmt % value}\n"

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

        lim = "#" + "=" * 79
        content = ""

        if parameters["xc"] == "LDA":
            parameters["xc"] = "pw-lda"

        spins = np.array([getattr(sp, "spin", 0) for sp in structure.species])
        magmom = structure.site_properties.get("magmom", spins)
        if (
            parameters.get("spin", "") == "collinear"
            and np.allclose(magmom, 0.0)
            and ("default_initial_moment" not in parameters)
        ):
            warn(
                "Removing spin from parameters since no spin information is in the structure",
                RuntimeWarning,
                stacklevel=2,
            )
            parameters.pop("spin")

        cubes = parameters.pop("cubes", None)

        if verbose_header:
            content += "# \n# List of parameters used to initialize the calculator:"
            for param, val in parameters.items():
                content += f"#     {param}:{val}\n"
        content += f"{lim}\n"

        if "smearing" in parameters and "occupation_type" in parameters:
            raise ValueError(f'both "smearing" and "occupation_type" in {parameters=}')

        for key, value in parameters.items():
            if key in ["species_dir", "plus_u"]:
                continue
            if key == "smearing":
                name = parameters["smearing"][0].lower()
                if name == "fermi-dirac":
                    name = "fermi"
                width = parameters["smearing"][1]
                if name == "methfessel-paxton":
                    order = parameters["smearing"][2]
                    order = f" {order:d}"
                else:
                    order = ""

                content += self.get_aims_control_parameter_str("occupation_type", (name, width, order), "%s %f%s")
            elif key == "output":
                for output_type in value:
                    content += self.get_aims_control_parameter_str(key, output_type, "%s")
            elif key == "vdw_correction_hirshfeld" and value:
                content += self.get_aims_control_parameter_str(key, "", "%s")
            elif key == "xc":
                if "libxc" in value:
                    content += self.get_aims_control_parameter_str("override_warning_libxc", ".true.", "%s")
                content += self.get_aims_control_parameter_str(key, value, "%s")
            elif isinstance(value, bool):
                content += self.get_aims_control_parameter_str(key, str(value).lower(), ".%s.")
            elif isinstance(value, tuple | list):
                content += self.get_aims_control_parameter_str(key, " ".join(map(str, value)), "%s")
            elif isinstance(value, str):
                content += self.get_aims_control_parameter_str(key, value, "%s")
            else:
                content += self.get_aims_control_parameter_str(key, value, "%r")

        if cubes:
            for cube in cubes:
                content += cube.control_block

        content += f"{lim}\n\n"
        species_defaults = self._parameters.get("species_dir", SETTINGS.get("AIMS_SPECIES_DIR", ""))
        if not species_defaults:
            raise KeyError("Species' defaults not specified in the parameters")

        species_dir = None
        if isinstance(species_defaults, str):
            species_defaults = Path(species_defaults)
            if species_defaults.is_absolute():
                species_dir = species_defaults.parent
                basis_set = species_defaults.name
            else:
                basis_set = str(species_defaults)
        else:
            basis_set = species_defaults
        content += self.get_species_block(structure, basis_set, species_dir=species_dir)

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

        with open(f"{directory}/control.in", mode="w", encoding="utf-8") as file:
            file.write(f"#{'=' * 72}\n")
            file.write(f"# FHI-aims geometry file: {directory}/geometry.in\n")
            file.write("# File generated from pymatgen\n")
            file.write(f"# {time.asctime()}\n")
            file.write(f"#{'=' * 72}\n")

            file.write(content)

    def get_species_block(
        self,
        structure: Structure | Molecule,
        basis_set: str | dict[str, str],
        species_dir: str | Path | None = None,
    ) -> str:
        """Get the basis set information for a structure

        Args:
            structure (Molecule or Structure): The structure to get the basis set information for
            basis_set (str | dict[str, str]):
                a name of a basis set (`light`, `tight`...) or a mapping from site labels to basis set names.
                The name of a basis set can either correspond to the subfolder in `defaults_2020` folder
                or be a full path from the `FHI-aims/species_defaults` directory.
            species_dir (str | Path | None): The base species directory

        Returns:
            The block to add to the control.in file for the species

        Raises:
            ValueError: If a file for the species is not found
        """
        species_defaults = SpeciesDefaults.from_structure(structure, basis_set, species_dir)
        return str(species_defaults)

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


@dataclass
class AimsSpeciesFile:
    """An FHI-aims single species' defaults file.

    Attributes:
        data (str): A string of the complete species defaults file
        label (str): A string representing the name of species
    """

    data: str = ""
    label: str | None = None

    def __post_init__(self) -> None:
        """Set default label"""
        if self.label is None:
            for line in self.data.splitlines():
                if "species" in line:
                    self.label = line.split()[1]

    @classmethod
    def from_file(cls, filename: str, label: str | None = None) -> Self:
        """Initialize from file.

        Args:
            filename (str): The filename of the species' defaults file
            label (str): A string representing the name of species

        Returns:
            AimsSpeciesFile
        """
        with zopen(filename, mode="rt", encoding="utf-8") as file:
            return cls(data=file.read(), label=label)  # type:ignore[arg-type]

    @classmethod
    def from_element_and_basis_name(
        cls,
        element: str,
        basis: str,
        *,
        species_dir: str | Path | None = None,
        label: str | None = None,
    ) -> Self:
        """Initialize from element and basis names.

        Args:
            element (str): the element name (not to confuse with the species)
            basis (str): the directory in which the species' defaults file is located relative to the
                root `species_defaults` (or `species_defaults/defaults_2020`) directory.`.
            label (str): A string representing the name of species. If not specified,
                then equal to element

        Returns:
            AimsSpeciesFile
        """
        # check if element is in the Periodic Table (+ Emptium)
        if element != "Emptium":
            if not hasattr(Element, element):
                raise ValueError(f"{element} is not a valid element name.")
            el_obj = Element(element)
            species_file_name = f"{el_obj.Z:02}_{element}_default"
        else:
            species_file_name = "00_Emptium_default"

        aims_species_dir = species_dir or SETTINGS.get("AIMS_SPECIES_DIR")

        if aims_species_dir is None:
            raise ValueError(
                "No AIMS_SPECIES_DIR variable found in the config file. "
                "Please set the variable in ~/.config/.pmgrc.yaml to the root of `species_defaults` "
                "folder in FHIaims/ directory."
            )
        paths_to_try = [
            (Path(aims_species_dir) / basis / species_file_name).expanduser().as_posix(),
            (Path(aims_species_dir) / "defaults_2020" / basis / species_file_name).expanduser().as_posix(),
        ]
        for path in paths_to_try:
            path = zpath(path)
            if os.path.isfile(path):
                return cls.from_file(path, label)

        raise RuntimeError(
            f"Can't find the species' defaults file for {element} in {basis} basis set. Paths tried: {paths_to_try}"
        )

    def __str__(self) -> str:
        """String representation of the species' defaults file"""
        return re.sub(
            r"^ *species +\w+",
            f"  species        {self.label}",
            self.data,
            flags=re.MULTILINE,
        )

    @property
    def element(self) -> str:
        if match := re.search(r"^ *species +(\w+)", self.data, flags=re.MULTILINE):
            return match[1]
        raise ValueError("Can't find element in species' defaults file")

    def as_dict(self) -> dict[str, Any]:
        """Dictionary representation of the species' defaults file."""
        return {
            "label": self.label,
            "data": self.data,
            "@module": type(self).__module__,
            "@class": type(self).__name__,
        }

    @classmethod
    def from_dict(cls, dct: dict[str, Any]) -> Self:
        """Deserialization of the AimsSpeciesFile object"""
        return cls(**dct)


class SpeciesDefaults(list, MSONable):
    """A list containing a set of species' defaults objects with
    methods to read and write them to files
    """

    def __init__(
        self,
        labels: Sequence[str],
        basis_set: str | dict[str, str],
        *,
        species_dir: str | Path | None = None,
        elements: dict[str, str] | None = None,
    ) -> None:
        """
        Args:
            labels (list[str]): a list of labels, for which to build species' defaults
            basis_set (str | dict[str, str]):
                a name of a basis set (`light`, `tight`...) or a mapping from site labels to basis set names.
                The name of a basis set can either correspond to the subfolder in `defaults_2020` folder
                or be a full path from the `FHI-aims/species_defaults` directory.
            species_dir (str | Path | None): The base species directory
            elements (dict[str, str] | None):
                a mapping from site labels to elements. If some label is not in this mapping,
                it coincides with an element.
        """
        super().__init__()
        self.labels = labels
        self.basis_set = basis_set
        self.species_dir = species_dir

        if elements is None:
            elements = {}

        self.elements = {}
        for label in self.labels:
            label = re.sub(r",\s*spin\s*=\s*[+-]?([0-9]*[.])?[0-9]+", "", label)
            self.elements[label] = elements.get(label, label)
        self._set_species()

    def _set_species(self) -> None:
        """Initialize species defaults from the instance data"""
        del self[:]

        for label, el in self.elements.items():
            if isinstance(self.basis_set, dict):
                basis_set = self.basis_set.get(label, None)
                if basis_set is None:
                    raise ValueError(f"Basis set not found for specie {label} (represented by element {el})")
            else:
                basis_set = self.basis_set
            self.append(
                AimsSpeciesFile.from_element_and_basis_name(el, basis_set, species_dir=self.species_dir, label=label)
            )

    def __str__(self):
        """String representation of the species' defaults"""
        return "".join([str(x) for x in self])

    @classmethod
    def from_structure(
        cls,
        struct: Structure | Molecule,
        basis_set: str | dict[str, str],
        species_dir: str | Path | None = None,
    ):
        """Initialize species defaults from a structure."""
        labels = []
        elements = {}
        for site in struct:
            el = site.specie
            if site.species_string not in labels:
                labels.append(site.species_string)
                elements[site.species_string] = el.name
        return SpeciesDefaults(sorted(labels), basis_set, species_dir=species_dir, elements=elements)

    def as_dict(self):
        """Dictionary representation of the species' defaults"""
        return {
            "labels": self.labels,
            "elements": self.elements,
            "basis_set": self.basis_set,
            "@module": type(self).__module__,
            "@class": type(self).__name__,
        }

    @deprecated(replacement=as_dict, deadline=(2026, 4, 4))
    def to_dict(self):
        return self.as_dict()

    @classmethod
    def from_dict(cls, dct: dict[str, Any]) -> SpeciesDefaults:
        """Deserialization of the SpeciesDefaults object"""
        return SpeciesDefaults(dct["labels"], dct["basis_set"], elements=dct["elements"])
