"""Classes for reading/manipulating/writing FHI-aims input files.

Works for aims cube objects, geometry.in and control.in
"""

from __future__ import annotations

import gzip
import os
import time
from copy import deepcopy
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from monty.io import zopen
from monty.json import MontyDecoder, MSONable

from pymatgen.core import Lattice, Molecule, Structure

if TYPE_CHECKING:
    from collections.abc import Sequence

    from typing_extensions import Self

__author__ = "Thomas A. R. Purcell"
__version__ = "1.0"
__email__ = "purcellt@arizona.edu"
__date__ = "November 2023"


@dataclass
class AimsGeometryIn(MSONable):
    """Representation of an aims geometry.in file

    Attributes:
        _content (str): The content of the input file
        _structure (Structure | Molecule): The structure or molecule
            representation of the file
    """

    _content: str
    _structure: Structure | Molecule

    @classmethod
    def from_str(cls, contents: str) -> Self:
        """Create an input from the content of an input file

        Args:
            contents (str): The content of the file

        Returns:
            The AimsGeometryIn file for the string contents
        """
        content_lines = [
            line.strip() for line in contents.split("\n") if len(line.strip()) > 0 and line.strip()[0] != "#"
        ]

        species, coords, is_frac, lattice_vectors = [], [], [], []
        charges_dct, moments_dct = {}, {}

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

        charge = np.zeros(len(coords))
        for key, val in charges_dct.items():
            charge[key] = val

        magmom = np.zeros(len(coords))
        for key, val in moments_dct.items():
            magmom[key] = val

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
        if lattice is None:
            structure = Molecule(species, coords, np.sum(charge), site_properties=site_props)
        else:
            structure = Structure(
                lattice, species, coords, np.sum(charge), coords_are_cartesian=True, site_properties=site_props
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
        content_lines: list[str] = []

        if isinstance(structure, Structure):
            for lv in structure.lattice.matrix:
                content_lines.append(f"lattice_vector {lv[0]: .12e} {lv[1]: .12e} {lv[2]: .12e}")

        charges = structure.site_properties.get("charge", np.zeros(len(structure.species)))
        magmoms = structure.site_properties.get("magmom", np.zeros(len(structure.species)))
        for species, coord, charge, magmom in zip(structure.species, structure.cart_coords, charges, magmoms):
            content_lines.append(f"atom {coord[0]: .12e} {coord[1]: .12e} {coord[2]: .12e} {species}")
            if charge != 0:
                content_lines.append(f"     initial_charge {charge:.12e}")
            if magmom != 0:
                content_lines.append(f"     initial_moment {magmom:.12e}")

        return cls(_content="\n".join(content_lines), _structure=structure)

    @property
    def structure(self) -> Structure | Molecule:
        """Access structure for the file"""
        return self._structure

    @property
    def content(self) -> str:
        """Access the contents of the file"""
        return self._content

    def write_file(self, directory: str | Path | None = None, overwrite: bool = False) -> None:
        """Write the geometry.in file

        Args:
            directory (str | Path | None): The directory to write the geometry.in file
            overwrite (bool): If True allow to overwrite existing files
        """
        directory = directory or Path.cwd()

        if not overwrite and (Path(directory) / "geometry.in").exists():
            raise ValueError(f"geometry.in file exists in {directory}")

        with open(f"{directory}/geometry.in", mode="w") as file:
            file.write(f"#{'=' * 72}\n")
            file.write(f"# FHI-aims geometry file: {directory}/geometry.in\n")
            file.write("# File generated from pymatgen\n")
            file.write(f"# {time.asctime()}\n")
            file.write(f"#{'=' * 72}\n")
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
    """Class representing the FHI-aims cubes

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
    edges: Sequence[Sequence[float]] = field(default_factory=lambda: 0.1 * np.eye(3))
    points: Sequence[int] | tuple[int, int, int] = field(default_factory=lambda: [0, 0, 0])
    format: str = "cube"
    spin_state: int | None = None
    kpoint: int | None = None
    filename: str | None = None
    elf_type: int | None = None

    def __eq__(self, other: object) -> bool:
        """Check if two cubes are equal to each other"""
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

        if self.elf_type != other.elf_type:
            return False

        return True

    def __post_init__(self) -> None:
        """Check the inputted variables to make sure they are correct

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
        """Get the block of text for the control.in file of the Cube"""
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
        attrs = ("type", "origin", "edges", "points", "format", "spin_state", "kpoint", "filename", "elf_type")
        decoded = {key: MontyDecoder().process_decoded(dct[key]) for key in attrs}

        return cls(**decoded)


@dataclass
class AimsControlIn(MSONable):
    """Class representing and FHI-aims control.in file

    Attributes:
        _parameters (dict[str, Any]): The parameters dictionary containing all input
            flags (key) and values for the control.in file
    """

    _parameters: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Initialize the output list of _parameters"""
        self._parameters.setdefault("output", [])

    def __getitem__(self, key: str) -> Any:
        """Get an input parameter

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
        """Set an attribute of the class

        Args:
            key (str): The parameter to get
            value (Any): The value for that parameter
        """
        if key == "output":
            if isinstance(value, str):
                value = [value]
            self._parameters[key] += value
        else:
            self._parameters[key] = value

    def __delitem__(self, key: str) -> Any:
        """Delete a parameter from the input object

        Args:
        key (str): The key in the parameter to remove

        Returns:
            Either the value of the deleted parameter or None if key is
            not in self._parameters
        """
        return self._parameters.pop(key, None)

    @property
    def parameters(self) -> dict[str, Any]:
        """The dictionary of input parameters for control.in"""
        return self._parameters

    @parameters.setter
    def parameters(self, parameters: dict[str, Any]) -> None:
        """Reset a control.in inputs from a parameters dictionary

        Args:
            parameters (dict[str, Any]): The new set of parameters to use
        """
        self._parameters = parameters
        self._parameters.setdefault("output", [])

    def get_aims_control_parameter_str(self, key: str, value: Any, fmt: str) -> str:
        """Get the string needed to add a parameter to the control.in file

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
        self, structure: Structure | Molecule, verbose_header: bool = False, directory: str | Path | None = None
    ) -> str:
        """Get the content of the file

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

        cubes = parameters.pop("cubes", None)

        if verbose_header:
            content += "# \n# List of parameters used to initialize the calculator:"
            for param, val in parameters.items():
                content += f"#     {param}:{val}\n"
        content += f"{lim}\n"

        assert ("smearing" in parameters and "occupation_type" in parameters) is False

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
                    order = " %d" % order
                else:
                    order = ""

                content += self.get_aims_control_parameter_str("occupation_type", (name, width, order), "%s %f%s")
            elif key == "output":
                for output_type in value:
                    content += self.get_aims_control_parameter_str(key, output_type, "%s")
            elif key == "vdw_correction_hirshfeld" and value:
                content += self.get_aims_control_parameter_str(key, "", "%s")
            elif isinstance(value, bool):
                content += self.get_aims_control_parameter_str(key, str(value).lower(), ".%s.")
            elif isinstance(value, (tuple, list)):
                content += self.get_aims_control_parameter_str(key, " ".join(map(str, value)), "%s")
            elif isinstance(value, str):
                content += self.get_aims_control_parameter_str(key, value, "%s")
            else:
                content += self.get_aims_control_parameter_str(key, value, "%r")

        if cubes:
            for cube in cubes:
                content += cube.control_block

        content += f"{lim}\n\n"
        species_dir = self._parameters.get("species_dir", os.environ.get("AIMS_SPECIES_DIR"))
        content += self.get_species_block(structure, species_dir)

        return content

    def write_file(
        self,
        structure: Structure | Molecule,
        directory: str | Path | None = None,
        verbose_header: bool = False,
        overwrite: bool = False,
    ) -> None:
        """Writes the control.in file

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
            file.write(f"#{'=' * 72}\n")
            file.write(f"# FHI-aims geometry file: {directory}/geometry.in\n")
            file.write("# File generated from pymatgen\n")
            file.write(f"# {time.asctime()}\n")
            file.write(f"#{'=' * 72}\n")

            file.write(content)

    def get_species_block(self, structure: Structure | Molecule, species_dir: str | Path) -> str:
        """Get the basis set information for a structure

        Args:
            structure (Molecule or Structure): The structure to get the basis set information for
            species_dir (str or Pat:): The directory to find the species files in

        Returns:
            The block to add to the control.in file for the species

        Raises:
            ValueError: If a file for the species is not found
        """
        block = ""
        species = np.unique(structure.species)
        for sp in species:
            filename = f"{species_dir}/{sp.Z:02d}_{sp.symbol}_default"
            if Path(filename).exists():
                with open(filename) as sf:
                    block += "".join(sf.readlines())
            elif Path(f"{filename}.gz").exists():
                with gzip.open(f"{filename}.gz", mode="rt") as sf:
                    block += "".join(sf.readlines())
            else:
                raise ValueError(f"Species file for {sp.symbol} not found.")

        return block

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
