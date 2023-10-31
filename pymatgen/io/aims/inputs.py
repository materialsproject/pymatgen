"""Classes for reading/manipulating/writing FHI-aims input files.

Works for geometry.in and control.in
"""

from __future__ import annotations

import gzip
import time
from copy import deepcopy
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from monty.json import MontyDecoder, MSONable

from pymatgen.core import Lattice, Molecule, Structure

if TYPE_CHECKING:
    from collections.abc import Sequence


@dataclass
class AimsGeometryIn(MSONable):
    """Class representing an aims geometry.in file"""

    _content: str
    _structure: Structure | Molecule

    @classmethod
    def from_str(cls, contents: str):
        """The contents of the input file

        self.__dict__
        ----------
        contents: str
            The content of the string
        """
        content_lines = [
            line.strip() for line in contents.split("\n") if len(line.strip()) > 0 and line.strip()[0] != "#"
        ]

        species = []
        coords = []
        is_frac = []
        lattice_vectors = []
        charges_dct = {}
        moments_dct = {}

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

        if lattice is None:
            structure = Molecule(
                species,
                coords,
                np.sum(charge),
                site_properties={"magmom": magmom, "charge": charge},
            )
        else:
            structure = Structure(
                lattice,
                species,
                coords,
                np.sum(charge),
                coords_are_cartesian=True,
                site_properties={"magmom": magmom, "charge": charge},
            )

        return cls(_content="\n".join(content_lines), _structure=structure)

    @classmethod
    def from_file(cls, filepath: str | Path):
        if str(filepath).endswith(".gz"):
            with gzip.open(filepath, "rt") as infile:
                content = infile.read()
        else:
            with open(filepath) as infile:
                content = infile.read()
        return cls.from_str(content)

    @classmethod
    def from_structure(cls, structure: Structure):
        """The contents of the input file

        self.__dict__
        ----------
        structure: Structure
            The structure for the file
        """
        content_lines = []

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
        """Accses structure for the file"""
        return self._structure

    @property
    def content(self) -> str:
        """Accses structure for the file"""
        return self._content

    def write_file(self, directory: str | Path | None = None, overwrite: bool = False):
        """Writes the geometry.in file

        self.__dict__
        ----------
        directory: str or Path
            The directory to write the geometry.in file
        overwrite: bool
            If True allow to overwrite existing files
        """
        if directory is None:
            directory = Path.cwd()

        if not overwrite and (Path(directory) / "geometry.in").exists():
            raise ValueError(f"geometry.in file exists in {directory}")

        with open(f"{directory}/geometry.in", "w") as fd:
            fd.write("#" + "=" * 72 + "\n")
            fd.write(f"# FHI-aims geometry file: {directory}/geometry.in\n")
            fd.write("# File generated from pymatgen\n")
            fd.write(f"# {time.asctime()}\n")
            fd.write("#" + "=" * 72 + "\n")
            fd.write(self.content)
            fd.write("\n")

    def as_dict(self) -> dict[str, Any]:
        """Get a dictionary representation of the geometry.in file.

        Returns
        -------
        The dictionary representation of the input file
        """
        dct = {}
        dct["@module"] = type(self).__module__
        dct["@class"] = type(self).__name__
        dct["content"] = self.content
        dct["structure"] = self.structure
        return dct

    @classmethod
    def from_dict(cls, d: dict[str, Any]):
        """Initialize from dictionary.

        self.__dict__
        ----------
        d: dict[str, Any]
            The MontyEncoded dictionary
        """
        decoded = {k: MontyDecoder().process_decoded(v) for k, v in d.items() if not k.startswith("@")}

        return cls(
            _content=decoded["content"],
            _structure=decoded["structure"],
        )


ALLOWED_AIMS_CUBE_TYPES = [
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
]

ALLOWED_AIMS_CUBE_TYPES_STATE = [
    "first_order_density",
    "eigenstate",
    "eigenstate_imag",
    "eigenstate_density",
]

ALLOWED_AIMS_CUBE_FORMATS = [
    "cube",
    "gOpenMol",
    "xsf",
]


@dataclass
class AimsCube(MSONable):
    """Class representing the FHI-aims cubes

    Parameters
    ----------
    type: str
        The value to be outputted as a cube file
    origin: Sequence[float] | tuple[float, float, float]
        The origin of the cube
    edges: Sequence[Sequence[float]]
        Specifies the edges of a cube: dx, dy, dz
        dx (float): The length of the step in the x direction
        dy (float): The length of the step in the y direction
        dx (float): The length of the step in the x direction
    points: Sequence[int] | tuple[int, int, int]
        The number of points along each edge
    spinstate: int
        The spin-channel to use either 1 or 2
    kpoint: int
        The k-point to use (the index of the list printed from `output k_point_list`)
    filename: str
        The filename to use
    format: str
        The format to output the cube file in: cube, gOpenMol, or xsf
    elf_type: int
        The type of electron localization function to use (see FHI-aims manual)
    """

    type: str = field(default_factory=str)
    origin: Sequence[float] | tuple[float, float, float] = field(default_factory=lambda: [0.0, 0.0, 0.0])
    edges: Sequence[Sequence[float]] = field(
        default_factory=lambda: [
            [0.1, 0.0, 0.0],
            [0.0, 0.1, 0.0],
            [0.0, 0.0, 0.1],
        ]
    )
    points: Sequence[int] | tuple[int, int, int] = field(default_factory=lambda: [0, 0, 0])
    format: str = "cube"
    spinstate: int | None = None
    kpoint: int | None = None
    filename: str | None = None
    elf_type: int | None = None

    def __post_init__(self):
        """Check the inputted variables to make sure they are correct"""
        split_type = self.type.split()
        if split_type[0] in ALLOWED_AIMS_CUBE_TYPES:
            if len(split_type) > 1:
                msg = f"Cube of type {split_type[0]} can not have a state associated with it"
                raise ValueError(msg)
        elif split_type[0] in ALLOWED_AIMS_CUBE_TYPES_STATE:
            if len(split_type) != 2:
                msg = f"Cube of type {split_type[0]} must have a state associated with it"
                raise ValueError(msg)
        else:
            raise ValueError("Cube type undefined")

        if self.format not in ALLOWED_AIMS_CUBE_FORMATS:
            raise ValueError("Cube file must have a format of cube, gOpenMol, or xsf")

        if self.spinstate is not None and (self.spinstate not in [1, 2]):
            raise ValueError("Spin state must be 1 or 2")

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
            raise ValueError("elf_type only used when the cube type is elf")

    @property
    def control_block(self):
        cb = f"output cube {self.type}\n"
        cb += f"    cube origin {self.origin[0]: .12e} {self.origin[1]: .12e} {self.origin[2]: .12e}\n"
        for ii in range(3):
            cb += f"    cube edge {self.points[ii]} "
            cb += f"{self.edges[ii][0]: .12e} "
            cb += f"{self.edges[ii][1]: .12e} "
            cb += f"{self.edges[ii][2]: .12e}\n"
        cb += f"    cube format {self.format}\n"
        if self.spinstate is not None:
            cb += f"    cube spinstate {self.spinstate}\n"
        if self.kpoint is not None:
            cb += f"    cube kpoint {self.kpoint}\n"
        if self.filename is not None:
            cb += f"    cube filename {self.filename}\n"
        if self.elf_type is not None:
            cb += f"    cube elf_type {self.elf_type}\n"

        return cb

    def as_dict(self) -> dict[str, Any]:
        """Get a dictionary representation of the geometry.in file.

        Returns
        -------
        The dictionary representation of the input file
        """
        dct: dict[str, Any] = {}
        dct["@module"] = type(self).__module__
        dct["@class"] = type(self).__name__
        dct["type"] = self.type
        dct["origin"] = self.origin
        dct["edges"] = self.edges
        dct["points"] = self.points
        dct["format"] = self.format
        dct["spinstate"] = self.spinstate
        dct["kpoint"] = self.kpoint
        dct["filename"] = self.filename
        dct["elf_type"] = self.elf_type
        return dct

    @classmethod
    def from_dict(cls, d: dict[str, Any]):
        """Initialize from dictionary.

        self.__dict__
        ----------
        d: dict[str, Any]
            The MontyEncoded dictionary
        """
        decoded = {k: MontyDecoder().process_decoded(v) for k, v in d.items() if not k.startswith("@")}

        return cls(
            type=decoded["type"],
            origin=decoded["origin"],
            edges=decoded["edges"],
            points=decoded["points"],
            format=decoded["format"],
            spinstate=decoded["spinstate"],
            kpoint=decoded["kpoint"],
            filename=decoded["filename"],
            elf_type=decoded["elf_type"],
        )


@dataclass
class AimsControlIn(MSONable):
    """Class representing and FHI-aims control.in file"""

    _parameters: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        if "output" not in self._parameters:
            self._parameters["output"] = []

    def __getitem__(self, key: str) -> Any:
        """Get an attribute of the class"""
        if key not in self._parameters:
            raise AttributeError(f"{key} not set in AimsControlIn")
        return self._parameters[key]

    def __setitem__(self, key: str, value: Any):
        """set an attribute of the class"""
        if key == "output":
            if isinstance(value, str):
                value = [value]
            self._parameters[key] += value
        else:
            self._parameters[key] = value

    def __delitem__(self, key: str) -> Any:
        """Delete a parameter from the input object

        Parameters
        ----------
        key: str
            The key in the parameter to remove
        """
        return self._parameters.pop(key, None)

    @property
    def parameters(self):
        return self._parameters

    @parameters.setter
    def parameters(self, parameters: dict[str, Any]):
        """reload a control.in inputs from a parameters dictionary"""
        self._parameters = parameters
        if "output" not in self._parameters:
            self._parameters["output"] = []

    def get_aims_control_parameter_str(self, key, value, format):
        return f"{key :35s}" + (format % value) + "\n"

    def write_file(
        self,
        structure: Structure,
        directory: str | Path | None = None,
        verbose_header: bool = False,
        overwrite: bool = False,
    ):
        """Writes the geometry.in file

        Parameters
        ----------
        structure: Structure
            The structure to write the input file for
        directory: str or Path
            The directory to write the geometry.in file
        verbose_header: bool
            If True print the input option dictionary
        overwrite: bool
            If True allow to overwrite existing files
        """
        if directory is None:
            directory = Path.cwd()

        if (Path(directory) / "control.in").exists() and not overwrite:
            raise ValueError(f"control.in file already in {directory}")

        lim = "#" + "=" * 79

        if isinstance(structure, Structure) and (
            "k_grid" not in self._parameters and "k_grid_density" not in self._parameters
        ):
            raise ValueError("k-grid must be defined for periodic systems")

        parameters = deepcopy(self._parameters)

        with open(f"{directory}/control.in", "w") as fd:
            fd.write("#" + "=" * 72 + "\n")
            fd.write(f"# FHI-aims geometry file: {directory}/geometry.in\n")
            fd.write("# File generated from pymatgen\n")
            fd.write(f"# {time.asctime()}\n")
            fd.write("#" + "=" * 72 + "\n")

            if parameters["xc"] == "LDA":
                parameters["xc"] = "pw-lda"

            cubes = parameters.pop("cubes", None)

            if verbose_header:
                fd.write("# \n# List of parameters used to initialize the calculator:")
                for p, v in parameters.items():
                    s = f"#     {p}:{v}\n"
                    fd.write(s)
            fd.write(lim + "\n")

            assert not ("smearing" in parameters and "occupation_type" in parameters)

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

                    fd.write(self.get_aims_control_parameter_str("occupation_type", (name, width, order), "%s %f%s"))
                elif key == "output":
                    for output_type in value:
                        fd.write(self.get_aims_control_parameter_str(key, output_type, "%s"))
                elif key == "vdw_correction_hirshfeld" and value:
                    fd.write(self.get_aims_control_parameter_str(key, "", "%s"))
                elif isinstance(value, bool):
                    fd.write(self.get_aims_control_parameter_str(key, str(value).lower(), ".%s."))
                elif isinstance(value, (tuple, list)):
                    fd.write(self.get_aims_control_parameter_str(key, " ".join([str(x) for x in value]), "%s"))
                elif isinstance(value, str):
                    fd.write(self.get_aims_control_parameter_str(key, value, "%s"))
                else:
                    fd.write(self.get_aims_control_parameter_str(key, value, "%r"))

            if cubes:
                for cube in cubes:
                    fd.write(cube.control_block)

            fd.write(lim + "\n\n")
            fd.write(self.get_species_block(structure, self._parameters["species_dir"]))

    def get_species_block(self, structure, species_dir):
        """Get the basis set information for a structure"""
        sb = ""
        species = np.unique(structure.species)
        for sp in species:
            filename = f"{species_dir}/{sp.Z:02d}_{sp.symbol}_default"
            if Path(filename).exists():
                with open(filename) as sf:
                    sb += "".join(sf.readlines())
            elif Path(f"{filename}.gz").exists():
                with gzip.open(f"{filename}.gz", "rt") as sf:
                    sb += "".join(sf.readlines())
            else:
                raise ValueError("Species file for {sp.symbol} not found.")

        return sb

    def as_dict(self) -> dict[str, Any]:
        """Get a dictionary representation of the geometry.in file.

        Returns
        -------
        The dictionary representation of the input file
        """
        dct = {}
        dct["@module"] = type(self).__module__
        dct["@class"] = type(self).__name__
        dct["parameters"] = self.parameters
        return dct

    @classmethod
    def from_dict(cls, d: dict[str, Any]):
        """Initialize from dictionary.

        self.__dict__
        ----------
        d: dict[str, Any]
            The MontyEncoded dictionary
        """
        decoded = {k: MontyDecoder().process_decoded(v) for k, v in d.items() if not k.startswith("@")}

        return cls(
            _parameters=decoded["parameters"],
        )
