"""Module for defining common data used and produced by atomistic simulation packages."""

from __future__ import annotations

import collections
import importlib
import itertools
import os
import typing
import warnings
from copy import deepcopy
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import orjson
from monty.io import zopen
from monty.json import MSONable
from scipy.interpolate import RegularGridInterpolator

from pymatgen.core import Element, Site, Structure
from pymatgen.core.units import ang_to_bohr, bohr_to_angstrom
from pymatgen.electronic_structure.core import Spin

if TYPE_CHECKING:
    from numpy.typing import NDArray
    from typing_extensions import Any, Self

    from pymatgen.core.structure import IStructure


class VolumetricData(MSONable):
    """
    A representation of volumetric data commonly used in atomistic simulation outputs,
    such as LOCPOT/CHGCAR files from VASP or cube files from other codes.

    Attributes:
        structure (Structure):
            The crystal structure associated with the volumetric data.
            Represents the lattice and atomic coordinates using the `Structure` class.

        is_spin_polarized (bool):
            Indicates if the simulation is spin-polarized. True for spin-polarized data
            (contains both total and spin-difference densities), False otherwise.

        dim (tuple[int, int, int]):
            The dimensions of the 3D volumetric grid along each axis in the format
            (nx, ny, nz), where nx, ny, and nz represent the number of grid points
            in the x, y, and z directions, respectively.

        data (dict[str, np.ndarray]):
            A dictionary containing the volumetric data. Keys include:
            - `"total"`: A 3D NumPy array representing the total spin density.
            - `"diff"` (optional): A 3D NumPy array representing the spin-difference
              density (spin up - spin down). Typically present in spin-polarized calculations.

        ngridpts (int):
            The total number of grid points in the volumetric data, calculated as
            `nx * ny * nz` using the grid dimensions.
    """

    def __init__(
        self,
        structure: Structure | IStructure,
        data: dict[str, np.ndarray],
        distance_matrix: dict | None = None,
        data_aug: dict[str, NDArray] | None = None,
    ) -> None:
        """
        Typically, this constructor is not used directly and the static
        from_file constructor is used. This constructor is designed to allow
        summation and other operations between VolumetricData objects.

        Args:
            structure (Structure): associated with the volumetric data
            data (dict[str, NDArray]): Actual volumetric data.
            distance_matrix (NDArray): A pre-computed distance matrix if available.
                Useful so pass distance_matrices between sums,
                short-circuiting an otherwise expensive operation.
            data_aug (NDArray): Any extra information associated with volumetric data
                (typically augmentation charges)
        """
        self.structure = structure
        self.is_spin_polarized = len(data) >= 2
        self.is_soc = len(data) >= 4
        # convert data to numpy arrays in case they were jsanitized as lists
        self.data = {k: np.array(v) for k, v in data.items()}
        self.dim = self.data["total"].shape
        self.data_aug = data_aug or {}
        self.ngridpts = self.dim[0] * self.dim[1] * self.dim[2]
        # lazy init the spin data since this is not always needed.
        self._spin_data: dict[Spin, float] = {}
        self._distance_matrix = distance_matrix if distance_matrix is not None else {}
        self.xpoints = np.linspace(0.0, 1.0, num=self.dim[0])
        self.ypoints = np.linspace(0.0, 1.0, num=self.dim[1])
        self.zpoints = np.linspace(0.0, 1.0, num=self.dim[2])
        self.interpolator = RegularGridInterpolator(
            (self.xpoints, self.ypoints, self.zpoints),
            self.data["total"],
            bounds_error=True,
        )
        self.name = "VolumetricData"

    @property
    def spin_data(self):
        """The data decomposed into actual spin data as {spin: data}.
        Essentially, this provides the actual Spin.up and Spin.down data
        instead of the total and diff. Note that by definition, a
        non-spin-polarized run would have Spin.up data == Spin.down data.
        """
        if not self._spin_data:
            spin_data = {}
            spin_data[Spin.up] = 0.5 * (self.data["total"] + self.data.get("diff", 0))
            spin_data[Spin.down] = 0.5 * (self.data["total"] - self.data.get("diff", 0))
            self._spin_data = spin_data
        return self._spin_data

    def get_axis_grid(self, ind):
        """Get the grid for a particular axis.

        Args:
            ind (int): Axis index.
        """
        ng = self.dim
        num_pts = ng[ind]
        lengths = self.structure.lattice.abc
        return [i / num_pts * lengths[ind] for i in range(num_pts)]

    def __add__(self, other):
        return self.linear_add(other, 1.0)

    def __radd__(self, other):
        if other == 0 or other is None:
            # sum() calls 0 + self first; we treat 0 as the identity element
            return self
        if isinstance(other, self.__class__):
            return self.__add__(other)

        raise TypeError(f"Unsupported operand type(s) for +: '{type(other).__name__}' and '{type(self).__name__}'")

    def __sub__(self, other):
        return self.linear_add(other, -1.0)

    def copy(self) -> VolumetricData:
        """Make a copy of VolumetricData object."""
        return VolumetricData(
            self.structure,
            {k: v.copy() for k, v in self.data.items()},
            distance_matrix=self._distance_matrix,  # type:ignore[arg-type]
            data_aug=self.data_aug,  # type:ignore[arg-type]
        )

    def linear_add(self, other, scale_factor=1.0) -> VolumetricData:
        """
        Method to do a linear sum of volumetric objects. Used by + and -
        operators as well. Returns a VolumetricData object containing the
        linear sum.

        Args:
            other (VolumetricData): Another VolumetricData object
            scale_factor (float): Factor to scale the other data by.

        Returns:
            VolumetricData corresponding to self + scale_factor * other.
        """
        if self.structure != other.structure:
            warnings.warn(
                "Structures are different. Make sure you know what you are doing...",
                stacklevel=2,
            )
        if list(self.data) != list(other.data):
            raise ValueError("Data have different keys! Maybe one is spin-polarized and the other is not?")

        # To add checks
        data = {}
        for k in self.data:
            data[k] = self.data[k] + scale_factor * other.data[k]

        new = deepcopy(self)
        new.data = data
        new.data_aug = {}
        return new

    def scale(self, factor):
        """Scale the data in place by a factor."""
        for k in self.data:
            self.data[k] = np.multiply(self.data[k], factor)

    def value_at(self, x, y, z):
        """Get a data value from self.data at a given point (x, y, z) in terms
        of fractional lattice parameters. Will be interpolated using a
        RegularGridInterpolator on self.data if (x, y, z) is not in the original
        set of data points.

        Args:
            x (float): Fraction of lattice vector a.
            y (float): Fraction of lattice vector b.
            z (float): Fraction of lattice vector c.

        Returns:
            Value from self.data (potentially interpolated) corresponding to
            the point (x, y, z).
        """
        return self.interpolator([x, y, z])[0]

    def linear_slice(self, p1, p2, n=100):
        """Get a linear slice of the volumetric data with n data points from
        point p1 to point p2, in the form of a list.

        Args:
            p1 (list): 3-element list containing fractional coordinates of the first point.
            p2 (list): 3-element list containing fractional coordinates of the second point.
            n (int): Number of data points to collect, defaults to 100.

        Returns:
            List of n data points (mostly interpolated) representing a linear slice of the
            data from point p1 to point p2.
        """
        if type(p1) not in {list, np.ndarray}:
            raise TypeError(f"type of p1 should be list or np.ndarray, got {type(p1).__name__}")
        if len(p1) != 3:
            raise ValueError(f"length of p1 should be 3, got {len(p1)}")
        if type(p2) not in {list, np.ndarray}:
            raise TypeError(f"type of p2 should be list or np.ndarray, got {type(p2).__name__}")
        if len(p2) != 3:
            raise ValueError(f"length of p2 should be 3, got {len(p2)}")

        x_pts = np.linspace(p1[0], p2[0], num=n)
        y_pts = np.linspace(p1[1], p2[1], num=n)
        z_pts = np.linspace(p1[2], p2[2], num=n)
        return [self.value_at(x_pts[i], y_pts[i], z_pts[i]) for i in range(n)]

    def get_integrated_diff(self, ind, radius, nbins=1):
        """Get integrated difference of atom index ind up to radius. This can be
        an extremely computationally intensive process, depending on how many
        grid points are in the VolumetricData.

        Args:
            ind (int): Index of atom.
            radius (float): Radius of integration.
            nbins (int): Number of bins. Defaults to 1. This allows one to
                obtain the charge integration up to a list of the cumulative
                charge integration values for radii for [radius/nbins,
                2 * radius/nbins, ....].

        Returns:
            Differential integrated charge as a np array of [[radius, value],
            ...]. Format is for ease of plotting. e.g. plt.plot(data[:,0],
            data[:,1])
        """
        # For non-spin-polarized runs, this is zero by definition.
        if not self.is_spin_polarized:
            radii = [radius / nbins * (i + 1) for i in range(nbins)]
            data = np.zeros((nbins, 2))
            data[:, 0] = radii
            return data

        struct = self.structure
        a = self.dim
        self._distance_matrix = {} if self._distance_matrix is None else self._distance_matrix
        if ind not in self._distance_matrix or self._distance_matrix[ind]["max_radius"] < radius:
            coords = []
            for x, y, z in itertools.product(*(list(range(i)) for i in a)):
                coords.append([x / a[0], y / a[1], z / a[2]])
            sites_dist = struct.lattice.get_points_in_sphere(coords, struct[ind].coords, radius)
            self._distance_matrix[ind] = {
                "max_radius": radius,
                "data": np.array(sites_dist, dtype=object),
            }

        data = self._distance_matrix[ind]["data"]

        # Use boolean indexing to find all charges within the desired distance.
        inds = data[:, 1] <= radius
        dists = data[inds, 1]
        data_inds = np.rint(np.mod(list(data[inds, 0]), 1) * np.tile(a, (len(dists), 1))).astype(int)
        vals = [self.data["diff"][x, y, z] for x, y, z in data_inds]

        hist, edges = np.histogram(dists, bins=nbins, range=[0, radius], weights=vals)
        data = np.zeros((nbins, 2))
        data[:, 0] = edges[1:]
        data[:, 1] = [sum(hist[0 : i + 1]) / self.ngridpts for i in range(nbins)]
        return data

    def get_average_along_axis(self, ind):
        """Get the averaged total of the volumetric data a certain axis direction.
        For example, useful for visualizing Hartree Potentials from a LOCPOT
        file.

        Args:
            ind (int): Index of axis.

        Returns:
            Average total along axis
        """
        total_spin_dens = self.data["total"]
        ng = self.dim
        if ind == 0:
            total = np.sum(np.sum(total_spin_dens, axis=1), 1)
        elif ind == 1:
            total = np.sum(np.sum(total_spin_dens, axis=0), 1)
        else:
            total = np.sum(np.sum(total_spin_dens, axis=0), 0)
        return total / ng[(ind + 1) % 3] / ng[(ind + 2) % 3]

    def to_hdf5(self, filename):
        """Write the VolumetricData to a HDF5 format, which is a highly optimized
        format for reading storing large data. The mapping of the VolumetricData
        to this file format is as follows:

        VolumetricData.data -> f["vdata"]
        VolumetricData.structure ->
            f["Z"]: Sequence of atomic numbers
            f["fcoords"]: Fractional coords
            f["lattice"]: Lattice in the pymatgen.core.Lattice matrix
                format
            f.attrs["structure_json"]: String of JSON representation

        Args:
            filename (str): Filename to output to.
        """
        import h5py

        with h5py.File(filename, mode="w") as file:
            ds = file.create_dataset("lattice", (3, 3), dtype="float")
            ds[...] = self.structure.lattice.matrix
            ds = file.create_dataset("Z", (len(self.structure.species),), dtype="i")
            ds[...] = np.array([sp.Z for sp in self.structure.species])
            ds = file.create_dataset("fcoords", self.structure.frac_coords.shape, dtype="float")
            ds[...] = self.structure.frac_coords
            dt = h5py.special_dtype(vlen=str)
            ds = file.create_dataset("species", (len(self.structure.species),), dtype=dt)
            ds[...] = [str(sp) for sp in self.structure.species]
            grp = file.create_group("vdata")
            for k in self.data:
                ds = grp.create_dataset(k, self.data[k].shape, dtype="float")
                ds[...] = self.data[k]
            file.attrs["name"] = self.name
            file.attrs["structure_json"] = orjson.dumps(self.structure.as_dict()).decode()

    @classmethod
    def from_hdf5(cls, filename: str, **kwargs) -> VolumetricData:
        """
        Reads VolumetricData from HDF5 file.

        Args:
            filename: Filename

        Returns:
            VolumetricData
        """
        import h5py

        with h5py.File(filename, mode="r") as file:
            data = {k: np.array(v) for k, v in file["vdata"].items()}
            data_aug = None
            if "vdata_aug" in file:
                data_aug = {k: np.array(v) for k, v in file["vdata_aug"].items()}
            structure = Structure.from_dict(orjson.loads(file.attrs["structure_json"]))
            return cls(structure, data=data, data_aug=data_aug, **kwargs)  # type:ignore[arg-type]

    def to_cube(self, filename, comment: str = ""):
        """Write the total volumetric data to a cube file format, which consists of two comment lines,
        a header section defining the structure IN BOHR, and the data.

        Args:
            filename (str): Name of the cube file to be written.
            comment (str): If provided, this will be added to the second comment line
        """
        with zopen(filename, mode="wt", encoding="utf-8") as file:
            file.write(f"# Cube file for {self.structure.formula} generated by Pymatgen\n")  # type:ignore[arg-type]
            file.write(f"# {comment}\n")  # type:ignore[arg-type]  # type:ignore[arg-type]
            file.write(f"\t {len(self.structure)} 0.000000 0.000000 0.000000\n")  # type:ignore[arg-type]

            for idx in range(3):
                lattice_matrix = self.structure.lattice.matrix[idx] / self.dim[idx] * ang_to_bohr
                file.write(
                    f"\t {self.dim[idx]} {lattice_matrix[0]:.6f} {lattice_matrix[1]:.6f} {lattice_matrix[2]:.6f}\n"  # type:ignore[arg-type]
                )

            for site in self.structure:
                file.write(
                    f"\t {Element(site.species_string).Z} 0.000000 "  # type:ignore[arg-type]
                    f"{ang_to_bohr * site.coords[0]} "  # type:ignore[arg-type]
                    f"{ang_to_bohr * site.coords[1]} "  # type:ignore[arg-type]
                    f"{ang_to_bohr * site.coords[2]} \n"  # type:ignore[arg-type]
                )

            for idx, dat in enumerate(self.data["total"].flatten(), start=1):
                file.write(f"{' ' if dat > 0 else ''}{dat:.6e} ")  # type:ignore[arg-type]
                if idx % 6 == 0:
                    file.write("\n")  # type:ignore[arg-type]

    @classmethod
    def from_cube(cls, filename: str | Path) -> Self:
        """
        Initialize the cube object and store the data as pymatgen objects.

        Args:
            filename (str): of the cube to read
        """
        # Limit the number of I/O operations by reading the entire file at once
        with zopen(filename, mode="rt", encoding="utf-8") as f:
            lines = f.readlines()

        # Parse the number of atoms from line 3 (first two lines are headers)
        n_atoms = int(lines[2].split()[0])

        # Pre-parse voxel data into arrays
        def parse_voxel(line):
            parts = line.split()
            return int(parts[0]), np.array(parts[1:], dtype=float) * bohr_to_angstrom

        # Extract voxel information from lines 4-6
        num_x_voxels, voxel_x = parse_voxel(lines[3])
        num_y_voxels, voxel_y = parse_voxel(lines[4])
        num_z_voxels, voxel_z = parse_voxel(lines[5])

        # Parse atomic site data (starts at line 7 and continues for n_atoms lines)
        # Each line contains: atomic number, charge, x, y, z coordinates (in Bohr units)
        atom_data_start = 6
        atom_data_end = atom_data_start + n_atoms
        sites = [
            Site(line.split()[0], np.array(line.split()[2:], dtype=float) * bohr_to_angstrom)  # type:ignore[arg-type]
            for line in lines[atom_data_start:atom_data_end]
        ]

        # Generate a Pymatgen Structure object from extracted data
        structure = Structure(
            lattice=[voxel_x * num_x_voxels, voxel_y * num_y_voxels, voxel_z * num_z_voxels],
            species=[site.specie for site in sites],
            coords=[site.coords for site in sites],
            coords_are_cartesian=True,
        )

        # Extract volumetric data (starts after atomic site data)
        volumetric_data_start = atom_data_end
        volumetric_data_lines = lines[volumetric_data_start:]

        # Convert the text-based volumetric data into a NumPy array (much faster)
        volumetric_data = np.fromstring(" ".join(volumetric_data_lines), sep=" ", dtype=float)  # type:ignore[arg-type, type-var]

        # Reshape 1D data array into a 3D grid matching the voxel dimensions
        data = volumetric_data.reshape((num_x_voxels, num_y_voxels, num_z_voxels))

        return cls(structure=structure, data={"total": data})


class PMGDir(collections.abc.Mapping):
    """
    User-friendly class to access all files in a directory as pymatgen objects in a dict. For now, only VASP files are
    implemented but there is no reason why this cannot be extended to other types of files.
    Note that the files are lazily parsed to minimize initialization costs since not all files will be needed by all
    users.

    Example:

    ```
    d = PMGDir(".")
    print(d["INCAR"]["NELM"])
    print(d["vasprun.xml"].parameters)
    ```
    """

    FILE_MAPPINGS: typing.ClassVar = {
        n: f"pymatgen.io.vasp.{n.capitalize()}"
        for n in [
            "INCAR",
            "POSCAR",
            "KPOINTS",
            "POTCAR",
            "vasprun",
            "OUTCAR",
            "OSZICAR",
            "CHGCAR",
            "WAVECAR",
            "WAVEDER",
            "LOCPOT",
            "XDATCAR",
            "EIGENVAL",
            "PROCAR",
            "ELFCAR",
            "DYNMAT",
        ]
    } | {
        "CONTCAR": "pymatgen.io.vasp.Poscar",
        "IBZKPT": "pymatgen.io.vasp.Kpoints",
        "WSWQ": "pymatgen.io.vasp.WSWQ",
    }

    def __init__(self, dirname: str | Path):
        """
        Args:
            dirname: The directory containing the VASP calculation as a string or Path.
        """
        self.path = Path(dirname).absolute()
        self.reset()

    def reset(self):
        """
        Reset all loaded files and recheck the directory for files. Use this when the contents of the directory has
        changed.
        """
        # Note that py3.12 has Path.walk(). But we need to use os.walk to ensure backwards compatibility for now.
        self._files: dict[str, Any] = {
            str((Path(d) / f).relative_to(self.path)): None for d, _, fnames in os.walk(self.path) for f in fnames
        }

    def __contains__(self, item):
        return item in self._files

    def __len__(self):
        return len(self._files)

    def __iter__(self):
        return iter(self._files)

    def __getitem__(self, item):
        if self._files.get(item):
            return self._files.get(item)
        fpath = self.path / item

        if not (self.path / item).exists():
            raise ValueError(f"{item} not found in {self.path}. List of files are {self._files.keys()}.")

        for k, cls_ in PMGDir.FILE_MAPPINGS.items():
            if k in item:
                modname, classname = cls_.rsplit(".", 1)
                module = importlib.import_module(modname)
                class_ = getattr(module, classname)
                try:
                    self._files[item] = class_.from_file(fpath)
                except AttributeError:
                    self._files[item] = class_(fpath)

                return self._files[item]

        warnings.warn(
            f"No parser defined for {item}. Contents are returned as a string.",
            stacklevel=2,
        )
        with zopen(fpath, mode="rt", encoding="utf-8") as f:
            return f.read()

    def get_files_by_name(self, name: str) -> dict[str, Any]:
        """
        Returns all files with a given name. E.g., if you want all the OUTCAR files, set name="OUTCAR".

        Returns:
            {filename: object from PMGDir[filename]}
        """
        return {f: self[f] for f in self._files if name in f}

    def __repr__(self):
        return f"PMGDir({self.path})"
