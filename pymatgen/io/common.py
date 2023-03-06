"""
Module for defining common data used and produced by atomistic simulation packages.
"""

from __future__ import annotations

import itertools
import json
import warnings
from copy import deepcopy

import numpy as np
from monty.io import zopen
from monty.json import MSONable
from scipy.interpolate import RegularGridInterpolator

from pymatgen.core import Element, Site, Structure
from pymatgen.core.units import ang_to_bohr, bohr_to_angstrom
from pymatgen.electronic_structure.core import Spin


class VolumetricData(MSONable):
    """
    Simple volumetric object. Used to read LOCPOT/CHGCAR files produced by
    vasp as well as cube files produced by other codes.

    .. attribute:: structure

        Structure associated with the Volumetric Data object

    ..attribute:: is_spin_polarized

        True if run is spin polarized

    ..attribute:: dim

        Tuple of dimensions of volumetric grid in each direction (nx, ny, nz).

    ..attribute:: data

        Actual data as a dict of {string: np.array}. The string are "total"
        and "diff", in accordance to the output format of Vasp LOCPOT and
        CHGCAR files where the total spin density is written first, followed
        by the difference spin density.

    .. attribute:: ngridpts

        Total number of grid points in volumetric data.
    """

    def __init__(self, structure: Structure, data, distance_matrix=None, data_aug=None):
        """
        Typically, this constructor is not used directly and the static
        from_file constructor is used. This constructor is designed to allow
        summation and other operations between VolumetricData objects.

        Args:
            structure: Structure associated with the volumetric data
            data: Actual volumetric data. If the data is provided as in list format,
                it will be converted into an np.array automatically
            data_aug: Any extra information associated with volumetric data
                (typically augmentation charges)
            distance_matrix: A pre-computed distance matrix if available.
                Useful so pass distance_matrices between sums,
                short-circuiting an otherwise expensive operation.
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
        self._distance_matrix = distance_matrix if distance_matrix else {}
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
        """
        The data decomposed into actual spin data as {spin: data}.
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
        """
        Returns the grid for a particular axis.

        Args:
            ind (int): Axis index.
        """
        ng = self.dim
        num_pts = ng[ind]
        lengths = self.structure.lattice.abc
        return [i / num_pts * lengths[ind] for i in range(num_pts)]

    def __add__(self, other):
        return self.linear_add(other, 1.0)

    def __sub__(self, other):
        return self.linear_add(other, -1.0)

    def copy(self):
        """
        :return: Copy of Volumetric object
        """
        return VolumetricData(
            self.structure,
            {k: v.copy() for k, v in self.data.items()},
            distance_matrix=self._distance_matrix,
            data_aug=self.data_aug,
        )

    def linear_add(self, other, scale_factor=1.0):
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
            warnings.warn("Structures are different. Make sure you know what you are doing...")
        if list(self.data) != list(other.data):
            raise ValueError("Data have different keys! Maybe one is spin-" "polarized and the other is not?")

        # To add checks
        data = {}
        for k in self.data:
            data[k] = self.data[k] + scale_factor * other.data[k]

        new = deepcopy(self)
        new.data = data
        new.data_aug = {}
        return new

    def scale(self, factor):
        """
        Scale the data in place by a factor.
        """
        for k in self.data:
            self.data[k] = np.multiply(self.data[k], factor)

    def value_at(self, x, y, z):
        """
        Get a data value from self.data at a given point (x, y, z) in terms
        of fractional lattice parameters. Will be interpolated using a
        RegularGridInterpolator on self.data if (x, y, z) is not in the original
        set of data points.

        Args:
            x (float): Fraction of lattice vector a.
            y (float): Fraction of lattice vector b.
            z (float): Fraction of lattice vector c.

        Returns:
            Value from self.data (potentially interpolated) correspondisng to
            the point (x, y, z).
        """
        return self.interpolator([x, y, z])[0]

    def linear_slice(self, p1, p2, n=100):
        """
        Get a linear slice of the volumetric data with n data points from
        point p1 to point p2, in the form of a list.

        Args:
            p1 (list): 3-element list containing fractional coordinates of the first point.
            p2 (list): 3-element list containing fractional coordinates of the second point.
            n (int): Number of data points to collect, defaults to 100.

        Returns:
            List of n data points (mostly interpolated) representing a linear slice of the
            data from point p1 to point p2.
        """
        assert type(p1) in [list, np.ndarray]
        assert type(p2) in [list, np.ndarray]
        assert len(p1) == 3
        assert len(p2) == 3
        xpts = np.linspace(p1[0], p2[0], num=n)
        ypts = np.linspace(p1[1], p2[1], num=n)
        zpts = np.linspace(p1[2], p2[2], num=n)
        return [self.value_at(xpts[i], ypts[i], zpts[i]) for i in range(n)]

    def get_integrated_diff(self, ind, radius, nbins=1):
        """
        Get integrated difference of atom index ind up to radius. This can be
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
            ...]. Format is for ease of plotting. E.g., plt.plot(data[:,0],
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
        """
        Get the averaged total of the volumetric data a certain axis direction.
        For example, useful for visualizing Hartree Potentials from a LOCPOT
        file.

        Args:
            ind (int): Index of axis.

        Returns:
            Average total along axis
        """
        m = self.data["total"]
        ng = self.dim
        if ind == 0:
            total = np.sum(np.sum(m, axis=1), 1)
        elif ind == 1:
            total = np.sum(np.sum(m, axis=0), 1)
        else:
            total = np.sum(np.sum(m, axis=0), 0)
        return total / ng[(ind + 1) % 3] / ng[(ind + 2) % 3]

    def to_hdf5(self, filename):
        """
        Writes the VolumetricData to a HDF5 format, which is a highly optimized
        format for reading storing large data. The mapping of the VolumetricData
        to this file format is as follows:

        VolumetricData.data -> f["vdata"]
        VolumetricData.structure ->
            f["Z"]: Sequence of atomic numbers
            f["fcoords"]: Fractional coords
            f["lattice"]: Lattice in the pymatgen.core.lattice.Lattice matrix
                format
            f.attrs["structure_json"]: String of json representation

        Args:
            filename (str): Filename to output to.
        """
        import h5py

        with h5py.File(filename, "w") as f:
            ds = f.create_dataset("lattice", (3, 3), dtype="float")
            ds[...] = self.structure.lattice.matrix
            ds = f.create_dataset("Z", (len(self.structure.species),), dtype="i")
            ds[...] = np.array([sp.Z for sp in self.structure.species])
            ds = f.create_dataset("fcoords", self.structure.frac_coords.shape, dtype="float")
            ds[...] = self.structure.frac_coords
            dt = h5py.special_dtype(vlen=str)
            ds = f.create_dataset("species", (len(self.structure.species),), dtype=dt)
            ds[...] = [str(sp) for sp in self.structure.species]
            grp = f.create_group("vdata")
            for k in self.data:
                ds = grp.create_dataset(k, self.data[k].shape, dtype="float")
                ds[...] = self.data[k]
            f.attrs["name"] = self.name
            f.attrs["structure_json"] = json.dumps(self.structure.as_dict())

    @classmethod
    def from_hdf5(cls, filename, **kwargs):
        """
        Reads VolumetricData from HDF5 file.

        :param filename: Filename
        :return: VolumetricData
        """
        import h5py

        with h5py.File(filename, "r") as f:
            data = {k: np.array(v) for k, v in f["vdata"].items()}
            data_aug = None
            if "vdata_aug" in f:
                data_aug = {k: np.array(v) for k, v in f["vdata_aug"].items()}
            structure = Structure.from_dict(json.loads(f.attrs["structure_json"]))
            return cls(structure, data=data, data_aug=data_aug, **kwargs)

    def to_cube(self, filename, comment=None):
        """
        Write the total volumetric data to a cube file format, which consists of two comment lines,
        a header section defining the structure IN BOHR, and the data.

        Args:
            filename (str): Name of the cube file to be written.
            comment (str): If provided, this will be added to the second comment line

        """
        with zopen(filename, "wt") as file:
            file.write(f"# Cube file for {self.structure.formula} generated by Pymatgen" + "\n")
            file.write(f"# {comment if comment else ''}" + "\n")
            file.write(f"\t {self.structure.num_sites} 0.000000 0.000000 0.000000" + "\n")
            file.write(
                f"\t {self.dim[0]} "
                f"{self.structure.lattice.matrix[0][0] * ang_to_bohr :.6f} "
                f"{self.structure.lattice.matrix[0][1] * ang_to_bohr :.6f} "
                f"{self.structure.lattice.matrix[0][2] * ang_to_bohr :.6f} \n"
            )
            file.write(
                f"\t {self.dim[1]} "
                f"{self.structure.lattice.matrix[1][0] * ang_to_bohr :.6f} "
                f"{self.structure.lattice.matrix[1][1] * ang_to_bohr :.6f} "
                f"{self.structure.lattice.matrix[1][2] * ang_to_bohr :.6f} \n"
            )
            file.write(
                f"\t {self.dim[2]} "
                f"{self.structure.lattice.matrix[2][0] * ang_to_bohr :.6f} "
                f"{self.structure.lattice.matrix[2][1] * ang_to_bohr :.6f} "
                f"{self.structure.lattice.matrix[2][2] * ang_to_bohr :.6f} \n"
            )

            for site in self.structure.sites:
                file.write(
                    f"\t {Element(site.species_string).Z} 0.000000 "
                    f"{ang_to_bohr * site.coords[0]} "
                    f"{ang_to_bohr * site.coords[1]} "
                    f"{ang_to_bohr * site.coords[2]} \n"
                )

            for i, dat in enumerate(self.data["total"].flatten()):
                file.write(f"{' ' if dat > 0 else ''}{dat:.6e} ")
                if (i + 1) % 6 == 0:
                    file.write("\n")

    @classmethod
    def from_cube(cls, filename):
        """
        Initialize the cube object and store the data as data

        Args:
            filename (str): of the cube to read
        """
        file = zopen(filename, "rt")

        # skip header lines
        file.readline()
        file.readline()

        # number of atoms followed by the position of the origin of the volumetric data
        line = file.readline().split()
        natoms = int(line[0])

        # The number of voxels along each axis (x, y, z) followed by the axis vector.
        line = file.readline().split()
        num_x_voxels = int(line[0])
        voxel_x = np.array([bohr_to_angstrom * float(l) for l in line[1:]])

        line = file.readline().split()
        num_y_voxels = int(line[0])
        voxel_y = np.array([bohr_to_angstrom * float(l) for l in line[1:]])

        line = file.readline().split()
        num_z_voxels = int(line[0])
        voxel_z = np.array([bohr_to_angstrom * float(l) for l in line[1:]])

        # The last section in the header is one line for each atom consisting of 5 numbers,
        # the first is the atom number, second is charge,
        # the last three are the x,y,z coordinates of the atom center.
        sites = []
        for _ in range(natoms):
            line = file.readline().split()
            sites.append(Site(line[0], np.multiply(bohr_to_angstrom, list(map(float, line[2:])))))

        structure = Structure(
            lattice=[voxel_x * num_x_voxels, voxel_y * num_y_voxels, voxel_z * num_z_voxels],
            species=[s.specie for s in sites],
            coords=[s.coords for s in sites],
            coords_are_cartesian=True,
        )

        # Volumetric data
        data = np.reshape(np.array(file.read().split()).astype(float), (num_x_voxels, num_y_voxels, num_z_voxels))
        return cls(structure=structure, data={"total": data})
