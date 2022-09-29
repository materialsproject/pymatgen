"""
Module for reading Gaussian cube files, which have become one of the standard file formats
for volumetric data in quantum chemistry and solid state physics software packages
(VASP being an exception).

Some basic info about cube files
(abridged info from http://paulbourke.net/dataformats/cube/ by Paul Bourke)

The file consists of a header which includes the atom information and the size as well
as orientation of the volumetric data. The first two lines of the header are comments. The
third line has the number of atoms included in the file followed by the position of the
origin of the volumetric data. The next three lines give the number of voxels along each axis
(x, y, z) followed by the axis vector. The last section in the header is one line for each
atom consisting of 5 numbers, the first is the atom number, the second is the charge, and
the last three are the x,y,z coordinates of the atom center. The volumetric data is straightforward,
one floating point number for each volumetric element.


Example
In the following example the volumetric data is a 40 by 40 by 40 grid, each voxel is 0.283459 units
wide and the volume is aligned with the coordinate axis. There are three atoms.

 CPMD CUBE FILE.
 OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z
    3    0.000000    0.000000    0.000000
   40    0.283459    0.000000    0.000000
   40    0.000000    0.283459    0.000000
   40    0.000000    0.000000    0.283459
    8    0.000000    5.570575    5.669178    5.593517
    1    0.000000    5.562867    5.669178    7.428055
    1    0.000000    7.340606    5.669178    5.111259
 -0.25568E-04  0.59213E-05  0.81068E-05  0.10868E-04  0.11313E-04  0.35999E-05
      :             :             :           :            :            :
      :             :             :           :            :            :
      :             :             :           :            :            :
        In this case there will be 40 x 40 x 40 floating point values
      :             :             :           :            :            :
      :             :             :           :            :            :
      :             :             :           :            :            :

"""

import numpy as np
from monty.io import zopen

from pymatgen.core.sites import Site
from pymatgen.core.structure import Structure
from pymatgen.core.units import bohr_to_angstrom


# TODO: can multiprocessing be incorporated without causing issues during drone assimilation?
class Cube:
    """
    Class to read Gaussian cube file formats for volumetric data.

    Cube files are, by default, written in atomic units, and this
    class assumes that convention.
    """

    def __init__(self, fname):
        """
        Initialize the cube object and store the data as self.data

        Args:
            fname (str): filename of the cube to read
        """
        f = zopen(fname, "rt")

        # skip header lines
        for _ in range(2):
            f.readline()

        # number of atoms followed by the position of the origin of the volumetric data
        line = f.readline().split()
        self.natoms = int(line[0])
        self.origin = np.array(list(map(float, line[1:])))

        # The number of voxels along each axis (x, y, z) followed by the axis vector.
        line = f.readline().split()
        self.NX = int(line[0])
        self.X = np.array([bohr_to_angstrom * float(l) for l in line[1:]])
        self.dX = np.linalg.norm(self.X)

        line = f.readline().split()
        self.NY = int(line[0])
        self.Y = np.array([bohr_to_angstrom * float(l) for l in line[1:]])
        self.dY = np.linalg.norm(self.Y)

        line = f.readline().split()
        self.NZ = int(line[0])
        self.Z = np.array([bohr_to_angstrom * float(l) for l in line[1:]])
        self.dZ = np.linalg.norm(self.Z)

        self.voxel_volume = abs(np.dot(np.cross(self.X, self.Y), self.Z))
        self.volume = abs(np.dot(np.cross(self.X.dot(self.NZ), self.Y.dot(self.NY)), self.Z.dot(self.NZ)))

        # The last section in the header is one line for each atom consisting of 5 numbers,
        # the first is the atom number, second is charge,
        # the last three are the x,y,z coordinates of the atom center.
        self.sites = []
        for _ in range(self.natoms):
            line = f.readline().split()
            self.sites.append(Site(line[0], np.multiply(bohr_to_angstrom, list(map(float, line[2:])))))

        self.structure = Structure(
            lattice=[self.X * self.NX, self.Y * self.NY, self.Z * self.NZ],
            species=[s.specie for s in self.sites],
            coords=[s.coords for s in self.sites],
            coords_are_cartesian=True,
        )

        # Volumetric data
        self.data = np.reshape(np.array(f.read().split()).astype(float), (self.NX, self.NY, self.NZ))

    def mask_sphere(self, radius, cx, cy, cz):
        """
        Create a mask for a sphere with radius=radius, centered at cx, cy, cz.

        Args:
            radius: (float) of the mask (in Angstroms)
            cx, cy, cz: (float) the fractional coordinates of the center of the sphere
        """
        dx, dy, dz = (
            np.floor(radius / np.linalg.norm(self.X)).astype(int),
            np.floor(radius / np.linalg.norm(self.Y)).astype(int),
            np.floor(radius / np.linalg.norm(self.Z)).astype(int),
        )
        gcd = max(np.gcd(dx, dy), np.gcd(dy, dz), np.gcd(dx, dz))
        sx, sy, sz = dx // gcd, dy // gcd, dz // gcd
        r = min(dx, dy, dz)

        x0, y0, z0 = int(np.round(self.NX * cx)), int(np.round(self.NY * cy)), int(np.round(self.NZ * cz))

        centerx, centery, centerz = self.NX // 2, self.NY // 2, self.NZ // 2
        a = np.roll(self.data, (centerx - x0, centery - y0, centerz - z0))

        i, j, k = np.indices(a.shape, sparse=True)
        a = np.sqrt((sx * i - sx * centerx) ** 2 + (sy * j - sy * centery) ** 2 + (sz * k - sz * centerz) ** 2)

        indices = a > r
        a[indices] = 0

        return a

    def get_atomic_site_averages(self, atomic_site_radii):
        """
        Get the average value around each atomic site.

        Args:
            atomic_site_radii (dict): dictionary determining the cutoff radius (in Angstroms)
                for averaging around atomic sites (e.g. {'Li': 0.97, 'B': 0.77, ...}. If
                not provided, then the
        returns:
            Array of site averages, [Average around site 1, Average around site 2, ...]
        """
        return [self._get_atomic_site_average(s, atomic_site_radii[s.species_string]) for s in self.structure.sites]

    def _get_atomic_site_average(self, site, radius):
        """
        Helper function for get_atomic_site_averages.

        Args:
            site: Site in the structure around which to get the average
            radius: (float) the atomic_site_radius (in Angstroms) for given atomic species

        returns:
            Average around the atomic site
        """
        mask = self.mask_sphere(radius, *site.frac_coords)
        return np.sum(self.data * mask) / np.count_nonzero(mask)

    def get_atomic_site_totals(self, atomic_site_radii):
        """
        Get the integrated total in a sphere around each atomic site.

        Args:
            atomic_site_radii (dict): dictionary determining the cutoff radius (in Angstroms)
                for averaging around atomic sites (e.g. {'Li': 0.97, 'B': 0.77, ...}. If
                not provided, then the
        returns:
            Array of site averages, [Average around site 1, Average around site 2, ...]
        """
        return [self._get_atomic_site_total(s, atomic_site_radii[s.species_string]) for s in self.structure.sites]

    def _get_atomic_site_total(self, site, radius):
        """
        Helper function for get_atomic_site_averages.

        Args:
            site: Site in the structure around which to get the total
            radius: (float) the atomic_site_radius (in Angstroms) for given atomic species

        returns:
            Average around the atomic site
        """
        mask = self.mask_sphere(radius, *site.frac_coords)
        return np.sum(self.data * mask)

    def get_axis_grid(self, ind):
        """
        Modified from pymatgen.io.vasp.outputs

        Returns the grid for a particular axis.

        Args:
            ind (int): Axis index.
        """
        ng = self.data.shape
        num_pts = ng[ind]
        lengths = self.structure.lattice.abc
        return [i / num_pts * lengths[ind] for i in range(num_pts)]

    def get_average_along_axis(self, ind):
        """
        Modified from pymatgen.io.vasp.outputs

        Get the averaged total of the volumetric data a certain axis direction.
        For example, useful for visualizing Hartree Potentials.

        Args:
            ind (int): Index of axis.

        Returns:
            Average total along axis
        """
        ng = self.data.shape
        m = self.data
        if ind == 0:
            total = np.sum(np.sum(m, axis=1), 1)
        elif ind == 1:
            total = np.sum(np.sum(m, axis=0), 1)
        else:
            total = np.sum(np.sum(m, axis=0), 0)
        return total / ng[(ind + 1) % 3] / ng[(ind + 2) % 3]
