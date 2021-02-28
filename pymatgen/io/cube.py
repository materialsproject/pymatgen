"""
Module for reading Gaussian cube files, which have become one of the standard file formats for volumetric data
in quantum chemistry and solid state physics software packages (VASP being an exception).

Some basic info about cube files (abridged info from http://paulbourke.net/dataformats/cube/ by Paul Bourke)

The file consists of a header which includes the atom information and the size as well as orientation of the
volumetric data. The first two lines of the header are comments. The third line has the number of atoms
included in the file followed by the position of the origin of the volumetric data. The next three lines
give the number of voxels along each axis (x, y, z) followed by the axis vector. The last section in the
header is one line for each atom consisting of 5 numbers, the first is the atom number, the second is the
charge, and the last three are the x,y,z coordinates of the atom center. The volumetric data is straightforward,
one floating point number for each volumetric element.


Example
In the following example the volumetric data is a 40 by 40 by 40 grid,
each voxel is 0.283459 units wide and the volume is aligned with the coordinate axis. There are three atoms.

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
from multiprocessing import Pool
from monty.io import zopen
from pymatgen import Site, Structure
from pymatgen.core.units import bohr_to_angstrom


class Cube:
    """
    Class to read Gaussian cube file formats for volumetric data.

    Cube files are, by default, written in units of 1/Bohr^3, and this file will convert
    the data into the pymatgen default of 1/Angstroms^3
    """

    def __init__(self, fname):
        """
        Initialize the cube object and store the data as self.data
        Args:
            fname (str): filename of the cube to read
        """
        f = zopen(fname, "rt")

        # skip header lines
        for i in range(2):
            f.readline()

        # number of atoms included in the file followed by the position of the origin of the volumetric data
        line = f.readline().split()
        self.natoms = int(line[0])
        self.origin = np.array(np.array(list(map(float, line[1:]))))

        # The next three lines give the number of voxels along each axis (x, y, z) followed by the axis vector.
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

        self.voxelVolume = abs(np.dot(np.cross(self.X, self.Y), self.Z))
        self.volume = abs(np.dot(np.cross(self.X.dot(self.NZ), self.Y.dot(self.NY)), self.Z.dot(self.NZ)))

        # The last section in the header is one line for each atom consisting of 5 numbers,
        # the first is the atom number, second is charge, the last three are the x,y,z coordinates of the atom center.
        self.sites = []
        for i in range(self.natoms):
            line = f.readline().split()
            self.sites.append(Site(line[0], np.multiply(bohr_to_angstrom, list(map(float, line[2:])))))

        self.structure = Structure(lattice=[self.X*self.NX, self.Y*self.NY, self.Z*self.NZ],
                                   species=[s.specie for s in self.sites],
                                   coords=[s.coords for s in self.sites], coords_are_cartesian=True)

        # Volumetric data
        self.data = np.zeros((self.NX, self.NY, self.NZ))
        i = 0
        for s in f:
            for v in s.split():
                self.data[
                    int(i / (self.NY * self.NZ)),
                    int((i / self.NZ) % self.NY),
                    int(i % self.NZ),
                ] = float(v)
                i += 1

    def mask_sphere(self, radius, cx, cy, cz):
        """
        Create a mask for a sphere with radius=radius, centered at cx, cy, cz.

        Args:
            radius: (flaot) of the mask (in Angstroms)
            Cx, Cy, Cz: (float) the fractional coordinates of the center of the sphere
        """

        r = np.floor(radius / np.linalg.norm(self.X)).astype(int), \
            np.floor(radius / np.linalg.norm(self.Y)).astype(int), \
            np.floor(radius / np.linalg.norm(self.Z)).astype(int)

        a = np.zeros((self.NX, self.NY, self.NZ))

        x0, y0, z0 = int(np.round(self.NX * cx)), \
            int(np.round(self.NY * cy)), \
            int(np.round(self.NZ * cz))

        for x in range(x0 - r[0], x0 + r[0] + 1):
            for y in range(y0 - r[1], y0 + r[1] + 1):
                for z in range(z0 - r[2], z0 + r[2] + 1):
                    dist = np.subtract(r, [abs(x0 - x), abs(y0 - y), abs(z0 - z)])

                    if all([_ > 0 for _ in dist]):
                        a[x % a.shape[0], y % a.shape[1], z % a.shape[2]] = 1
        return a

    def get_atomic_site_averages(self, atomic_site_radii=None, nproc=4):
        """
        Given a cube (pymatgen.io.cube.Cube), get the average value around each atomic site.

        Args:
            cube (Cube): pymatgen cube object
            atomic_site_radii (dict): dictionary determining the cutoff radius (in Angstroms)
                for averaging around atomic sites (e.g. {'Li': 0.97, 'B': 0.77, ...}. If
                not provided, then the
            nproc (int or None): number of processes to use in calculating atomic site averages. If
                nproc=None, then Pool will use os.cpu_count() to set nproc as the number of available
                cpus.

        returns:
            Array of site averages, [Average around site 1, Average around site 2, ...]
        """
        pool = Pool(nproc)
        results = pool.map(self, [(s, atomic_site_radii[s.species_string]) for s in self.structure.sites])
        pool.close()
        pool.join()
        return results

    def _get_atomic_site_averages(self, args):
        """
        Helper function for get_atomic_site_averages.

        Args:
            args: (tuple) Contains the site and the atomic_site_radius for given atomic species

        returns:
            Average around the atomic site
        """
        _s, _r = args[0], args[1]
        mask = self.mask_sphere(_r, *_s.frac_coords)
        return np.sum(self.data * mask) / np.count_nonzero(mask)

    def __call__(self, args):
        """
        Call function used by site averaging Pool
        """
        return self._get_atomic_site_averages(args)

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
        For example, useful for visualizing Hartree Potentials from a LOCPOT
        file.

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
