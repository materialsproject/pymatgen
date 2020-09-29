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
from monty.io import zopen
from pymatgen import Site, Structure
from pymatgen.core.units import bohr_to_angstrom


class Cube:
    """
    Class to read Gaussian cube file formats for volumetric data.
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

        line = f.readline().split()
        self.NY = int(line[0])
        self.Y = np.array([bohr_to_angstrom * float(l) for l in line[1:]])

        line = f.readline().split()
        self.NZ = int(line[0])
        self.Z = np.array([bohr_to_angstrom * float(l) for l in line[1:]])

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
