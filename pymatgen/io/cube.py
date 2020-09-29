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
        # the first is the atom number, second (?), the last three are the x,y,z coordinates of the atom center.
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
