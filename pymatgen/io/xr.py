# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module provides input and output mechanisms
for the xr file format, which is a modified CSSR
file format and, for example, used in GULP.
In particular, the module makes it easy
to remove shell positions from relaxations
that employed core-shell models.
"""

import re
from math import fabs

import numpy as np
from monty.io import zopen

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure

__author__ = "Nils Edvin Richard Zimmermann"
__copyright__ = "Copyright 2016, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Nils Edvin Richard Zimmermann"
__email__ = "nils.e.r.zimmermann@gmail.com"
__date__ = "June 23, 2016"


class Xr:
    """
    Basic object for working with xr files.
    """

    def __init__(self, structure):
        """
        Args:
            structure (Structure/IStructure): Structure object to create the
                    Xr object.
        """
        if not structure.is_ordered:
            raise ValueError("Xr file can only be constructed from ordered structure")
        self.structure = structure

    def __str__(self):
        a, b, c = self.structure.lattice.abc
        alpha, beta, gamma = self.structure.lattice.angles
        output = [
            f"pymatgen   {a:.4f} {b:.4f} {c:.4f}",
            f"{alpha:.3f} {beta:.3f} {gamma:.3f}",
            f"{len(self.structure)} 0",
            f"0 {self.structure.formula}",
        ]
        # There are actually 10 more fields per site
        # in a typical xr file from GULP, for example.
        for idx, site in enumerate(self.structure.sites):
            output.append(f"{idx + 1} {site.specie} {site.x:.4f} {site.y:.4f} {site.z:.4f}")
        mat = self.structure.lattice.matrix
        for _ in range(2):
            for j in range(3):
                output.append(f"{mat[j][0]:.4f} {mat[j][1]:.4f} {mat[j][2]:.4f}")
        return "\n".join(output)

    def write_file(self, filename):
        """
        Write out an xr file.

        Args:
            filename (str): name of the file to write to.
        """
        with zopen(filename, "wt") as f:
            f.write(str(self) + "\n")

    @staticmethod
    def from_string(string, use_cores=True, thresh=1.0e-4):
        """
        Creates an Xr object from a string representation.

        Args:
            string (str): string representation of an Xr object.
            use_cores (bool): use core positions and discard shell
                    positions if set to True (default). Otherwise,
                    use shell positions and discard core positions.
            thresh (float): relative threshold for consistency check
                    between cell parameters (lengths and angles) from
                    header information and cell vectors, respectively.

        Returns:
            xr (Xr): Xr object corresponding to the input
                    string representation.
        """
        lines = string.split("\n")
        toks = lines[0].split()
        lengths = [float(toks[i]) for i in range(1, len(toks))]
        toks = lines[1].split()
        angles = [float(i) for i in toks[0:3]]
        toks = lines[2].split()
        nsites = int(toks[0])
        mat = np.zeros((3, 3), dtype=float)
        for i in range(3):
            toks = lines[4 + nsites + i].split()
            toks2 = lines[4 + nsites + i + 3].split()
            for j, item in enumerate(toks):
                if item != toks2[j]:
                    raise RuntimeError("expected both matrices to be the same in xr file")
            mat[i] = np.array([float(w) for w in toks])
        lat = Lattice(mat)
        if (
            fabs(lat.a - lengths[0]) / fabs(lat.a) > thresh
            or fabs(lat.b - lengths[1]) / fabs(lat.b) > thresh
            or fabs(lat.c - lengths[2]) / fabs(lat.c) > thresh
            or fabs(lat.alpha - angles[0]) / fabs(lat.alpha) > thresh
            or fabs(lat.beta - angles[1]) / fabs(lat.beta) > thresh
            or fabs(lat.gamma - angles[2]) / fabs(lat.gamma) > thresh
        ):
            raise RuntimeError(
                "cell parameters in header ("
                + str(lengths)
                + ", "
                + str(angles)
                + ") are not consistent with Cartesian"
                + " lattice vectors ("
                + str(lat.abc)
                + ", "
                + str(lat.angles)
                + ")"
            )
        # Ignore line w/ index 3.
        sp = []
        coords = []
        for j in range(nsites):
            m = re.match(
                r"\d+\s+(\w+)\s+([0-9\-\.]+)\s+([0-9\-\.]+)\s+" + r"([0-9\-\.]+)",
                lines[4 + j].strip(),
            )
            if m:
                tmp_sp = m.group(1)
                if use_cores and tmp_sp[len(tmp_sp) - 2 :] == "_s":
                    continue
                if not use_cores and tmp_sp[len(tmp_sp) - 2 :] == "_c":
                    continue
                if tmp_sp[len(tmp_sp) - 2] == "_":
                    sp.append(tmp_sp[0 : len(tmp_sp) - 2])
                else:
                    sp.append(tmp_sp)
                coords.append([float(m.group(i)) for i in range(2, 5)])
        return Xr(Structure(lat, sp, coords, coords_are_cartesian=True))

    @staticmethod
    def from_file(filename, use_cores=True, thresh=1.0e-4):
        """
        Reads an xr-formatted file to create an Xr object.

        Args:
            filename (str): name of file to read from.
            use_cores (bool): use core positions and discard shell
                    positions if set to True (default). Otherwise,
                    use shell positions and discard core positions.
            thresh (float): relative threshold for consistency check
                    between cell parameters (lengths and angles) from
                    header information and cell vectors, respectively.

        Returns:
            xr (Xr): Xr object corresponding to the input
                    file.
        """
        with zopen(filename, "rt") as f:
            return Xr.from_string(f.read(), use_cores=use_cores, thresh=thresh)
