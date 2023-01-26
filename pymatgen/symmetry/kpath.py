# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Provides classes for generating high-symmetry k-paths using different conventions.
"""

from __future__ import annotations

import abc
import itertools
import operator
from math import ceil, cos, e, pi, sin, tan
from typing import Any
from warnings import warn

import networkx as nx
import numpy as np
import spglib
from monty.dev import requires
from scipy.linalg import sqrtm

from pymatgen.core.lattice import Lattice
from pymatgen.core.operations import MagSymmOp, SymmOp
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.typing import SpeciesLike

try:
    from seekpath import get_path
except ImportError:
    get_path = None

__author__ = "Geoffroy Hautier, Katherine Latimer, Jason Munro"
__copyright__ = "Copyright 2020, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Jason Munro"
__email__ = "jmunro@lbl.gov"
__status__ = "Development"
__date__ = "March 2020"


class KPathBase(metaclass=abc.ABCMeta):
    """
    This is the base class for classes used to generate high-symmetry
    paths in reciprocal space (k-paths) for band structure calculations.
    """

    @abc.abstractmethod
    def __init__(self, structure: Structure, symprec: float = 0.01, angle_tolerance=5, atol=1e-5, *args, **kwargs):
        """
        Args:
        structure (Structure): Structure object.
        symprec (float): Tolerance for symmetry finding.
        angle_tolerance (float): Angle tolerance for symmetry finding.
        atol (float): Absolute tolerance used to compare structures
            and determine symmetric equivalence of points and lines
            in the BZ.
        """
        self._structure = structure
        self._latt = self._structure.lattice
        self._rec_lattice = self._structure.lattice.reciprocal_lattice
        self._kpath: dict[str, Any] | None = None

        self._symprec = symprec
        self._atol = atol
        self._angle_tolerance = angle_tolerance

    @property
    def structure(self):
        """
        Returns:
        The input structure
        """
        return self._structure

    @property
    def lattice(self):
        """
        Returns:
        The real space lattice
        """
        return self._latt

    @property
    def rec_lattice(self):
        """
        Returns:
        The reciprocal space lattice
        """
        return self._rec_lattice

    @property
    def kpath(self):
        """
        Returns:
        The symmetry line path in reciprocal space
        """
        return self._kpath

    def get_kpoints(self, line_density=20, coords_are_cartesian=True):
        """
        Returns:
        kpoints along the path in Cartesian coordinates
        together with the critical-point labels.
        """
        list_k_points = []
        sym_point_labels = []
        for b in self.kpath["path"]:
            for i in range(1, len(b)):
                start = np.array(self.kpath["kpoints"][b[i - 1]])
                end = np.array(self.kpath["kpoints"][b[i]])
                distance = np.linalg.norm(
                    self._rec_lattice.get_cartesian_coords(start) - self._rec_lattice.get_cartesian_coords(end)
                )
                nb = int(ceil(distance * line_density))
                if nb == 0:
                    continue
                sym_point_labels.extend([b[i - 1]] + [""] * (nb - 1) + [b[i]])
                list_k_points.extend(
                    [
                        self._rec_lattice.get_cartesian_coords(start)
                        + float(i)
                        / float(nb)
                        * (self._rec_lattice.get_cartesian_coords(end) - self._rec_lattice.get_cartesian_coords(start))
                        for i in range(0, nb + 1)
                    ]
                )
        if coords_are_cartesian:
            return list_k_points, sym_point_labels

        frac_k_points = [self._rec_lattice.get_fractional_coords(k) for k in list_k_points]
        return frac_k_points, sym_point_labels


class KPathSetyawanCurtarolo(KPathBase):
    """
    This class looks for a path along high-symmetry lines in
    the Brillouin zone.
    It is based on Setyawan, W., & Curtarolo, S. (2010).
    High-throughput electronic band structure calculations:
    Challenges and tools. Computational Materials Science,
    49(2), 299-312. doi:10.1016/j.commatsci.2010.05.010
    It should be used with primitive structures that
    comply with the definition given in the paper.
    The symmetry is determined by spglib using the
    SpacegroupAnalyzer class. The analyzer can be used to
    produce the correct primitive structure with the method
    get_primitive_standard_structure(international_monoclinic=False).
    A warning will signal possible compatibility problems
    with the given structure. k-points generated using the get_kpoints() method
    are returned for the reciprocal cell basis defined in the paper.
    """

    def __init__(self, structure: Structure, symprec: float = 0.01, angle_tolerance=5, atol=1e-5):
        """
        Args:
        structure (Structure): Structure object.
        symprec (float): Tolerance for symmetry finding.
        angle_tolerance (float): Angle tolerance for symmetry finding.
        atol (float): Absolute tolerance used to compare the input
            structure with the one expected as primitive standard.
            A warning will be issued if the cells don't match.
        """
        if "magmom" in structure.site_properties:
            warn(
                "'magmom' entry found in site properties but will be ignored \
                  for the Setyawan and Curtarolo convention."
            )

        super().__init__(structure, symprec=symprec, angle_tolerance=angle_tolerance, atol=atol)

        self._sym = SpacegroupAnalyzer(structure, symprec=symprec, angle_tolerance=angle_tolerance)
        self._prim = self._sym.get_primitive_standard_structure(international_monoclinic=False)
        self._conv = self._sym.get_conventional_standard_structure(international_monoclinic=False)
        self._rec_lattice = self._prim.lattice.reciprocal_lattice

        # Note: this warning will be issued for space groups 38-41, since the primitive cell must be
        # reformatted to match the Setyawan/Curtarolo convention in order to work with the current k-path
        # generation scheme.
        if not np.allclose(self._structure.lattice.matrix, self._prim.lattice.matrix, atol=atol):
            warn(
                "The input structure does not match the expected standard primitive! "
                "The path may be incorrect. Use at your own risk."
            )

        lattice_type = self._sym.get_lattice_type()
        spg_symbol = self._sym.get_space_group_symbol()

        if lattice_type == "cubic":
            if "P" in spg_symbol:
                self._kpath = self.cubic()
            elif "F" in spg_symbol:
                self._kpath = self.fcc()
            elif "I" in spg_symbol:
                self._kpath = self.bcc()
            else:
                warn(f"Unexpected value for spg_symbol: {spg_symbol}")

        elif lattice_type == "tetragonal":
            if "P" in spg_symbol:
                self._kpath = self.tet()
            elif "I" in spg_symbol:
                a = self._conv.lattice.abc[0]
                c = self._conv.lattice.abc[2]
                if c < a:
                    self._kpath = self.bctet1(c, a)
                else:
                    self._kpath = self.bctet2(c, a)
            else:
                warn(f"Unexpected value for spg_symbol: {spg_symbol}")

        elif lattice_type == "orthorhombic":
            a = self._conv.lattice.abc[0]
            b = self._conv.lattice.abc[1]
            c = self._conv.lattice.abc[2]

            if "P" in spg_symbol:
                self._kpath = self.orc()

            elif "F" in spg_symbol:
                if 1 / a**2 > 1 / b**2 + 1 / c**2:
                    self._kpath = self.orcf1(a, b, c)
                elif 1 / a**2 < 1 / b**2 + 1 / c**2:
                    self._kpath = self.orcf2(a, b, c)
                else:
                    self._kpath = self.orcf3(a, b, c)

            elif "I" in spg_symbol:
                self._kpath = self.orci(a, b, c)

            elif "C" in spg_symbol or "A" in spg_symbol:
                self._kpath = self.orcc(a, b, c)
            else:
                warn(f"Unexpected value for spg_symbol: {spg_symbol}")

        elif lattice_type == "hexagonal":
            self._kpath = self.hex()

        elif lattice_type == "rhombohedral":
            alpha = self._prim.lattice.parameters[3]
            if alpha < 90:
                self._kpath = self.rhl1(alpha * pi / 180)
            else:
                self._kpath = self.rhl2(alpha * pi / 180)

        elif lattice_type == "monoclinic":
            a, b, c = self._conv.lattice.abc
            alpha = self._conv.lattice.parameters[3]
            # beta = self._conv.lattice.parameters[4]

            if "P" in spg_symbol:
                self._kpath = self.mcl(b, c, alpha * pi / 180)

            elif "C" in spg_symbol:
                kgamma = self._rec_lattice.parameters[5]
                if kgamma > 90:
                    self._kpath = self.mclc1(a, b, c, alpha * pi / 180)
                if kgamma == 90:
                    self._kpath = self.mclc2(a, b, c, alpha * pi / 180)
                if kgamma < 90:
                    if b * cos(alpha * pi / 180) / c + b**2 * sin(alpha * pi / 180) ** 2 / a**2 < 1:
                        self._kpath = self.mclc3(a, b, c, alpha * pi / 180)
                    if b * cos(alpha * pi / 180) / c + b**2 * sin(alpha * pi / 180) ** 2 / a**2 == 1:
                        self._kpath = self.mclc4(a, b, c, alpha * pi / 180)
                    if b * cos(alpha * pi / 180) / c + b**2 * sin(alpha * pi / 180) ** 2 / a**2 > 1:
                        self._kpath = self.mclc5(a, b, c, alpha * pi / 180)
            else:
                warn(f"Unexpected value for spg_symbol: {spg_symbol}")

        elif lattice_type == "triclinic":
            kalpha = self._rec_lattice.parameters[3]
            kbeta = self._rec_lattice.parameters[4]
            kgamma = self._rec_lattice.parameters[5]
            if kalpha > 90 and kbeta > 90 and kgamma > 90:
                self._kpath = self.tria()
            if kalpha < 90 and kbeta < 90 and kgamma < 90:
                self._kpath = self.trib()
            if kalpha > 90 and kbeta > 90 and kgamma == 90:
                self._kpath = self.tria()
            if kalpha < 90 and kbeta < 90 and kgamma == 90:
                self._kpath = self.trib()

        else:
            warn(f"Unknown lattice type {lattice_type}")

    @property
    def conventional(self):
        """
        Returns:
            The conventional cell structure
        """
        return self._conv

    @property
    def prim(self):
        """
        Returns:
            The primitive cell structure
        """
        return self._prim

    @property
    def prim_rec(self):
        """
        Returns:
            The primitive reciprocal cell structure
        """
        return self._rec_lattice

    def cubic(self):
        """
        CUB Path
        """
        self.name = "CUB"
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "X": np.array([0.0, 0.5, 0.0]),
            "R": np.array([0.5, 0.5, 0.5]),
            "M": np.array([0.5, 0.5, 0.0]),
        }
        path = [["\\Gamma", "X", "M", "\\Gamma", "R", "X"], ["M", "R"]]
        return {"kpoints": kpoints, "path": path}

    def fcc(self):
        """
        FCC Path
        """
        self.name = "FCC"
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "K": np.array([3.0 / 8.0, 3.0 / 8.0, 3.0 / 4.0]),
            "L": np.array([0.5, 0.5, 0.5]),
            "U": np.array([5.0 / 8.0, 1.0 / 4.0, 5.0 / 8.0]),
            "W": np.array([0.5, 1.0 / 4.0, 3.0 / 4.0]),
            "X": np.array([0.5, 0.0, 0.5]),
        }
        path = [
            ["\\Gamma", "X", "W", "K", "\\Gamma", "L", "U", "W", "L", "K"],
            ["U", "X"],
        ]
        return {"kpoints": kpoints, "path": path}

    def bcc(self):
        """
        BCC Path
        """
        self.name = "BCC"
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "H": np.array([0.5, -0.5, 0.5]),
            "P": np.array([0.25, 0.25, 0.25]),
            "N": np.array([0.0, 0.0, 0.5]),
        }
        path = [["\\Gamma", "H", "N", "\\Gamma", "P", "H"], ["P", "N"]]
        return {"kpoints": kpoints, "path": path}

    def tet(self):
        """
        TET Path
        """
        self.name = "TET"
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "A": np.array([0.5, 0.5, 0.5]),
            "M": np.array([0.5, 0.5, 0.0]),
            "R": np.array([0.0, 0.5, 0.5]),
            "X": np.array([0.0, 0.5, 0.0]),
            "Z": np.array([0.0, 0.0, 0.5]),
        }
        path = [
            ["\\Gamma", "X", "M", "\\Gamma", "Z", "R", "A", "Z"],
            ["X", "R"],
            ["M", "A"],
        ]
        return {"kpoints": kpoints, "path": path}

    def bctet1(self, c, a):
        """
        BCT1 Path
        """
        self.name = "BCT1"
        eta = (1 + c**2 / a**2) / 4.0
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "M": np.array([-0.5, 0.5, 0.5]),
            "N": np.array([0.0, 0.5, 0.0]),
            "P": np.array([0.25, 0.25, 0.25]),
            "X": np.array([0.0, 0.0, 0.5]),
            "Z": np.array([eta, eta, -eta]),
            "Z_1": np.array([-eta, 1 - eta, eta]),
        }
        path = [["\\Gamma", "X", "M", "\\Gamma", "Z", "P", "N", "Z_1", "M"], ["X", "P"]]
        return {"kpoints": kpoints, "path": path}

    def bctet2(self, c, a):
        """
        BCT2 Path
        """
        self.name = "BCT2"
        eta = (1 + a**2 / c**2) / 4.0
        zeta = a**2 / (2 * c**2)
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "N": np.array([0.0, 0.5, 0.0]),
            "P": np.array([0.25, 0.25, 0.25]),
            "\\Sigma": np.array([-eta, eta, eta]),
            "\\Sigma_1": np.array([eta, 1 - eta, -eta]),
            "X": np.array([0.0, 0.0, 0.5]),
            "Y": np.array([-zeta, zeta, 0.5]),
            "Y_1": np.array([0.5, 0.5, -zeta]),
            "Z": np.array([0.5, 0.5, -0.5]),
        }
        path = [
            [
                "\\Gamma",
                "X",
                "Y",
                "\\Sigma",
                "\\Gamma",
                "Z",
                "\\Sigma_1",
                "N",
                "P",
                "Y_1",
                "Z",
            ],
            ["X", "P"],
        ]
        return {"kpoints": kpoints, "path": path}

    def orc(self):
        """
        ORC Path
        """
        self.name = "ORC"
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "R": np.array([0.5, 0.5, 0.5]),
            "S": np.array([0.5, 0.5, 0.0]),
            "T": np.array([0.0, 0.5, 0.5]),
            "U": np.array([0.5, 0.0, 0.5]),
            "X": np.array([0.5, 0.0, 0.0]),
            "Y": np.array([0.0, 0.5, 0.0]),
            "Z": np.array([0.0, 0.0, 0.5]),
        }
        path = [
            ["\\Gamma", "X", "S", "Y", "\\Gamma", "Z", "U", "R", "T", "Z"],
            ["Y", "T"],
            ["U", "X"],
            ["S", "R"],
        ]
        return {"kpoints": kpoints, "path": path}

    def orcf1(self, a, b, c):
        """
        ORFC1 Path
        """
        self.name = "ORCF1"
        zeta = (1 + a**2 / b**2 - a**2 / c**2) / 4
        eta = (1 + a**2 / b**2 + a**2 / c**2) / 4

        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "A": np.array([0.5, 0.5 + zeta, zeta]),
            "A_1": np.array([0.5, 0.5 - zeta, 1 - zeta]),
            "L": np.array([0.5, 0.5, 0.5]),
            "T": np.array([1, 0.5, 0.5]),
            "X": np.array([0.0, eta, eta]),
            "X_1": np.array([1, 1 - eta, 1 - eta]),
            "Y": np.array([0.5, 0.0, 0.5]),
            "Z": np.array([0.5, 0.5, 0.0]),
        }
        path = [
            ["\\Gamma", "Y", "T", "Z", "\\Gamma", "X", "A_1", "Y"],
            ["T", "X_1"],
            ["X", "A", "Z"],
            ["L", "\\Gamma"],
        ]
        return {"kpoints": kpoints, "path": path}

    def orcf2(self, a, b, c):
        """
        ORFC2 Path
        """
        self.name = "ORCF2"
        phi = (1 + c**2 / b**2 - c**2 / a**2) / 4
        eta = (1 + a**2 / b**2 - a**2 / c**2) / 4
        delta = (1 + b**2 / a**2 - b**2 / c**2) / 4
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "C": np.array([0.5, 0.5 - eta, 1 - eta]),
            "C_1": np.array([0.5, 0.5 + eta, eta]),
            "D": np.array([0.5 - delta, 0.5, 1 - delta]),
            "D_1": np.array([0.5 + delta, 0.5, delta]),
            "L": np.array([0.5, 0.5, 0.5]),
            "H": np.array([1 - phi, 0.5 - phi, 0.5]),
            "H_1": np.array([phi, 0.5 + phi, 0.5]),
            "X": np.array([0.0, 0.5, 0.5]),
            "Y": np.array([0.5, 0.0, 0.5]),
            "Z": np.array([0.5, 0.5, 0.0]),
        }
        path = [
            ["\\Gamma", "Y", "C", "D", "X", "\\Gamma", "Z", "D_1", "H", "C"],
            ["C_1", "Z"],
            ["X", "H_1"],
            ["H", "Y"],
            ["L", "\\Gamma"],
        ]
        return {"kpoints": kpoints, "path": path}

    def orcf3(self, a, b, c):
        """
        ORFC3 Path
        """
        self.name = "ORCF3"
        zeta = (1 + a**2 / b**2 - a**2 / c**2) / 4
        eta = (1 + a**2 / b**2 + a**2 / c**2) / 4
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "A": np.array([0.5, 0.5 + zeta, zeta]),
            "A_1": np.array([0.5, 0.5 - zeta, 1 - zeta]),
            "L": np.array([0.5, 0.5, 0.5]),
            "T": np.array([1, 0.5, 0.5]),
            "X": np.array([0.0, eta, eta]),
            "X_1": np.array([1, 1 - eta, 1 - eta]),
            "Y": np.array([0.5, 0.0, 0.5]),
            "Z": np.array([0.5, 0.5, 0.0]),
        }
        path = [
            ["\\Gamma", "Y", "T", "Z", "\\Gamma", "X", "A_1", "Y"],
            ["X", "A", "Z"],
            ["L", "\\Gamma"],
        ]
        return {"kpoints": kpoints, "path": path}

    def orci(self, a, b, c):
        """
        ORCI Path
        """
        self.name = "ORCI"
        zeta = (1 + a**2 / c**2) / 4
        eta = (1 + b**2 / c**2) / 4
        delta = (b**2 - a**2) / (4 * c**2)
        mu = (a**2 + b**2) / (4 * c**2)
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "L": np.array([-mu, mu, 0.5 - delta]),
            "L_1": np.array([mu, -mu, 0.5 + delta]),
            "L_2": np.array([0.5 - delta, 0.5 + delta, -mu]),
            "R": np.array([0.0, 0.5, 0.0]),
            "S": np.array([0.5, 0.0, 0.0]),
            "T": np.array([0.0, 0.0, 0.5]),
            "W": np.array([0.25, 0.25, 0.25]),
            "X": np.array([-zeta, zeta, zeta]),
            "X_1": np.array([zeta, 1 - zeta, -zeta]),
            "Y": np.array([eta, -eta, eta]),
            "Y_1": np.array([1 - eta, eta, -eta]),
            "Z": np.array([0.5, 0.5, -0.5]),
        }
        path = [
            ["\\Gamma", "X", "L", "T", "W", "R", "X_1", "Z", "\\Gamma", "Y", "S", "W"],
            ["L_1", "Y"],
            ["Y_1", "Z"],
        ]
        return {"kpoints": kpoints, "path": path}

    def orcc(self, a, b, c):
        """
        ORCC Path
        """
        self.name = "ORCC"
        zeta = (1 + a**2 / b**2) / 4
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "A": np.array([zeta, zeta, 0.5]),
            "A_1": np.array([-zeta, 1 - zeta, 0.5]),
            "R": np.array([0.0, 0.5, 0.5]),
            "S": np.array([0.0, 0.5, 0.0]),
            "T": np.array([-0.5, 0.5, 0.5]),
            "X": np.array([zeta, zeta, 0.0]),
            "X_1": np.array([-zeta, 1 - zeta, 0.0]),
            "Y": np.array([-0.5, 0.5, 0]),
            "Z": np.array([0.0, 0.0, 0.5]),
        }
        path = [
            [
                "\\Gamma",
                "X",
                "S",
                "R",
                "A",
                "Z",
                "\\Gamma",
                "Y",
                "X_1",
                "A_1",
                "T",
                "Y",
            ],
            ["Z", "T"],
        ]
        return {"kpoints": kpoints, "path": path}

    def hex(self):
        """
        HEX Path
        """
        self.name = "HEX"
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "A": np.array([0.0, 0.0, 0.5]),
            "H": np.array([1.0 / 3.0, 1.0 / 3.0, 0.5]),
            "K": np.array([1.0 / 3.0, 1.0 / 3.0, 0.0]),
            "L": np.array([0.5, 0.0, 0.5]),
            "M": np.array([0.5, 0.0, 0.0]),
        }
        path = [
            ["\\Gamma", "M", "K", "\\Gamma", "A", "L", "H", "A"],
            ["L", "M"],
            ["K", "H"],
        ]
        return {"kpoints": kpoints, "path": path}

    def rhl1(self, alpha):
        """
        RHL1 Path
        """
        self.name = "RHL1"
        eta = (1 + 4 * cos(alpha)) / (2 + 4 * cos(alpha))
        nu = 3.0 / 4.0 - eta / 2.0
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "B": np.array([eta, 0.5, 1.0 - eta]),
            "B_1": np.array([1.0 / 2.0, 1.0 - eta, eta - 1.0]),
            "F": np.array([0.5, 0.5, 0.0]),
            "L": np.array([0.5, 0.0, 0.0]),
            "L_1": np.array([0.0, 0.0, -0.5]),
            "P": np.array([eta, nu, nu]),
            "P_1": np.array([1.0 - nu, 1.0 - nu, 1.0 - eta]),
            "P_2": np.array([nu, nu, eta - 1.0]),
            "Q": np.array([1.0 - nu, nu, 0.0]),
            "X": np.array([nu, 0.0, -nu]),
            "Z": np.array([0.5, 0.5, 0.5]),
        }
        path = [
            ["\\Gamma", "L", "B_1"],
            ["B", "Z", "\\Gamma", "X"],
            ["Q", "F", "P_1", "Z"],
            ["L", "P"],
        ]
        return {"kpoints": kpoints, "path": path}

    def rhl2(self, alpha):
        """
        RHL2 Path
        """
        self.name = "RHL2"
        eta = 1 / (2 * tan(alpha / 2.0) ** 2)
        nu = 3.0 / 4.0 - eta / 2.0
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "F": np.array([0.5, -0.5, 0.0]),
            "L": np.array([0.5, 0.0, 0.0]),
            "P": np.array([1 - nu, -nu, 1 - nu]),
            "P_1": np.array([nu, nu - 1.0, nu - 1.0]),
            "Q": np.array([eta, eta, eta]),
            "Q_1": np.array([1.0 - eta, -eta, -eta]),
            "Z": np.array([0.5, -0.5, 0.5]),
        }
        path = [["\\Gamma", "P", "Z", "Q", "\\Gamma", "F", "P_1", "Q_1", "L", "Z"]]
        return {"kpoints": kpoints, "path": path}

    def mcl(self, b, c, beta):
        """
        MCL Path
        """
        self.name = "MCL"
        eta = (1 - b * cos(beta) / c) / (2 * sin(beta) ** 2)
        nu = 0.5 - eta * c * cos(beta) / b
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "A": np.array([0.5, 0.5, 0.0]),
            "C": np.array([0.0, 0.5, 0.5]),
            "D": np.array([0.5, 0.0, 0.5]),
            "D_1": np.array([0.5, 0.5, -0.5]),
            "E": np.array([0.5, 0.5, 0.5]),
            "H": np.array([0.0, eta, 1.0 - nu]),
            "H_1": np.array([0.0, 1.0 - eta, nu]),
            "H_2": np.array([0.0, eta, -nu]),
            "M": np.array([0.5, eta, 1.0 - nu]),
            "M_1": np.array([0.5, 1 - eta, nu]),
            "M_2": np.array([0.5, 1 - eta, nu]),
            "X": np.array([0.0, 0.5, 0.0]),
            "Y": np.array([0.0, 0.0, 0.5]),
            "Y_1": np.array([0.0, 0.0, -0.5]),
            "Z": np.array([0.5, 0.0, 0.0]),
        }
        path = [
            ["\\Gamma", "Y", "H", "C", "E", "M_1", "A", "X", "H_1"],
            ["M", "D", "Z"],
            ["Y", "D"],
        ]
        return {"kpoints": kpoints, "path": path}

    def mclc1(self, a, b, c, alpha):
        """
        MCLC1 Path
        """
        self.name = "MCLC1"
        zeta = (2 - b * cos(alpha) / c) / (4 * sin(alpha) ** 2)
        eta = 0.5 + 2 * zeta * c * cos(alpha) / b
        psi = 0.75 - a**2 / (4 * b**2 * sin(alpha) ** 2)
        phi = psi + (0.75 - psi) * b * cos(alpha) / c
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "N": np.array([0.5, 0.0, 0.0]),
            "N_1": np.array([0.0, -0.5, 0.0]),
            "F": np.array([1 - zeta, 1 - zeta, 1 - eta]),
            "F_1": np.array([zeta, zeta, eta]),
            "F_2": np.array([-zeta, -zeta, 1 - eta]),
            "I": np.array([phi, 1 - phi, 0.5]),
            "I_1": np.array([1 - phi, phi - 1, 0.5]),
            "L": np.array([0.5, 0.5, 0.5]),
            "M": np.array([0.5, 0.0, 0.5]),
            "X": np.array([1 - psi, psi - 1, 0.0]),
            "X_1": np.array([psi, 1 - psi, 0.0]),
            "X_2": np.array([psi - 1, -psi, 0.0]),
            "Y": np.array([0.5, 0.5, 0.0]),
            "Y_1": np.array([-0.5, -0.5, 0.0]),
            "Z": np.array([0.0, 0.0, 0.5]),
        }
        path = [
            ["\\Gamma", "Y", "F", "L", "I"],
            ["I_1", "Z", "F_1"],
            ["Y", "X_1"],
            ["X", "\\Gamma", "N"],
            ["M", "\\Gamma"],
        ]
        return {"kpoints": kpoints, "path": path}

    def mclc2(self, a, b, c, alpha):
        """
        MCLC2 Path
        """
        self.name = "MCLC2"
        zeta = (2 - b * cos(alpha) / c) / (4 * sin(alpha) ** 2)
        eta = 0.5 + 2 * zeta * c * cos(alpha) / b
        psi = 0.75 - a**2 / (4 * b**2 * sin(alpha) ** 2)
        phi = psi + (0.75 - psi) * b * cos(alpha) / c
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "N": np.array([0.5, 0.0, 0.0]),
            "N_1": np.array([0.0, -0.5, 0.0]),
            "F": np.array([1 - zeta, 1 - zeta, 1 - eta]),
            "F_1": np.array([zeta, zeta, eta]),
            "F_2": np.array([-zeta, -zeta, 1 - eta]),
            "F_3": np.array([1 - zeta, -zeta, 1 - eta]),
            "I": np.array([phi, 1 - phi, 0.5]),
            "I_1": np.array([1 - phi, phi - 1, 0.5]),
            "L": np.array([0.5, 0.5, 0.5]),
            "M": np.array([0.5, 0.0, 0.5]),
            "X": np.array([1 - psi, psi - 1, 0.0]),
            "X_1": np.array([psi, 1 - psi, 0.0]),
            "X_2": np.array([psi - 1, -psi, 0.0]),
            "Y": np.array([0.5, 0.5, 0.0]),
            "Y_1": np.array([-0.5, -0.5, 0.0]),
            "Z": np.array([0.0, 0.0, 0.5]),
        }
        path = [
            ["\\Gamma", "Y", "F", "L", "I"],
            ["I_1", "Z", "F_1"],
            ["N", "\\Gamma", "M"],
        ]
        return {"kpoints": kpoints, "path": path}

    def mclc3(self, a, b, c, alpha):
        """
        MCLC3 Path
        """
        self.name = "MCLC3"
        mu = (1 + b**2 / a**2) / 4.0
        delta = b * c * cos(alpha) / (2 * a**2)
        zeta = mu - 0.25 + (1 - b * cos(alpha) / c) / (4 * sin(alpha) ** 2)
        eta = 0.5 + 2 * zeta * c * cos(alpha) / b
        phi = 1 + zeta - 2 * mu
        psi = eta - 2 * delta
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "F": np.array([1 - phi, 1 - phi, 1 - psi]),
            "F_1": np.array([phi, phi - 1, psi]),
            "F_2": np.array([1 - phi, -phi, 1 - psi]),
            "H": np.array([zeta, zeta, eta]),
            "H_1": np.array([1 - zeta, -zeta, 1 - eta]),
            "H_2": np.array([-zeta, -zeta, 1 - eta]),
            "I": np.array([0.5, -0.5, 0.5]),
            "M": np.array([0.5, 0.0, 0.5]),
            "N": np.array([0.5, 0.0, 0.0]),
            "N_1": np.array([0.0, -0.5, 0.0]),
            "X": np.array([0.5, -0.5, 0.0]),
            "Y": np.array([mu, mu, delta]),
            "Y_1": np.array([1 - mu, -mu, -delta]),
            "Y_2": np.array([-mu, -mu, -delta]),
            "Y_3": np.array([mu, mu - 1, delta]),
            "Z": np.array([0.0, 0.0, 0.5]),
        }
        path = [
            ["\\Gamma", "Y", "F", "H", "Z", "I", "F_1"],
            ["H_1", "Y_1", "X", "\\Gamma", "N"],
            ["M", "\\Gamma"],
        ]
        return {"kpoints": kpoints, "path": path}

    def mclc4(self, a, b, c, alpha):
        """
        MCLC4 Path
        """
        self.name = "MCLC4"
        mu = (1 + b**2 / a**2) / 4.0
        delta = b * c * cos(alpha) / (2 * a**2)
        zeta = mu - 0.25 + (1 - b * cos(alpha) / c) / (4 * sin(alpha) ** 2)
        eta = 0.5 + 2 * zeta * c * cos(alpha) / b
        phi = 1 + zeta - 2 * mu
        psi = eta - 2 * delta
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "F": np.array([1 - phi, 1 - phi, 1 - psi]),
            "F_1": np.array([phi, phi - 1, psi]),
            "F_2": np.array([1 - phi, -phi, 1 - psi]),
            "H": np.array([zeta, zeta, eta]),
            "H_1": np.array([1 - zeta, -zeta, 1 - eta]),
            "H_2": np.array([-zeta, -zeta, 1 - eta]),
            "I": np.array([0.5, -0.5, 0.5]),
            "M": np.array([0.5, 0.0, 0.5]),
            "N": np.array([0.5, 0.0, 0.0]),
            "N_1": np.array([0.0, -0.5, 0.0]),
            "X": np.array([0.5, -0.5, 0.0]),
            "Y": np.array([mu, mu, delta]),
            "Y_1": np.array([1 - mu, -mu, -delta]),
            "Y_2": np.array([-mu, -mu, -delta]),
            "Y_3": np.array([mu, mu - 1, delta]),
            "Z": np.array([0.0, 0.0, 0.5]),
        }
        path = [
            ["\\Gamma", "Y", "F", "H", "Z", "I"],
            ["H_1", "Y_1", "X", "\\Gamma", "N"],
            ["M", "\\Gamma"],
        ]
        return {"kpoints": kpoints, "path": path}

    def mclc5(self, a, b, c, alpha):
        """
        MCLC5 Path
        """
        self.name = "MCLC5"
        zeta = (b**2 / a**2 + (1 - b * cos(alpha) / c) / sin(alpha) ** 2) / 4
        eta = 0.5 + 2 * zeta * c * cos(alpha) / b
        mu = eta / 2 + b**2 / (4 * a**2) - b * c * cos(alpha) / (2 * a**2)
        nu = 2 * mu - zeta
        rho = 1 - zeta * a**2 / b**2
        omega = (4 * nu - 1 - b**2 * sin(alpha) ** 2 / a**2) * c / (2 * b * cos(alpha))
        delta = zeta * c * cos(alpha) / b + omega / 2 - 0.25
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "F": np.array([nu, nu, omega]),
            "F_1": np.array([1 - nu, 1 - nu, 1 - omega]),
            "F_2": np.array([nu, nu - 1, omega]),
            "H": np.array([zeta, zeta, eta]),
            "H_1": np.array([1 - zeta, -zeta, 1 - eta]),
            "H_2": np.array([-zeta, -zeta, 1 - eta]),
            "I": np.array([rho, 1 - rho, 0.5]),
            "I_1": np.array([1 - rho, rho - 1, 0.5]),
            "L": np.array([0.5, 0.5, 0.5]),
            "M": np.array([0.5, 0.0, 0.5]),
            "N": np.array([0.5, 0.0, 0.0]),
            "N_1": np.array([0.0, -0.5, 0.0]),
            "X": np.array([0.5, -0.5, 0.0]),
            "Y": np.array([mu, mu, delta]),
            "Y_1": np.array([1 - mu, -mu, -delta]),
            "Y_2": np.array([-mu, -mu, -delta]),
            "Y_3": np.array([mu, mu - 1, delta]),
            "Z": np.array([0.0, 0.0, 0.5]),
        }
        path = [
            ["\\Gamma", "Y", "F", "L", "I"],
            ["I_1", "Z", "H", "F_1"],
            ["H_1", "Y_1", "X", "\\Gamma", "N"],
            ["M", "\\Gamma"],
        ]
        return {"kpoints": kpoints, "path": path}

    def tria(self):
        """
        TRI1a Path
        """
        self.name = "TRI1a"
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "L": np.array([0.5, 0.5, 0.0]),
            "M": np.array([0.0, 0.5, 0.5]),
            "N": np.array([0.5, 0.0, 0.5]),
            "R": np.array([0.5, 0.5, 0.5]),
            "X": np.array([0.5, 0.0, 0.0]),
            "Y": np.array([0.0, 0.5, 0.0]),
            "Z": np.array([0.0, 0.0, 0.5]),
        }
        path = [
            ["X", "\\Gamma", "Y"],
            ["L", "\\Gamma", "Z"],
            ["N", "\\Gamma", "M"],
            ["R", "\\Gamma"],
        ]
        return {"kpoints": kpoints, "path": path}

    def trib(self):
        """
        TRI1b Path
        """
        self.name = "TRI1b"
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "L": np.array([0.5, -0.5, 0.0]),
            "M": np.array([0.0, 0.0, 0.5]),
            "N": np.array([-0.5, -0.5, 0.5]),
            "R": np.array([0.0, -0.5, 0.5]),
            "X": np.array([0.0, -0.5, 0.0]),
            "Y": np.array([0.5, 0.0, 0.0]),
            "Z": np.array([-0.5, 0.0, 0.5]),
        }
        path = [
            ["X", "\\Gamma", "Y"],
            ["L", "\\Gamma", "Z"],
            ["N", "\\Gamma", "M"],
            ["R", "\\Gamma"],
        ]
        return {"kpoints": kpoints, "path": path}


class KPathSeek(KPathBase):
    """
    This class looks for a path along high-symmetry lines in the Brillouin zone. It is based on
    Hinuma, Y., Pizzi, G., Kumagai, Y., Oba, F., & Tanaka, I. (2017). Band structure diagram paths
    based on crystallography. Computational Materials Science, 128, 140-184.
    https://doi.org/10.1016/j.commatsci.2016.10.015. It should be used with primitive structures that
    comply with the definition given in the paper. The symmetry is determined by spglib using the
    SpacegroupAnalyzer class. k-points are generated using the get_kpoints() method for the
    reciprocal cell basis defined in the paper.
    """

    @requires(
        get_path is not None,
        "SeeK-path needs to be installed to use the convention of Hinuma et al. (2015)",
    )
    def __init__(self, structure: Structure, symprec: float = 0.01, angle_tolerance=5, atol=1e-5, system_is_tri=True):
        """
        Args:
            structure (Structure): Structure object
            symprec (float): Tolerance for symmetry finding
            angle_tolerance (float): Angle tolerance for symmetry finding.
            atol (float): Absolute tolerance used to determine edge cases
                for settings of structures.
            system_is_tri (bool): Indicates if the system is time-reversal
                invariant.
        """
        super().__init__(structure, symprec=symprec, angle_tolerance=angle_tolerance, atol=atol)

        positions = structure.frac_coords

        sp = structure.site_properties
        species = [site.species for site in structure]
        site_data = species

        if not system_is_tri:
            warn("Non-zero 'magmom' data will be used to define unique atoms in the cell.")
            site_data = zip(species, [tuple(vec) for vec in sp["magmom"]])  # type: ignore

        unique_species: list[SpeciesLike] = []
        numbers = []

        for species, g in itertools.groupby(site_data):
            if species in unique_species:
                ind = unique_species.index(species)
                numbers.extend([ind + 1] * len(tuple(g)))
            else:
                unique_species.append(species)
                numbers.extend([len(unique_species)] * len(tuple(g)))

        cell = (self._latt.matrix, positions, numbers)

        lattice, scale_pos, atom_num = spglib.standardize_cell(
            cell, to_primitive=False, no_idealize=True, symprec=symprec
        )

        spg_struct = (lattice, scale_pos, atom_num)
        spath_dat = get_path(spg_struct, system_is_tri, "hpkot", atol, symprec, angle_tolerance)

        self._tmat = self._trans_sc_to_Hin(spath_dat["bravais_lattice_extended"])
        self._rec_lattice = Lattice(spath_dat["reciprocal_primitive_lattice"])

        spath_data_formatted = [[spath_dat["path"][0][0]]]
        count = 0
        for pnum in range(len(spath_dat["path"]) - 1):
            if spath_dat["path"][pnum][1] == spath_dat["path"][pnum + 1][0]:
                spath_data_formatted[count].append(spath_dat["path"][pnum][1])
            else:
                spath_data_formatted[count].append(spath_dat["path"][pnum][1])
                spath_data_formatted.append([])
                count += 1
                spath_data_formatted[count].append(spath_dat["path"][pnum + 1][0])

        spath_data_formatted[-1].append(spath_dat["path"][-1][1])

        self._kpath = {
            "kpoints": spath_dat["point_coords"],
            "path": spath_data_formatted,
        }

    @staticmethod
    def _trans_sc_to_Hin(sub_class):

        if sub_class in [
            "cP1",
            "cP2",
            "cF1",
            "cF2",
            "cI1",
            "tP1",
            "oP1",
            "hP1",
            "hP2",
            "tI1",
            "tI2",
            "oF1",
            "oF3",
            "oI1",
            "oI3",
            "oC1",
            "hR1",
            "hR2",
            "aP1",
            "aP2",
            "aP3",
            "oA1",
        ]:
            return np.eye(3)
        if sub_class == "oF2":
            return np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
        if sub_class == "oI2":
            return np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
        if sub_class == "oI3":
            return np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
        if sub_class == "oA2":
            return np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])
        if sub_class == "oC2":
            return np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])
        if sub_class in ["mP1", "mC1", "mC2", "mC3"]:
            return np.array([[0, 1, 0], [-1, 0, 0], [0, 0, 1]])
        raise RuntimeError("Sub-classification of crystal not found!")


class KPathLatimerMunro(KPathBase):
    """
    This class looks for a path along high-symmetry lines in the
    Brillouin zone. It is based on the method outlined in:
    npj Comput Mater 6, 112 (2020). 10.1038/s41524-020-00383-7
    The user should ensure that the unit cell of the input structure
    is as reduced as possible, i.e. that there is no linear
    combination of lattice vectors which can produce a vector of
    lesser magnitude than the given set (this is required to
    obtain the correct Brillouin zone within the current
    implementation). This is checked during initialization and a
    warning is issued if the condition is not fulfilled.
    In the case of magnetic structures, care must also be taken to
    provide the magnetic primitive cell (i.e. that which reproduces
    the entire crystal, including the correct magnetic ordering,
    upon application of lattice translations). There is no algorithm to
        check for this, so if the input structure is
    incorrect, the class will output the incorrect k-path without
    any warning being issued.
    """

    def __init__(
        self,
        structure,
        has_magmoms=False,
        magmom_axis=None,
        symprec=0.01,
        angle_tolerance=5,
        atol=1e-5,
    ):
        """
        Args:
            structure (Structure): Structure object
            has_magmoms (bool): Whether the input structure contains
                magnetic moments as site properties with the key 'magmom.'
                Values may be in the form of 3-component vectors given in
                the basis of the input lattice vectors, or as scalars, in
                which case the spin axis will default to a_3, the third
                real-space lattice vector (this triggers a warning).
            magmom_axis (list or numpy array): 3-component vector specifying
                direction along which magnetic moments given as scalars
                should point. If all magnetic moments are provided as
                vectors then this argument is not used.
            symprec (float): Tolerance for symmetry finding
            angle_tolerance (float): Angle tolerance for symmetry finding.
            atol (float): Absolute tolerance used to determine symmetric
                equivalence of points and lines in the BZ.
        """
        super().__init__(structure, symprec=symprec, angle_tolerance=angle_tolerance, atol=atol)

        # Check to see if input cell is reducible. Ref: B Gruber in Acta. Cryst. Vol. A29,
        # pp. 433-440 ('The Relationship between Reduced Cells in a General Bravais lattice').
        # The correct BZ will still be obtained if the lattice vectors are reducible by any
        # linear combination of themselves with coefficients of absolute value less than 2,
        # hence a missing factor of 2 as compared to the reference.
        reducible = []
        for i in range(3):
            for j in range(3):
                if i != j:
                    if (
                        np.absolute(np.dot(self._latt.matrix[i], self._latt.matrix[j]))
                        > np.dot(self._latt.matrix[i], self._latt.matrix[i])
                        and np.absolute(
                            np.dot(self._latt.matrix[i], self._latt.matrix[j])
                            - np.dot(self._latt.matrix[i], self._latt.matrix[i])
                        )
                        > atol
                    ):
                        reducible.append(True)
                    else:
                        reducible.append(False)
        if np.any(reducible):
            print("reducible")
            warn(
                "The unit cell of the input structure is not fully reduced!"
                "The path may be incorrect. Use at your own risk."
            )

        if magmom_axis is None:
            magmom_axis = np.array([0, 0, 1])
            axis_specified = False
        else:
            axis_specified = True

        self._kpath = self._get_ksymm_kpath(has_magmoms, magmom_axis, axis_specified, symprec, angle_tolerance, atol)

    @property
    def mag_type(self):
        """
        Returns:
        The type of magnetic space group as a string.
        Current implementation does not distinguish
        between types 3 and 4, so return value is '3/4'.
        If has_magmoms is False, returns '0'.
        """
        return self._mag_type

    def _get_ksymm_kpath(self, has_magmoms, magmom_axis, axis_specified, symprec, angle_tolerance, atol):
        ID = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        # parity, aka the inversion operation (not calling it
        PAR = np.array([[-1, 0, 0], [0, -1, 0], [0, 0, -1]])
        # INV to avoid confusion with np.linalg.inv() function)

        # 1: Get lattices of real and reciprocal structures, and reciprocal
        # point group, and Brillouin zone (BZ)

        V = self._latt.matrix.T  # fractional real space to Cartesian real space
        # fractional reciprocal space to Cartesian reciprocal space
        W = self._rec_lattice.matrix.T
        # fractional real space to fractional reciprocal space
        A = np.dot(np.linalg.inv(W), V)

        if has_magmoms:
            grey_struct = self._structure.copy()
            grey_struct.remove_site_property("magmom")
            sga = SpacegroupAnalyzer(grey_struct, symprec=symprec, angle_tolerance=angle_tolerance)
            grey_ops = sga.get_symmetry_operations()
            self._structure = self._convert_all_magmoms_to_vectors(magmom_axis, axis_specified)
            mag_ops = self._get_magnetic_symmetry_operations(self._structure, grey_ops, atol)

            D = [
                SymmOp.from_rotation_and_translation(
                    rotation_matrix=op.rotation_matrix,
                    translation_vec=op.translation_vector,
                )
                for op in mag_ops
                if op.time_reversal == 1
            ]

            fD = [
                SymmOp.from_rotation_and_translation(
                    rotation_matrix=op.rotation_matrix,
                    translation_vec=op.translation_vector,
                )
                for op in mag_ops
                if op.time_reversal == -1
            ]

            if np.array([m == np.array([0, 0, 0]) for m in self._structure.site_properties["magmom"]]).all():
                fD = D
                D = []

            if len(fD) == 0:  # no operations contain time reversal; type 1
                self._mag_type = "1"
                isomorphic_point_group = [d.rotation_matrix for d in D]
                recip_point_group = self._get_reciprocal_point_group(isomorphic_point_group, ID, A)
            elif len(D) == 0:  # all operations contain time reversal / all magmoms zero; type 2
                self._mag_type = "2"
                isomorphic_point_group = [d.rotation_matrix for d in fD]
                recip_point_group = self._get_reciprocal_point_group(isomorphic_point_group, PAR, A)
            else:  # half and half; type 3 or 4
                self._mag_type = "3/4"
                f = self._get_coset_factor(D + fD, D)
                isomorphic_point_group = [d.rotation_matrix for d in D]
                recip_point_group = self._get_reciprocal_point_group(
                    isomorphic_point_group, np.dot(PAR, f.rotation_matrix), A
                )
        else:
            self._mag_type = "0"
            if "magmom" in self._structure.site_properties:
                warn(
                    "The parameter has_magmoms is False, but site_properties contains the key magmom."
                    "This property will be removed and could result in different symmetry operations."
                )
                self._structure.remove_site_property("magmom")
            sga = SpacegroupAnalyzer(self._structure)
            ops = sga.get_symmetry_operations()
            isomorphic_point_group = [op.rotation_matrix for op in ops]

            recip_point_group = self._get_reciprocal_point_group(isomorphic_point_group, PAR, A)

        self._rpg = recip_point_group

        # 2: Get all vertices, edge- and face- center points of BZ ("key points")

        key_points, bz_as_key_point_inds, face_center_inds = self._get_key_points()

        # 3: Find symmetry-equivalent points, which can be mapped to each other by a combination of point group
        # operations and integer translations by lattice vectors. The integers will only be -1, 0, or 1, since
        # we are restricted to the BZ.

        key_points_inds_orbits = self._get_key_point_orbits(key_points=key_points)

        # 4: Get all lines in BZ between adjacent key points and between gamma
        # and key points ("key lines")

        key_lines = self._get_key_lines(key_points=key_points, bz_as_key_point_inds=bz_as_key_point_inds)

        # 5: Find symmetry-equivalent key lines, defined as end points of first line being equivalent
        # to end points of second line, and a random point in between being equivalent to the mapped
        # random point.

        key_lines_inds_orbits = self._get_key_line_orbits(
            key_points=key_points,
            key_lines=key_lines,
            key_points_inds_orbits=key_points_inds_orbits,
        )

        # 6 & 7: Get little groups for key points (group of symmetry elements present at that point).
        # Get little groups for key lines (group of symmetry elements present at every point
        # along the line). This is implemented by testing the symmetry at a point e/pi of the
        # way between the two endpoints.

        little_groups_points, little_groups_lines = self._get_little_groups(
            key_points=key_points,
            key_points_inds_orbits=key_points_inds_orbits,
            key_lines_inds_orbits=key_lines_inds_orbits,
        )

        # 8: Choose key lines for k-path. Loose criteria set: choose any points / segments
        # with spatial symmetry greater than the general position (AKA more symmetry operations
        # than just the identity or identity * TR in the little group).
        # This function can be edited to alter high-symmetry criteria for choosing points and lines

        point_orbits_in_path, line_orbits_in_path = self._choose_path(
            key_points=key_points,
            key_points_inds_orbits=key_points_inds_orbits,
            key_lines_inds_orbits=key_lines_inds_orbits,
            little_groups_points=little_groups_points,
            little_groups_lines=little_groups_lines,
        )

        # 10: Consolidate selected segments into a single irreducible section of the BZ (as determined
        # by the reciprocal point and lattice symmetries). This is accomplished by identifying the boundary
        # planes of the IRBZ. Also, get labels for points according to distance from axes.

        IRBZ_points_inds = self._get_IRBZ(recip_point_group, W, key_points, face_center_inds, atol)
        lines_in_path_inds = []
        for ind in line_orbits_in_path:
            for tup in key_lines_inds_orbits[ind]:
                if tup[0] in IRBZ_points_inds and tup[1] in IRBZ_points_inds:
                    lines_in_path_inds.append(tup)
                    break
        G = nx.Graph(lines_in_path_inds)
        lines_in_path_inds = list(nx.edge_dfs(G))
        points_in_path_inds = [ind for tup in lines_in_path_inds for ind in tup]
        points_in_path_inds_unique = list(set(points_in_path_inds))

        orbit_cosines = []
        for orbit in key_points_inds_orbits[:-1]:

            orbit_cosines.append(
                sorted(
                    sorted(
                        (
                            (
                                j,
                                np.round(
                                    np.dot(key_points[k], self.LabelPoints(j))
                                    / (np.linalg.norm(key_points[k]) * np.linalg.norm(self.LabelPoints(j))),
                                    decimals=3,
                                ),
                            )
                            for k in orbit
                            for j in range(26)
                        ),
                        key=operator.itemgetter(0),
                    ),
                    key=operator.itemgetter(1),
                    reverse=True,
                )
            )

        orbit_labels = self._get_orbit_labels(orbit_cosines, key_points_inds_orbits, atol)
        key_points_labels = ["" for i in range(len(key_points))]
        for i, orbit in enumerate(key_points_inds_orbits):
            for point_ind in orbit:
                key_points_labels[point_ind] = self.LabelSymbol(int(orbit_labels[i]))

        kpoints = {}
        reverse_kpoints = {}
        for point_ind in points_in_path_inds_unique:
            point_label = key_points_labels[point_ind]
            if point_label not in kpoints:
                kpoints[point_label] = key_points[point_ind]
                reverse_kpoints[point_ind] = point_label
            else:
                existing_labels = [key for key in kpoints if point_label in key]
                if "'" not in point_label:
                    existing_labels[:] = [label for label in existing_labels if "'" not in label]
                if len(existing_labels) == 1:
                    max_occurence = 0
                else:
                    if "'" not in point_label:
                        max_occurence = max(int(label[3:-1]) for label in existing_labels[1:])
                    else:
                        max_occurence = max(int(label[4:-1]) for label in existing_labels[1:])
                kpoints[point_label + "_{" + str(max_occurence + 1) + "}"] = key_points[point_ind]
                reverse_kpoints[point_ind] = point_label + "_{" + str(max_occurence + 1) + "}"

        path = []
        i = 0
        start_of_subpath = True
        while i < len(points_in_path_inds):
            if start_of_subpath:
                path.append([reverse_kpoints[points_in_path_inds[i]]])
                i += 1
                start_of_subpath = False
            elif points_in_path_inds[i] == points_in_path_inds[i + 1]:
                path[-1].append(reverse_kpoints[points_in_path_inds[i]])
                i += 2
            else:
                path[-1].append(reverse_kpoints[points_in_path_inds[i]])
                i += 1
                start_of_subpath = True
            if i == len(points_in_path_inds) - 1:
                path[-1].append(reverse_kpoints[points_in_path_inds[i]])
                i += 1

        return {"kpoints": kpoints, "path": path}

    def _choose_path(
        self,
        key_points,
        key_points_inds_orbits,
        key_lines_inds_orbits,
        little_groups_points,
        little_groups_lines,
    ):
        #
        # This function can be edited to alter high-symmetry criteria for choosing points and lines
        #

        ID = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        PAR = np.array([[-1, 0, 0], [0, -1, 0], [0, 0, -1]])

        gamma_ind = len(key_points) - 1

        line_orbits_in_path = []
        point_orbits_in_path = []
        for (i, little_group) in enumerate(little_groups_lines):
            add_rep = False
            nC2 = 0
            nC3 = 0
            nsig = 0
            for opind in little_group:
                op = self._rpg[opind]
                if not (op == ID).all():
                    if (np.dot(op, op) == ID).all():
                        if np.linalg.det(op) == 1:
                            nC2 += 1
                            break
                        if not (op == PAR).all():
                            nsig += 1
                            break
                    elif (np.dot(op, np.dot(op, op)) == ID).all():
                        nC3 += 1
                        break
            if nC2 > 0 or nC3 > 0 or nsig > 0:
                add_rep = True

            if add_rep:
                line_orbits_in_path.append(i)
                l = key_lines_inds_orbits[i][0]
                ind0 = l[0]
                ind1 = l[1]
                found0 = False
                found1 = False
                for (j, orbit) in enumerate(key_points_inds_orbits):
                    if ind0 in orbit:
                        point_orbits_in_path.append(j)
                        found0 = True
                    if ind1 in orbit:
                        point_orbits_in_path.append(j)
                        found1 = True
                    if found0 and found1:
                        break
        point_orbits_in_path = list(set(point_orbits_in_path))

        # Choose remaining unconnected key points for k-path. The ones that remain are
        # those with inversion symmetry. Connect them to gamma.

        unconnected = []

        for i in range(len(key_points_inds_orbits)):
            if i not in point_orbits_in_path:
                unconnected.append(i)

        for ind in unconnected:
            connect = False
            for op_ind in little_groups_points[ind]:
                op = self._rpg[op_ind]
                if (op == ID).all():
                    pass
                elif (op == PAR).all():
                    connect = True
                    break
                elif np.linalg.det(op) == 1:
                    if (np.dot(op, np.dot(op, op)) == ID).all():
                        pass
                    else:
                        connect = True
                        break
                else:
                    pass
            if connect:
                l = (key_points_inds_orbits[ind][0], gamma_ind)
                for (j, orbit) in enumerate(key_lines_inds_orbits):
                    if l in orbit:
                        line_orbits_in_path.append(j)
                        break
                if gamma_ind not in point_orbits_in_path:
                    point_orbits_in_path.append(gamma_ind)
                point_orbits_in_path.append(ind)

        return point_orbits_in_path, line_orbits_in_path

    def _get_key_points(self):
        decimals = ceil(-1 * np.log10(self._atol)) - 1
        bz = self._rec_lattice.get_wigner_seitz_cell()

        key_points = []
        face_center_inds = []
        bz_as_key_point_inds = []

        # pymatgen gives BZ in Cartesian coordinates; convert to fractional in
        # the primitive basis for reciprocal space
        for (i, facet) in enumerate(bz):
            for (j, vert) in enumerate(facet):
                vert = self._rec_lattice.get_fractional_coords(vert)
                bz[i][j] = vert
        pop = []
        for i, facet in enumerate(bz):
            rounded_facet = np.around(facet, decimals=decimals)
            u, indices = np.unique(rounded_facet, axis=0, return_index=True)
            if len(u) in [1, 2]:
                pop.append(i)
            else:
                bz[i] = [facet[j] for j in np.sort(indices)]
        bz = [bz[i] for i in range(len(bz)) if i not in pop]

        # use vertex points to calculate edge- and face- centers
        for (i, facet) in enumerate(bz):
            bz_as_key_point_inds.append([])
            for (j, vert) in enumerate(facet):
                edge_center = (vert + facet[j + 1]) / 2.0 if j != len(facet) - 1 else (vert + facet[0]) / 2.0
                duplicatevert = False
                duplicateedge = False
                for (k, point) in enumerate(key_points):
                    if np.allclose(vert, point, atol=self._atol):
                        bz_as_key_point_inds[i].append(k)
                        duplicatevert = True
                        break
                for (k, point) in enumerate(key_points):
                    if np.allclose(edge_center, point, atol=self._atol):
                        bz_as_key_point_inds[i].append(k)
                        duplicateedge = True
                        break
                if not duplicatevert:
                    key_points.append(vert)
                    bz_as_key_point_inds[i].append(len(key_points) - 1)
                if not duplicateedge:
                    key_points.append(edge_center)
                    bz_as_key_point_inds[i].append(len(key_points) - 1)
            if len(facet) == 4:  # parallelogram facet
                face_center = (facet[0] + facet[1] + facet[2] + facet[3]) / 4.0
                key_points.append(face_center)
                face_center_inds.append(len(key_points) - 1)
                bz_as_key_point_inds[i].append(len(key_points) - 1)
            else:  # hexagonal facet
                face_center = (facet[0] + facet[1] + facet[2] + facet[3] + facet[4] + facet[5]) / 6.0
                key_points.append(face_center)
                face_center_inds.append(len(key_points) - 1)
                bz_as_key_point_inds[i].append(len(key_points) - 1)

        # add gamma point
        key_points.append(np.array([0, 0, 0]))
        return key_points, bz_as_key_point_inds, face_center_inds

    def _get_key_point_orbits(self, key_points):
        key_points_copy = dict(zip(range(len(key_points) - 1), key_points[0 : len(key_points) - 1]))
        # gamma not equivalent to any in BZ and is last point added to
        # key_points
        key_points_inds_orbits = []

        i = 0
        while len(key_points_copy) > 0:
            key_points_inds_orbits.append([])
            k0ind = list(key_points_copy)[0]
            k0 = key_points_copy[k0ind]
            key_points_inds_orbits[i].append(k0ind)
            key_points_copy.pop(k0ind)

            for op in self._rpg:
                to_pop = []
                k1 = np.dot(op, k0)
                for ind_key in key_points_copy:
                    diff = k1 - key_points_copy[ind_key]
                    if self._all_ints(diff, atol=self._atol):
                        key_points_inds_orbits[i].append(ind_key)
                        to_pop.append(ind_key)

                for key in to_pop:
                    key_points_copy.pop(key)
            i += 1

        key_points_inds_orbits.append([len(key_points) - 1])

        return key_points_inds_orbits

    @staticmethod
    def _get_key_lines(key_points, bz_as_key_point_inds):
        key_lines = []
        gamma_ind = len(key_points) - 1

        for facet_as_key_point_inds in bz_as_key_point_inds:
            facet_as_key_point_inds_bndy = facet_as_key_point_inds[: len(facet_as_key_point_inds) - 1]
            # not the face center point (don't need to check it since it's not
            # shared with other facets)
            face_center_ind = facet_as_key_point_inds[-1]
            for (j, ind) in enumerate(facet_as_key_point_inds_bndy):
                if (
                    min(ind, facet_as_key_point_inds_bndy[j - 1]),
                    max(ind, facet_as_key_point_inds_bndy[j - 1]),
                ) not in key_lines:
                    key_lines.append(
                        (
                            min(ind, facet_as_key_point_inds_bndy[j - 1]),
                            max(ind, facet_as_key_point_inds_bndy[j - 1]),
                        )
                    )
                k = j + 1 if j != len(facet_as_key_point_inds_bndy) - 1 else 0
                if (
                    min(ind, facet_as_key_point_inds_bndy[k]),
                    max(ind, facet_as_key_point_inds_bndy[k]),
                ) not in key_lines:
                    key_lines.append(
                        (
                            min(ind, facet_as_key_point_inds_bndy[k]),
                            max(ind, facet_as_key_point_inds_bndy[k]),
                        )
                    )
                if (ind, gamma_ind) not in key_lines:
                    key_lines.append((ind, gamma_ind))
                key_lines.append((min(ind, face_center_ind), max(ind, face_center_ind)))
            key_lines.append((face_center_ind, gamma_ind))

        return key_lines

    def _get_key_line_orbits(self, key_points, key_lines, key_points_inds_orbits):
        key_lines_copy = dict(zip(range(len(key_lines)), key_lines))
        key_lines_inds_orbits = []

        i = 0
        while len(key_lines_copy) > 0:
            key_lines_inds_orbits.append([])
            l0ind = list(key_lines_copy)[0]
            l0 = key_lines_copy[l0ind]
            key_lines_inds_orbits[i].append(l0)
            key_lines_copy.pop(l0ind)
            to_pop = []
            p00 = key_points[l0[0]]
            p01 = key_points[l0[1]]
            pmid0 = p00 + e / pi * (p01 - p00)
            for ind_key in key_lines_copy:

                l1 = key_lines_copy[ind_key]
                p10 = key_points[l1[0]]
                p11 = key_points[l1[1]]
                equivptspar = False
                equivptsperp = False
                equivline = False

                if (
                    np.array([l0[0] in orbit and l1[0] in orbit for orbit in key_points_inds_orbits]).any()
                    and np.array([l0[1] in orbit and l1[1] in orbit for orbit in key_points_inds_orbits]).any()
                ):
                    equivptspar = True
                elif (
                    np.array([l0[1] in orbit and l1[0] in orbit for orbit in key_points_inds_orbits]).any()
                    and np.array([l0[0] in orbit and l1[1] in orbit for orbit in key_points_inds_orbits]).any()
                ):
                    equivptsperp = True

                if equivptspar:
                    pmid1 = p10 + e / pi * (p11 - p10)
                    for op in self._rpg:
                        if not equivline:
                            p00pr = np.dot(op, p00)
                            diff0 = p10 - p00pr
                            if self._all_ints(diff0, atol=self._atol):
                                pmid0pr = np.dot(op, pmid0) + diff0
                                p01pr = np.dot(op, p01) + diff0
                                if np.allclose(p11, p01pr, atol=self._atol) and np.allclose(
                                    pmid1, pmid0pr, atol=self._atol
                                ):
                                    equivline = True

                elif equivptsperp:
                    pmid1 = p11 + e / pi * (p10 - p11)
                    for op in self._rpg:
                        if not equivline:
                            p00pr = np.dot(op, p00)
                            diff0 = p11 - p00pr
                            if self._all_ints(diff0, atol=self._atol):
                                pmid0pr = np.dot(op, pmid0) + diff0
                                p01pr = np.dot(op, p01) + diff0
                                if np.allclose(p10, p01pr, atol=self._atol) and np.allclose(
                                    pmid1, pmid0pr, atol=self._atol
                                ):
                                    equivline = True

                if equivline:
                    key_lines_inds_orbits[i].append(l1)
                    to_pop.append(ind_key)

            for key in to_pop:
                key_lines_copy.pop(key)
            i += 1

        return key_lines_inds_orbits

    def _get_little_groups(self, key_points, key_points_inds_orbits, key_lines_inds_orbits):

        little_groups_points = []  # elements are lists of indices of recip_point_group. the
        # list little_groups_points[i] is the little group for the
        # orbit key_points_inds_orbits[i]
        for (i, orbit) in enumerate(key_points_inds_orbits):
            k0 = key_points[orbit[0]]
            little_groups_points.append([])
            for (j, op) in enumerate(self._rpg):
                gamma_to = np.dot(op, -1 * k0) + k0
                check_gamma = True
                if not self._all_ints(gamma_to, atol=self._atol):
                    check_gamma = False
                if check_gamma:
                    little_groups_points[i].append(j)

        # elements are lists of indices of recip_point_group. the list
        # little_groups_lines[i] is
        little_groups_lines = []
        # the little group for the orbit key_points_inds_lines[i]

        for (i, orbit) in enumerate(key_lines_inds_orbits):
            l0 = orbit[0]
            v = key_points[l0[1]] - key_points[l0[0]]
            k0 = key_points[l0[0]] + np.e / pi * v
            little_groups_lines.append([])
            for (j, op) in enumerate(self._rpg):
                gamma_to = np.dot(op, -1 * k0) + k0
                check_gamma = True
                if not self._all_ints(gamma_to, atol=self._atol):
                    check_gamma = False
                if check_gamma:
                    little_groups_lines[i].append(j)

        return little_groups_points, little_groups_lines

    def _convert_all_magmoms_to_vectors(self, magmom_axis, axis_specified):
        struct = self._structure.copy()
        magmom_axis = np.array(magmom_axis)
        if "magmom" not in struct.site_properties:
            warn(
                "The 'magmom' property is not set in the structure's site properties."
                "All magnetic moments are being set to zero."
            )
            struct.add_site_property("magmom", [np.array([0, 0, 0]) for i in range(len(struct.sites))])

            return struct

        old_magmoms = struct.site_properties["magmom"]
        new_magmoms = []
        found_scalar = False

        for magmom in old_magmoms:
            if isinstance(magmom, np.ndarray):
                new_magmoms.append(magmom)
            elif isinstance(magmom, list):
                new_magmoms.append(np.array(magmom))
            else:
                found_scalar = True
                new_magmoms.append(magmom * magmom_axis)

        if found_scalar and not axis_specified:
            warn("At least one magmom had a scalar value and magmom_axis was not specified. Defaulted to z+ spinor.")

        struct.remove_site_property("magmom")
        struct.add_site_property("magmom", new_magmoms)
        return struct

    def _get_magnetic_symmetry_operations(self, struct, grey_ops, atol):
        mag_ops = []
        magmoms = struct.site_properties["magmom"]
        nonzero_magmom_inds = [i for i in range(len(struct.sites)) if not (magmoms[i] == np.array([0, 0, 0])).all()]
        init_magmoms = [site.properties["magmom"] for (i, site) in enumerate(struct.sites) if i in nonzero_magmom_inds]
        sites = [site for (i, site) in enumerate(struct.sites) if i in nonzero_magmom_inds]
        init_site_coords = [site.frac_coords for site in sites]
        for op in grey_ops:
            r = op.rotation_matrix
            t = op.translation_vector
            xformed_magmoms = [self._apply_op_to_magmom(r, magmom) for magmom in init_magmoms]
            xformed_site_coords = [np.dot(r, site.frac_coords) + t for site in sites]
            permutation = ["a" for i in range(len(sites))]
            not_found = list(range(len(sites)))
            for i in range(len(sites)):
                xformed = xformed_site_coords[i]
                for k, j in enumerate(not_found):
                    init = init_site_coords[j]
                    diff = xformed - init
                    if self._all_ints(diff, atol=atol):
                        permutation[i] = j
                        not_found.pop(k)
                        break

            same = np.zeros(len(sites))
            flipped = np.zeros(len(sites))
            for i, magmom in enumerate(xformed_magmoms):
                if (magmom == init_magmoms[permutation[i]]).all():
                    same[i] = 1
                elif (magmom == -1 * init_magmoms[permutation[i]]).all():
                    flipped[i] = 1

            if same.all():  # add symm op without tr
                mag_ops.append(
                    MagSymmOp.from_rotation_and_translation_and_time_reversal(
                        rotation_matrix=op.rotation_matrix,
                        translation_vec=op.translation_vector,
                        time_reversal=1,
                    )
                )
            if flipped.all():  # add symm op with tr
                mag_ops.append(
                    MagSymmOp.from_rotation_and_translation_and_time_reversal(
                        rotation_matrix=op.rotation_matrix,
                        translation_vec=op.translation_vector,
                        time_reversal=-1,
                    )
                )

        return mag_ops

    @staticmethod
    def _get_reciprocal_point_group(ops, R, A):
        Ainv = np.linalg.inv(A)
        # convert to reciprocal primitive basis
        recip_point_group = [np.around(np.dot(A, np.dot(R, Ainv)), decimals=2)]
        for op in ops:
            op = np.around(np.dot(A, np.dot(op, Ainv)), decimals=2)
            new = True
            new_coset = True
            for thing in recip_point_group:
                if (thing == op).all():
                    new = False
                if (thing == np.dot(R, op)).all():
                    new_coset = False

            if new:
                recip_point_group.append(op)
            if new_coset:
                recip_point_group.append(np.dot(R, op))

        return recip_point_group

    @staticmethod
    def _closewrapped(pos1, pos2, tolerance):
        pos1 = pos1 % 1.0
        pos2 = pos2 % 1.0

        if len(pos1) != len(pos2):
            return False
        for idx, p1 in enumerate(pos1):
            if abs(p1 - pos2[idx]) > tolerance[idx] and abs(p1 - pos2[idx]) < 1.0 - tolerance[idx]:
                return False
        return True

    def _get_coset_factor(self, G, H):
        # finds g for left coset decomposition G = H + gH (H must be subgroup of G with index two.)
        # in this implementation, G and H are lists of objects of type
        # SymmOp
        gH = []
        for op1 in G:
            in_H = False
            for op2 in H:
                if np.allclose(op1.rotation_matrix, op2.rotation_matrix, atol=self._atol) and self._closewrapped(
                    op1.translation_vector,
                    op2.translation_vector,
                    np.ones(3) * self._atol,
                ):
                    in_H = True
                    break
            if not in_H:
                gH.append(op1)

        for op in gH:
            opH = [op * h for h in H]
            is_coset_factor = True
            for op1 in opH:
                for op2 in H:
                    if np.allclose(op1.rotation_matrix, op2.rotation_matrix, atol=self._atol) and self._closewrapped(
                        op1.translation_vector,
                        op2.translation_vector,
                        np.ones(3) * self._atol,
                    ):

                        is_coset_factor = False
                        break
                if not is_coset_factor:
                    break
            if is_coset_factor:
                return op

        return "No coset factor found."

    @staticmethod
    def _apply_op_to_magmom(r, magmom):
        if np.linalg.det(r) == 1:
            return np.dot(r, magmom)
        return -1 * np.dot(r, magmom)

    @staticmethod
    def _all_ints(arr, atol):
        rounded_arr = np.around(arr, decimals=0)
        return np.allclose(rounded_arr, arr, atol=atol)

    def _get_IRBZ(self, recip_point_group, W, key_points, face_center_inds, atol):
        rpgdict = self._get_reciprocal_point_group_dict(recip_point_group, atol)

        g = np.dot(W.T, W)  # just using change of basis matrix rather than
        # Lattice.get_cartesian_coordinates for conciseness
        ginv = np.linalg.inv(g)
        D = np.linalg.det(W)

        primary_orientation = None
        secondary_orientation = None
        tertiary_orientation = None

        planar_boundaries = []
        IRBZ_points = list(enumerate(key_points))

        for sigma in rpgdict["reflections"]:
            norm = sigma["normal"]
            if primary_orientation is None:
                primary_orientation = norm
                planar_boundaries.append(norm)
            elif np.isclose(np.dot(primary_orientation, np.dot(g, norm)), 0, atol=atol):
                if secondary_orientation is None:
                    secondary_orientation = norm
                    planar_boundaries.append(norm)
                elif np.isclose(np.dot(secondary_orientation, np.dot(g, norm)), 0, atol=atol):
                    if tertiary_orientation is None:
                        tertiary_orientation = norm
                        planar_boundaries.append(norm)
                    elif np.allclose(norm, -1 * tertiary_orientation, atol=atol):
                        pass
                elif np.dot(secondary_orientation, np.dot(g, norm)) < 0:
                    planar_boundaries.append(-1 * norm)
                else:
                    planar_boundaries.append(norm)
            elif np.dot(primary_orientation, np.dot(g, norm)) < 0:
                planar_boundaries.append(-1 * norm)
            else:
                planar_boundaries.append(norm)

        IRBZ_points = self._reduce_IRBZ(IRBZ_points, planar_boundaries, g, atol)

        used_axes = []

        # six-fold rotoinversion always comes with horizontal mirror so don't
        # need to check
        for rotn in rpgdict["rotations"]["six-fold"]:
            ax = rotn["axis"]
            op = rotn["op"]
            if not np.any([np.allclose(ax, usedax, atol) for usedax in used_axes]):
                if self._op_maps_IRBZ_to_self(op, IRBZ_points, atol):
                    face_center_found = False
                    for point in IRBZ_points:
                        if point[0] in face_center_inds:
                            cross = D * np.dot(ginv, np.cross(ax, point[1]))
                            if not np.allclose(cross, 0, atol=atol):
                                rot_boundaries = [cross, -1 * np.dot(op, cross)]
                                face_center_found = True
                                used_axes.append(ax)
                                break
                    if not face_center_found:
                        print("face center not found")
                        for point in IRBZ_points:
                            cross = D * np.dot(ginv, np.cross(ax, point[1]))
                            if not np.allclose(cross, 0, atol=atol):
                                rot_boundaries = [cross, -1 * np.dot(op, cross)]
                                used_axes.append(ax)
                                break
                    IRBZ_points = self._reduce_IRBZ(IRBZ_points, rot_boundaries, g, atol)

        for rotn in rpgdict["rotations"]["rotoinv-four-fold"]:
            ax = rotn["axis"]
            op = rotn["op"]
            if not np.any([np.allclose(ax, usedax, atol) for usedax in used_axes]):
                if self._op_maps_IRBZ_to_self(op, IRBZ_points, atol):
                    face_center_found = False
                    for point in IRBZ_points:
                        if point[0] in face_center_inds:
                            cross = D * np.dot(ginv, np.cross(ax, point[1]))
                            if not np.allclose(cross, 0, atol=atol):
                                rot_boundaries = [cross, np.dot(op, cross)]
                                face_center_found = True
                                used_axes.append(ax)
                                break
                    if not face_center_found:
                        print("face center not found")
                        for point in IRBZ_points:
                            cross = D * np.dot(ginv, np.cross(ax, point[1]))
                            if not np.allclose(cross, 0, atol=atol):
                                rot_boundaries = [cross, -1 * np.dot(op, cross)]
                                used_axes.append(ax)
                                break
                    IRBZ_points = self._reduce_IRBZ(IRBZ_points, rot_boundaries, g, atol)

        for rotn in rpgdict["rotations"]["four-fold"]:
            ax = rotn["axis"]
            op = rotn["op"]
            if not np.any([np.allclose(ax, usedax, atol) for usedax in used_axes]):
                if self._op_maps_IRBZ_to_self(op, IRBZ_points, atol):
                    face_center_found = False
                    for point in IRBZ_points:
                        if point[0] in face_center_inds:
                            cross = D * np.dot(ginv, np.cross(ax, point[1]))
                            if not np.allclose(cross, 0, atol=atol):
                                rot_boundaries = [cross, -1 * np.dot(op, cross)]
                                face_center_found = True
                                used_axes.append(ax)
                                break
                    if not face_center_found:
                        print("face center not found")
                        for point in IRBZ_points:
                            cross = D * np.dot(ginv, np.cross(ax, point[1]))
                            if not np.allclose(cross, 0, atol=atol):
                                rot_boundaries = [cross, -1 * np.dot(op, cross)]
                                used_axes.append(ax)
                                break
                    IRBZ_points = self._reduce_IRBZ(IRBZ_points, rot_boundaries, g, atol)

        for rotn in rpgdict["rotations"]["rotoinv-three-fold"]:
            ax = rotn["axis"]
            op = rotn["op"]
            if not np.any([np.allclose(ax, usedax, atol) for usedax in used_axes]):
                if self._op_maps_IRBZ_to_self(op, IRBZ_points, atol):
                    face_center_found = False
                    for point in IRBZ_points:
                        if point[0] in face_center_inds:
                            cross = D * np.dot(ginv, np.cross(ax, point[1]))
                            if not np.allclose(cross, 0, atol=atol):
                                rot_boundaries = [
                                    cross,
                                    -1 * np.dot(sqrtm(-1 * op), cross),
                                ]
                                face_center_found = True
                                used_axes.append(ax)
                                break
                    if not face_center_found:
                        print("face center not found")
                        for point in IRBZ_points:
                            cross = D * np.dot(ginv, np.cross(ax, point[1]))
                            if not np.allclose(cross, 0, atol=atol):
                                rot_boundaries = [cross, -1 * np.dot(op, cross)]
                                used_axes.append(ax)
                                break
                    IRBZ_points = self._reduce_IRBZ(IRBZ_points, rot_boundaries, g, atol)

        for rotn in rpgdict["rotations"]["three-fold"]:
            ax = rotn["axis"]
            op = rotn["op"]
            if not np.any([np.allclose(ax, usedax, atol) for usedax in used_axes]):
                if self._op_maps_IRBZ_to_self(op, IRBZ_points, atol):
                    face_center_found = False
                    for point in IRBZ_points:
                        if point[0] in face_center_inds:
                            cross = D * np.dot(ginv, np.cross(ax, point[1]))
                            if not np.allclose(cross, 0, atol=atol):
                                rot_boundaries = [cross, -1 * np.dot(op, cross)]
                                face_center_found = True
                                used_axes.append(ax)
                                break
                    if not face_center_found:
                        print("face center not found")
                        for point in IRBZ_points:
                            cross = D * np.dot(ginv, np.cross(ax, point[1]))
                            if not np.allclose(cross, 0, atol=atol):
                                rot_boundaries = [cross, -1 * np.dot(op, cross)]
                                used_axes.append(ax)
                                break
                    IRBZ_points = self._reduce_IRBZ(IRBZ_points, rot_boundaries, g, atol)

        for rotn in rpgdict["rotations"]["two-fold"]:
            ax = rotn["axis"]
            op = rotn["op"]
            if not np.any([np.allclose(ax, usedax, atol) for usedax in used_axes]):
                if self._op_maps_IRBZ_to_self(op, IRBZ_points, atol):
                    face_center_found = False
                    for point in IRBZ_points:
                        if point[0] in face_center_inds:
                            cross = D * np.dot(ginv, np.cross(ax, point[1]))
                            if not np.allclose(cross, 0, atol=atol):
                                rot_boundaries = [cross, -1 * np.dot(op, cross)]
                                face_center_found = True
                                used_axes.append(ax)
                                break
                    if not face_center_found:
                        print("face center not found")
                        for point in IRBZ_points:
                            cross = D * np.dot(ginv, np.cross(ax, point[1]))
                            if not np.allclose(cross, 0, atol=atol):
                                rot_boundaries = [cross, -1 * np.dot(op, cross)]
                                used_axes.append(ax)
                                break
                    IRBZ_points = self._reduce_IRBZ(IRBZ_points, rot_boundaries, g, atol)

        return [point[0] for point in IRBZ_points]

    @staticmethod
    def _get_reciprocal_point_group_dict(recip_point_group, atol):
        PAR = np.array([[-1, 0, 0], [0, -1, 0], [0, 0, -1]])

        d = {
            "reflections": [],
            "rotations": {
                "two-fold": [],
                "three-fold": [],
                "four-fold": [],
                "six-fold": [],
                "rotoinv-three-fold": [],
                "rotoinv-four-fold": [],
                "rotoinv-six-fold": [],
            },
            "inversion": [],
        }

        for i, op in enumerate(recip_point_group):
            evals, evects = np.linalg.eig(op)

            tr = np.trace(op)
            det = np.linalg.det(op)

            # Proper rotations
            if np.isclose(det, 1, atol=atol):
                if np.isclose(tr, 3, atol=atol):
                    continue
                if np.isclose(tr, -1, atol=atol):  # two-fold rotation
                    for j in range(3):
                        if np.isclose(evals[j], 1, atol=atol):
                            ax = evects[:, j]
                    d["rotations"]["two-fold"].append({"ind": i, "axis": ax, "op": op})
                elif np.isclose(tr, 0, atol=atol):  # three-fold rotation
                    for j in range(3):
                        if np.isreal(evals[j]) and np.isclose(np.absolute(evals[j]), 1, atol=atol):
                            ax = evects[:, j]
                    d["rotations"]["three-fold"].append({"ind": i, "axis": ax, "op": op})
                # four-fold rotation
                elif np.isclose(tr, 1, atol=atol):
                    for j in range(3):
                        if np.isreal(evals[j]) and np.isclose(np.absolute(evals[j]), 1, atol=atol):
                            ax = evects[:, j]
                    d["rotations"]["four-fold"].append({"ind": i, "axis": ax, "op": op})
                elif np.isclose(tr, 2, atol=atol):  # six-fold rotation
                    for j in range(3):
                        if np.isreal(evals[j]) and np.isclose(np.absolute(evals[j]), 1, atol=atol):
                            ax = evects[:, j]
                    d["rotations"]["six-fold"].append({"ind": i, "axis": ax, "op": op})

            # Improper rotations
            if np.isclose(det, -1, atol=atol):
                if np.isclose(tr, -3, atol=atol):
                    d["inversion"].append({"ind": i, "op": PAR})
                elif np.isclose(tr, 1, atol=atol):  # two-fold rotation
                    for j in range(3):
                        if np.isclose(evals[j], -1, atol=atol):
                            norm = evects[:, j]
                    d["reflections"].append({"ind": i, "normal": norm, "op": op})
                elif np.isclose(tr, 0, atol=atol):  # three-fold rotoinversion
                    for j in range(3):
                        if np.isreal(evals[j]) and np.isclose(np.absolute(evals[j]), 1, atol=atol):
                            ax = evects[:, j]
                    d["rotations"]["rotoinv-three-fold"].append({"ind": i, "axis": ax, "op": op})
                # four-fold rotoinversion
                elif np.isclose(tr, -1, atol=atol):
                    for j in range(3):
                        if np.isreal(evals[j]) and np.isclose(np.absolute(evals[j]), 1, atol=atol):
                            ax = evects[:, j]
                    d["rotations"]["rotoinv-four-fold"].append({"ind": i, "axis": ax, "op": op})
                # six-fold rotoinversion
                elif np.isclose(tr, -2, atol=atol):
                    for j in range(3):
                        if np.isreal(evals[j]) and np.isclose(np.absolute(evals[j]), 1, atol=atol):
                            ax = evects[:, j]
                    d["rotations"]["rotoinv-six-fold"].append({"ind": i, "axis": ax, "op": op})

        return d

    @staticmethod
    def _op_maps_IRBZ_to_self(op, IRBZ_points, atol):
        point_coords = [point[1] for point in IRBZ_points]
        for point in point_coords:
            point_prime = np.dot(op, point)
            mapped_back = False
            for checkpoint in point_coords:
                if np.allclose(point_prime, checkpoint, atol):
                    mapped_back = True
                    break
            if not mapped_back:
                return False

        return True

    @staticmethod
    def _reduce_IRBZ(IRBZ_points, boundaries, g, atol):
        in_reduced_section = []
        for point in IRBZ_points:
            in_reduced_section.append(
                np.all(
                    [
                        (
                            np.dot(point[1], np.dot(g, boundary)) >= 0
                            or np.isclose(np.dot(point[1], np.dot(g, boundary)), 0, atol=atol)
                        )
                        for boundary in boundaries
                    ]
                )
            )

        return [IRBZ_points[i] for i in range(len(IRBZ_points)) if in_reduced_section[i]]

    def _get_orbit_labels(self, orbit_cosines_orig, key_points_inds_orbits, atol):

        orbit_cosines_copy = orbit_cosines_orig.copy()
        orbit_labels_unsorted = [(len(key_points_inds_orbits) - 1, 26)]
        orbit_inds_remaining = range(len(key_points_inds_orbits) - 1)
        pop_orbits = []
        pop_labels = []

        for i, orb_cos in enumerate(orbit_cosines_copy):
            if np.isclose(orb_cos[0][1], 1.0, atol=atol):
                # (point orbit index, label index)
                orbit_labels_unsorted.append((i, orb_cos[0][0]))
                pop_orbits.append(i)
                pop_labels.append(orb_cos[0][0])

        orbit_cosines_copy = self._reduce_cosines_array(orbit_cosines_copy, pop_orbits, pop_labels)
        orbit_inds_remaining = [i for i in orbit_inds_remaining if i not in pop_orbits]

        # orbit_labels_unsorted already contains gamma orbit
        while len(orbit_labels_unsorted) < len(orbit_cosines_orig) + 1:
            pop_orbits = []
            pop_labels = []
            max_cosine_value = max(orb_cos[0][1] for orb_cos in orbit_cosines_copy)
            max_cosine_value_inds = [
                j for j in range(len(orbit_cosines_copy)) if orbit_cosines_copy[j][0][1] == max_cosine_value
            ]
            max_cosine_label_inds = self._get_max_cosine_labels(
                [orbit_cosines_copy[j] for j in max_cosine_value_inds],
                key_points_inds_orbits,
                atol,
            )

            for j, label_ind in enumerate(max_cosine_label_inds):
                orbit_labels_unsorted.append((orbit_inds_remaining[max_cosine_value_inds[j]], label_ind))
                pop_orbits.append(max_cosine_value_inds[j])
                pop_labels.append(label_ind)
            orbit_cosines_copy = self._reduce_cosines_array(orbit_cosines_copy, pop_orbits, pop_labels)
            orbit_inds_remaining = [
                orbit_inds_remaining[j] for j in range(len(orbit_inds_remaining)) if j not in pop_orbits
            ]

        orbit_labels = np.zeros(len(key_points_inds_orbits))
        for tup in orbit_labels_unsorted:
            orbit_labels[tup[0]] = tup[1]

        return orbit_labels

    @staticmethod
    def _reduce_cosines_array(orbit_cosines, pop_orbits, pop_labels):
        return [
            [orb_cos[i] for i in range(len(orb_cos)) if orb_cos[i][0] not in pop_labels]
            for j, orb_cos in enumerate(orbit_cosines)
            if j not in pop_orbits
        ]

    def _get_max_cosine_labels(self, max_cosine_orbits_orig, key_points_inds_orbits, atol):
        max_cosine_orbits_copy = max_cosine_orbits_orig.copy()
        max_cosine_label_inds = np.zeros(len(max_cosine_orbits_copy))
        initial_max_cosine_label_inds = [max_cos_orb[0][0] for max_cos_orb in max_cosine_orbits_copy]
        u, inds, counts = np.unique(initial_max_cosine_label_inds, return_index=True, return_counts=True)
        grouped_inds = [
            [
                i
                for i in range(len(initial_max_cosine_label_inds))
                if max_cosine_orbits_copy[i][0][0] == max_cosine_orbits_copy[ind][0][0]
            ]
            for ind in inds
        ]
        pop_orbits = []
        pop_labels = []
        unassigned_orbits = []
        for i, ind in enumerate(inds):
            if counts[i] == 1:
                max_cosine_label_inds[ind] = initial_max_cosine_label_inds[ind]
                pop_orbits.append(ind)
                pop_labels.append(initial_max_cosine_label_inds[ind])
            else:
                next_choices = []
                for grouped_ind in grouped_inds[i]:
                    j = 1
                    while True:
                        if max_cosine_orbits_copy[grouped_ind][j][0] not in initial_max_cosine_label_inds:
                            next_choices.append(max_cosine_orbits_copy[grouped_ind][j][1])
                            break
                        j += 1
                worst_next_choice = next_choices.index(min(next_choices))
                for grouped_ind in grouped_inds[i]:
                    if grouped_ind != worst_next_choice:
                        unassigned_orbits.append(grouped_ind)
                max_cosine_label_inds[grouped_inds[i][worst_next_choice]] = initial_max_cosine_label_inds[
                    grouped_inds[i][worst_next_choice]
                ]
                pop_orbits.append(grouped_inds[i][worst_next_choice])
                pop_labels.append(initial_max_cosine_label_inds[grouped_inds[i][worst_next_choice]])

        if len(unassigned_orbits) != 0:
            max_cosine_orbits_copy = self._reduce_cosines_array(max_cosine_orbits_copy, pop_orbits, pop_labels)
            unassigned_orbits_labels = self._get_orbit_labels(max_cosine_orbits_copy, key_points_inds_orbits, atol)
            for i, unassigned_orbit in enumerate(unassigned_orbits):
                max_cosine_label_inds[unassigned_orbit] = unassigned_orbits_labels[i]

        return max_cosine_label_inds

    @staticmethod
    def LabelPoints(index):
        """
        Axes used in generating labels for Latimer-Munro convention
        """
        points = [
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 1, 0],
            [1, 0, 1],
            [0, 1, 1],
            [1, 1, 1],
            [1, 2, 0],
            [1, 0, 2],
            [1, 2, 2],
            [2, 1, 0],
            [0, 1, 2],
            [2, 1, 2],
            [2, 0, 1],
            [0, 2, 1],
            [2, 2, 1],
            [1, 1, 2],
            [1, 2, 1],
            [2, 1, 1],
            [3, 3, 2],
            [3, 2, 3],
            [2, 3, 3],
            [2, 2, 2],
            [3, 2, 2],
            [2, 3, 2],
            [1e-10, 1e-10, 1e-10],
        ]

        return points[index]

    @staticmethod
    def LabelSymbol(index):
        """
        Letters used in generating labels for the Latimer-Munro convention
        """
        symbols = [
            "a",
            "b",
            "c",
            "d",
            "e",
            "f",
            "g",
            "h",
            "i",
            "j",
            "k",
            "l",
            "m",
            "n",
            "o",
            "p",
            "q",
            "r",
            "s",
            "t",
            "u",
            "v",
            "w",
            "x",
            "y",
            "z",
            "Γ",
        ]
        return symbols[index]
