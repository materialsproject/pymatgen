# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import numpy as np
import warnings
from math import ceil
from math import cos
from math import sin
from math import tan
from math import pi
from warnings import warn
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

"""
Created on March 25, 2013

@author: geoffroy
"""

class HighSymmKpath(object):
    """
    This class looks for path along high symmetry lines in
    the Brillouin Zone.
    It is based on Setyawan, W., & Curtarolo, S. (2010).
    High-throughput electronic band structure calculations:
    Challenges and tools. Computational Materials Science,
    49(2), 299-312. doi:10.1016/j.commatsci.2010.05.010
    It should be used with primitive structures that
    comply with the definition from the paper.
    The symmetry is determined by spglib through the
    SpacegroupAnalyzer class. The analyzer can be used to
    produce the correct primitive structure (method
    get_primitive_standard_structure(international_monoclinic=False)).
    A warning will signal possible compatibility problems
    with the given structure.

    Args:
        structure (Structure): Structure object
        symprec (float): Tolerance for symmetry finding
        angle_tolerance (float): Angle tolerance for symmetry finding.
        atol (float): Absolute tolerance used to compare the input
            structure with the one expected as primitive standard.
            A warning will be issued if the lattices don't match.
    """

    def __init__(self, structure, symprec=0.01, angle_tolerance=5, atol=1e-8):
        self._structure = structure
        self._sym = SpacegroupAnalyzer(structure, symprec=symprec,
                                   angle_tolerance=angle_tolerance)
        self._prim = self._sym\
            .get_primitive_standard_structure(international_monoclinic=False)
        self._conv = self._sym.get_conventional_standard_structure(international_monoclinic=False)
        self._prim_rec = self._prim.lattice.reciprocal_lattice
        self._kpath = None

        #Note: this warning will be issued for space groups 38-41, since the primitive cell must be 
        #reformatted to match Setyawan/Curtarolo convention in order to work with the current k-path 
        #generation scheme.
        if not np.allclose(self._structure.lattice.matrix, self._prim.lattice.matrix, atol=atol):
            warnings.warn("The input structure does not match the expected standard primitive! "
                          "The path can be incorrect. Use at your own risk.")

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
                warn("Unexpected value for spg_symbol: %s" % spg_symbol)

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
                warn("Unexpected value for spg_symbol: %s" % spg_symbol)

        elif lattice_type == "orthorhombic":
            a = self._conv.lattice.abc[0]
            b = self._conv.lattice.abc[1]
            c = self._conv.lattice.abc[2]

            if "P" in spg_symbol:
                self._kpath = self.orc()

            elif "F" in spg_symbol:
                if 1 / a ** 2 > 1 / b ** 2 + 1 / c ** 2:
                    self._kpath = self.orcf1(a, b, c)
                elif 1 / a ** 2 < 1 / b ** 2 + 1 / c ** 2:
                    self._kpath = self.orcf2(a, b, c)
                else:
                    self._kpath = self.orcf3(a, b, c)

            elif "I" in spg_symbol:
                self._kpath = self.orci(a, b, c)

            elif "C" in spg_symbol or "A" in spg_symbol:
                self._kpath = self.orcc(a, b, c)
            else:
                warn("Unexpected value for spg_symbol: %s" % spg_symbol)

        elif lattice_type == "hexagonal":
            self._kpath = self.hex()

        elif lattice_type == "rhombohedral":
            alpha = self._prim.lattice.lengths_and_angles[1][0]
            if alpha < 90:
                self._kpath = self.rhl1(alpha * pi / 180)
            else:
                self._kpath = self.rhl2(alpha * pi / 180)

        elif lattice_type == "monoclinic":
            a, b, c = self._conv.lattice.abc
            alpha = self._conv.lattice.lengths_and_angles[1][0]
            #beta = self._conv.lattice.lengths_and_angles[1][1]

            if "P" in spg_symbol:
                self._kpath = self.mcl(b, c, alpha * pi / 180)

            elif "C" in spg_symbol:
                kgamma = self._prim_rec.lengths_and_angles[1][2]
                if kgamma > 90:
                    self._kpath = self.mclc1(a, b, c, alpha * pi / 180)
                if kgamma == 90:
                    self._kpath = self.mclc2(a, b, c, alpha * pi / 180)
                if kgamma < 90:
                    if b * cos(alpha * pi / 180) / c\
                            + b ** 2 * sin(alpha * pi / 180) ** 2 / a ** 2 < 1:
                        self._kpath = self.mclc3(a, b, c, alpha * pi / 180)
                    if b * cos(alpha * pi / 180) / c \
                            + b ** 2 * sin(alpha * pi / 180) ** 2 / a ** 2 == 1:
                        self._kpath = self.mclc4(a, b, c, alpha * pi / 180)
                    if b * cos(alpha * pi / 180) / c \
                            + b ** 2 * sin(alpha * pi / 180) ** 2 / a ** 2 > 1:
                        self._kpath = self.mclc5(a, b, c, alpha * pi / 180)
            else:
                warn("Unexpected value for spg_symbol: %s" % spg_symbol)

        elif lattice_type == "triclinic":
            kalpha = self._prim_rec.lengths_and_angles[1][0]
            kbeta = self._prim_rec.lengths_and_angles[1][1]
            kgamma = self._prim_rec.lengths_and_angles[1][2]
            if kalpha > 90 and kbeta > 90 and kgamma > 90:
                self._kpath = self.tria()
            if kalpha < 90 and kbeta < 90 and kgamma < 90:
                self._kpath = self.trib()
            if kalpha > 90 and kbeta > 90 and kgamma == 90:
                self._kpath = self.tria()
            if kalpha < 90 and kbeta < 90 and kgamma == 90:
                self._kpath = self.trib()

        else:
            warn("Unknown lattice type %s" % lattice_type)

    @property
    def structure(self):
        """
        Returns:
            The standardized primitive structure
        """
        return self._prim

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
        return self._prim_rec

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
            the kpoints along the paths in cartesian coordinates
            together with the labels for symmetry points -Wei
        """
        list_k_points = []
        sym_point_labels = []
        for b in self.kpath['path']:
            for i in range(1, len(b)):
                start = np.array(self.kpath['kpoints'][b[i - 1]])
                end = np.array(self.kpath['kpoints'][b[i]])
                distance = np.linalg.norm(
                    self._prim_rec.get_cartesian_coords(start) -
                    self._prim_rec.get_cartesian_coords(end))
                nb = int(ceil(distance * line_density))
                sym_point_labels.extend([b[i - 1]] + [''] * (nb - 1) + [b[i]])
                list_k_points.extend(
                    [self._prim_rec.get_cartesian_coords(start)
                     + float(i) / float(nb) *
                     (self._prim_rec.get_cartesian_coords(end)
                      - self._prim_rec.get_cartesian_coords(start))
                     for i in range(0, nb + 1)])
        if coords_are_cartesian:
            return list_k_points, sym_point_labels
        else:
            frac_k_points = [self._prim_rec.get_fractional_coords(k)
                             for k in list_k_points]
            return frac_k_points, sym_point_labels

    def cubic(self):
        self.name = "CUB"
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'X': np.array([0.0, 0.5, 0.0]),
                   'R': np.array([0.5, 0.5, 0.5]),
                   'M': np.array([0.5, 0.5, 0.0])}
        path = [["\\Gamma", "X", "M", "\\Gamma", "R", "X"], ["M", "R"]]
        return {'kpoints': kpoints, 'path': path}

    def fcc(self):
        self.name = "FCC"
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'K': np.array([3.0 / 8.0, 3.0 / 8.0, 3.0 / 4.0]),
                   'L': np.array([0.5, 0.5, 0.5]),
                   'U': np.array([5.0 / 8.0, 1.0 / 4.0, 5.0 / 8.0]),
                   'W': np.array([0.5, 1.0 / 4.0, 3.0 / 4.0]),
                   'X': np.array([0.5, 0.0, 0.5])}
        path = [["\\Gamma", "X", "W", "K",
                 "\\Gamma", "L", "U", "W", "L", "K"], ["U", "X"]]
        return {'kpoints': kpoints, 'path': path}

    def bcc(self):
        self.name = "BCC"
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'H': np.array([0.5, -0.5, 0.5]),
                   'P': np.array([0.25, 0.25, 0.25]),
                   'N': np.array([0.0, 0.0, 0.5])}
        path = [["\\Gamma", "H", "N", "\\Gamma", "P", "H"], ["P", "N"]]
        return {'kpoints': kpoints, 'path': path}

    def tet(self):
        self.name = "TET"
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'A': np.array([0.5, 0.5, 0.5]),
                   'M': np.array([0.5, 0.5, 0.0]),
                   'R': np.array([0.0, 0.5, 0.5]),
                   'X': np.array([0.0, 0.5, 0.0]),
                   'Z': np.array([0.0, 0.0, 0.5])}
        path = [["\\Gamma", "X", "M", "\\Gamma", "Z", "R", "A", "Z"], ["X", "R"],
                ["M", "A"]]
        return {'kpoints': kpoints, 'path': path}

    def bctet1(self, c, a):
        self.name = "BCT1"
        eta = (1 + c ** 2 / a ** 2) / 4.0
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'M': np.array([-0.5, 0.5, 0.5]),
                   'N': np.array([0.0, 0.5, 0.0]),
                   'P': np.array([0.25, 0.25, 0.25]),
                   'X': np.array([0.0, 0.0, 0.5]),
                   'Z': np.array([eta, eta, -eta]),
                   'Z_1': np.array([-eta, 1 - eta, eta])}
        path = [["\\Gamma", "X", "M", "\\Gamma", "Z", "P", "N", "Z_1", "M"],
                ["X", "P"]]
        return {'kpoints': kpoints, 'path': path}

    def bctet2(self, c, a):
        self.name = "BCT2"
        eta = (1 + a ** 2 / c ** 2) / 4.0
        zeta = a ** 2 / (2 * c ** 2)
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'N': np.array([0.0, 0.5, 0.0]),
                   'P': np.array([0.25, 0.25, 0.25]),
                   '\\Sigma': np.array([-eta, eta, eta]),
                   '\\Sigma_1': np.array([eta, 1 - eta, -eta]),
                   'X': np.array([0.0, 0.0, 0.5]),
                   'Y': np.array([-zeta, zeta, 0.5]),
                   'Y_1': np.array([0.5, 0.5, -zeta]),
                   'Z': np.array([0.5, 0.5, -0.5])}
        path = [["\\Gamma", "X", "Y", "\\Sigma", "\\Gamma", "Z",
                 "\\Sigma_1", "N", "P", "Y_1", "Z"], ["X", "P"]]
        return {'kpoints': kpoints, 'path': path}

    def orc(self):
        self.name = "ORC"
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'R': np.array([0.5, 0.5, 0.5]),
                   'S': np.array([0.5, 0.5, 0.0]),
                   'T': np.array([0.0, 0.5, 0.5]),
                   'U': np.array([0.5, 0.0, 0.5]),
                   'X': np.array([0.5, 0.0, 0.0]),
                   'Y': np.array([0.0, 0.5, 0.0]),
                   'Z': np.array([0.0, 0.0, 0.5])}
        path = [["\\Gamma", "X", "S", "Y", "\\Gamma",
                 "Z", "U", "R", "T", "Z"], ["Y", "T"], ["U", "X"], ["S", "R"]]
        return {'kpoints': kpoints, 'path': path}

    def orcf1(self, a, b, c):
        self.name = "ORCF1"
        zeta = (1 + a ** 2 / b ** 2 - a ** 2 / c ** 2) / 4
        eta = (1 + a ** 2 / b ** 2 + a ** 2 / c ** 2) / 4

        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'A': np.array([0.5, 0.5 + zeta, zeta]),
                   'A_1': np.array([0.5, 0.5 - zeta, 1 - zeta]),
                   'L': np.array([0.5, 0.5, 0.5]),
                   'T': np.array([1, 0.5, 0.5]),
                   'X': np.array([0.0, eta, eta]),
                   'X_1': np.array([1, 1 - eta, 1 - eta]),
                   'Y': np.array([0.5, 0.0, 0.5]),
                   'Z': np.array([0.5, 0.5, 0.0])}
        path = [["\\Gamma", "Y", "T", "Z", "\\Gamma", "X", "A_1", "Y"],
                ["T", "X_1"], ["X", "A", "Z"], ["L", "\\Gamma"]]
        return {'kpoints': kpoints, 'path': path}

    def orcf2(self, a, b, c):
        self.name = "ORCF2"
        phi = (1 + c ** 2 / b ** 2 - c ** 2 / a ** 2) / 4
        eta = (1 + a ** 2 / b ** 2 - a ** 2 / c ** 2) / 4
        delta = (1 + b ** 2 / a ** 2 - b ** 2 / c ** 2) / 4
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'C': np.array([0.5, 0.5 - eta, 1 - eta]),
                   'C_1': np.array([0.5, 0.5 + eta, eta]),
                   'D': np.array([0.5 - delta, 0.5, 1 - delta]),
                   'D_1': np.array([0.5 + delta, 0.5, delta]),
                   'L': np.array([0.5, 0.5, 0.5]),
                   'H': np.array([1 - phi, 0.5 - phi, 0.5]),
                   'H_1': np.array([phi, 0.5 + phi, 0.5]),
                   'X': np.array([0.0, 0.5, 0.5]),
                   'Y': np.array([0.5, 0.0, 0.5]),
                   'Z': np.array([0.5, 0.5, 0.0])}
        path = [["\\Gamma", "Y", "C", "D", "X", "\\Gamma",
                 "Z", "D_1", "H", "C"], ["C_1", "Z"], ["X", "H_1"], ["H", "Y"],
                ["L", "\\Gamma"]]
        return {'kpoints': kpoints, 'path': path}

    def orcf3(self, a, b, c):
        self.name = "ORCF3"
        zeta = (1 + a ** 2 / b ** 2 - a ** 2 / c ** 2) / 4
        eta = (1 + a ** 2 / b ** 2 + a ** 2 / c ** 2) / 4
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'A': np.array([0.5, 0.5 + zeta, zeta]),
                   'A_1': np.array([0.5, 0.5 - zeta, 1 - zeta]),
                   'L': np.array([0.5, 0.5, 0.5]),
                   'T': np.array([1, 0.5, 0.5]),
                   'X': np.array([0.0, eta, eta]),
                   'X_1': np.array([1, 1 - eta, 1 - eta]),
                   'Y': np.array([0.5, 0.0, 0.5]),
                   'Z': np.array([0.5, 0.5, 0.0])}
        path = [["\\Gamma", "Y", "T", "Z", "\\Gamma", "X", "A_1", "Y"],
                ["X", "A", "Z"], ["L", "\\Gamma"]]
        return {'kpoints': kpoints, 'path': path}

    def orci(self, a, b, c):
        self.name = "ORCI"
        zeta = (1 + a ** 2 / c ** 2) / 4
        eta = (1 + b ** 2 / c ** 2) / 4
        delta = (b ** 2 - a ** 2) / (4 * c ** 2)
        mu = (a ** 2 + b ** 2) / (4 * c ** 2)
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'L': np.array([-mu, mu, 0.5 - delta]),
                   'L_1': np.array([mu, -mu, 0.5 + delta]),
                   'L_2': np.array([0.5 - delta, 0.5 + delta, -mu]),
                   'R': np.array([0.0, 0.5, 0.0]),
                   'S': np.array([0.5, 0.0, 0.0]),
                   'T': np.array([0.0, 0.0, 0.5]),
                   'W': np.array([0.25, 0.25, 0.25]),
                   'X': np.array([-zeta, zeta, zeta]),
                   'X_1': np.array([zeta, 1 - zeta, -zeta]),
                   'Y': np.array([eta, -eta, eta]),
                   'Y_1': np.array([1 - eta, eta, -eta]),
                   'Z': np.array([0.5, 0.5, -0.5])}
        path = [["\\Gamma", "X", "L", "T", "W", "R", "X_1", "Z",
                 "\\Gamma", "Y", "S", "W"], ["L_1", "Y"], ["Y_1", "Z"]]
        return {'kpoints': kpoints, 'path': path}

    def orcc(self, a, b, c):
        self.name = "ORCC"
        zeta = (1 + a ** 2 / b ** 2) / 4
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'A': np.array([zeta, zeta, 0.5]),
                   'A_1': np.array([-zeta, 1 - zeta, 0.5]),
                   'R': np.array([0.0, 0.5, 0.5]),
                   'S': np.array([0.0, 0.5, 0.0]),
                   'T': np.array([-0.5, 0.5, 0.5]),
                   'X': np.array([zeta, zeta, 0.0]),
                   'X_1': np.array([-zeta, 1 - zeta, 0.0]),
                   'Y': np.array([-0.5, 0.5, 0]),
                   'Z': np.array([0.0, 0.0, 0.5])}
        path = [["\\Gamma", "X", "S", "R", "A", "Z",
                 "\\Gamma", "Y", "X_1", "A_1", "T", "Y"], ["Z", "T"]]
        return {'kpoints': kpoints, 'path': path}

    def hex(self):
        self.name = "HEX"
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'A': np.array([0.0, 0.0, 0.5]),
                   'H': np.array([1.0 / 3.0, 1.0 / 3.0, 0.5]),
                   'K': np.array([1.0 / 3.0, 1.0 / 3.0, 0.0]),
                   'L': np.array([0.5, 0.0, 0.5]),
                   'M': np.array([0.5, 0.0, 0.0])}
        path = [["\\Gamma", "M", "K", "\\Gamma", "A", "L", "H", "A"], ["L", "M"],
                ["K", "H"]]
        return {'kpoints': kpoints, 'path': path}

    def rhl1(self, alpha):
        self.name = "RHL1"
        eta = (1 + 4 * cos(alpha)) / (2 + 4 * cos(alpha))
        nu = 3.0 / 4.0 - eta / 2.0
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'B': np.array([eta, 0.5, 1.0 - eta]),
                   'B_1': np.array([1.0 / 2.0, 1.0 - eta, eta - 1.0]),
                   'F': np.array([0.5, 0.5, 0.0]),
                   'L': np.array([0.5, 0.0, 0.0]),
                   'L_1': np.array([0.0, 0.0, -0.5]),
                   'P': np.array([eta, nu, nu]),
                   'P_1': np.array([1.0 - nu, 1.0 - nu, 1.0 - eta]),
                   'P_2': np.array([nu, nu, eta - 1.0]),
                   'Q': np.array([1.0 - nu, nu, 0.0]),
                   'X': np.array([nu, 0.0, -nu]),
                   'Z': np.array([0.5, 0.5, 0.5])}
        path = [["\\Gamma", "L", "B_1"], ["B", "Z", "\\Gamma", "X"],
                ["Q", "F", "P_1", "Z"], ["L", "P"]]
        return {'kpoints': kpoints, 'path': path}

    def rhl2(self, alpha):
        self.name = "RHL2"
        eta = 1 / (2 * tan(alpha / 2.0) ** 2)
        nu = 3.0 / 4.0 - eta / 2.0
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'F': np.array([0.5, -0.5, 0.0]),
                   'L': np.array([0.5, 0.0, 0.0]),
                   'P': np.array([1 - nu, -nu, 1 - nu]),
                   'P_1': np.array([nu, nu - 1.0, nu - 1.0]),
                   'Q': np.array([eta, eta, eta]),
                   'Q_1': np.array([1.0 - eta, -eta, -eta]),
                   'Z': np.array([0.5, -0.5, 0.5])}
        path = [["\\Gamma", "P", "Z", "Q", "\\Gamma",
                 "F", "P_1", "Q_1", "L", "Z"]]
        return {'kpoints': kpoints, 'path': path}

    def mcl(self, b, c, beta):
        self.name = "MCL"
        eta = (1 - b * cos(beta) / c) / (2 * sin(beta) ** 2)
        nu = 0.5 - eta * c * cos(beta) / b
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'A': np.array([0.5, 0.5, 0.0]),
                   'C': np.array([0.0, 0.5, 0.5]),
                   'D': np.array([0.5, 0.0, 0.5]),
                   'D_1': np.array([0.5, 0.5, -0.5]),
                   'E': np.array([0.5, 0.5, 0.5]),
                   'H': np.array([0.0, eta, 1.0 - nu]),
                   'H_1': np.array([0.0, 1.0 - eta, nu]),
                   'H_2': np.array([0.0, eta, -nu]),
                   'M': np.array([0.5, eta, 1.0 - nu]),
                   'M_1': np.array([0.5, 1 - eta, nu]),
                   'M_2': np.array([0.5, 1 - eta, nu]),
                   'X': np.array([0.0, 0.5, 0.0]),
                   'Y': np.array([0.0, 0.0, 0.5]),
                   'Y_1': np.array([0.0, 0.0, -0.5]),
                   'Z': np.array([0.5, 0.0, 0.0])}
        path = [["\\Gamma", "Y", "H", "C", "E", "M_1", "A", "X", "H_1"],
                ["M", "D", "Z"], ["Y", "D"]]
        return {'kpoints': kpoints, 'path': path}

    def mclc1(self, a, b, c, alpha):
        self.name = "MCLC1"
        zeta = (2 - b * cos(alpha) / c) / (4 * sin(alpha) ** 2)
        eta = 0.5 + 2 * zeta * c * cos(alpha) / b
        psi = 0.75 - a ** 2 / (4 * b ** 2 * sin(alpha) ** 2)
        phi = psi + (0.75 - psi) * b * cos(alpha) / c
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'N': np.array([0.5, 0.0, 0.0]),
                   'N_1': np.array([0.0, -0.5, 0.0]),
                   'F': np.array([1 - zeta, 1 - zeta, 1 - eta]),
                   'F_1': np.array([zeta, zeta, eta]),
                   'F_2': np.array([-zeta, -zeta, 1 - eta]),
                   #'F_3': np.array([1 - zeta, -zeta, 1 - eta]),
                   'I': np.array([phi, 1 - phi, 0.5]),
                   'I_1': np.array([1 - phi, phi - 1, 0.5]),
                   'L': np.array([0.5, 0.5, 0.5]),
                   'M': np.array([0.5, 0.0, 0.5]),
                   'X': np.array([1 - psi, psi - 1, 0.0]),
                   'X_1': np.array([psi, 1 - psi, 0.0]),
                   'X_2': np.array([psi - 1, -psi, 0.0]),
                   'Y': np.array([0.5, 0.5, 0.0]),
                   'Y_1': np.array([-0.5, -0.5, 0.0]),
                   'Z': np.array([0.0, 0.0, 0.5])}
        path = [["\\Gamma", "Y", "F", "L", "I"], ["I_1", "Z", "F_1"],
                ["Y", "X_1"], ["X", "\\Gamma", "N"], ["M", "\\Gamma"]]
        return {'kpoints': kpoints, 'path': path}

    def mclc2(self, a, b, c, alpha):
        self.name = "MCLC2"
        zeta = (2 - b * cos(alpha) / c) / (4 * sin(alpha) ** 2)
        eta = 0.5 + 2 * zeta * c * cos(alpha) / b
        psi = 0.75 - a ** 2 / (4 * b ** 2 * sin(alpha) ** 2)
        phi = psi + (0.75 - psi) * b * cos(alpha) / c
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'N': np.array([0.5, 0.0, 0.0]),
                   'N_1': np.array([0.0, -0.5, 0.0]),
                   'F': np.array([1 - zeta, 1 - zeta, 1 - eta]),
                   'F_1': np.array([zeta, zeta, eta]),
                   'F_2': np.array([-zeta, -zeta, 1 - eta]),
                   'F_3': np.array([1 - zeta, -zeta, 1 - eta]),
                   'I': np.array([phi, 1 - phi, 0.5]),
                   'I_1': np.array([1 - phi, phi - 1, 0.5]),
                   'L': np.array([0.5, 0.5, 0.5]),
                   'M': np.array([0.5, 0.0, 0.5]),
                   'X': np.array([1 - psi, psi - 1, 0.0]),
                   'X_1': np.array([psi, 1 - psi, 0.0]),
                   'X_2': np.array([psi - 1, -psi, 0.0]),
                   'Y': np.array([0.5, 0.5, 0.0]),
                   'Y_1': np.array([-0.5, -0.5, 0.0]),
                   'Z': np.array([0.0, 0.0, 0.5])}
        path = [["\\Gamma", "Y", "F", "L", "I"], ["I_1", "Z", "F_1"],
                ["N", "\\Gamma", "M"]]
        return {'kpoints': kpoints, 'path': path}

    def mclc3(self, a, b, c, alpha):
        self.name = "MCLC3"
        mu = (1 + b ** 2 / a ** 2) / 4.0
        delta = b * c * cos(alpha) / (2 * a ** 2)
        zeta = mu - 0.25 + (1 - b * cos(alpha) / c)\
            / (4 * sin(alpha) ** 2)
        eta = 0.5 + 2 * zeta * c * cos(alpha) / b
        phi = 1 + zeta - 2 * mu
        psi = eta - 2 * delta
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'F': np.array([1 - phi, 1 - phi, 1 - psi]),
                   'F_1': np.array([phi, phi - 1, psi]),
                   'F_2': np.array([1 - phi, -phi, 1 - psi]),
                   'H': np.array([zeta, zeta, eta]),
                   'H_1': np.array([1 - zeta, -zeta, 1 - eta]),
                   'H_2': np.array([-zeta, -zeta, 1 - eta]),
                   'I': np.array([0.5, -0.5, 0.5]),
                   'M': np.array([0.5, 0.0, 0.5]),
                   'N': np.array([0.5, 0.0, 0.0]),
                   'N_1': np.array([0.0, -0.5, 0.0]),
                   'X': np.array([0.5, -0.5, 0.0]),
                   'Y': np.array([mu, mu, delta]),
                   'Y_1': np.array([1 - mu, -mu, -delta]),
                   'Y_2': np.array([-mu, -mu, -delta]),
                   'Y_3': np.array([mu, mu - 1, delta]),
                   'Z': np.array([0.0, 0.0, 0.5])}
        path = [["\\Gamma", "Y", "F", "H", "Z", "I", "F_1"],
                ["H_1", "Y_1", "X", "\\Gamma", "N"], ["M", "\\Gamma"]]
        return {'kpoints': kpoints, 'path': path}

    def mclc4(self, a, b, c, alpha):
        self.name = "MCLC4"
        mu = (1 + b ** 2 / a ** 2) / 4.0
        delta = b * c * cos(alpha) / (2 * a ** 2)
        zeta = mu - 0.25 + (1 - b * cos(alpha) / c)\
            / (4 * sin(alpha) ** 2)
        eta = 0.5 + 2 * zeta * c * cos(alpha) / b
        phi = 1 + zeta - 2 * mu
        psi = eta - 2 * delta
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'F': np.array([1 - phi, 1 - phi, 1 - psi]),
                   'F_1': np.array([phi, phi - 1, psi]),
                   'F_2': np.array([1 - phi, -phi, 1 - psi]),
                   'H': np.array([zeta, zeta, eta]),
                   'H_1': np.array([1 - zeta, -zeta, 1 - eta]),
                   'H_2': np.array([-zeta, -zeta, 1 - eta]),
                   'I': np.array([0.5, -0.5, 0.5]),
                   'M': np.array([0.5, 0.0, 0.5]),
                   'N': np.array([0.5, 0.0, 0.0]),
                   'N_1': np.array([0.0, -0.5, 0.0]),
                   'X': np.array([0.5, -0.5, 0.0]),
                   'Y': np.array([mu, mu, delta]),
                   'Y_1': np.array([1 - mu, -mu, -delta]),
                   'Y_2': np.array([-mu, -mu, -delta]),
                   'Y_3': np.array([mu, mu - 1, delta]),
                   'Z': np.array([0.0, 0.0, 0.5])}
        path = [["\\Gamma", "Y", "F", "H", "Z", "I"],
                ["H_1", "Y_1", "X", "\\Gamma", "N"], ["M", "\\Gamma"]]
        return {'kpoints': kpoints, 'path': path}

    def mclc5(self, a, b, c, alpha):
        self.name = "MCLC5"
        zeta = (b ** 2 / a ** 2 + (1 - b * cos(alpha) / c)
                / sin(alpha) ** 2) / 4
        eta = 0.5 + 2 * zeta * c * cos(alpha) / b
        mu = eta / 2 + b ** 2 / (4 * a ** 2) \
            - b * c * cos(alpha) / (2 * a ** 2)
        nu = 2 * mu - zeta
        rho = 1 - zeta * a ** 2 / b ** 2
        omega = (4 * nu - 1 - b ** 2 * sin(alpha) ** 2 / a ** 2)\
            * c / (2 * b * cos(alpha))
        delta = zeta * c * cos(alpha) / b + omega / 2 - 0.25
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'F': np.array([nu, nu, omega]),
                   'F_1': np.array([1 - nu, 1 - nu, 1 - omega]),
                   'F_2': np.array([nu, nu - 1, omega]),
                   'H': np.array([zeta, zeta, eta]),
                   'H_1': np.array([1 - zeta, -zeta, 1 - eta]),
                   'H_2': np.array([-zeta, -zeta, 1 - eta]),
                   'I': np.array([rho, 1 - rho, 0.5]),
                   'I_1': np.array([1 - rho, rho - 1, 0.5]),
                   'L': np.array([0.5, 0.5, 0.5]),
                   'M': np.array([0.5, 0.0, 0.5]),
                   'N': np.array([0.5, 0.0, 0.0]),
                   'N_1': np.array([0.0, -0.5, 0.0]),
                   'X': np.array([0.5, -0.5, 0.0]),
                   'Y': np.array([mu, mu, delta]),
                   'Y_1': np.array([1 - mu, -mu, -delta]),
                   'Y_2': np.array([-mu, -mu, -delta]),
                   'Y_3': np.array([mu, mu - 1, delta]),
                   'Z': np.array([0.0, 0.0, 0.5])}
        path = [["\\Gamma", "Y", "F", "L", "I"], ["I_1", "Z", "H", "F_1"],
                ["H_1", "Y_1", "X", "\\Gamma", "N"], ["M", "\\Gamma"]]
        return {'kpoints': kpoints, 'path': path}

    def tria(self):
        self.name = "TRI1a"
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'L': np.array([0.5, 0.5, 0.0]),
                   'M': np.array([0.0, 0.5, 0.5]),
                   'N': np.array([0.5, 0.0, 0.5]),
                   'R': np.array([0.5, 0.5, 0.5]),
                   'X': np.array([0.5, 0.0, 0.0]),
                   'Y': np.array([0.0, 0.5, 0.0]),
                   'Z': np.array([0.0, 0.0, 0.5])}
        path = [["X", "\\Gamma", "Y"], ["L", "\\Gamma", "Z"],
                ["N", "\\Gamma", "M"], ["R", "\\Gamma"]]
        return {'kpoints': kpoints, 'path': path}

    def trib(self):
        self.name = "TRI1b"
        kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                   'L': np.array([0.5, -0.5, 0.0]),
                   'M': np.array([0.0, 0.0, 0.5]),
                   'N': np.array([-0.5, -0.5, 0.5]),
                   'R': np.array([0.0, -0.5, 0.5]),
                   'X': np.array([0.0, -0.5, 0.0]),
                   'Y': np.array([0.5, 0.0, 0.0]),
                   'Z': np.array([-0.5, 0.0, 0.5])}
        path = [["X", "\\Gamma", "Y"], ["L", "\\Gamma", "Z"],
                ["N", "\\Gamma", "M"], ["R", "\\Gamma"]]
        return {'kpoints': kpoints, 'path': path}
