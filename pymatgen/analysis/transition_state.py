# coding: utf-8

from __future__ import division, unicode_literals

"""
Some reimplementation of Henkelman's Transition State Analysis utilities,
which are originally in Perl. Additional features beyond those offered by
Henkelman's utilities will be added.

This allows the usage and customization in Python.
"""

__author__ = 'Shyue Ping Ong'
__copyright__ = 'Copyright 2013, The Materials Virtual Lab'
__version__ = '0.1'
__maintainer__ = 'Shyue Ping Ong'
__email__ = 'ongsp@ucsd.edu'
__date__ = '6/1/15'

from pymatgen.io.vaspio import Poscar, Outcar
import numpy as np
import os
import glob
import collections
from pymatgen.util.plotting_utils import get_publication_quality_plot


def is_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def get_cubic_spline(r1, e1, f1, r2, e2, f2):
    r = (r2 - r1)
    a = (2 * e2 - e1) / r ** 3 - (f2 + f1) / r ** 2
    b = (3 * e2 - e1) / r ** 2 + (2 * f2 + f1) / r
    c = -f1
    d = e1
    return lambda x: a * x ** 3 + b * x ** 2 + c * x + d


class CINEBAnalysis(object):
    """
    Note that this only works for CINEB calculations as it depends on certain
    lines present in the OUTCAR. Also, even the start and end point directories
    must contain an OUTCAR.
    """

    def __init__(self, outcars, structures):
        """
        Initializes a CINEBAnalysis from Outcar and Structure objects. Use
        the static constructors if you prefer to have these automatically
        generated from a directory of NEB calculations.

        Args:
            outcars ([Outcar]): List of Outcar objects. Note that these have
                to be ordered from start to end along reaction coordinates.
            structures ([Structure]): List of Structures along reaction
                coordinate. Must be same length as outcar.
        """
        if len(outcars) != len(structures):
            raise ValueError("# of Outcars must be same as # of Structures")
        dists = [0]
        prev = structures[0]
        for s in structures[1:]:
            ss = 0
            for site1, site2 in zip(prev, s):
                dist = site2.distance(site1)
                ss += dist ** 2
            prev_dist = dists[-1]
            dists.append(np.sqrt(ss) + prev_dist)
            prev = s

        energies = np.array([o.data["energy"] for o in outcars])
        energies -= energies[0]
        forces = np.array(
            [o.data.get("tangent_force", 0) for o in outcars])
        self.r = np.array(dists)
        self.energies = energies
        self.forces = forces
        self._fit_splines()

    def _fit_splines(self):
        # Internal spline fitting procedure based on boundary conditions of
        # tangent force and energies.
        spline_params = {}
        for i in range(0, len(self.r) - 1):
            dr = self.r[i + 1] - self.r[i]
            f1 = self.forces[i] * dr
            f2 = self.forces[i + 1] * dr
            u1 = self.energies[i]
            u2 = self.energies[i + 1]
            fs = f1 + f2
            ud = u2 - u1
            d = u1
            c = -f1
            b = 3 * ud + f1 + fs
            a = -2 * ud - fs
            spline_params[i] = a, b, c, d
        self.spline_params = spline_params

    def get_extrema(self, normalize_rxn_coordinate=True):
        x = []
        y = []
        scale = 1 if not normalize_rxn_coordinate else 1 / self.r[-1]
        for i in range(0, len(self.r) - 1):
            dr = self.r[i + 1] - self.r[i]
            a, b, c, d = self.spline_params[i]
            f = np.arange(0, 1.0, 0.01)
            x.extend(self.r[i] + f * dr)
            y.extend(a * f ** 3 + b * f ** 2 + c * f + d)
        for i in range(1, len(x) - 1):
            if y[i] < y[i-1] and y[i] < y[i+1]:
                print "Minimum at r = %s, E = %s" % (x[i] * scale, y[i])
            elif y[i] > y[i-1] and y[i] > y[i+1]:
                print "Maximum at r = %s, E = %s" % (x[i] * scale, y[i])

    def get_plot(self, normalize_rxn_coordinate=True):
        """
        Returns the NEB plot. Uses Henkelman's approach of spline fitting
        each section of the reaction path based on tangent force and energies.

        Args:
            normalize_rxn_coordinate (bool): Whether to normalize the
                reaction coordinate to between 0 and 1. Defaults to True.

        Returns:
            matplotlib.pyplot object.
        """
        plt = get_publication_quality_plot(12, 8)
        scale = 1 if not normalize_rxn_coordinate else 1 / self.r[-1]
        all_y = []
        for i in range(0, len(self.r) - 1):
            dr = self.r[i + 1] - self.r[i]
            a, b, c, d = self.spline_params[i]
            f = np.arange(0, 1.01, 0.01)
            x = self.r[i] + f * dr
            y = (a * f ** 3 + b * f ** 2 + c * f + d) * 1000
            all_y.extend(y)
            plt.plot(x * scale, y, 'k-', linewidth=2)
        plt.plot(self.r * scale, self.energies * 1000, 'ro', markersize=10)
        plt.xlabel("Reaction coordinate")
        plt.ylabel("Energy (meV)")
        plt.ylim((np.min(all_y) - 10, np.max(all_y) + 10))
        self.get_extrema(normalize_rxn_coordinate=normalize_rxn_coordinate)
        return plt

    @staticmethod
    def from_dir(root_dir):
        """
        Initializes a CINEBAnalysis object from a directory of a NEB run.
        Note that the directory must have OUTCARs in all image directories,
        including the terminal points. These are usually obtained from
        relaxation calculations. For the non-terminal points, the CONTCAR is
        read to obtain structures. For terminal points, the POSCAR is used.
        The image directories are assumed to be the only directories that can
        be resolved to integers. E.g., "00", "01", "02", "03", "04", "05", "06".

        Args:
            root_dir (str): Path to the root directory of the NEB calculation.

        Returns:
            CINEBAnalysis object.
        """
        neb_dirs = []
        for d in os.listdir(root_dir):
            if os.path.isdir(d) and is_int(d):
                i = int(os.path.basename(d))
                neb_dirs.append((i, d))
        neb_dirs = sorted(neb_dirs, key=lambda d: d[0])
        outcars = []
        structures = []
        for i, d in neb_dirs:
            outcar = glob.glob(os.path.join(d, "OUTCAR*"))
            contcar = glob.glob(os.path.join(d, "CONTCAR*"))
            poscar = glob.glob(os.path.join(d, "POSCAR*"))
            terminal = i == 0 or i == neb_dirs[-1][0]
            if terminal:
                outcars.append(Outcar(outcar[0], force=False))
                structures.append(Poscar.from_file(contcar[0]).structure)
            else:
                outcars.append(Outcar(outcar[0], force=True))
                structures.append(Poscar.from_file(poscar[0]).structure)
        return CINEBAnalysis(outcars, structures)