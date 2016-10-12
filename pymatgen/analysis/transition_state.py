# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import os
import glob

import numpy as np
from monty.json import jsanitize
from monty.json import MSONable
scipy_old_piecewisepolynomial = True
try:
    from scipy.interpolate import PiecewisePolynomial
except ImportError:
    from scipy.interpolate import CubicSpline
    scipy_old_piecewisepolynomial = False

from pymatgen.util.plotting_utils import get_publication_quality_plot
from pymatgen.io.vasp import Poscar, Outcar

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


class NEBAnalysis(MSONable):
    """
    An NEBAnalysis class.
    """

    def __init__(self, r, energies, forces, structures, **kwargs):
        """
        Initializes an NEBAnalysis from the cumulative root mean squared distances
        between structures, the energies, the forces, the structures and the
        interpolation_order for the analysis.

        Args:
            r: Root mean square distances between structures
            energies: Energies of each structure along reaction coordinate
            forces: Tangent forces along the reaction coordinate.
            structures ([Structure]): List of Structures along reaction
                coordinate.
        """
        self.r = np.array(r)
        self.energies = np.array(energies)
        self.forces = np.array(forces)
        self.structures = structures

        # We do a piecewise interpolation between the points. Each spline (
        # cubic by default) is constrained by the boundary conditions of the
        # energies and the tangent force, i.e., the derivative of
        # the energy at each pair of points.
        if scipy_old_piecewisepolynomial:
            self.spline = PiecewisePolynomial(
                self.r, np.array([self.energies, -self.forces]).T,
                orders=3)
        else:
            # New scipy implementation for scipy > 0.18.0
            self.spline = CubicSpline(x=self.r, y=self.energies, bc_type=((1, 0.0), (1, 0.0)))

    @classmethod
    def from_outcars(cls, outcars, structures, **kwargs):
        """
        Initializes an NEBAnalysis from Outcar and Structure objects. Use
        the static constructors, e.g., :class:`from_dir` instead if you
        prefer to have these automatically generated from a directory of NEB
        calculations.

        Args:
            outcars ([Outcar]): List of Outcar objects. Note that these have
                to be ordered from start to end along reaction coordinates.
            structures ([Structure]): List of Structures along reaction
                coordinate. Must be same length as outcar.
            interpolation_order (int): Order of polynomial to use to
                interpolate between images. Same format as order parameter in
                scipy.interplotate.PiecewisePolynomial.
        """
        if len(outcars) != len(structures):
            raise ValueError("# of Outcars must be same as # of Structures")

        # Calculate cumulative root mean square distance between structures,
        # which serves as the reaction coordinate. Note that these are
        # calculated from the final relaxed structures as the coordinates may
        # have changed from the initial interpolation.
        r = [0]
        prev = structures[0]
        for st in structures[1:]:
            dists = np.array([s2.distance(s1) for s1, s2 in zip(prev, st)])
            r.append(np.sqrt(np.sum(dists ** 2)))
            prev = st
        r = np.cumsum(r)

        energies = []
        forces = []
        for i, o in enumerate(outcars):
            o.read_neb()
            energies.append(o.data["energy"])
            if i in [0, len(outcars) - 1]:
                forces.append(0)
            else:
                forces.append(o.data["tangent_force"])
        energies = np.array(energies)
        energies -= energies[0]
        forces = np.array(forces)
        r = np.array(r)
        return cls(r=r, energies=energies, forces=forces, structures=structures, **kwargs)

    def get_extrema(self, normalize_rxn_coordinate=True):
        """
        Returns the positions of the extrema along the MEP. Both local
        minimums and maximums are returned.

        Args:
            normalize_rxn_coordinate (bool): Whether to normalize the
                reaction coordinate to between 0 and 1. Defaults to True.

        Returns:
            (min_extrema, max_extrema), where the extrema are given as
            [(x1, y1), (x2, y2), ...].
        """
        x = np.arange(0, np.max(self.r), 0.01)
        y = self.spline(x) * 1000

        scale = 1 if not normalize_rxn_coordinate else 1 / self.r[-1]
        min_extrema = []
        max_extrema = []
        for i in range(1, len(x) - 1):
            if y[i] < y[i-1] and y[i] < y[i+1]:
                min_extrema.append((x[i] * scale, y[i]))
            elif y[i] > y[i-1] and y[i] > y[i+1]:
                max_extrema.append((x[i] * scale, y[i]))
        return min_extrema, max_extrema

    def get_plot(self, normalize_rxn_coordinate=True, label_barrier=True):
        """
        Returns the NEB plot. Uses Henkelman's approach of spline fitting
        each section of the reaction path based on tangent force and energies.

        Args:
            normalize_rxn_coordinate (bool): Whether to normalize the
                reaction coordinate to between 0 and 1. Defaults to True.
            label_barrier (bool): Whether to label the maximum barrier.

        Returns:
            matplotlib.pyplot object.
        """
        plt = get_publication_quality_plot(12, 8)
        scale = 1 if not normalize_rxn_coordinate else 1 / self.r[-1]
        x = np.arange(0, np.max(self.r), 0.01)
        y = self.spline(x) * 1000
        plt.plot(self.r * scale, self.energies * 1000, 'ro',
                 x * scale, y, 'k-', linewidth=2, markersize=10)
        plt.xlabel("Reaction coordinate")
        plt.ylabel("Energy (meV)")
        plt.ylim((np.min(y) - 10, np.max(y) * 1.02 + 20))
        if label_barrier:
            data = zip(x * scale, y)
            barrier = max(data, key=lambda d: d[1])
            plt.plot([0, barrier[0]], [barrier[1], barrier[1]], 'k--')
            plt.annotate('%.0f meV' % barrier[1],
                         xy=(barrier[0] / 2, barrier[1] * 1.02),
                         xytext=(barrier[0] / 2, barrier[1] * 1.02),
                         horizontalalignment='center')
        plt.tight_layout()
        return plt

    @classmethod
    def from_dir(cls, root_dir, relaxation_dirs=None, **kwargs):
        """
        Initializes a NEBAnalysis object from a directory of a NEB run.
        Note that OUTCARs must be present in all image directories. For the
        terminal OUTCARs from relaxation calculations, you can specify the
        locations using relaxation_dir. If these are not specified, the code
        will attempt to look for the OUTCARs in 00 and 0n directories,
        followed by subdirs "start", "end" or "initial", "final" in the
        root_dir. These are just some typical conventions used
        preferentially in Shyue Ping's MAVRL research group. For the
        non-terminal points, the CONTCAR is read to obtain structures. For
        terminal points, the POSCAR is used. The image directories are
        assumed to be the only directories that can be resolved to integers.
        E.g., "00", "01", "02", "03", "04", "05", "06". The minimum
        sub-directory structure that can be parsed is of the following form (
        a 5-image example is shown):

        00:
        - POSCAR
        - OUTCAR
        01, 02, 03, 04, 05:
        - CONTCAR
        - OUTCAR
        06:
        - POSCAR
        - OUTCAR

        Args:
            root_dir (str): Path to the root directory of the NEB calculation.
            relaxation_dirs (tuple): This specifies the starting and ending
                relaxation directories from which the OUTCARs are read for the
                terminal points for the energies.

        Returns:
            NEBAnalysis object.
        """
        neb_dirs = []

        for d in os.listdir(root_dir):
            pth = os.path.join(root_dir, d)
            if os.path.isdir(pth) and d.isdigit():
                i = int(d)
                neb_dirs.append((i, pth))
        neb_dirs = sorted(neb_dirs, key=lambda d: d[0])
        outcars = []
        structures = []

        # Setup the search sequence for the OUTCARs for the terminal
        # directories.
        terminal_dirs = []
        if relaxation_dirs is not None:
            terminal_dirs.append(relaxation_dirs)
        terminal_dirs.append((neb_dirs[0][1], neb_dirs[-1][1]))
        terminal_dirs.append([os.path.join(root_dir, d)
                              for d in ["start", "end"]])
        terminal_dirs.append([os.path.join(root_dir, d)
                              for d in ["initial", "final"]])

        for i, d in neb_dirs:
            outcar = glob.glob(os.path.join(d, "OUTCAR*"))
            contcar = glob.glob(os.path.join(d, "CONTCAR*"))
            poscar = glob.glob(os.path.join(d, "POSCAR*"))
            terminal = i == 0 or i == neb_dirs[-1][0]
            if terminal:
                found = False
                for ds in terminal_dirs:
                    od = ds[0] if i == 0 else ds[1]
                    outcar = glob.glob(os.path.join(od, "OUTCAR*"))
                    if outcar:
                        outcar = sorted(outcar)
                        outcars.append(Outcar(outcar[-1]))
                        found = True
                        break
                if not found:
                    raise ValueError("OUTCAR cannot be found for terminal "
                                     "point %s" % d)
                structures.append(Poscar.from_file(poscar[0]).structure)
            else:
                outcars.append(Outcar(outcar[0]))
                structures.append(Poscar.from_file(contcar[0]).structure)
        return NEBAnalysis.from_outcars(outcars, structures, **kwargs)

    def as_dict(self):
        """
        Dict representation of NEBAnalysis.

        Returns:
            JSON serializable dict representation.
        """
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                'r': jsanitize(self.r),
                'energies': jsanitize(self.energies),
                'forces': jsanitize(self.forces),
                'structures': [s.as_dict() for s in self.structures]}