"""
Some reimplementation of Henkelman's Transition State Analysis utilities,
which are originally in Perl. Additional features beyond those offered by
Henkelman's utilities will be added.

This allows the usage and customization in Python.
"""

from __future__ import annotations

import os
from glob import glob

import numpy as np
from monty.json import MSONable, jsanitize
from scipy.interpolate import CubicSpline

from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.io.vasp import Outcar, Poscar
from pymatgen.util.plotting import pretty_plot


class NEBAnalysis(MSONable):
    """An NEBAnalysis class."""

    def __init__(self, r, energies, forces, structures, spline_options=None):
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
            spline_options (dict): Options for cubic spline. For example,
                {"saddle_point": "zero_slope"} forces the slope at the saddle to
                be zero.
        """
        self.r = np.array(r)
        self.energies = np.array(energies)
        self.forces = np.array(forces)
        self.structures = structures
        self.spline_options = spline_options if spline_options is not None else {}

        # We do a piecewise interpolation between the points. Each spline (
        # cubic by default) is constrained by the boundary conditions of the
        # energies and the tangent force, i.e., the derivative of
        # the energy at each pair of points.

        self.setup_spline(spline_options=self.spline_options)

    def setup_spline(self, spline_options=None):
        """
        Setup of the options for the spline interpolation.

        Args:
            spline_options (dict): Options for cubic spline. For example,
                {"saddle_point": "zero_slope"} forces the slope at the saddle to
                be zero.
        """
        self.spline_options = spline_options
        relative_energies = self.energies - self.energies[0]
        if self.spline_options.get("saddle_point", "") == "zero_slope":
            imax = np.argmax(relative_energies)
            self.spline = CubicSpline(
                x=self.r[: imax + 1],
                y=relative_energies[: imax + 1],
                bc_type=((1, 0.0), (1, 0.0)),
            )
            cspline2 = CubicSpline(
                x=self.r[imax:],
                y=relative_energies[imax:],
                bc_type=((1, 0.0), (1, 0.0)),
            )
            self.spline.extend(c=cspline2.c, x=cspline2.x[1:])
        else:
            self.spline = CubicSpline(x=self.r, y=relative_energies, bc_type=((1, 0.0), (1, 0.0)))

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
            r.append(np.sqrt(np.sum(dists**2)))
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
            if y[i] < y[i - 1] and y[i] < y[i + 1]:
                min_extrema.append((x[i] * scale, y[i]))
            elif y[i] > y[i - 1] and y[i] > y[i + 1]:
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
        plt = pretty_plot(12, 8)
        scale = 1 if not normalize_rxn_coordinate else 1 / self.r[-1]
        x = np.arange(0, np.max(self.r), 0.01)
        y = self.spline(x) * 1000
        relative_energies = self.energies - self.energies[0]
        plt.plot(
            self.r * scale,
            relative_energies * 1000,
            "ro",
            x * scale,
            y,
            "k-",
            linewidth=2,
            markersize=10,
        )
        plt.xlabel("Reaction coordinate")
        plt.ylabel("Energy (meV)")
        plt.ylim((np.min(y) - 10, np.max(y) * 1.02 + 20))
        if label_barrier:
            data = zip(x * scale, y)
            barrier = max(data, key=lambda d: d[1])
            plt.plot([0, barrier[0]], [barrier[1], barrier[1]], "k--")
            plt.annotate(
                f"{np.max(y) - np.min(y):.0f} meV",
                xy=(barrier[0] / 2, barrier[1] * 1.02),
                xytext=(barrier[0] / 2, barrier[1] * 1.02),
                horizontalalignment="center",
            )
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
        terminal_dirs.append([os.path.join(root_dir, d) for d in ["start", "end"]])
        terminal_dirs.append([os.path.join(root_dir, d) for d in ["initial", "final"]])

        for i, d in neb_dirs:
            outcar = glob(f"{d}/OUTCAR*")
            contcar = glob(f"{d}/CONTCAR*")
            poscar = glob(f"{d}/POSCAR*")
            terminal = i in [0, neb_dirs[-1][0]]
            if terminal:
                for ds in terminal_dirs:
                    od = ds[0] if i == 0 else ds[1]
                    outcar = glob(f"{od}/OUTCAR*")
                    if outcar:
                        outcar = sorted(outcar)
                        outcars.append(Outcar(outcar[-1]))
                        break
                else:
                    raise ValueError(f"OUTCAR cannot be found for terminal point {d}")
                structures.append(Poscar.from_file(poscar[0]).structure)
            else:
                outcars.append(Outcar(outcar[0]))
                structures.append(Poscar.from_file(contcar[0]).structure)
        return NEBAnalysis.from_outcars(outcars, structures, **kwargs)

    def as_dict(self):
        """
        Dict representation of NEBAnalysis.

        Returns:
            JSON-serializable dict representation.
        """
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "r": jsanitize(self.r),
            "energies": jsanitize(self.energies),
            "forces": jsanitize(self.forces),
            "structures": [s.as_dict() for s in self.structures],
        }


def combine_neb_plots(neb_analyses, arranged_neb_analyses=False, reverse_plot=False):
    """
    neb_analyses: a list of NEBAnalysis objects.

    arranged_neb_analyses: The code connects two end points with the
    smallest-energy difference. If all end points have very close energies, it's
    likely to result in an inaccurate connection. Manually arrange neb_analyses
    if the combined plot is not as expected compared with all individual plots.
    E.g., if there are two NEBAnalysis objects to combine, arrange in such a
    way that the end-point energy of the first NEBAnalysis object is the
    start-point energy of the second NEBAnalysis object.
    Note that the barrier labeled in y-axis in the combined plot might be
    different from that in the individual plot due to the reference energy used.
    reverse_plot: reverse the plot or percolation direction.
    return: a NEBAnalysis object
    """
    x = StructureMatcher()
    for neb_index, neb in enumerate(neb_analyses):
        if neb_index == 0:
            neb1 = neb
            neb1_energies = list(neb1.energies)
            neb1_structures = neb1.structures
            neb1_forces = neb1.forces
            neb1_r = neb1.r
            continue

        neb2 = neb
        neb2_energies = list(neb2.energies)

        matching = 0
        for neb1_s in [neb1_structures[0], neb1_structures[-1]]:
            if x.fit(neb1_s, neb2.structures[0]) or x.fit(neb1_s, neb2.structures[-1]):
                matching += 1
                break
        if matching == 0:
            raise ValueError("no matched structures for connection!")

        neb1_start_e, neb1_end_e = neb1_energies[0], neb1_energies[-1]
        neb2_start_e, neb2_end_e = neb2_energies[0], neb2_energies[-1]
        min_e_diff = min(
            [
                abs(neb1_start_e - neb2_start_e),
                abs(neb1_start_e - neb2_end_e),
                abs(neb1_end_e - neb2_start_e),
                abs(neb1_end_e - neb2_end_e),
            ]
        )

        if arranged_neb_analyses:
            neb1_energies = (
                neb1_energies[0 : len(neb1_energies) - 1]
                + [(neb1_energies[-1] + neb2_energies[0]) / 2]
                + neb2_energies[1:]
            )
            neb1_structures = neb1_structures + neb2.structures[1:]
            neb1_forces = list(neb1_forces) + list(neb2.forces)[1:]
            neb1_r = list(neb1_r) + [i + neb1_r[-1] for i in list(neb2.r)[1:]]

        elif abs(neb1_start_e - neb2_start_e) == min_e_diff:
            neb1_energies = list(reversed(neb1_energies[1:])) + neb2_energies
            neb1_structures = list(reversed(neb1_structures[1:])) + neb2.structures
            neb1_forces = list(reversed(list(neb1_forces)[1:])) + list(neb2.forces)
            neb1_r = list(reversed([i * -1 - neb1_r[-1] * -1 for i in list(neb1_r)[1:]])) + [
                i + neb1_r[-1] for i in list(neb2.r)
            ]

        elif abs(neb1_start_e - neb2_end_e) == min_e_diff:
            neb1_energies = neb2_energies + neb1_energies[1:]
            neb1_structures = neb2.structures + neb1_structures[1:]
            neb1_forces = list(neb2.forces) + list(neb1_forces)[1:]
            neb1_r = list(neb2.r) + [i + list(neb2.r)[-1] for i in list(neb1_r)[1:]]

        elif abs(neb1_end_e - neb2_start_e) == min_e_diff:
            neb1_energies = neb1_energies + neb2_energies[1:]
            neb1_structures = neb1_structures + neb2.structures[1:]
            neb1_forces = list(neb1_forces) + list(neb2.forces)[1:]
            neb1_r = list(neb1_r) + [i + neb1_r[-1] for i in list(neb2.r)[1:]]

        else:
            neb1_energies = neb1_energies + list(reversed(neb2_energies))[1:]
            neb1_structures = neb1_structures + list(reversed(neb2.structures))[1:]
            neb1_forces = list(neb1_forces) + list(reversed(list(neb2.forces)))[1:]
            neb1_r = list(neb1_r) + list(
                reversed([i * -1 - list(neb2.r)[-1] * -1 + list(neb1_r)[-1] for i in list(neb2.r)[:-1]])
            )

    if reverse_plot:
        na = NEBAnalysis(
            list(reversed([i * -1 - neb1_r[-1] * -1 for i in list(neb1_r)])),
            list(reversed(neb1_energies)),
            list(reversed(neb1_forces)),
            list(reversed(neb1_structures)),
        )
    else:
        na = NEBAnalysis(neb1_r, neb1_energies, neb1_forces, neb1_structures)
    return na
