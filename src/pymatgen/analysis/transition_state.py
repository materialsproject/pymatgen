"""
Some reimplementation of Henkelman's Transition State Analysis utilities,
which are originally in Perl. Additional features beyond those offered by
Henkelman's utilities will be added.
"""

from __future__ import annotations

import os
import warnings
from glob import glob
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
from monty.json import MSONable, jsanitize
from scipy.interpolate import CubicSpline

from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Structure
from pymatgen.io.vasp import Outcar
from pymatgen.util.plotting import pretty_plot

if TYPE_CHECKING:
    from typing import Any

    from numpy.typing import ArrayLike, NDArray
    from typing_extensions import Self

    from pymatgen.util.typing import PathLike


class NEBAnalysis(MSONable):
    """An NEBAnalysis class."""

    def __init__(
        self,
        r: ArrayLike,
        energies: ArrayLike,
        forces: ArrayLike,
        structures: list[Structure],
        spline_options: dict | None = None,
        zero_slope_saddle: bool = False,
    ) -> None:
        """Initialize an NEBAnalysis from the cumulative root mean squared distances
        between structures, the energies, the forces and the structures for the analysis.

        Args:
            r (ArrayLike): Root mean square distances between structures.
            energies (ArrayLike): Energies of each structure along reaction coordinate.
            forces (ArrayLike): Tangent forces along the reaction coordinate.
            structures (list[Structure]): Structures along reaction coordinate.
            spline_options (dict, optional): [Deprecated] Use `zero_slope_saddle` instead.
            zero_slope_saddle (bool): If True, enforces zero slope at saddle point.
        """
        self.zero_slope_saddle = zero_slope_saddle
        if spline_options is not None:
            warnings.warn(
                "`spline_options` is deprecated and will be removed in 2026-08-01. "
                "Use `zero_slope_saddle=True` instead.",
                DeprecationWarning,
                stacklevel=2,
            )
            self.zero_slope_saddle = spline_options.get("saddle_point") == "zero_slope"

        self.r = np.asarray(r)
        self.energies = np.asarray(energies)
        self.relative_energies: NDArray = self.energies - self.energies[0]
        self.forces = np.asarray(forces)
        self.structures = structures

        # We do a piecewise interpolation between the points. Each spline (
        # cubic by default) is constrained by the boundary conditions of the
        # energies and the tangent force, i.e., the derivative of
        # the energy at each pair of points.
        self.setup_spline(zero_slope_saddle=self.zero_slope_saddle)

    def setup_spline(
        self,
        spline_options: dict | None = None,
        zero_slope_saddle: bool = False,
    ) -> None:
        """
        Set up the cubic spline interpolation of the MEP.

        Args:
            spline_options (dict, optional): [Deprecated] Dictionary of options.
                Set to {"saddle_point": "zero_slope"} to enforce zero slope at saddle.
            zero_slope_saddle (bool): New preferred argument.
                If True, enforces zero slope at the saddle point.
        """
        self.zero_slope_saddle: bool = zero_slope_saddle
        if spline_options is not None:
            warnings.warn(
                "`spline_options` is deprecated and will be removed in 2026-08-01. "
                "Use `zero_slope_saddle=True` instead.",
                DeprecationWarning,
                stacklevel=2,
            )
            self.zero_slope_saddle = spline_options.get("saddle_point") == "zero_slope"

        if self.zero_slope_saddle:
            idx_max: int = np.argmax(self.relative_energies)
            self.spline = CubicSpline(
                x=self.r[: idx_max + 1],
                y=self.relative_energies[: idx_max + 1],
                bc_type=((1, 0.0), (1, 0.0)),
            )
            cspline2 = CubicSpline(
                x=self.r[idx_max:],
                y=self.relative_energies[idx_max:],
                bc_type=((1, 0.0), (1, 0.0)),
            )
            self.spline.extend(c=cspline2.c, x=cspline2.x[1:])
        else:
            self.spline = CubicSpline(x=self.r, y=self.relative_energies, bc_type=((1, 0.0), (1, 0.0)))

    @classmethod
    def from_outcars(
        cls,
        outcars: list[Outcar],
        structures: list[Structure],
        **kwargs,
    ) -> Self:
        """Initialize an NEBAnalysis from Outcar and Structure objects. Use
        the static constructors, e.g. `from_dir` instead if you prefer to
        have these automatically generated from a directory of NEB calculations.

        Args:
            outcars (list[Outcar]): Outcar objects. Note that these have
                to be ordered from start to end along reaction coordinates.
            structures (list[Structure]): Structures along reaction
                coordinate. Must be same length as outcars.
        """
        if len(outcars) != len(structures):
            raise ValueError("# of Outcars must be same as # of Structures")

        # Calculate cumulative root mean square distance between structures,
        # which serves as the reaction coordinate. Note that these are
        # calculated from the final relaxed structures as the coordinates may
        # have changed from the initial interpolation.
        _rms_dist: list[float] = [0]
        prev_struct: Structure = structures[0]
        for struct in structures[1:]:
            dists = np.array([s2.distance(s1) for s1, s2 in zip(prev_struct, struct, strict=True)])
            _rms_dist.append(np.sqrt(np.sum(dists**2)))
            prev_struct = struct
        rms_dist: NDArray = np.cumsum(_rms_dist)

        energies, forces = [], []
        for idx, outcar in enumerate(outcars):
            outcar.read_neb()
            energies.append(outcar.data["energy"])
            if idx in {0, len(outcars) - 1}:
                forces.append(0)
            else:
                forces.append(outcar.data["tangent_force"])

        return cls(
            r=rms_dist,
            energies=np.array(energies),
            forces=np.array(forces),
            structures=structures,
            **kwargs,
        )

    def get_extrema(
        self, normalize_rxn_coordinate: bool = True
    ) -> tuple[list[tuple[float, float]], list[tuple[float, float]]]:
        """Get the positions of the extrema in meV along the minimum energy path (MEP).
        Both local minimums and maximums are returned.

        Args:
            normalize_rxn_coordinate (bool): Whether to normalize the
                reaction coordinate to between 0 and 1. Defaults to True.

        Returns:
            tuple[min_extrema, max_extrema]: where the extrema are given as [(x1, y1), (x2, y2), ...].
        """
        x: NDArray = np.arange(0, np.max(self.r), 0.01)
        y: NDArray = self.spline(x) * 1000

        scale: float = 1 if not normalize_rxn_coordinate else 1 / self.r[-1]
        min_extrema: list[tuple[float, float]] = []
        max_extrema: list[tuple[float, float]] = []
        for i in range(1, len(x) - 1):
            if y[i] < y[i - 1] and y[i] < y[i + 1]:
                min_extrema.append((x[i] * scale, y[i]))
            elif y[i] > y[i - 1] and y[i] > y[i + 1]:
                max_extrema.append((x[i] * scale, y[i]))

        return min_extrema, max_extrema

    def get_plot(
        self,
        normalize_rxn_coordinate: bool = True,
        label_barrier: bool = True,
    ) -> plt.Axes:
        """Get an NEB plot. Uses Henkelman's approach of spline fitting
        each section of the reaction path based on tangent force and energies.

        Args:
            normalize_rxn_coordinate (bool): Whether to normalize the
                reaction coordinate to between 0 and 1. Defaults to True.
            label_barrier (bool): Whether to label the maximum barrier. Defaults to True.

        Returns:
            plt.Axes: matplotlib axes object.
        """
        scale: float = 1 / self.r[-1] if normalize_rxn_coordinate else 1
        xs: NDArray = np.arange(0, np.max(self.r), 0.01)
        ys_mev: NDArray = self.spline(xs) * 1000
        relative_energies_mev: NDArray = self.relative_energies * 1000

        ax = pretty_plot(12, 8)
        ax.plot(
            self.r * scale,
            relative_energies_mev,
            "ro",
            xs * scale,
            ys_mev,
            "k-",
            linewidth=2,
            markersize=10,
        )

        ax.set_xlabel("Reaction Coordinate")
        ax.set_ylabel("Energy (meV)")
        ax.set_ylim((np.min(ys_mev) - 10, np.max(ys_mev) * 1.02 + 20))

        if label_barrier:
            data = zip(xs * scale, ys_mev, strict=True)
            barrier: tuple[float, float] = max(data, key=lambda d: d[1])
            ax.plot([0, barrier[0]], [barrier[1], barrier[1]], "k--", linewidth=0.5)
            ax.annotate(
                f"{max(relative_energies_mev):.0f} meV",
                xy=(barrier[0] / 2, barrier[1] * 1.02),
                xytext=(barrier[0] / 2, barrier[1] * 1.02),
                horizontalalignment="center",
                fontsize=18,
            )

        plt.tight_layout()
        return ax

    @classmethod
    def from_dir(
        cls,
        root_dir: PathLike,
        relaxation_dirs: tuple[PathLike, PathLike] | None = None,
        **kwargs,
    ) -> Self:
        """Initialize an NEBAnalysis object from a directory of an NEB run.

        Note that OUTCARs must be present in all image directories. For the
        terminal OUTCARs from relaxation calculations, you can specify the
        locations using `relaxation_dir`. If these are not specified, the code
        will attempt to look for the OUTCARs in 00 and 0n directories,
        followed by subdirs ("start", "end") or ("initial", "final") in the
        `root_dir`. These are just some typical conventions used
        preferentially in Shyue Ping's MAVRL research group.

        For the non-terminal points, the CONTCAR is read to obtain structures.
        For terminal points, the POSCAR is used. The image directories are
        assumed to be the only directories that can be resolved to integers.
        e.g. "00", "01", ..., "06". The minimum sub-directory structure that
        can be parsed is of the following form (a 5-image example is shown):
            00:
                - POSCAR
                - OUTCAR
            01-05:
                - CONTCAR
                - OUTCAR
            06:
                - POSCAR
                - OUTCAR

        Args:
            root_dir (PathLike): Path to the root directory of the NEB calculation.
            relaxation_dirs (tuple): The (start, end) directories from which
                the OUTCARs are read for the terminal points for the energies.

        Returns:
            NEBAnalysis object.
        """
        neb_dirs = []

        for digit in os.listdir(root_dir):
            pth = os.path.join(root_dir, digit)
            if os.path.isdir(pth) and digit.isdigit():
                idx = int(digit)
                neb_dirs.append((idx, pth))
        neb_dirs = sorted(neb_dirs, key=lambda d: d[0])

        # Setup the search sequence for the OUTCARs for the terminal directories.
        terminal_dirs = [] if relaxation_dirs is None else [relaxation_dirs]
        terminal_dirs += [(neb_dirs[0][1], neb_dirs[-1][1])]
        terminal_dirs += [[os.path.join(root_dir, folder) for folder in ("start", "end")]]
        terminal_dirs += [[os.path.join(root_dir, folder) for folder in ("initial", "final")]]

        outcars = []
        structures = []
        for idx, neb_dir in neb_dirs:
            outcar = glob(f"{neb_dir}/OUTCAR") or glob(f"{neb_dir}/OUTCAR*")
            contcar = glob(f"{neb_dir}/CONTCAR") or glob(f"{neb_dir}/CONTCAR*")
            poscar = glob(f"{neb_dir}/POSCAR") or glob(f"{neb_dir}/POSCAR*")
            terminal: bool = idx in {0, neb_dirs[-1][0]}
            if terminal:
                for dirs in terminal_dirs:
                    od = dirs[0] if idx == 0 else dirs[1]
                    if outcar := (glob(f"{od}/OUTCAR") or glob(f"{od}/OUTCAR*")):
                        outcar = sorted(outcar)
                        outcars.append(Outcar(outcar[-1]))
                        break
                else:
                    raise ValueError(f"OUTCAR cannot be found for terminal point {neb_dir}")
                structures.append(Structure.from_file(poscar[0]))
            else:
                outcars.append(Outcar(outcar[0]))
                structures.append(Structure.from_file(contcar[0]))

        return NEBAnalysis.from_outcars(outcars, structures, **kwargs)

    def as_dict(self) -> dict[str, Any]:
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
            "structures": [struct.as_dict() for struct in self.structures],
        }


def combine_neb_plots(
    neb_analyses: list[NEBAnalysis],
    arranged_neb_analyses: bool = False,
    reverse_plot: bool = False,
) -> NEBAnalysis:
    """
    Combine NEB plots.

    Args:
        neb_analyses: a list of NEBAnalysis objects.
        arranged_neb_analyses: The code connects two end points with the
            smallest-energy difference. If all end points have very close energies, it's
            likely to result in an inaccurate connection. Manually arrange neb_analyses
            if the combined plot is not as expected compared with all individual plots.
            e.g. if there are two NEBAnalysis objects to combine, arrange in such a
            way that the end-point energy of the first NEBAnalysis object is the
            start-point energy of the second NEBAnalysis object.
            Note that the barrier labeled in y-axis in the combined plot might be
            different from that in the individual plot due to the reference energy used.
            reverse_plot: reverse the plot or percolation direction.

    Returns:
        NEBAnalysis object
    """
    x = StructureMatcher()
    neb1_structures = []
    neb1_energies = []
    neb1_forces = []
    neb1_r = []
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
            neb1_energies = [
                *neb1_energies[0 : len(neb1_energies) - 1],
                (neb1_energies[-1] + neb2_energies[0]) / 2,
                *neb2_energies[1:],
            ]
            neb1_structures += neb2.structures[1:]
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
            neb1_energies += neb2_energies[1:]
            neb1_structures += neb2.structures[1:]
            neb1_forces = list(neb1_forces) + list(neb2.forces)[1:]
            neb1_r = list(neb1_r) + [i + neb1_r[-1] for i in list(neb2.r)[1:]]

        else:
            neb1_energies += list(reversed(neb2_energies))[1:]
            neb1_structures += list(reversed(neb2.structures))[1:]
            neb1_forces = list(neb1_forces) + list(reversed(list(neb2.forces)))[1:]
            neb1_r = list(neb1_r) + list(
                reversed([i * -1 - list(neb2.r)[-1] * -1 + list(neb1_r)[-1] for i in list(neb2.r)[:-1]])
            )

    if reverse_plot:
        return NEBAnalysis(
            list(reversed([i * -1 - neb1_r[-1] * -1 for i in list(neb1_r)])),
            list(reversed(neb1_energies)),
            list(reversed(neb1_forces)),
            list(reversed(neb1_structures)),
        )

    return NEBAnalysis(neb1_r, neb1_energies, neb1_forces, neb1_structures)
