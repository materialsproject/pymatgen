"""This module implements core classes for calculation of diffraction patterns."""

from __future__ import annotations

import abc
from collections import defaultdict
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np

from pymatgen.core.spectrum import Spectrum
from pymatgen.util.plotting import add_fig_kwargs, pretty_plot

if TYPE_CHECKING:
    from pymatgen.core import Structure


class DiffractionPattern(Spectrum):
    """A representation of a diffraction pattern."""

    XLABEL = "$2\\Theta$"
    YLABEL = "Intensity"

    def __init__(self, x, y, hkls, d_hkls):
        """
        Args:
            x: Two theta angles.
            y: Intensities
            hkls: [{"hkl": (h, k, l), "multiplicity": mult}],
                where {"hkl": (h, k, l), "multiplicity": mult}
                is a dict of Miller
                indices for all diffracted lattice facets contributing to each
                intensity.
            d_hkls: List of interplanar spacings.
        """
        super().__init__(x, y, hkls, d_hkls)
        self.hkls = hkls
        self.d_hkls = d_hkls


class AbstractDiffractionPatternCalculator(abc.ABC):
    """Abstract base class for computing the diffraction pattern of a crystal."""

    # Tolerance in which to treat two peaks as having the same two theta.
    TWO_THETA_TOL = 1e-5

    # Tolerance in which to treat a peak as effectively 0 if the scaled
    # intensity is less than this number. Since the max intensity is 100,
    # this means the peak must be less than 1e-5 of the peak intensity to be
    # considered as zero. This deals with numerical issues where systematic
    # absences do not cancel exactly to zero.
    SCALED_INTENSITY_TOL = 1e-3

    @abc.abstractmethod
    def get_pattern(self, structure: Structure, scaled=True, two_theta_range=(0, 90)):
        """
        Calculates the diffraction pattern for a structure.

        Args:
            structure (Structure): Input structure
            scaled (bool): Whether to return scaled intensities. The maximum
                peak is set to a value of 100. Defaults to True. Use False if
                you need the absolute values to combine XRD plots.
            two_theta_range ([float of length 2]): Tuple for range of
                two_thetas to calculate in degrees. Defaults to (0, 90). Set to
                None if you want all diffracted beams within the limiting
                sphere of radius 2 / wavelength.

        Returns:
            DiffractionPattern
        """
        raise NotImplementedError

    def get_plot(
        self,
        structure: Structure,
        two_theta_range: tuple[float, float] = (0, 90),
        annotate_peaks="compact",
        ax: plt.Axes = None,
        with_labels=True,
        fontsize=16,
    ) -> plt.Axes:
        """Get the diffraction plot as a matplotlib Axes.

        Args:
            structure: Input structure
            two_theta_range (tuple[float, float]): Range of two_thetas to calculate in degrees.
                Defaults to (0, 90). Set to None if you want all diffracted beams within the limiting
                sphere of radius 2 / wavelength.
            annotate_peaks (str | None): Whether and how to annotate the peaks
                with hkl indices. Default is 'compact', i.e. show short
                version (oriented vertically), e.g. 100. If 'full', show
                long version, e.g. (1, 0, 0). If None, do not show anything.
            ax: matplotlib Axes or None if a new figure should be
                created.
            with_labels: True to add xlabels and ylabels to the plot.
            fontsize: (int) fontsize for peak labels.

        Returns:
            plt.Axes: matplotlib Axes object
        """
        ax = ax or pretty_plot(16, 10)

        xrd = self.get_pattern(structure, two_theta_range=two_theta_range)
        imax = max(xrd.y)

        for two_theta, i, hkls in zip(xrd.x, xrd.y, xrd.hkls, strict=True):
            if two_theta_range[0] <= two_theta <= two_theta_range[1]:
                hkl_tuples = [hkl["hkl"] for hkl in hkls]
                label = ", ".join(map(str, hkl_tuples))  # 'full' label
                ax.plot([two_theta, two_theta], [0, i], color="k", linewidth=3, label=label)

                if annotate_peaks == "full":
                    ax.annotate(
                        label,
                        xy=[two_theta, i],
                        xytext=[two_theta, i],
                        fontsize=fontsize,
                    )
                elif annotate_peaks == "compact":
                    if all(all(i < 10 for i in hkl_tuple) for hkl_tuple in hkl_tuples):
                        label = ",".join("".join(map(str, hkl_tuple)) for hkl_tuple in hkl_tuples)
                        # 'compact' label. Would be unclear for indices >= 10
                        # It would have more than 3 figures, e.g. 1031

                    if i / imax > 0.5:  # Big peak: annotation on the side
                        xytext = [-fontsize / 4, 0]
                        ha = "right"
                        va = "top"
                    else:  # Small peak: annotation on top
                        xytext = [0, 10]
                        ha = "center"
                        va = "bottom"

                    ax.annotate(
                        label,
                        xy=[two_theta, i],
                        xytext=xytext,
                        textcoords="offset points",
                        va=va,
                        ha=ha,
                        rotation=90,
                        fontsize=fontsize,
                    )

        if with_labels:
            ax.set_xlabel(r"$2\theta$ ($^\circ$)")
            ax.set_ylabel("Intensities (scaled)")

        plt.tight_layout()

        return ax

    def show_plot(self, structure: Structure, **kwargs):
        """Show the diffraction plot.

        Args:
            structure (Structure): Input structure
            two_theta_range ([float of length 2]): Tuple for range of
                two_thetas to calculate in degrees. Defaults to (0, 90). Set to
                None if you want all diffracted beams within the limiting
                sphere of radius 2 / wavelength.
            annotate_peaks (str | None): Whether and how to annotate the peaks
                with hkl indices. Default is 'compact', i.e. show short
                version (oriented vertically), e.g. 100. If 'full', show
                long version, e.g. (1, 0, 0). If None, do not show anything.
        """
        self.get_plot(structure, **kwargs).get_figure().show()

    @add_fig_kwargs
    def plot_structures(self, structures, fontsize=6, **kwargs):
        """
        Plot diffraction patterns for multiple structures on the same figure.

        Args:
            structures (Structure): List of structures
            two_theta_range ([float of length 2]): Tuple for range of
                two_thetas to calculate in degrees. Defaults to (0, 90). Set to
                None if you want all diffracted beams within the limiting
                sphere of radius 2 / wavelength.
            annotate_peaks (str | None): Whether and how to annotate the peaks
                with hkl indices. Default is 'compact', i.e. show short
                version (oriented vertically), e.g. 100. If 'full', show
                long version, e.g. (1, 0, 0). If None, do not show anything.
            fontsize: (int) fontsize for peak labels.
        """
        n_rows = len(structures)
        fig, axes = plt.subplots(nrows=n_rows, ncols=1, sharex=True, squeeze=False)

        for idx, (ax, structure) in enumerate(zip(axes.ravel(), structures, strict=True)):
            self.get_plot(structure, fontsize=fontsize, ax=ax, with_labels=idx == n_rows - 1, **kwargs)
            spg_symbol, spg_number = structure.get_space_group_info()
            ax.set_title(f"{structure.formula} {spg_symbol} ({spg_number}) ")

        return fig


def get_unique_families(hkls):
    """Get unique families of Miller indices. Families must be permutations
    of each other.

    Args:
        hkls ([h, k, l]): List of Miller indices.

    Returns:
        {hkl: multiplicity}: A dict with unique hkl and multiplicity.
    """

    # TODO can definitely be sped up
    def is_perm(hkl1, hkl2) -> bool:
        h1 = np.abs(hkl1)
        h2 = np.abs(hkl2)
        return np.all(np.sort(h1) == np.sort(h2))

    unique = defaultdict(list)
    for hkl1 in hkls:
        found = False
        for hkl2, v2 in unique.items():
            if is_perm(hkl1, hkl2):
                found = True
                v2.append(hkl1)
                break
        if not found:
            unique[hkl1].append(hkl1)

    pretty_unique = {}
    for val in unique.values():
        pretty_unique[max(val)] = len(val)

    return pretty_unique
