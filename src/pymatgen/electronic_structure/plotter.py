"""This module implements plotter for DOS and band structure."""

from __future__ import annotations

import copy
import itertools
import logging
import math
import warnings
from collections import Counter
from typing import TYPE_CHECKING, cast

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import palettable
import scipy.interpolate as scint
from matplotlib.collections import LineCollection
from matplotlib.gridspec import GridSpec
from monty.dev import requires
from monty.json import jsanitize
from numpy.typing import ArrayLike

from pymatgen.core import Element
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from pymatgen.electronic_structure.boltztrap import BoltztrapError
from pymatgen.electronic_structure.core import OrbitalType, Spin
from pymatgen.util.plotting import add_fig_kwargs, get_ax3d_fig, pretty_plot

try:
    from mayavi import mlab
except ImportError:
    mlab = None

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Literal

    from pymatgen.electronic_structure.dos import CompleteDos, Dos

logger = logging.getLogger(__name__)

__author__ = "Shyue Ping Ong, Geoffroy Hautier, Anubhav Jain"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "May 1, 2012"


class DosPlotter:
    """Plot DOS. The interface is extremely flexible given there are many
    different ways in which people want to view DOS.
    Typical usage is:
        # Initialize plotter with some optional args. Defaults are usually fine
        plotter = PhononDosPlotter().

        # Add DOS with a label
        plotter.add_dos("Total DOS", dos)

        # Alternatively, you can add a dict of DOS. This is the typical form
        # returned by CompletePhononDos.get_element_dos().
        plotter.add_dos_dict({"dos1": dos1, "dos2": dos2})
        plotter.add_dos_dict(complete_dos.get_spd_dos())
    """

    def __init__(
        self,
        zero_at_efermi: bool = True,
        stack: bool = False,
        sigma: float | None = None,
    ) -> None:
        """
        Args:
            zero_at_efermi (bool): Whether to shift all Dos to have zero energy at the
                fermi energy. Defaults to True.
            stack (bool): Whether to plot the DOS as a stacked area graph
            sigma (float): Specify a standard deviation for Gaussian smearing
                the DOS for nicer looking plots. Defaults to None for no
                smearing.
        """
        self.zero_at_efermi = zero_at_efermi
        self.stack = stack
        self.sigma = sigma
        self._norm_val = True
        self._doses: dict[
            str,
            dict[
                Literal["energies", "densities", "efermi"],
                float | ArrayLike | dict[Spin, ArrayLike],
            ],
        ] = {}

    def add_dos(self, label: str, dos: Dos) -> None:
        """Add a DOS for plotting.

        Args:
            label: a unique label for the DOS.
            dos: Dos object
        """
        if dos.norm_vol is None:
            self._norm_val = False
        energies = dos.energies - dos.efermi if self.zero_at_efermi else dos.energies
        densities = dos.get_smeared_densities(self.sigma) if self.sigma else dos.densities
        efermi = dos.efermi
        self._doses[label] = {
            "energies": energies,
            "densities": densities,
            "efermi": efermi,
        }

    def add_dos_dict(self, dos_dict, key_sort_func=None) -> None:
        """Add a dictionary of DOSs, with an optional sorting
        function for the keys.

        Args:
            dos_dict: dict of {label: Dos}
            key_sort_func: function used to sort the dos_dict keys.
        """
        keys = sorted(dos_dict, key=key_sort_func) if key_sort_func else list(dos_dict)
        for label in keys:
            self.add_dos(label, dos_dict[label])

    def get_dos_dict(self):
        """Get the added doses as a json-serializable dict. Note that if you
        have specified smearing for the DOS plot, the densities returned will
        be the smeared densities, not the original densities.

        Returns:
            dict: Dict of dos data. Generally of the form
            {label: {'energies':..., 'densities': {'up':...}, 'efermi':efermi}}
        """
        return jsanitize(self._doses)

    def get_plot(
        self,
        xlim: tuple[float, float] | None = None,
        ylim: tuple[float, float] | None = None,
        invert_axes: bool = False,
        beta_dashed: bool = False,
    ) -> plt.Axes:
        """Get a matplotlib plot showing the DOS.

        Args:
            xlim (tuple[float, float]): The energy axis limits. Defaults to None for automatic
                determination.
            ylim (tuple[float, float]): The y-axis limits. Defaults to None for automatic determination.
            invert_axes (bool): Whether to invert the x and y axes. Enables chemist style DOS plotting.
                Defaults to False.
            beta_dashed (bool): Plots the beta spin channel with a dashed line. Defaults to False.

        Returns:
            plt.Axes: matplotlib Axes object.
        """
        n_colors = min(9, max(3, len(self._doses)))

        colors = palettable.colorbrewer.qualitative.Set1_9.mpl_colors

        ys = None
        all_densities = []
        all_energies = []
        ax = pretty_plot(12, 8)

        # Note that this complicated processing of energies is to allow for
        # stacked plots in matplotlib.
        for dos in self._doses.values():
            energies = cast(ArrayLike, dos["energies"])
            densities = cast(ArrayLike, dos["densities"])
            if not ys:
                ys = {
                    Spin.up: np.zeros(energies.shape),
                    Spin.down: np.zeros(energies.shape),
                }
            new_dens = {}
            for spin in [Spin.up, Spin.down]:
                if spin in densities:
                    if self.stack:
                        ys[spin] += densities[spin]
                        new_dens[spin] = ys[spin].copy()
                    else:
                        new_dens[spin] = densities[spin]
            all_energies.append(energies)
            all_densities.append(new_dens)

        keys = list(reversed(self._doses))
        all_densities.reverse()
        all_energies.reverse()
        all_pts = []

        for idx, key in enumerate(keys):
            for spin in [Spin.up, Spin.down]:
                if spin in all_densities[idx]:
                    energy = all_energies[idx]
                    densities = list(int(spin) * all_densities[idx][spin])
                    if invert_axes:
                        x = densities
                        y = energy
                    else:
                        x = energy
                        y = densities
                    all_pts.extend(list(zip(x, y, strict=True)))
                    if self.stack:
                        ax.fill(x, y, color=colors[idx % n_colors], label=str(key))
                    elif spin == Spin.down and beta_dashed:
                        ax.plot(
                            x,
                            y,
                            color=colors[idx % n_colors],
                            label=str(key),
                            linestyle="--",
                            linewidth=3,
                        )
                    else:
                        ax.plot(
                            x,
                            y,
                            color=colors[idx % n_colors],
                            label=str(key),
                            linewidth=3,
                        )

        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)
        elif not invert_axes:
            xlim = ax.get_xlim()
            if xlim is None:
                raise RuntimeError("xlim cannot be None.")
            relevant_y = [p[1] for p in all_pts if xlim[0] < p[0] < xlim[1]]
            ax.set_ylim((min(relevant_y), max(relevant_y)))
        if not xlim and invert_axes:
            ylim = ax.get_ylim()
            if ylim is None:
                raise RuntimeError("ylim cannot be None.")
            relevant_y = [p[0] for p in all_pts if ylim[0] < p[1] < ylim[1]]
            ax.set_xlim((min(relevant_y), max(relevant_y)))

        if self.zero_at_efermi:
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            (ax.plot(xlim, [0, 0], "k--", linewidth=2) if invert_axes else ax.plot([0, 0], ylim, "k--", linewidth=2))

        if invert_axes:
            ax.set_ylabel("Energies (eV)")
            ax.set_xlabel(f"Density of states (states/eV{'/Å³' if self._norm_val else ''})")
            ax.axvline(x=0, color="k", linestyle="--", linewidth=2)
        else:
            ax.set_xlabel("Energies (eV)")
            if self._norm_val:
                ax.set_ylabel("Density of states (states/eV/Å³)")
            else:
                ax.set_ylabel("Density of states (states/eV)")
            ax.axhline(y=0, color="k", linestyle="--", linewidth=2)

        # Remove duplicate labels with a dictionary
        handles, labels = ax.get_legend_handles_labels()
        label_dict = dict(zip(labels, handles, strict=True))
        ax.legend(label_dict.values(), label_dict)
        legend_text = ax.get_legend().get_texts()  # all the text.Text instance in the legend
        plt.setp(legend_text, fontsize=30)
        plt.tight_layout()
        return ax

    def save_plot(self, filename: str, xlim=None, ylim=None, invert_axes=False, beta_dashed=False) -> None:
        """Save matplotlib plot to a file.

        Args:
            filename (str): Filename to write to. Must include extension to specify image format.
            xlim: Specifies the x-axis limits. Set to None for automatic
                determination.
            ylim: Specifies the y-axis limits.
            invert_axes (bool): Whether to invert the x and y axes. Enables chemist style DOS plotting.
                Defaults to False.
            beta_dashed (bool): Plots the beta spin channel with a dashed line. Defaults to False.
        """
        self.get_plot(xlim, ylim, invert_axes, beta_dashed)
        plt.savefig(filename)

    def show(self, xlim=None, ylim=None, invert_axes=False, beta_dashed=False) -> None:
        """Show the plot using matplotlib.

        Args:
            xlim: Specifies the x-axis limits. Set to None for automatic
                determination.
            ylim: Specifies the y-axis limits.
            invert_axes (bool): Whether to invert the x and y axes. Enables chemist style DOS plotting.
                Defaults to False.
            beta_dashed (bool): Plots the beta spin channel with a dashed line. Defaults to False.
        """
        self.get_plot(xlim, ylim, invert_axes, beta_dashed)
        plt.show()


class BSPlotter:
    """Plot or get data to facilitate the plotting of band structure."""

    def __init__(self, bs: BandStructureSymmLine) -> None:
        """
        Args:
            bs: A BandStructureSymmLine object.
        """
        self._bs: list[BandStructureSymmLine] = []
        self._nb_bands: list[int] = []

        self.add_bs(bs)

    def _check_bs_kpath(self, band_structs: list[BandStructureSymmLine]) -> Literal[True]:
        """Helper method that check all the band objs in bs_list are
        BandStructureSymmLine objs and they all have the same kpath.
        """
        # check obj type
        for bs in band_structs:
            if not isinstance(bs, BandStructureSymmLine):
                raise TypeError(
                    "BSPlotter only works with BandStructureSymmLine objects. "
                    "A BandStructure object (on a uniform grid for instance and "
                    "not along symmetry lines won't work)"
                )

        # check the kpath
        if len(band_structs) == 1 and not self._bs:
            return True

        if self._bs:
            kpath_ref = [br["name"] for br in self._bs[0].branches]
        else:
            kpath_ref = [br["name"] for br in band_structs[0].branches]

        for bs in band_structs:
            if kpath_ref != [br["name"] for br in bs.branches]:
                msg = (
                    f"BSPlotter only works with BandStructureSymmLine "
                    f"which have the same kpath. \n{bs} has a different kpath!"
                )
                raise ValueError(msg)

        return True

    def add_bs(self, bs: BandStructureSymmLine | list[BandStructureSymmLine]) -> None:
        """Add bands objects to the BSPlotter."""
        if not isinstance(bs, list):
            bs = [bs]

        if self._check_bs_kpath(bs):
            self._bs.extend(bs)
            # TODO: come with an intelligent way to cut the highest unconverged
            # bands
            self._nb_bands.extend([b.nb_bands for b in bs])

    def _make_ticks(self, ax: plt.Axes) -> plt.Axes:
        """Utility private method to add ticks to a band structure."""
        ticks = self.get_ticks()
        # Sanitize only plot the uniq values
        uniq_d = []
        uniq_l = []
        temp_ticks = list(zip(ticks["distance"], ticks["label"], strict=True))
        for idx, t in enumerate(temp_ticks):
            if idx == 0:
                uniq_d.append(t[0])
                uniq_l.append(t[1])
                logger.debug(f"Adding label {t[0]} at {t[1]}")
            elif t[1] == temp_ticks[idx - 1][1]:
                logger.debug(f"Skipping label {t[1]}")
            else:
                logger.debug(f"Adding label {t[0]} at {t[1]}")
                uniq_d.append(t[0])
                uniq_l.append(t[1])

        logger.debug(f"Unique labels are {list(zip(uniq_d, uniq_l, strict=True))}")
        ax.set_xticks(uniq_d)
        ax.set_xticklabels(uniq_l)

        for idx, label in enumerate(ticks["label"]):
            if label is not None:
                # don't print the same label twice
                if idx != 0:
                    if label == ticks["label"][idx - 1]:
                        logger.debug(f"already print label... skipping label {ticks['label'][idx]}")
                    else:
                        logger.debug(f"Adding a line at {ticks['distance'][idx]} for label {ticks['label'][idx]}")
                        ax.axvline(ticks["distance"][idx], color="k")
                else:
                    logger.debug(f"Adding a line at {ticks['distance'][idx]} for label {ticks['label'][idx]}")
                    ax.axvline(ticks["distance"][idx], color="k")
        return ax

    @staticmethod
    def _get_branch_steps(branches):
        """Find discontinuous branches."""
        steps = [0]
        for b1, b2 in itertools.pairwise(branches):
            if b2["name"].split("-")[0] != b1["name"].split("-")[-1]:
                steps.append(b2["start_index"])
        steps.append(branches[-1]["end_index"] + 1)
        return steps

    @staticmethod
    def _rescale_distances(bs_ref, bs):
        """Rescale distances of bs to distances in bs_ref.
        This is used for plotting two bandstructures (same k-path)
        of different materials.
        """
        scaled_distances = []

        for br, br2 in zip(bs_ref.branches, bs.branches, strict=True):
            start = br["start_index"]
            end = br["end_index"]
            max_d = bs_ref.distance[end]
            min_d = bs_ref.distance[start]
            s2 = br2["start_index"]
            e2 = br2["end_index"]
            np = e2 - s2
            if np == 0:
                # it deals with single point branches
                scaled_distances.extend([min_d])
            else:
                scaled_distances.extend([(max_d - min_d) / np * i + min_d for i in range(np + 1)])

        return scaled_distances

    def bs_plot_data(self, zero_to_efermi=True, bs=None, bs_ref=None, split_branches=True):
        """Get the data nicely formatted for a plot.

        Args:
            zero_to_efermi: Automatically set the Fermi level as the plot's origin (i.e. subtract E_f).
                Defaults to True.
            bs: the bandstructure to get the data from. If not provided, the first
                one in the self._bs list will be used.
            bs_ref: is the bandstructure of reference when a rescale of the distances
                is need to plot multiple bands
            split_branches: if True distances and energies are split according to the
                branches. If False distances and energies are split only where branches
                are discontinuous (reducing the number of lines to plot).

        Returns:
            dict: A dictionary of the following format:
            ticks: A dict with the 'distances' at which there is a kpoint (the
            x axis) and the labels (None if no label).
            energy: A dict storing bands for spin up and spin down data
            {Spin:[np.array(nb_bands,kpoints),...]} as a list of discontinuous kpath
            of energies. The energy of multiple continuous branches are stored together.
            vbm: A list of tuples (distance,energy) marking the vbms. The
            energies are shifted with respect to the Fermi level is the
            option has been selected.
            cbm: A list of tuples (distance,energy) marking the cbms. The
            energies are shifted with respect to the Fermi level is the
            option has been selected.
            lattice: The reciprocal lattice.
            zero_energy: This is the energy used as zero for the plot.
            band_gap:A string indicating the band gap and its nature (empty if
            it's a metal).
            is_metal: True if the band structure is metallic (i.e., there is at
            least one band crossing the Fermi level).
        """
        if bs is None:
            # if: BSPlotter, else: BSPlotterProjected
            bs = self._bs[0] if isinstance(self._bs, list) else self._bs

        energies = {str(sp): [] for sp in bs.bands}

        bs_is_metal = bs.is_metal()

        if not bs_is_metal:
            vbm = bs.get_vbm()
            cbm = bs.get_cbm()
        else:
            vbm = {}
            cbm = {}

        zero_energy = 0.0
        if zero_to_efermi:
            zero_energy = bs.efermi if bs_is_metal else vbm["energy"]

        # rescale distances when a bs_ref is given as reference,
        # and when bs and bs_ref have different points in branches.
        # Usually bs_ref is the first one in self._bs list is bs_ref
        distances = bs.distance
        if bs_ref is not None and bs_ref.branches != bs.branches:
            distances = self._rescale_distances(bs_ref, bs)

        if split_branches:
            steps = [br["end_index"] + 1 for br in bs.branches][:-1]
        else:
            # join all the continuous branches
            # to reduce the total number of branches to plot
            steps = self._get_branch_steps(bs.branches)[1:-1]

        distances = np.split(distances, steps)
        for sp in bs.bands:
            energies[str(sp)] = np.hsplit(bs.bands[sp] - zero_energy, steps)

        ticks = self.get_ticks()

        vbm_plot = []
        cbm_plot = []
        bg_str = ""

        if not bs_is_metal:
            for index in cbm["kpoint_index"]:
                cbm_plot.append(
                    (
                        bs.distance[index],
                        (cbm["energy"] - zero_energy if zero_to_efermi else cbm["energy"]),
                    )
                )

            for index in vbm["kpoint_index"]:
                vbm_plot.append(
                    (
                        bs.distance[index],
                        (vbm["energy"] - zero_energy if zero_to_efermi else vbm["energy"]),
                    )
                )

            bg = bs.get_band_gap()
            direct = "Indirect"
            if bg["direct"]:
                direct = "Direct"

            bg_str = f"{direct} {bg['transition']} bandgap = {bg['energy']}"

        return {
            "ticks": ticks,
            "distances": distances,
            "energy": energies,
            "vbm": vbm_plot,
            "cbm": cbm_plot,
            "lattice": bs.lattice_rec.as_dict(),
            "zero_energy": zero_energy,
            "is_metal": bs_is_metal,
            "band_gap": bg_str,
        }

    @staticmethod
    def _interpolate_bands(distances, energies, smooth_tol=0, smooth_k=3, smooth_np=100):
        """Interpolate the provided energies using B-splines as implemented in scipy.interpolate.
        Distances and energies has to provided already split into pieces (branches work good,
        for longer segments the interpolation may fail).

        Interpolation failure can be caused by trying to fit an entire
        band with one spline rather than fitting with piecewise splines
        (splines are ill-suited to fit discontinuities).

        The number of splines used to fit a band is determined by the
        number of branches (high symmetry lines) defined in the
        BandStructureSymmLine object (see BandStructureSymmLine._branches).
        """
        int_energies, int_distances = [], []
        smooth_k_orig = smooth_k

        for dist, ene in zip(distances, energies, strict=True):
            br_en = []
            warning_nan = (
                f"WARNING! Distance / branch, band cannot be "
                f"interpolated. See full warning in source. "
                f"If this is not a mistake, try increasing "
                f"smooth_tol. Current {smooth_tol=}."
            )

            warning_m_fewer_k = (
                f"The number of points (m) has to be higher then "
                f"the order (k) of the splines. In this branch {len(dist)} "
                f"points are found, while k is set to {smooth_k}. "
                f"Smooth_k will be reduced to {smooth_k - 1} for this branch."
            )

            # skip single point branches
            if len(dist) in (2, 3):
                # reducing smooth_k when the number
                # of points are fewer then k
                smooth_k = len(dist) - 1
                warnings.warn(warning_m_fewer_k)
            elif len(dist) == 1:
                warnings.warn("Skipping single point branch")
                continue

            int_distances.append(np.linspace(dist[0], dist[-1], smooth_np))

            for ien in ene:
                tck = scint.splrep(dist, ien, s=smooth_tol, k=smooth_k)

                br_en.append(scint.splev(int_distances[-1], tck))

            smooth_k = smooth_k_orig

            int_energies.append(np.vstack(br_en))

            if np.any(np.isnan(int_energies[-1])):
                warnings.warn(warning_nan)

        return int_distances, int_energies

    def get_plot(
        self,
        zero_to_efermi=True,
        ylim=None,
        smooth=False,
        vbm_cbm_marker=False,
        smooth_tol=0,
        smooth_k=3,
        smooth_np=100,
        bs_labels=None,
    ):
        """Get a matplotlib object for the bandstructures plot.
        Multiple bandstructure objs are plotted together if they have the
        same high symm path.

        Args:
            zero_to_efermi: Automatically set the Fermi level as the plot's origin (i.e. subtract E_f).
                Defaults to True.
            ylim: Specify the y-axis (energy) limits; by default None let
                the code choose. It is vbm-4 and cbm+4 if insulator
                efermi-10 and efermi+10 if metal
            smooth (bool or list(bools)): interpolates the bands by a spline cubic.
                A single bool values means to interpolate all the bandstructure objs.
                A list of bools allows to select the bandstructure obs to interpolate.
            vbm_cbm_marker (bool): if True, a marker is added to the vbm and cbm.
            smooth_tol (float) : tolerance for fitting spline to band data.
                Default is None such that no tolerance will be used.
            smooth_k (int): degree of splines 1<k<5
            smooth_np (int): number of interpolated points per each branch.
            bs_labels: labels for each band for the plot legend.
        """
        ax = pretty_plot(12, 8)

        if isinstance(smooth, bool):
            smooth = [smooth] * len(self._bs)

        handles = []
        vbm_min, cbm_max = [], []
        data = []
        one_is_metal = False

        colors = next(iter(plt.rcParams["axes.prop_cycle"].by_key().values()))
        for ibs, bs in enumerate(self._bs):
            # set first bs in the list as ref for rescaling the distances of the other bands
            bs_ref = self._bs[0] if len(self._bs) > 1 and ibs > 0 else None

            if smooth[ibs]:
                # interpolation works good on short segments like branches
                data = self.bs_plot_data(zero_to_efermi, bs, bs_ref, split_branches=True)
            else:
                data = self.bs_plot_data(zero_to_efermi, bs, bs_ref, split_branches=False)

            # remember if one bs is a metal for setting the ylim later
            one_is_metal = data["is_metal"]

            # remember all the cbm and vbm for setting the ylim later
            if not data["is_metal"]:
                cbm_max.append(data["cbm"][0][1])
                vbm_min.append(data["vbm"][0][1])
            else:
                cbm_max.append(bs.efermi)
                vbm_min.append(bs.efermi)

            for sp in bs.bands:
                ls = "-" if str(sp) == "1" else "--"

                # else case assumes bs_labels is Sequence[str]
                bs_label = f"Band {ibs} {sp.name}" if bs_labels is None else f"{bs_labels[ibs]} {sp.name}"

                handles.append(mlines.Line2D([], [], lw=2, ls=ls, color=colors[ibs], label=bs_label))

                distances, energies = data["distances"], data["energy"][str(sp)]

                if smooth[ibs]:
                    distances, energies = self._interpolate_bands(
                        distances,
                        energies,
                        smooth_tol=smooth_tol,
                        smooth_k=smooth_k,
                        smooth_np=smooth_np,
                    )
                    # join all branches together
                    distances = np.hstack(distances)
                    energies = np.hstack(energies)
                    # split only discontinuous branches
                    steps = self._get_branch_steps(bs.branches)[1:-1]
                    distances = np.split(distances, steps)
                    energies = np.hsplit(energies, steps)

                for dist, ene in zip(distances, energies, strict=True):
                    ax.plot(dist, ene.T, c=colors[ibs], ls=ls)

            # plot markers for vbm and cbm
            if vbm_cbm_marker:
                for cbm in data["cbm"]:
                    ax.scatter(cbm[0], cbm[1], color="r", marker="o", s=100)
                for vbm in data["vbm"]:
                    ax.scatter(vbm[0], vbm[1], color="g", marker="o", s=100)

            # Draw Fermi energy, only if not the zero
            if not zero_to_efermi:
                ef = bs.efermi
                ax.axhline(ef, lw=2, ls="-.", color=colors[ibs])

        # defaults for ylim
        e_min = -4
        e_max = 4
        if one_is_metal:
            e_min = -10
            e_max = 10

        if ylim is None:
            if zero_to_efermi:
                ax.set_ylim(e_min, e_max if one_is_metal else max(cbm_max) + e_max)
            else:
                all_efermi = [band_struct.efermi for band_struct in self._bs]
                ll = min([min(vbm_min), min(all_efermi)])
                hh = max([max(cbm_max), max(all_efermi)])
                ax.set_ylim(ll + e_min, hh + e_max)
        else:
            ax.set_ylim(ylim)

        self._make_ticks(ax)

        # Main X and Y Labels
        ax.set_xlabel(r"$\mathrm{Wave\ Vector}$", fontsize=30)
        ylabel = r"$\mathrm{E\ -\ E_f\ (eV)}$" if zero_to_efermi else r"$\mathrm{Energy\ (eV)}$"
        ax.set_ylabel(ylabel, fontsize=30)

        # X range (K)
        # last distance point
        x_max = data["distances"][-1][-1]
        ax.set_xlim(0, x_max)

        ax.legend(handles=handles)

        plt.tight_layout()

        # auto tight_layout when resizing or pressing t
        def fix_layout(event) -> None:
            if (event.name == "key_press_event" and event.key == "t") or event.name == "resize_event":
                plt.tight_layout()
                plt.gcf().canvas.draw()

        ax.figure.canvas.mpl_connect("key_press_event", fix_layout)
        ax.figure.canvas.mpl_connect("resize_event", fix_layout)

        return ax

    def show(self, zero_to_efermi=True, ylim=None, smooth=False, smooth_tol=None) -> None:
        """Show the plot using matplotlib.

        Args:
            zero_to_efermi: Set the Fermi level as the plot's origin
                (i.e. subtract E_f). Defaults to True.
            ylim: Specify the y-axis (energy) limits; by default None let
                the code choose. It is vbm-4 and cbm+4 if insulator
                efermi-10 and efermi+10 if metal
            smooth: interpolates the bands by a spline cubic
            smooth_tol (float) : tolerance for fitting spline to band data.
                Default is None such that no tolerance will be used.
        """
        self.get_plot(zero_to_efermi, ylim, smooth)
        plt.show()

    def save_plot(self, filename: str, ylim=None, zero_to_efermi=True, smooth=False) -> None:
        """Save matplotlib plot to a file.

        Args:
            filename (str): Filename to write to. Must include extension to specify image format.
            ylim: Specifies the y-axis limits.
            zero_to_efermi: Automatically set the Fermi level as the plot's origin (i.e. subtract E - E_f).
                Defaults to True.
            smooth: Cubic spline interpolation of the bands.
        """
        self.get_plot(ylim=ylim, zero_to_efermi=zero_to_efermi, smooth=smooth)
        plt.savefig(filename)
        plt.close()

    def get_ticks(self):
        """Get all ticks and labels for a band structure plot.

        Returns:
            dict: A dictionary with 'distance': a list of distance at which
            ticks should be set and 'label': a list of label for each of those
            ticks.
        """
        bs = self._bs[0] if isinstance(self._bs, list) else self._bs
        ticks, distance = [], []
        for br in bs.branches:
            start, end = br["start_index"], br["end_index"]
            labels = br["name"].split("-")

            # skip those branches with only one point
            if labels[0] == labels[1]:
                continue

            # add latex $$
            for idx, label in enumerate(labels):
                if label.startswith("\\") or "_" in label:
                    labels[idx] = f"${label}$"

            # If next branch is not continuous,
            # join the first lbl to the previous tick label
            # and add the second lbl to ticks list
            # otherwise add to ticks list both new labels.
            # Similar for distances.
            if ticks and labels[0] != ticks[-1]:
                ticks[-1] += f"$\\mid${labels[0]}"
                ticks.append(labels[1])
                distance.append(bs.distance[end])
            else:
                ticks.extend(labels)
                distance.extend([bs.distance[start], bs.distance[end]])

        return {"distance": distance, "label": ticks}

    def get_ticks_old(self):
        """Get all ticks and labels for a band structure plot.

        Returns:
            dict: A dictionary with 'distance': a list of distance at which
            ticks should be set and 'label': a list of label for each of those
            ticks.
        """
        bs = self._bs[0]
        tick_distance = []
        tick_labels = []
        previous_label = bs.kpoints[0].label
        previous_branch = bs.branches[0]["name"]
        for idx, kpt in enumerate(bs.kpoints):
            if kpt.label is not None:
                tick_distance.append(bs.distance[idx])
                this_branch = None
                for b in bs.branches:
                    if b["start_index"] <= idx <= b["end_index"]:
                        this_branch = b["name"]
                        break
                if kpt.label != previous_label and previous_branch != this_branch:
                    label1 = kpt.label
                    if label1.startswith("\\") or label1.find("_") != -1:
                        label1 = f"${label1}$"
                    label0 = previous_label
                    if label0.startswith("\\") or label0.find("_") != -1:
                        label0 = f"${label0}$"
                    tick_labels.pop()
                    tick_distance.pop()
                    tick_labels.append(f"{label0}$\\mid${label1}")
                elif kpt.label.startswith("\\") or kpt.label.find("_") != -1:
                    tick_labels.append(f"${kpt.label}$")
                else:
                    tick_labels.append(kpt.label)
                previous_label = kpt.label
                previous_branch = this_branch
        return {"distance": tick_distance, "label": tick_labels}

    def plot_compare(self, other_plotter, legend=True) -> plt.Axes:
        """Plot two band structure for comparison. One is in red the other in blue
        (no difference in spins). The two band structures need to be defined
        on the same symmetry lines! and the distance between symmetry lines is
        the one of the band structure used to build the BSPlotter.

        Args:
            other_plotter: Another band structure object defined along the same symmetry lines
            legend: True to add a legend to the plot

        Returns:
            plt.Axes: matplotlib Axes object with both band structures
        """
        warnings.warn("Deprecated method. Use BSPlotter([sbs1,sbs2,...]).get_plot() instead.")

        # TODO: add exception if the band structures are not compatible
        ax = self.get_plot()
        data_orig = self.bs_plot_data()
        data = other_plotter.bs_plot_data()
        band_linewidth = 1
        for i in range(other_plotter._nb_bands):
            for d in range(len(data_orig["distances"])):
                ax.plot(
                    data_orig["distances"][d],
                    [e[str(Spin.up)][i] for e in data["energy"]][d],
                    "c-",
                    linewidth=band_linewidth,
                )
                if other_plotter._bs.is_spin_polarized:
                    ax.plot(
                        data_orig["distances"][d],
                        [e[str(Spin.down)][i] for e in data["energy"]][d],
                        "m--",
                        linewidth=band_linewidth,
                    )
        if legend:
            handles = [
                mlines.Line2D([], [], linewidth=2, color="b", label="bs 1 up"),
                mlines.Line2D([], [], linewidth=2, color="r", label="bs 1 down", linestyle="--"),
                mlines.Line2D([], [], linewidth=2, color="c", label="bs 2 up"),
                mlines.Line2D([], [], linewidth=2, color="m", linestyle="--", label="bs 2 down"),
            ]

            ax.legend(handles=handles)
        return ax

    def plot_brillouin(self):
        """Plot the Brillouin zone.

        Returns:
            plt.Figure: A matplotlib figure object with the Brillouin zone.
        """
        # make labels and lines
        labels = {}
        for k in self._bs[0].kpoints:
            if k.label:
                labels[k.label] = k.frac_coords

        lines = []
        for branch in self._bs[0].branches:
            kpts = self._bs[0].kpoints
            start_idx, end_idx = branch["start_index"], branch["end_index"]
            lines.append([kpts[start_idx].frac_coords, kpts[end_idx].frac_coords])

        return plot_brillouin_zone(self._bs[0].lattice_rec, lines=lines, labels=labels)


class BSPlotterProjected(BSPlotter):
    """Plot or get data to facilitate plotting of projected
    band structure along orbitals, elements or sites.
    """

    def __init__(self, bs: BandStructureSymmLine) -> None:
        """
        Args:
            bs: A BandStructureSymmLine object with projections
            e.g. from a VASP calculation.
        """
        if isinstance(bs, list):
            warnings.warn(
                "Multiple band structures are not handled by BSPlotterProjected. "
                "Only the first in the list will be considered"
            )
            bs = bs[0]

        if len(bs.projections) == 0:
            raise ValueError("Can't plot projections on a band structure without projections data")

        self._bs: BandStructureSymmLine = bs
        self._nb_bands: int = bs.nb_bands  # type: ignore[assignment]

    def _get_projections_by_branches(self, project_onto):
        proj = self._bs.get_projections_on_elements_and_orbitals(project_onto)
        proj_br = []
        for branch in self._bs.branches:
            if self._bs.is_spin_polarized:
                proj_br.append(
                    {
                        str(Spin.up): [[] for _ in range(self._nb_bands)],
                        str(Spin.down): [[] for _ in range(self._nb_bands)],
                    }
                )
            else:
                proj_br.append({str(Spin.up): [[] for _ in range(self._nb_bands)]})

            for band_idx in range(self._nb_bands):
                for j in range(branch["start_index"], branch["end_index"] + 1):
                    proj_br[-1][str(Spin.up)][band_idx].append(
                        {
                            e: {o: proj[Spin.up][band_idx][j][e][o] for o in proj[Spin.up][band_idx][j][e]}
                            for e in proj[Spin.up][band_idx][j]
                        }
                    )
            if self._bs.is_spin_polarized:
                for branch in self._bs.branches:
                    for band_idx in range(self._nb_bands):
                        for j in range(branch["start_index"], branch["end_index"] + 1):
                            proj_br[-1][str(Spin.down)][band_idx].append(
                                {
                                    e: {o: proj[Spin.down][band_idx][j][e][o] for o in proj[Spin.down][band_idx][j][e]}
                                    for e in proj[Spin.down][band_idx][j]
                                }
                            )
        return proj_br

    def get_projected_plots_dots(
        self,
        dictio: dict[str, list],
        zero_to_efermi: bool = True,
        ylim: tuple[float, float] | None = None,
        vbm_cbm_marker: bool = False,
        band_linewidth: float = 1.0,
        marker_size: float = 15.0,
    ) -> plt.Axes:
        """Generate a plot with subplots for each element-orbital pair.

        The orbitals are named as in the FATBAND file, e.g. "2p" or "2p_x".

        he blue and red colors are for spin up and spin down
        The size of the dot in the plot corresponds to the value
        for the specific point.

        Args:
            dictio: The element and orbitals you want a projection on. The
                format is {Element: [*Orbitals]} for instance
                {"Cu":["d", "s"], "O":["p"]} will yield projections for
                Cu on d and s orbitals and oxygen on p.
            zero_to_efermi: Set the Fermi level as the plot's origin
                (i.e. subtract E_f). Defaults to True.
            ylim: The y-axis limits. Defaults to None.
            vbm_cbm_marker (bool): Add markers for the VBM and CBM. Defaults to False.
            band_linewidth (float): The width of the lines. Defaults to 1.0.
            marker_size (float): The size of the markers. Defaults to 15.0.

        Returns:
            plt.Axes
        """
        n_rows = max(map(len, dictio.values()))  # largest number of orbitals for an element
        n_cols = len(dictio)  # number of elements
        proj = self._get_projections_by_branches(dictio)
        data = self.bs_plot_data(zero_to_efermi)
        e_min, e_max = (-10, 10) if self._bs.is_metal() else (-4, 4)

        fig, axs = plt.subplots(n_rows, n_cols, figsize=(12, 8), constrained_layout=True)

        for col_idx, element in enumerate(dictio):
            for row_idx in range(n_rows):
                if n_cols == 1 and n_cols == 1:
                    ax = axs
                elif n_rows == 1:
                    ax = axs[col_idx]
                else:
                    ax = axs[row_idx, col_idx]

                # Skip empty orbitals
                try:
                    orbital = dictio[element][row_idx]
                except IndexError:
                    ax.set_visible(False)
                    continue

                self._make_ticks(ax)

                # Set title (with orbital name as subscript)
                ax.set_title(rf"${{\mathrm{{{element}}}}}_{{\mathrm{{{orbital}}}}}$", fontsize=18)

                # Walk through high symmetry points of the band structure
                # (Gamma->X, X->M, ...)
                for k_path_idx in range(len(data["distances"])):
                    for band_idx in range(self._nb_bands):
                        ax.plot(
                            data["distances"][k_path_idx],
                            data["energy"][str(Spin.up)][k_path_idx][band_idx],
                            "b-",
                            linewidth=band_linewidth,
                        )
                        for j in range(len(data["energy"][str(Spin.up)][k_path_idx][band_idx])):
                            ax.plot(
                                data["distances"][k_path_idx][j],
                                data["energy"][str(Spin.up)][k_path_idx][band_idx][j],
                                "bo",
                                markersize=proj[k_path_idx][str(Spin.up)][band_idx][j][str(element)][orbital]
                                * marker_size,
                            )

                        # Plot spin-down if spin polarized
                        if self._bs.is_spin_polarized:
                            ax.plot(
                                data["distances"][k_path_idx],
                                data["energy"][str(Spin.down)][k_path_idx][band_idx],
                                "r--",
                                linewidth=band_linewidth,
                            )
                            for j in range(len(data["energy"][str(Spin.up)][k_path_idx][band_idx])):
                                ax.plot(
                                    data["distances"][k_path_idx][j],
                                    data["energy"][str(Spin.down)][k_path_idx][band_idx][j],
                                    "ro",
                                    markersize=proj[k_path_idx][str(Spin.down)][band_idx][j][str(element)][orbital]
                                    * marker_size,
                                )

                # Set x range
                ax.set_xlim(np.min(data["distances"]), np.max(data["distances"]))

                # Set y range
                if ylim is None:
                    if self._bs.is_metal():
                        if zero_to_efermi:
                            ax.set_ylim(e_min, e_max)
                        else:
                            ax.set_ylim(self._bs.efermi + e_min, self._bs.efermi + e_max)
                    else:
                        if vbm_cbm_marker:
                            for cbm in data["cbm"]:
                                ax.scatter(cbm[0], cbm[1], color="r", marker="o", s=100)

                            for vbm in data["vbm"]:
                                ax.scatter(vbm[0], vbm[1], color="g", marker="o", s=100)

                        ax.set_ylim(data["vbm"][0][1] + e_min, data["cbm"][0][1] + e_max)
                else:
                    ax.set_ylim(*ylim)

        return fig.axes

    def get_elt_projected_plots(
        self,
        zero_to_efermi: bool = True,
        ylim: tuple[float, float] | None = None,
        vbm_cbm_marker: bool = False,
        band_linewidth: float = 1.0,
    ) -> plt.Axes:
        """Generate a plot with subplots for different elements.

        The blue and red colors are for spin up and spin down
        The size of the dot in the plot corresponds to the value
        for the specific point.

        Returns:
            np.ndarray[plt.Axes]: 2x2 array of plt.Axes with different
                subplots for each projection.
        """
        if self._bs.structure is None:
            raise RuntimeError("Band structure cannot be None.")
        proj = self._get_projections_by_branches({el.symbol: ["s", "p", "d"] for el in self._bs.structure.elements})

        data = self.bs_plot_data(zero_to_efermi)
        _fig, axs = plt.subplots(2, 2, figsize=(12, 8))  # Adjust the layout as needed
        e_min, e_max = (-10, 10) if self._bs.is_metal() else (-4, 4)

        for idx, el in enumerate(self._bs.structure.elements):
            ax = pretty_plot(12, 8, ax=axs.flat[idx])

            self._make_ticks(ax)
            ax.set_title(str(el))

            for b in range(len(data["distances"])):
                for band_idx in range(self._nb_bands):
                    ax.plot(
                        data["distances"][b],
                        data["energy"][str(Spin.up)][b][band_idx],
                        "-",
                        color=[192 / 255, 192 / 255, 192 / 255],
                        linewidth=band_linewidth,
                    )
                    if self._bs.is_spin_polarized:
                        ax.plot(
                            data["distances"][b],
                            data["energy"][str(Spin.down)][b][band_idx],
                            "--",
                            color=[128 / 255, 128 / 255, 128 / 255],
                            linewidth=band_linewidth,
                        )
                        for j in range(len(data["energy"][str(Spin.up)][b][band_idx])):
                            markerscale = sum(
                                proj[b][str(Spin.down)][band_idx][j][str(el)][o]
                                for o in proj[b][str(Spin.down)][band_idx][j][str(el)]
                            )
                            ax.plot(
                                data["distances"][b][j],
                                data["energy"][str(Spin.down)][b][band_idx][j],
                                "bo",
                                markersize=markerscale * 15.0,
                                color=[
                                    markerscale,
                                    0.3 * markerscale,
                                    0.4 * markerscale,
                                ],
                            )
                    for j in range(len(data["energy"][str(Spin.up)][b][band_idx])):
                        markerscale = sum(
                            proj[b][str(Spin.up)][band_idx][j][str(el)][o]
                            for o in proj[b][str(Spin.up)][band_idx][j][str(el)]
                        )
                        ax.plot(
                            data["distances"][b][j],
                            data["energy"][str(Spin.up)][b][band_idx][j],
                            "o",
                            markersize=markerscale * 15.0,
                            color=[markerscale, 0.3 * markerscale, 0.4 * markerscale],
                        )

            # Set ylim
            if ylim is None:
                if self._bs.is_metal():
                    if zero_to_efermi:
                        ax.set_ylim(e_min, e_max)
                    else:
                        ax.set_ylim(self._bs.efermi + e_min, self._bs.efermi + e_max)
                else:
                    if vbm_cbm_marker:
                        for cbm in data["cbm"]:
                            ax.scatter(cbm[0], cbm[1], color="r", marker="o", s=100)

                        for vbm in data["vbm"]:
                            ax.scatter(vbm[0], vbm[1], color="g", marker="o", s=100)

                    ax.set_ylim(data["vbm"][0][1] + e_min, data["cbm"][0][1] + e_max)
            else:
                ax.set_ylim(ylim)

        return axs

    def get_elt_projected_plots_color(
        self,
        zero_to_efermi: bool = True,
        elt_ordered: list | None = None,
        band_linewidth: float = 3,
    ) -> plt.Axes:
        """Generate a pyplot plot where the band structure
        line color depends on the element of the band. where each
        element is associated with red, green or blue.

        The method can only deal with binary and ternary compounds.

        Spin up and spin down are differentiated by a '-' and a '--' line.

        Args:
            zero_to_efermi: set the Fermi level as the plot's origin
                (i.e. subtract E_f). Defaults to True.
            elt_ordered: A list of ordered Elements.
                The first one is red, second green, last blue.
            band_linewidth (float): width of the line.

        Raises:
            RuntimeError: if the band structure is None.
            ValueError: if the number of elements is not 2 or 3.

        Returns:
            a pyplot object
        """
        if self._bs.structure is None:
            raise RuntimeError("Band structure cannot be None.")

        n_elems = len(self._bs.structure.elements)
        if n_elems > 3:
            raise ValueError(f"Can only plot binary and ternary compounds, got {n_elems} elements")
        if elt_ordered is None:
            elt_ordered = self._bs.structure.elements
        proj = self._get_projections_by_branches({e.symbol: ["s", "p", "d"] for e in self._bs.structure.elements})
        data = self.bs_plot_data(zero_to_efermi)
        ax = pretty_plot(12, 8)

        spins = [Spin.up]
        if self._bs.is_spin_polarized:
            spins = [Spin.up, Spin.down]
        self._make_ticks(ax)
        for spin in spins:
            for b in range(len(data["distances"])):
                for band_idx in range(self._nb_bands):
                    for j in range(len(data["energy"][str(spin)][b][band_idx]) - 1):
                        sum_e = 0.0
                        for el in elt_ordered:
                            sum_e += sum(
                                proj[b][str(spin)][band_idx][j][str(el)][o]
                                for o in proj[b][str(spin)][band_idx][j][str(el)]
                            )
                        if sum_e == 0.0:
                            color = [0.0] * len(elt_ordered)
                        else:
                            color = [
                                sum(
                                    proj[b][str(spin)][band_idx][j][str(el)][o]
                                    for o in proj[b][str(spin)][band_idx][j][str(el)]
                                )
                                / sum_e
                                for el in elt_ordered
                            ]
                        if len(color) == 2:
                            color.append(0.0)
                            color[2] = color[1]
                            color[1] = 0.0

                        sign = "--" if spin == Spin.down else "-"

                        ax.plot(
                            [data["distances"][b][j], data["distances"][b][j + 1]],
                            [
                                data["energy"][str(spin)][b][band_idx][j],
                                data["energy"][str(spin)][b][band_idx][j + 1],
                            ],
                            sign,
                            color=color,
                            linewidth=band_linewidth,
                        )

        if self._bs.is_metal():
            if zero_to_efermi:
                e_min = -10
                e_max = 10
                ax.set_ylim(e_min, e_max)
                ax.set_ylim(self._bs.efermi + e_min, self._bs.efermi + e_max)
        else:
            ax.set_ylim(data["vbm"][0][1] - 4.0, data["cbm"][0][1] + 2.0)
        # https://github.com/materialsproject/pymatgen/issues/562
        x_max = data["distances"][-1][-1]
        ax.set_xlim(0, x_max)
        return ax

    def _get_projections_by_branches_patom_pmorb(self, dictio, dictpa, sum_atoms, sum_morbs, selected_branches):
        setos = {
            "s": 0,
            "py": 1,
            "pz": 2,
            "px": 3,
            "dxy": 4,
            "dyz": 5,
            "dz2": 6,
            "dxz": 7,
            "dx2": 8,
            "f_3": 9,
            "f_2": 10,
            "f_1": 11,
            "f0": 12,
            "f1": 13,
            "f2": 14,
            "f3": 15,
        }

        n_branches = len(self._bs.branches)
        if selected_branches is not None:
            indices = []
            if not isinstance(selected_branches, list):
                raise TypeError("You do not give a correct type of 'selected_branches'. It should be 'list' type.")
            if len(selected_branches) == 0:
                raise ValueError("The 'selected_branches' is empty. We cannot do anything.")
            for index in selected_branches:
                if not isinstance(index, int):
                    raise TypeError(
                        "You do not give a correct type of index of symmetry lines. It should be 'int' type"
                    )
                if index > n_branches or index < 1:
                    raise ValueError(
                        f"You give a incorrect index of symmetry lines: {index}. The index should be in range of "
                        f"[1, {n_branches}]."
                    )
                indices.append(index - 1)
        else:
            indices = range(n_branches)

        proj = self._bs.projections
        proj_br = []
        for index in indices:
            b = self._bs.branches[index]
            print(b)
            if self._bs.is_spin_polarized:
                proj_br.append(
                    {
                        str(Spin.up): [[] for _ in range(self._nb_bands)],
                        str(Spin.down): [[] for _ in range(self._nb_bands)],
                    }
                )
            else:
                proj_br.append({str(Spin.up): [[] for _ in range(self._nb_bands)]})

            for band_idx in range(self._nb_bands):
                for j in range(b["start_index"], b["end_index"] + 1):
                    edict = {}
                    for elt in dictpa:
                        for anum in dictpa[elt]:
                            edict[f"{elt}{anum}"] = {}
                            for morb in dictio[elt]:
                                edict[f"{elt}{anum}"][morb] = proj[Spin.up][band_idx][j][setos[morb]][anum - 1]
                    proj_br[-1][str(Spin.up)][band_idx].append(edict)

            if self._bs.is_spin_polarized:
                for band_idx in range(self._nb_bands):
                    for j in range(b["start_index"], b["end_index"] + 1):
                        edict = {}
                        for elt in dictpa:
                            for anum in dictpa[elt]:
                                edict[f"{elt}{anum}"] = {}
                                for morb in dictio[elt]:
                                    edict[f"{elt}{anum}"][morb] = proj[Spin.up][band_idx][j][setos[morb]][anum - 1]
                        proj_br[-1][str(Spin.down)][band_idx].append(edict)

        # Adjust projections for plot
        dictio_d, dictpa_d = self._summarize_keys_for_plot(dictio, dictpa, sum_atoms, sum_morbs)

        if (sum_atoms is None) and (sum_morbs is None):
            proj_br_d = copy.deepcopy(proj_br)
        else:
            proj_br_d = []
            branch = -1
            for index in indices:
                branch += 1
                br = self._bs.branches[index]
                if self._bs.is_spin_polarized:
                    proj_br_d.append(
                        {
                            str(Spin.up): [[] for _ in range(self._nb_bands)],
                            str(Spin.down): [[] for _ in range(self._nb_bands)],
                        }
                    )
                else:
                    proj_br_d.append({str(Spin.up): [[] for _ in range(self._nb_bands)]})

                if (sum_atoms is not None) and (sum_morbs is None):
                    for band_idx in range(self._nb_bands):
                        for j in range(br["end_index"] - br["start_index"] + 1):
                            atoms_morbs = copy.deepcopy(proj_br[branch][str(Spin.up)][band_idx][j])
                            edict = {}
                            for elt in dictpa:
                                if elt in sum_atoms:
                                    for anum in dictpa_d[elt][:-1]:
                                        edict[elt + anum] = copy.deepcopy(atoms_morbs[elt + anum])
                                    edict[elt + dictpa_d[elt][-1]] = {}
                                    for morb in dictio[elt]:
                                        sprojection = 0.0
                                        for anum in sum_atoms[elt]:
                                            sprojection += atoms_morbs[f"{elt}{anum}"][morb]
                                        edict[elt + dictpa_d[elt][-1]][morb] = sprojection
                                else:
                                    for anum in dictpa_d[elt]:
                                        edict[elt + anum] = copy.deepcopy(atoms_morbs[elt + anum])
                            proj_br_d[-1][str(Spin.up)][band_idx].append(edict)
                    if self._bs.is_spin_polarized:
                        for band_idx in range(self._nb_bands):
                            for j in range(br["end_index"] - br["start_index"] + 1):
                                atoms_morbs = copy.deepcopy(proj_br[branch][str(Spin.down)][band_idx][j])
                                edict = {}
                                for elt in dictpa:
                                    if elt in sum_atoms:
                                        for anum in dictpa_d[elt][:-1]:
                                            edict[elt + anum] = copy.deepcopy(atoms_morbs[elt + anum])
                                        edict[elt + dictpa_d[elt][-1]] = {}
                                        for morb in dictio[elt]:
                                            sprojection = 0.0
                                            for anum in sum_atoms[elt]:
                                                sprojection += atoms_morbs[f"{elt}{anum}"][morb]
                                            edict[elt + dictpa_d[elt][-1]][morb] = sprojection
                                    else:
                                        for anum in dictpa_d[elt]:
                                            edict[elt + anum] = copy.deepcopy(atoms_morbs[elt + anum])
                                proj_br_d[-1][str(Spin.down)][band_idx].append(edict)

                elif (sum_atoms is None) and (sum_morbs is not None):
                    for band_idx in range(self._nb_bands):
                        for j in range(br["end_index"] - br["start_index"] + 1):
                            atoms_morbs = copy.deepcopy(proj_br[branch][str(Spin.up)][band_idx][j])
                            edict = {}
                            for elt in dictpa:
                                if elt in sum_morbs:
                                    for anum in dictpa_d[elt]:
                                        edict[elt + anum] = {}
                                        for morb in dictio_d[elt][:-1]:
                                            edict[elt + anum][morb] = atoms_morbs[elt + anum][morb]
                                        sprojection = 0.0
                                        for morb in sum_morbs[elt]:
                                            sprojection += atoms_morbs[elt + anum][morb]
                                        edict[elt + anum][dictio_d[elt][-1]] = sprojection
                                else:
                                    for anum in dictpa_d[elt]:
                                        edict[elt + anum] = copy.deepcopy(atoms_morbs[elt + anum])
                            proj_br_d[-1][str(Spin.up)][band_idx].append(edict)
                    if self._bs.is_spin_polarized:
                        for band_idx in range(self._nb_bands):
                            for j in range(br["end_index"] - br["start_index"] + 1):
                                atoms_morbs = copy.deepcopy(proj_br[branch][str(Spin.down)][band_idx][j])
                                edict = {}
                                for elt in dictpa:
                                    if elt in sum_morbs:
                                        for anum in dictpa_d[elt]:
                                            edict[elt + anum] = {}
                                            for morb in dictio_d[elt][:-1]:
                                                edict[elt + anum][morb] = atoms_morbs[elt + anum][morb]
                                            sprojection = 0.0
                                            for morb in sum_morbs[elt]:
                                                sprojection += atoms_morbs[elt + anum][morb]
                                            edict[elt + anum][dictio_d[elt][-1]] = sprojection
                                    else:
                                        for anum in dictpa_d[elt]:
                                            edict[elt + anum] = copy.deepcopy(atoms_morbs[elt + anum])
                                proj_br_d[-1][str(Spin.down)][band_idx].append(edict)

                else:
                    for band_idx in range(self._nb_bands):
                        for j in range(br["end_index"] - br["start_index"] + 1):
                            atoms_morbs = copy.deepcopy(proj_br[branch][str(Spin.up)][band_idx][j])
                            edict = {}
                            for elt in dictpa:
                                if (elt in sum_atoms) and (elt in sum_morbs):
                                    for anum in dictpa_d[elt][:-1]:
                                        edict[elt + anum] = {}
                                        for morb in dictio_d[elt][:-1]:
                                            edict[elt + anum][morb] = atoms_morbs[elt + anum][morb]
                                        sprojection = 0.0
                                        for morb in sum_morbs[elt]:
                                            sprojection += atoms_morbs[elt + anum][morb]
                                        edict[elt + anum][dictio_d[elt][-1]] = sprojection

                                    edict[elt + dictpa_d[elt][-1]] = {}
                                    for morb in dictio_d[elt][:-1]:
                                        sprojection = 0.0
                                        for anum in sum_atoms[elt]:
                                            sprojection += atoms_morbs[f"{elt}{anum}"][morb]
                                        edict[elt + dictpa_d[elt][-1]][morb] = sprojection

                                    sprojection = 0.0
                                    for anum in sum_atoms[elt]:
                                        for morb in sum_morbs[elt]:
                                            sprojection += atoms_morbs[f"{elt}{anum}"][morb]
                                    edict[elt + dictpa_d[elt][-1]][dictio_d[elt][-1]] = sprojection

                                elif (elt in sum_atoms) and (elt not in sum_morbs):
                                    for anum in dictpa_d[elt][:-1]:
                                        edict[elt + anum] = copy.deepcopy(atoms_morbs[elt + anum])
                                    edict[elt + dictpa_d[elt][-1]] = {}
                                    for morb in dictio[elt]:
                                        sprojection = 0.0
                                        for anum in sum_atoms[elt]:
                                            sprojection += atoms_morbs[f"{elt}{anum}"][morb]
                                        edict[elt + dictpa_d[elt][-1]][morb] = sprojection

                                elif (elt not in sum_atoms) and (elt in sum_morbs):
                                    for anum in dictpa_d[elt]:
                                        edict[elt + anum] = {}
                                        for morb in dictio_d[elt][:-1]:
                                            edict[elt + anum][morb] = atoms_morbs[elt + anum][morb]
                                        sprojection = 0.0
                                        for morb in sum_morbs[elt]:
                                            sprojection += atoms_morbs[elt + anum][morb]
                                        edict[elt + anum][dictio_d[elt][-1]] = sprojection

                                else:
                                    for anum in dictpa_d[elt]:
                                        edict[elt + anum] = {}
                                        for morb in dictio_d[elt]:
                                            edict[elt + anum][morb] = atoms_morbs[elt + anum][morb]
                            proj_br_d[-1][str(Spin.up)][band_idx].append(edict)

                    if self._bs.is_spin_polarized:
                        for band_idx in range(self._nb_bands):
                            for j in range(br["end_index"] - br["start_index"] + 1):
                                atoms_morbs = copy.deepcopy(proj_br[branch][str(Spin.down)][band_idx][j])
                                edict = {}
                                for elt in dictpa:
                                    if (elt in sum_atoms) and (elt in sum_morbs):
                                        for anum in dictpa_d[elt][:-1]:
                                            edict[elt + anum] = {}
                                            for morb in dictio_d[elt][:-1]:
                                                edict[elt + anum][morb] = atoms_morbs[elt + anum][morb]
                                            sprojection = 0.0
                                            for morb in sum_morbs[elt]:
                                                sprojection += atoms_morbs[elt + anum][morb]
                                            edict[elt + anum][dictio_d[elt][-1]] = sprojection

                                        edict[elt + dictpa_d[elt][-1]] = {}
                                        for morb in dictio_d[elt][:-1]:
                                            sprojection = 0.0
                                            for anum in sum_atoms[elt]:
                                                sprojection += atoms_morbs[f"{elt}{anum}"][morb]
                                            edict[elt + dictpa_d[elt][-1]][morb] = sprojection

                                        sprojection = 0.0
                                        for anum in sum_atoms[elt]:
                                            for morb in sum_morbs[elt]:
                                                sprojection += atoms_morbs[f"{elt}{anum}"][morb]
                                        edict[elt + dictpa_d[elt][-1]][dictio_d[elt][-1]] = sprojection

                                    elif (elt in sum_atoms) and (elt not in sum_morbs):
                                        for anum in dictpa_d[elt][:-1]:
                                            edict[elt + anum] = copy.deepcopy(atoms_morbs[elt + anum])
                                        edict[elt + dictpa_d[elt][-1]] = {}
                                        for morb in dictio[elt]:
                                            sprojection = 0.0
                                            for anum in sum_atoms[elt]:
                                                sprojection += atoms_morbs[f"{elt}{anum}"][morb]
                                            edict[elt + dictpa_d[elt][-1]][morb] = sprojection

                                    elif (elt not in sum_atoms) and (elt in sum_morbs):
                                        for anum in dictpa_d[elt]:
                                            edict[elt + anum] = {}
                                            for morb in dictio_d[elt][:-1]:
                                                edict[elt + anum][morb] = atoms_morbs[elt + anum][morb]
                                            sprojection = 0.0
                                            for morb in sum_morbs[elt]:
                                                sprojection += atoms_morbs[elt + anum][morb]
                                            edict[elt + anum][dictio_d[elt][-1]] = sprojection

                                    else:
                                        for anum in dictpa_d[elt]:
                                            edict[elt + anum] = {}
                                            for morb in dictio_d[elt]:
                                                edict[elt + anum][morb] = atoms_morbs[elt + anum][morb]
                                proj_br_d[-1][str(Spin.down)][band_idx].append(edict)

        return proj_br_d, dictio_d, dictpa_d, indices

    def get_projected_plots_dots_patom_pmorb(
        self,
        dictio,
        dictpa,
        sum_atoms=None,
        sum_morbs=None,
        zero_to_efermi=True,
        ylim=None,
        vbm_cbm_marker=False,
        selected_branches=None,
        w_h_size=(12, 8),
        num_column=None,
    ):
        """Return a plot composed of subplots for different atoms and orbitals (subshell
        orbitals such as 's', 'p', 'd' and 'f' defined by azimuthal quantum numbers l = 0,
        1, 2 and 3, respectively or individual orbitals like 'px', 'py' and 'pz' defined
        by magnetic quantum numbers m = -1, 1 and 0, respectively). This is an extension
        of "get_projected_plots_dots" method.

        Args:
            dictio: The elements and the orbitals you need to project on. The
                format is {Element:[Orbitals]}, for instance:
                {'Cu':['dxy','s','px'],'O':['px','py','pz']} will give projections for Cu on
                orbitals dxy, s, px and for O on orbitals px, py, pz. If you want to sum over all
                individual orbitals of subshell orbitals, for example, 'px', 'py' and 'pz' of O,
                just simply set {'Cu':['dxy','s','px'],'O':['p']} and set sum_morbs (see
                explanations below) as {'O':[p],...}. Otherwise, you will get an error.
            dictpa: The elements and their sites (defined by site numbers) you
                need to project on. The format is {Element: [Site numbers]}, for instance:
                {'Cu':[1,5],'O':[3,4]} will give projections for Cu on site-1 and on site-5, O on
                site-3 and on site-4 in the cell. The correct site numbers of atoms are consistent
                with themselves in the structure computed. Normally, the structure should be totally
                similar with POSCAR file, however, sometimes VASP can rotate or translate the cell.
                Thus, it would be safe if using Vasprun class to get the final_structure and as a
                result, correct index numbers of atoms.
            sum_atoms: Sum projection of the similar atoms together (e.g.: Cu
                on site-1 and Cu on site-5). The format is {Element: [Site numbers]}, for instance:
                {'Cu': [1,5], 'O': [3,4]} means summing projections over Cu on site-1 and Cu on
                site-5 and O on site-3 and on site-4. If you do not want to use this functional,
                just turn it off by setting sum_atoms = None.
            sum_morbs: Sum projections of individual orbitals of similar atoms
                together (e.g.: 'dxy' and 'dxz'). The format is {Element: [individual orbitals]},
                for instance: {'Cu': ['dxy', 'dxz'], 'O': ['px', 'py']} means summing projections
                over 'dxy' and 'dxz' of Cu and 'px' and 'py' of O. If you do not want to use this
                functional, just turn it off by setting sum_morbs = None.
            zero_to_efermi: Automatically set the Fermi level as the plot's origin (i.e. subtract E_f).
                Defaults to True.
            ylim: The y-axis limit. Defaults to None.
            vbm_cbm_marker: Whether to plot points to indicate valence band maxima and conduction
                band minima positions. Defaults to False.
            selected_branches: The index of symmetry lines you chose for
                plotting. This can be useful when the number of symmetry lines (in KPOINTS file) are
                manny while you only want to show for certain ones. The format is [index of line],
                for instance: [1, 3, 4] means you just need to do projection along lines number 1, 3
                and 4 while neglecting lines number 2 and so on. By default, this is None type and
                all symmetry lines will be plotted.
            w_h_size: This variable help you to control the width and height
                of figure. By default, width = 12 and height = 8 (inches). The width/height ratio is
                kept the same for subfigures and the size of each depends on how many number of
                subfigures are plotted.
            num_column: This variable help you to manage how the subfigures are
                arranged in the figure by setting up the number of columns of subfigures. The value
                should be an int number. For example, num_column = 3 means you want to plot
                subfigures in 3 columns. By default, num_column = None and subfigures are aligned in
                2 columns.

        Returns:
            A pyplot object with different subfigures for different projections.
            The blue and red colors lines are bands
            for spin up and spin down. The green and cyan dots are projections
            for spin up and spin down. The bigger
            the green or cyan dots in the projected band structures, the higher
            character for the corresponding elements
            and orbitals. List of individual orbitals and their numbers (set up
            by VASP and no special meaning):
            s = 0; py = 1 pz = 2 px = 3; dxy = 4 dyz = 5 dz2 = 6 dxz = 7 dx2 = 8;
            f_3 = 9 f_2 = 10 f_1 = 11 f0 = 12 f1 = 13 f2 = 14 f3 = 15
        """
        dictio, sum_morbs = self._Orbitals_SumOrbitals(dictio, sum_morbs)
        dictpa, sum_atoms, n_figs = self._number_of_subfigures(dictio, dictpa, sum_atoms, sum_morbs)
        print(f"Number of subfigures: {n_figs}")
        if n_figs > 9:
            print(
                f"The number of subfigures {n_figs} might be too manny and the implementation might take a long "
                f"time.\n A smaller number or a plot with selected symmetry lines (selected_branches) might be better."
            )

        band_linewidth = 0.5
        ax = pretty_plot(w_h_size[0], w_h_size[1])
        proj_br_d, dictio_d, dictpa_d, branches = self._get_projections_by_branches_patom_pmorb(
            dictio, dictpa, sum_atoms, sum_morbs, selected_branches
        )
        data = self.bs_plot_data(zero_to_efermi)
        e_min = -4
        e_max = 4
        if self._bs.is_metal():
            e_min = -10
            e_max = 10

        count = 0
        for elt in dictpa_d:
            for numa in dictpa_d[elt]:
                for o in dictio_d[elt]:
                    count += 1
                    if num_column is None:
                        if n_figs == 1:
                            plt.subplot(1, 1, 1)
                        else:
                            row = n_figs // 2
                            if n_figs % 2 == 0:
                                plt.subplot(row, 2, count)
                            else:
                                plt.subplot(row + 1, 2, count)
                    elif isinstance(num_column, int):
                        row = n_figs / num_column
                        if n_figs % num_column == 0:
                            plt.subplot(row, num_column, count)
                        else:
                            plt.subplot(row + 1, num_column, count)
                    else:
                        raise ValueError("The invalid 'num_column' is assigned. It should be an integer.")

                    ax, shift = self._make_ticks_selected(ax, branches)
                    br = -1
                    for b in branches:
                        br += 1
                        for band_idx in range(self._nb_bands):
                            ax.plot(
                                [x - shift[br] for x in data["distances"][b]],
                                [
                                    data["energy"][str(Spin.up)][b][band_idx][j]
                                    for j in range(len(data["distances"][b]))
                                ],
                                "b-",
                                linewidth=band_linewidth,
                            )

                            if self._bs.is_spin_polarized:
                                ax.plot(
                                    [x - shift[br] for x in data["distances"][b]],
                                    [
                                        data["energy"][str(Spin.down)][b][band_idx][j]
                                        for j in range(len(data["distances"][b]))
                                    ],
                                    "r--",
                                    linewidth=band_linewidth,
                                )
                                for j in range(len(data["energy"][str(Spin.up)][b][band_idx])):
                                    ax.plot(
                                        data["distances"][b][j] - shift[br],
                                        data["energy"][str(Spin.down)][b][band_idx][j],
                                        "co",
                                        markersize=proj_br_d[br][str(Spin.down)][band_idx][j][elt + numa][o] * 15.0,
                                    )

                            for j in range(len(data["energy"][str(Spin.up)][b][band_idx])):
                                ax.plot(
                                    data["distances"][b][j] - shift[br],
                                    data["energy"][str(Spin.up)][b][band_idx][j],
                                    "go",
                                    markersize=proj_br_d[br][str(Spin.up)][band_idx][j][elt + numa][o] * 15.0,
                                )

                    if ylim is None:
                        if self._bs.is_metal():
                            if zero_to_efermi:
                                ax.set_ylim(e_min, e_max)
                            else:
                                ax.set_ylim(self._bs.efermi + e_min, self._bs._efermi + e_max)
                        else:
                            if vbm_cbm_marker:
                                for cbm in data["cbm"]:
                                    ax.scatter(cbm[0], cbm[1], color="r", marker="o", s=100)

                                for vbm in data["vbm"]:
                                    ax.scatter(vbm[0], vbm[1], color="g", marker="o", s=100)

                            ax.set_ylim(data["vbm"][0][1] + e_min, data["cbm"][0][1] + e_max)
                    else:
                        ax.set_ylim(ylim)
                    ax.set_title(f"{elt} {numa} {o}")

        return ax

    @classmethod
    def _Orbitals_SumOrbitals(cls, dictio, sum_morbs):
        all_orbitals = [
            "s",
            "p",
            "d",
            "f",
            "px",
            "py",
            "pz",
            "dxy",
            "dyz",
            "dxz",
            "dx2",
            "dz2",
            "f_3",
            "f_2",
            "f_1",
            "f0",
            "f1",
            "f2",
            "f3",
        ]
        individual_orbs = {
            "p": ["px", "py", "pz"],
            "d": ["dxy", "dyz", "dxz", "dx2", "dz2"],
            "f": ["f_3", "f_2", "f_1", "f0", "f1", "f2", "f3"],
        }

        if not isinstance(dictio, dict):
            raise TypeError("The invalid type of 'dictio' was bound. It should be dict type.")
        if len(dictio) == 0:
            raise KeyError("The 'dictio' is empty. We cannot do anything.")

        orb = None
        for elt in dictio:
            if Element.is_valid_symbol(elt):
                if isinstance(dictio[elt], list):
                    if len(dictio[elt]) == 0:
                        raise ValueError(f"The dictio[{elt}] is empty. We cannot do anything")
                    for orb in dictio[elt]:
                        if not isinstance(orb, str):
                            raise TypeError(
                                f"The invalid format of orbitals is in 'dictio[{elt}]': {orb}. They should be string."
                            )
                        if orb not in all_orbitals:
                            raise ValueError(f"The invalid name of orbital is given in 'dictio[{elt}]'.")
                        if orb in individual_orbs and len(set(dictio[elt]).intersection(individual_orbs[orb])) != 0:
                            raise ValueError(f"The 'dictio[{elt}]' contains orbitals repeated.")
                    nelems = Counter(dictio[elt]).values()
                    if sum(nelems) > len(nelems):
                        raise ValueError(f"You put in at least two similar orbitals in dictio[{elt}].")
                else:
                    raise TypeError(f"The invalid type of value was put into 'dictio[{elt}]'. It should be list type.")
            else:
                raise KeyError(f"The invalid element was put into 'dictio' as a key: {elt}")

        if sum_morbs is None:
            print("You do not want to sum projection over orbitals.")
        elif not isinstance(sum_morbs, dict):
            raise TypeError("The invalid type of 'sum_orbs' was bound. It should be dict or 'None' type.")
        elif len(sum_morbs) == 0:
            raise KeyError("The 'sum_morbs' is empty. We cannot do anything")
        else:
            for elt in sum_morbs:
                if Element.is_valid_symbol(elt):
                    if isinstance(sum_morbs[elt], list):
                        for orb in sum_morbs[elt]:
                            if not isinstance(orb, str):
                                raise TypeError(
                                    f"The invalid format of orbitals is in 'sum_morbs[{elt}]': {orb}. "
                                    "They should be string."
                                )
                            if orb not in all_orbitals:
                                raise ValueError(f"The invalid name of orbital in 'sum_morbs[{elt}]' is given.")
                            if orb in individual_orbs and len(set(sum_morbs[elt]) & set(individual_orbs[orb])) != 0:
                                raise ValueError(f"The 'sum_morbs[{elt}]' contains orbitals repeated.")
                        nelems = Counter(sum_morbs[elt]).values()
                        if sum(nelems) > len(nelems):
                            raise ValueError(f"You put in at least two similar orbitals in sum_morbs[{elt}].")
                    else:
                        raise TypeError(
                            f"The invalid type of value was put into 'sum_morbs[{elt}]'. It should be list type."
                        )
                    if elt not in dictio:
                        raise ValueError(
                            f"You cannot sum projection over orbitals of atoms {elt!r} because they are not "
                            "mentioned in 'dictio'."
                        )
                else:
                    raise KeyError(f"The invalid element was put into 'sum_morbs' as a key: {elt}")

        for elt in dictio:
            if len(dictio[elt]) == 1:
                if len(dictio[elt][0]) > 1:
                    if elt in sum_morbs:
                        raise ValueError(
                            f"You cannot sum projection over one individual orbital {dictio[elt][0]!r} of {elt!r}."
                        )
                elif sum_morbs is None:
                    pass
                elif elt not in sum_morbs:
                    print(f"You do not want to sum projection over orbitals of element: {elt}")
                else:
                    if len(sum_morbs[elt]) == 0:
                        raise ValueError(f"The empty list is an invalid value for sum_morbs[{elt}].")
                    if len(sum_morbs[elt]) > 1:
                        for orb in sum_morbs[elt]:
                            if dictio[elt][0] not in orb:
                                raise ValueError(f"The invalid orbital {orb!r} was put into 'sum_morbs[{elt}]'.")
                    else:
                        if orb == "s" or len(orb) > 1:
                            raise ValueError(f"The invalid orbital {orb!r} was put into sum_orbs[{elt!r}].")
                        sum_morbs[elt] = individual_orbs[dictio[elt][0]]
                        dictio[elt] = individual_orbs[dictio[elt][0]]
            else:
                duplicate = copy.deepcopy(dictio[elt])
                for orb in dictio[elt]:
                    if orb in individual_orbs:
                        duplicate.remove(orb)
                        duplicate += individual_orbs[orb]
                dictio[elt] = copy.deepcopy(duplicate)

                if sum_morbs is None:
                    pass
                elif elt not in sum_morbs:
                    print(f"You do not want to sum projection over orbitals of element: {elt}")
                else:
                    if len(sum_morbs[elt]) == 0:
                        raise ValueError(f"The empty list is an invalid value for sum_morbs[{elt}].")
                    if len(sum_morbs[elt]) == 1:
                        orb = sum_morbs[elt][0]
                        if orb == "s":
                            raise ValueError(
                                "We do not sum projection over only 's' orbital of the same type of element."
                            )
                        if orb in individual_orbs:
                            sum_morbs[elt].pop(0)
                            sum_morbs[elt] += individual_orbs[orb]
                        else:
                            raise ValueError(f"You never sum projection over one orbital in sum_morbs[{elt}]")
                    else:
                        duplicate = copy.deepcopy(sum_morbs[elt])
                        for orb in sum_morbs[elt]:
                            if orb in individual_orbs:
                                duplicate.remove(orb)
                                duplicate += individual_orbs[orb]
                        sum_morbs[elt] = copy.deepcopy(duplicate)

                    for orb in sum_morbs[elt]:
                        if orb not in dictio[elt]:
                            raise ValueError(f"The orbitals of sum_morbs[{elt}] conflict with those of dictio[{elt}].")

        return dictio, sum_morbs

    def _number_of_subfigures(self, dictio, dictpa, sum_atoms, sum_morbs):
        if not isinstance(dictpa, dict):
            raise TypeError("The invalid type of 'dictpa' was bound. It should be dict type.")
        if len(dictpa) == 0:
            raise KeyError("The 'dictpa' is empty. We cannot do anything.")
        for elt in dictpa:
            if Element.is_valid_symbol(elt):
                if isinstance(dictpa[elt], list):
                    if len(dictpa[elt]) == 0:
                        raise ValueError(f"The dictpa[{elt}] is empty. We cannot do anything")
                    _sites = self._bs.structure.sites
                    indices = []
                    for site_idx, site in enumerate(_sites):
                        if next(iter(site._species)) == Element(elt):
                            indices.append(site_idx + 1)
                    for number in dictpa[elt]:
                        if isinstance(number, str):
                            if number.lower() == "all":
                                dictpa[elt] = indices
                                print(f"You want to consider all {elt!r} atoms.")
                                break

                            raise ValueError(f"You put wrong site numbers in 'dictpa[{elt}]': {number}.")
                        if isinstance(number, int):
                            if number not in indices:
                                raise ValueError(f"You put wrong site numbers in 'dictpa[{elt}]': {number}.")
                        else:
                            raise TypeError(f"You put wrong site numbers in 'dictpa[{elt}]': {number}.")
                    nelems = Counter(dictpa[elt]).values()
                    if sum(nelems) > len(nelems):
                        raise ValueError(f"You put at least two similar site numbers into 'dictpa[{elt}]'.")
                else:
                    raise TypeError(f"The invalid type of value was put into 'dictpa[{elt}]'. It should be list type.")
            else:
                raise KeyError(f"The invalid element was put into 'dictpa' as a key: {elt}")

        if len(list(dictio)) != len(list(dictpa)):
            raise KeyError("The number of keys in 'dictio' and 'dictpa' are not the same.")
        for elt in dictio:
            if elt not in dictpa:
                raise KeyError(f"The element {elt!r} is not in both dictpa and dictio.")
        for elt in dictpa:
            if elt not in dictio:
                raise KeyError(f"The element {elt!r} in not in both dictpa and dictio.")

        if sum_atoms is None:
            print("You do not want to sum projection over atoms.")
        elif not isinstance(sum_atoms, dict):
            raise TypeError("The invalid type of 'sum_atoms' was bound. It should be dict type.")
        elif len(sum_atoms) == 0:
            raise KeyError("The 'sum_atoms' is empty. We cannot do anything.")
        else:
            for elt in sum_atoms:
                if Element.is_valid_symbol(elt):
                    if isinstance(sum_atoms[elt], list):
                        if len(sum_atoms[elt]) == 0:
                            raise ValueError(f"The sum_atoms[{elt}] is empty. We cannot do anything")
                        _sites = self._bs.structure.sites
                        indices = []
                        for site_idx, site in enumerate(_sites):
                            if next(iter(site._species)) == Element(elt):
                                indices.append(site_idx + 1)
                        for number in sum_atoms[elt]:
                            if isinstance(number, str):
                                if number.lower() == "all":
                                    sum_atoms[elt] = indices
                                    print(f"You want to sum projection over all {elt!r} atoms.")
                                    break
                                raise ValueError(f"You put wrong site numbers in 'sum_atoms[{elt}]'.")
                            if isinstance(number, int):
                                if number not in indices:
                                    raise ValueError(f"You put wrong site numbers in 'sum_atoms[{elt}]'.")
                                if number not in dictpa[elt]:
                                    raise ValueError(
                                        f"You cannot sum projection with atom number {number!r} because it is not "
                                        f"mentioned in dicpta[{elt}]"
                                    )
                            else:
                                raise TypeError(f"You put wrong site numbers in 'sum_atoms[{elt}]'.")

                        nelems = Counter(sum_atoms[elt]).values()
                        if sum(nelems) > len(nelems):
                            raise ValueError(f"You put at least two similar site numbers into 'sum_atoms[{elt}]'.")
                    else:
                        raise TypeError(
                            f"The invalid type of value was put into 'sum_atoms[{elt}]'. It should be list type."
                        )
                    if elt not in dictpa:
                        raise ValueError(
                            f"You cannot sum projection over atoms {elt!r} because it is not mentioned in 'dictio'."
                        )
                else:
                    raise KeyError(f"The invalid element was put into 'sum_atoms' as a key: {elt}")
                if len(sum_atoms[elt]) == 1:
                    raise ValueError(f"We do not sum projection over only one atom: {elt}")

        max_number_figs = decrease = 0
        for elt in dictio:
            max_number_figs += len(dictio[elt]) * len(dictpa[elt])

        if (sum_atoms is None) and (sum_morbs is None):
            number_figs = max_number_figs
        elif (sum_atoms is not None) and (sum_morbs is None):
            for elt in sum_atoms:
                decrease += (len(sum_atoms[elt]) - 1) * len(dictio[elt])
            number_figs = max_number_figs - decrease
        elif (sum_atoms is None) and (sum_morbs is not None):
            for elt in sum_morbs:
                decrease += (len(sum_morbs[elt]) - 1) * len(dictpa[elt])
            number_figs = max_number_figs - decrease
        elif (sum_atoms is not None) and (sum_morbs is not None):
            for elt in sum_atoms:
                decrease += (len(sum_atoms[elt]) - 1) * len(dictio[elt])
            for elt in sum_morbs:
                if elt in sum_atoms:
                    decrease += (len(sum_morbs[elt]) - 1) * (len(dictpa[elt]) - len(sum_atoms[elt]) + 1)
                else:
                    decrease += (len(sum_morbs[elt]) - 1) * len(dictpa[elt])
            number_figs = max_number_figs - decrease
        else:
            raise ValueError("Invalid format of 'sum_atoms' and 'sum_morbs'.")

        return dictpa, sum_atoms, number_figs

    def _summarize_keys_for_plot(self, dictio, dictpa, sum_atoms, sum_morbs):
        individual_orbs = {
            "p": ["px", "py", "pz"],
            "d": ["dxy", "dyz", "dxz", "dx2", "dz2"],
            "f": ["f_3", "f_2", "f_1", "f0", "f1", "f2", "f3"],
        }

        def number_label(list_numbers):
            list_numbers = sorted(list_numbers)
            divide = [[]]
            divide[0].append(list_numbers[0])
            group = 0
            for idx in range(1, len(list_numbers)):
                if list_numbers[idx] == list_numbers[idx - 1] + 1:
                    divide[group].append(list_numbers[idx])
                else:
                    group += 1
                    divide.append([list_numbers[idx]])
            label = ""
            for elem in divide:
                if len(elem) > 1:
                    label += f"{elem[0]}-{elem[-1]},"
                else:
                    label += f"{elem[0]},"
            return label[:-1]

        def orbital_label(list_orbitals):
            divide = {}
            for orb in list_orbitals:
                if orb[0] in divide:
                    divide[orb[0]].append(orb)
                else:
                    divide[orb[0]] = []
                    divide[orb[0]].append(orb)
            label = ""
            for elem, orbs in divide.items():
                if elem == "s":
                    label += "s,"
                elif len(orbs) == len(individual_orbs[elem]):
                    label += f"{elem},"
                else:
                    orb_label = [orb[1:] for orb in orbs]
                    label += f"{elem}{str(orb_label).replace('[', '').replace(']', '').replace(', ', '-')},"
            return label[:-1]

        if sum_atoms is None and sum_morbs is None:
            dictio_d = dictio
            dictpa_d = {elt: [str(anum) for anum in dictpa[elt]] for elt in dictpa}

        elif (sum_atoms is not None) and (sum_morbs is None):
            dictio_d = dictio
            dictpa_d = {}
            for elt in dictpa:
                dictpa_d[elt] = []
                if elt in sum_atoms:
                    _sites = self._bs.structure.sites
                    indices = []
                    for site_idx in range(len(_sites)):
                        if next(iter(_sites[site_idx]._species)) == Element(elt):
                            indices.append(site_idx + 1)
                    flag_1 = len(set(dictpa[elt]).intersection(indices))
                    flag_2 = len(set(sum_atoms[elt]).intersection(indices))
                    if flag_1 == len(indices) and flag_2 == len(indices):
                        dictpa_d[elt].append("all")
                    else:
                        for anum in dictpa[elt]:
                            if anum not in sum_atoms[elt]:
                                dictpa_d[elt].append(str(anum))
                        label = number_label(sum_atoms[elt])
                        dictpa_d[elt].append(label)
                else:
                    for anum in dictpa[elt]:
                        dictpa_d[elt].append(str(anum))

        elif (sum_atoms is None) and (sum_morbs is not None):
            dictio_d = {}
            for elt in dictio:
                dictio_d[elt] = []
                if elt in sum_morbs:
                    for morb in dictio[elt]:
                        if morb not in sum_morbs[elt]:
                            dictio_d[elt].append(morb)
                    label = orbital_label(sum_morbs[elt])
                    dictio_d[elt].append(label)
                else:
                    dictio_d[elt] = dictio[elt]
            dictpa_d = {elt: [str(anum) for anum in dictpa[elt]] for elt in dictpa}

        else:
            dictio_d = {}
            for elt in dictio:
                dictio_d[elt] = []
                if elt in sum_morbs:
                    for morb in dictio[elt]:
                        if morb not in sum_morbs[elt]:
                            dictio_d[elt].append(morb)
                    label = orbital_label(sum_morbs[elt])
                    dictio_d[elt].append(label)
                else:
                    dictio_d[elt] = dictio[elt]
            dictpa_d = {}
            for elt in dictpa:
                dictpa_d[elt] = []
                if elt in sum_atoms:
                    _sites = self._bs.structure.sites
                    indices = []
                    for site_idx in range(len(_sites)):
                        if next(iter(_sites[site_idx]._species)) == Element(elt):
                            indices.append(site_idx + 1)
                    flag_1 = len(set(dictpa[elt]).intersection(indices))
                    flag_2 = len(set(sum_atoms[elt]).intersection(indices))
                    if flag_1 == len(indices) and flag_2 == len(indices):
                        dictpa_d[elt].append("all")
                    else:
                        for anum in dictpa[elt]:
                            if anum not in sum_atoms[elt]:
                                dictpa_d[elt].append(str(anum))
                        label = number_label(sum_atoms[elt])
                        dictpa_d[elt].append(label)
                else:
                    for anum in dictpa[elt]:
                        dictpa_d[elt].append(str(anum))

        return dictio_d, dictpa_d

    def _make_ticks_selected(self, ax: plt.Axes, branches: list[int]) -> tuple[plt.Axes, list[float]]:
        """Utility private method to add ticks to a band structure with selected branches."""
        if not ax.figure:
            fig = plt.figure()  # Create a figure object
            ax.set_figure(fig)

        ticks = self.get_ticks()
        distance = []
        label = []
        rm_elems = []
        for idx in range(1, len(ticks["distance"])):
            if ticks["label"][idx] == ticks["label"][idx - 1]:
                rm_elems.append(idx)
        for idx in range(len(ticks["distance"])):
            if idx not in rm_elems:
                distance.append(ticks["distance"][idx])
                label.append(ticks["label"][idx])
        l_branches = [distance[i] - distance[i - 1] for i in range(1, len(distance))]
        n_distance = []
        n_label = []
        for branch in branches:
            n_distance.append(l_branches[branch])
            if ("$\\mid$" not in label[branch]) and ("$\\mid$" not in label[branch + 1]):
                n_label.append([label[branch], label[branch + 1]])
            elif ("$\\mid$" in label[branch]) and ("$\\mid$" not in label[branch + 1]):
                n_label.append([label[branch].split("$")[-1], label[branch + 1]])
            elif ("$\\mid$" not in label[branch]) and ("$\\mid$" in label[branch + 1]):
                n_label.append([label[branch], label[branch + 1].split("$")[0]])
            else:
                n_label.append([label[branch].split("$")[-1], label[branch + 1].split("$")[0]])

        f_distance: list[float] = []
        rf_distance: list[float] = []
        f_label: list[str] = []
        f_label.extend((n_label[0][0], n_label[0][1]))
        f_distance.extend((0.0, n_distance[0]))
        rf_distance.extend((0.0, n_distance[0]))
        length = n_distance[0]
        for idx in range(1, len(n_distance)):
            if n_label[idx][0] == n_label[idx - 1][1]:
                f_distance.extend((length, length + n_distance[idx]))
                f_label.extend((n_label[idx][0], n_label[idx][1]))
            else:
                f_distance.append(length + n_distance[idx])
                f_label[-1] = f"{n_label[idx - 1][1]}$\\mid${n_label[idx][0]}"
                f_label.append(n_label[idx][1])
            rf_distance.append(length + n_distance[idx])
            length += n_distance[idx]

        uniq_d = []
        uniq_l = []
        temp_ticks = list(zip(f_distance, f_label, strict=True))
        for idx, tick in enumerate(temp_ticks):
            if idx == 0:
                uniq_d.append(tick[0])
                uniq_l.append(tick[1])
                logger.debug(f"Adding label {tick[0]} at {tick[1]}")
            elif tick[1] == temp_ticks[idx - 1][1]:
                logger.debug(f"Skipping label {tick[1]}")
            else:
                logger.debug(f"Adding label {tick[0]} at {tick[1]}")
                uniq_d.append(tick[0])
                uniq_l.append(tick[1])

        logger.debug(f"Unique labels are {list(zip(uniq_d, uniq_l, strict=True))}")
        ax.set_xticks(uniq_d)
        ax.set_xticklabels(uniq_l)

        for idx in range(len(f_label)):
            if f_label[idx] is not None:
                # don't print the same label twice
                if idx != 0:
                    if f_label[idx] == f_label[idx - 1]:
                        logger.debug(f"already print label... skipping label {f_label[idx]}")
                    else:
                        logger.debug(f"Adding a line at {f_distance[idx]} for label {f_label[idx]}")
                        ax.axvline(f_distance[idx], color="k")
                else:
                    logger.debug(f"Adding a line at {f_distance[idx]} for label {f_label[idx]}")
                    ax.axvline(f_distance[idx], color="k")

        shift = []
        br = -1
        for branch in branches:
            br += 1
            shift.append(distance[branch] - rf_distance[br])

        return ax, shift


class BSDOSPlotter:
    """A joint, aligned band structure and density of states plot. Contributions
    from Jan Pohls as well as the online example from Germain Salvato-Vallverdu:
    https://gvallver.perso.univ-pau.fr/?p=587.
    """

    def __init__(
        self,
        bs_projection: Literal["elements"] | None = "elements",
        dos_projection: str = "elements",
        vb_energy_range: float = 4,
        cb_energy_range: float = 4,
        fixed_cb_energy: bool = False,
        egrid_interval: float = 1,
        font: str = "Times New Roman",
        axis_fontsize: float = 20,
        tick_fontsize: float = 15,
        legend_fontsize: float = 14,
        bs_legend: str = "best",
        dos_legend: str = "best",
        rgb_legend: bool = True,
        fig_size: tuple[float, float] = (11, 8.5),
    ) -> None:
        """Instantiate plotter settings.

        Args:
            bs_projection ('elements' | None): Whether to project the bands onto elements.
            dos_projection (str): "elements", "orbitals", or None
            vb_energy_range (float): energy in eV to show of valence bands
            cb_energy_range (float): energy in eV to show of conduction bands
            fixed_cb_energy (bool): If true, the cb_energy_range will be interpreted
                as constant (i.e., no gap correction for cb energy)
            egrid_interval (float): interval for grid marks
            font (str): font family
            axis_fontsize (float): font size for axis
            tick_fontsize (float): font size for axis tick labels
            legend_fontsize (float): font size for legends
            bs_legend (str): matplotlib string location for legend or None
            dos_legend (str): matplotlib string location for legend or None
            rgb_legend (bool): (T/F) whether to draw RGB triangle/bar for element proj.
            fig_size(tuple): dimensions of figure size (width, height)
        """
        self.bs_projection = bs_projection
        self.dos_projection = dos_projection
        self.vb_energy_range = vb_energy_range
        self.cb_energy_range = cb_energy_range
        self.fixed_cb_energy = fixed_cb_energy
        self.egrid_interval = egrid_interval
        self.font = font
        self.axis_fontsize = axis_fontsize
        self.tick_fontsize = tick_fontsize
        self.legend_fontsize = legend_fontsize
        self.bs_legend = bs_legend
        self.dos_legend = dos_legend
        self.rgb_legend = rgb_legend
        self.fig_size = fig_size

    def get_plot(
        self, bs: BandStructureSymmLine, dos: Dos | CompleteDos | None = None
    ) -> plt.Axes | tuple[plt.Axes, plt.Axes]:
        """Get a matplotlib plot object.

        Args:
            bs (BandStructureSymmLine): the bandstructure to plot. Projection
                data must exist for projected plots.
            dos (Dos): the Dos to plot. Projection data must exist (i.e.,
                CompleteDos) for projected plots.

        Returns:
            plt.Axes | tuple[plt.Axes, plt.Axes]: matplotlib axes for the band structure and DOS, resp.
        """
        # make sure the user-specified band structure projection is valid
        bs_projection = self.bs_projection
        if dos:
            elements = [e.symbol for e in dos.structure.elements]
        elif bs_projection and bs.structure:
            elements = [e.symbol for e in bs.structure.elements]
        else:
            elements = []

        rgb_legend = (
            self.rgb_legend and bs_projection and bs_projection.lower() == "elements" and len(elements) in [2, 3, 4]
        )

        if (
            bs_projection
            and bs_projection.lower() == "elements"
            and (len(elements) not in [2, 3, 4] or not bs.get_projection_on_elements())
        ):
            warnings.warn(
                "Cannot get element projected data; either the projection data "
                "doesn't exist, or you don't have a compound with exactly 2 "
                "or 3 or 4 unique elements."
            )
            bs_projection = None

        # specify energy range of plot
        emin = -self.vb_energy_range
        emax = self.cb_energy_range if self.fixed_cb_energy else self.cb_energy_range + bs.get_band_gap()["energy"]

        # initialize all the k-point labels and k-point x-distances for bs plot
        xlabels = []  # all symmetry point labels on x-axis
        xlabel_distances = []  # positions of symmetry point x-labels

        x_distances_list = []
        prev_right_klabel = None  # used to determine which branches require a midline separator

        for branch in bs.branches:
            x_distances = []

            # get left and right kpoint labels of this branch
            left_k, right_k = branch["name"].split("-")

            # add $ notation for LaTeX kpoint labels
            if left_k[0] == "\\" or "_" in left_k:
                left_k = f"${left_k}$"
            if right_k[0] == "\\" or "_" in right_k:
                right_k = f"${right_k}$"

            # add left k label to list of labels
            if prev_right_klabel is None:
                xlabels.append(left_k)
                xlabel_distances.append(0)
            elif prev_right_klabel != left_k:  # used for pipe separator
                xlabels[-1] = f"{xlabels[-1]}$\\mid$ {left_k}"

            # add right k label to list of labels
            xlabels.append(right_k)
            prev_right_klabel = right_k

            # add x-coordinates for labels
            left_kpoint = bs.kpoints[branch["start_index"]].cart_coords
            right_kpoint = bs.kpoints[branch["end_index"]].cart_coords
            distance = np.linalg.norm(right_kpoint - left_kpoint)
            xlabel_distances.append(xlabel_distances[-1] + distance)

            # add x-coordinates for kpoint data
            npts = branch["end_index"] - branch["start_index"]
            distance_interval = distance / npts
            x_distances.append(xlabel_distances[-2])
            for _ in range(npts):
                x_distances.append(x_distances[-1] + distance_interval)
            x_distances_list.append(x_distances)

        # set up bs and dos plot
        gs = GridSpec(1, 2, width_ratios=[2, 1]) if dos else GridSpec(1, 1)

        fig = plt.figure(figsize=self.fig_size)
        fig.patch.set_facecolor("white")
        bs_ax = plt.subplot(gs[0])
        dos_ax = plt.subplot(gs[1]) if dos else None

        # set basic axes limits for the plot
        bs_ax.set_xlim(0, x_distances_list[-1][-1])
        bs_ax.set_ylim(emin, emax)
        if dos and dos_ax is not None:
            dos_ax.set_ylim(emin, emax)

        # add BS xticks, labels, etc.
        bs_ax.set_xticks(xlabel_distances)
        bs_ax.set_xticklabels(xlabels, size=self.tick_fontsize)
        bs_ax.set_xlabel("Wavevector $k$", fontsize=self.axis_fontsize, family=self.font)
        bs_ax.set_ylabel("$E-E_F$ / eV", fontsize=self.axis_fontsize, family=self.font)

        # add BS Fermi level line at E=0 and gridlines
        bs_ax.hlines(y=0, xmin=0, xmax=x_distances_list[-1][-1], color="k", lw=2)
        bs_ax.set_yticks(np.arange(emin, emax + 1e-5, self.egrid_interval))
        bs_ax.set_yticklabels(np.arange(emin, emax + 1e-5, self.egrid_interval), size=self.tick_fontsize)
        bs_ax.set_axisbelow(b=True)
        bs_ax.grid(color=[0.5, 0.5, 0.5], linestyle="dotted", linewidth=1)
        if dos and dos_ax is not None:
            dos_ax.set_yticks(np.arange(emin, emax + 1e-5, self.egrid_interval))
            dos_ax.set_yticklabels([])
            dos_ax.grid(color=[0.5, 0.5, 0.5], linestyle="dotted", linewidth=1)

        # renormalize the band energy to the Fermi level
        band_energies: dict[Spin, list[float]] = {}
        for spin in (Spin.up, Spin.down):
            if spin in bs.bands:
                band_energies[spin] = []
                for band in bs.bands[spin]:
                    band = cast(list[float], band)
                    band_energies[spin].append([e - bs.efermi for e in band])  # type: ignore[arg-type]

        # renormalize the DOS energies to Fermi level
        dos_energies = [e - dos.efermi for e in dos.energies] if dos else []

        # get the projection data to set colors for the band structure
        colordata = self._get_colordata(bs, elements, bs_projection)

        # plot the colored band structure lines
        for spin in (Spin.up, Spin.down):
            if spin in band_energies:
                linestyles = "solid" if spin == Spin.up else "dotted"
                for band_idx, band in enumerate(band_energies[spin]):
                    current_pos = 0
                    for x_distances in x_distances_list:
                        sub_band = band[current_pos : current_pos + len(x_distances)]

                        self._rgbline(
                            bs_ax,
                            x_distances,
                            sub_band,
                            colordata[spin][band_idx, :, 0][current_pos : current_pos + len(x_distances)],
                            colordata[spin][band_idx, :, 1][current_pos : current_pos + len(x_distances)],
                            colordata[spin][band_idx, :, 2][current_pos : current_pos + len(x_distances)],
                            linestyles=linestyles,
                        )

                        current_pos += len(x_distances)

        if dos and dos_ax is not None:
            # Plot the DOS and projected DOS
            for spin in (Spin.up, Spin.down):
                if spin in dos.densities:
                    # plot the total DOS
                    dos_densities = dos.densities[spin] * int(spin)
                    label: str | None = "total" if spin == Spin.up else None
                    dos_ax.plot(dos_densities, dos_energies, color=(0.6, 0.6, 0.6), label=label)
                    dos_ax.fill_betweenx(
                        dos_energies,
                        0,
                        dos_densities,
                        color=(0.7, 0.7, 0.7),
                        facecolor=(0.7, 0.7, 0.7),
                    )

                    if self.dos_projection is None:
                        pass

                    elif self.dos_projection.lower() == "elements":
                        # plot the atom-projected DOS
                        colors = ["b", "r", "g", "m", "y", "c", "k", "w"]
                        el_dos = dos.get_element_dos()
                        for idx, el in enumerate(elements):
                            dos_densities = el_dos[Element(el)].densities[spin] * int(spin)
                            label = el if spin == Spin.up else None
                            dos_ax.plot(
                                dos_densities,
                                dos_energies,
                                color=colors[idx],
                                label=label,
                            )

                    elif self.dos_projection.lower() == "orbitals":
                        # plot each of the atomic projected DOS
                        colors = ["b", "r", "g", "m"]
                        spd_dos = dos.get_spd_dos()
                        for idx, orb in enumerate([OrbitalType.s, OrbitalType.p, OrbitalType.d, OrbitalType.f]):
                            if orb in spd_dos:
                                dos_densities = spd_dos[orb].densities[spin] * int(spin)
                                label = orb if spin == Spin.up else None  # type: ignore[assignment]
                                dos_ax.plot(
                                    dos_densities,
                                    dos_energies,
                                    color=colors[idx],
                                    label=label,
                                )

            # get index of lowest and highest energy being plotted, used to help auto-scale DOS x-axis
            emin_idx = next(x[0] for x in enumerate(dos_energies) if x[1] >= emin)
            emax_idx = len(dos_energies) - next(x[0] for x in enumerate(reversed(dos_energies)) if x[1] <= emax)

            # determine DOS x-axis range
            dos_xmin = (
                0 if Spin.down not in dos.densities else -max(dos.densities[Spin.down][emin_idx : emax_idx + 1] * 1.05)
            )
            dos_xmax = max([max(dos.densities[Spin.up][emin_idx:emax_idx]) * 1.05, abs(dos_xmin)])

            # set up the DOS x-axis and add Fermi level line
            dos_ax.set_xlim(dos_xmin, dos_xmax)
            dos_ax.set_xticklabels([])
            dos_ax.hlines(y=0, xmin=dos_xmin, xmax=dos_xmax, color="k", lw=2)
            dos_ax.set_xlabel("DOS", fontsize=self.axis_fontsize, family=self.font)

        # add legend for band structure
        if self.bs_legend and not rgb_legend:
            handles = []

            if bs_projection is None:
                handles = [
                    mlines.Line2D([], [], linewidth=2, color="k", label="spin up"),
                    mlines.Line2D(
                        [],
                        [],
                        linewidth=2,
                        color="b",
                        linestyle="dotted",
                        label="spin down",
                    ),
                ]

            elif bs_projection.lower() == "elements":
                colors = ["b", "r", "g"]
                for idx, el in enumerate(elements):
                    handles.append(mlines.Line2D([], [], linewidth=2, color=colors[idx], label=el))

            bs_ax.legend(
                handles=handles,
                fancybox=True,
                prop={"size": self.legend_fontsize, "family": self.font},
                loc=self.bs_legend,
            )

        elif self.bs_legend and rgb_legend:
            if len(elements) == 2:
                self._rb_line(bs_ax, elements[1], elements[0], loc=self.bs_legend)
            elif len(elements) == 3:
                self._rgb_triangle(bs_ax, elements[1], elements[2], elements[0], loc=self.bs_legend)
            elif len(elements) == 4:
                self._cmyk_triangle(
                    bs_ax,
                    elements[1],
                    elements[2],
                    elements[0],
                    elements[3],
                    loc=self.bs_legend,
                )
        # add legend for DOS
        if dos and self.dos_legend and dos_ax is not None:
            dos_ax.legend(
                fancybox=True,
                prop={"size": self.legend_fontsize, "family": self.font},
                loc=self.dos_legend,
            )

        plt.subplots_adjust(wspace=0.1)
        if dos and dos_ax is not None:
            return bs_ax, dos_ax
        return bs_ax

    @staticmethod
    def _rgbline(ax, k, e, red, green, blue, alpha=1, linestyles="solid") -> None:
        """An RGB colored line for plotting.
        creation of segments based on:
        http://nbviewer.ipython.org/urls/raw.github.com/dpsanders/matplotlib-examples/master/colorline.ipynb.

        Args:
            ax: matplotlib axis
            k: x-axis data (k-points)
            e: y-axis data (energies)
            red: red data
            green: green data
            blue: blue data
            alpha: alpha values data
            linestyles: linestyle for plot (e.g., "solid" or "dotted").
        """
        pts = np.array([k, e]).T.reshape(-1, 1, 2)
        seg = np.concatenate([pts[:-1], pts[1:]], axis=1)

        n_seg = len(k) - 1
        red = [0.5 * (red[i] + red[i + 1]) for i in range(n_seg)]
        green = [0.5 * (green[i] + green[i + 1]) for i in range(n_seg)]
        blue = [0.5 * (blue[i] + blue[i + 1]) for i in range(n_seg)]
        alpha = np.ones(n_seg, float) * alpha
        lc = LineCollection(
            seg,
            colors=list(zip(red, green, blue, alpha, strict=True)),
            linewidth=2,
            linestyles=linestyles,
        )
        ax.add_collection(lc)

    @staticmethod
    def _get_colordata(bs, elements, bs_projection):
        """Get color data, including projected band structures.

        Args:
            bs: Bandstructure object
            elements: elements (in desired order) for setting to blue, red, green
            bs_projection: None for no projection, "elements" for element projection

        Returns:
            Dictionary representation of color data.
        """
        contribs = {}
        projections = None
        if bs_projection and bs_projection.lower() == "elements":
            projections = bs.get_projection_on_elements()

        for spin in (Spin.up, Spin.down):
            if spin in bs.bands:
                contribs[spin] = []
                for band_idx in range(bs.nb_bands):
                    colors = []
                    for k_idx in range(len(bs.kpoints)):
                        if bs_projection and bs_projection.lower() == "elements":
                            c = [0, 0, 0, 0]
                            projs = projections[spin][band_idx][k_idx]
                            # note: squared color interpolations are smoother
                            # see: https://youtu.be/LKnqECcg6Gw
                            projs = {k: v**2 for k, v in projs.items()}
                            total = sum(projs.values())
                            if total > 0:
                                for idx, e in enumerate(elements):
                                    c[idx] = math.sqrt(projs[e] / total)  # min is to handle round errors

                            c = [
                                c[1],
                                c[2],
                                c[0],
                                c[3],
                            ]  # prefer blue, then red, then green or magenta, then yellow, then cyan, then black
                            if len(elements) == 4:
                                # convert cmyk to rgb
                                c = [
                                    (1 - c[0]) * (1 - c[3]),
                                    ((1 - c[1]) * (1 - c[3])),
                                    ((1 - c[2]) * (1 - c[3])),
                                ]
                            else:
                                c = [c[0], c[1], c[2]]

                        else:
                            c = [0, 0, 0] if spin == Spin.up else [0, 0, 1]  # black for spin up, blue for spin down

                        colors.append(c)

                    contribs[spin].append(colors)
                contribs[spin] = np.array(contribs[spin])

        return contribs

    @staticmethod
    def _cmyk_triangle(ax, c_label, m_label, y_label, k_label, loc) -> None:
        """Draw an RGB triangle legend on the desired axis."""
        if loc not in range(1, 11):
            loc = 2

        from mpl_toolkits.axes_grid1.inset_locator import inset_axes

        inset_ax = inset_axes(ax, width=1.5, height=1.5, loc=loc)
        mesh = 35
        x = []
        y = []
        color = []
        for c in range(mesh):
            for ye in range(mesh):
                for m in range(mesh):
                    if not (c == mesh - 1 and ye == mesh - 1 and m == mesh - 1) and not (c == 0 and ye == 0 and m == 0):
                        c1 = c / (c + ye + m)
                        ye1 = ye / (c + ye + m)
                        m1 = m / (c + ye + m)
                        x.append(0.33 * (2.0 * ye1 + c1) / (c1 + ye1 + m1))
                        y.append(0.33 * np.sqrt(3) * c1 / (c1 + ye1 + m1))
                        rc = 1 - c / (mesh - 1)
                        gc = 1 - m / (mesh - 1)
                        bc = 1 - ye / (mesh - 1)
                        color.append([rc, gc, bc])

        # x = [n + 0.25 for n in x]  # nudge x coordinates
        # y = [n + (max_y - 1) for n in y]  # shift y coordinates to top
        # plot the triangle
        inset_ax.scatter(x, y, s=7, marker=".", edgecolor=color)
        inset_ax.set_xlim([-0.35, 1.00])
        inset_ax.set_ylim([-0.35, 1.00])

        # add the labels
        common = {"fontsize": 13, "family": "Times New Roman"}
        inset_ax.text(0.70, -0.2, m_label, **common, color=(0, 0, 0), horizontalalignment="left")
        inset_ax.text(
            0.325,
            0.70,
            c_label,
            **common,
            color=(0, 0, 0),
            horizontalalignment="center",
        )
        inset_ax.text(-0.05, -0.2, y_label, **common, color=(0, 0, 0), horizontalalignment="right")
        inset_ax.text(
            0.325,
            0.22,
            k_label,
            **common,
            color=(1, 1, 1),
            horizontalalignment="center",
        )

        inset_ax.axis("off")

    @staticmethod
    def _rgb_triangle(ax, r_label, g_label, b_label, loc) -> None:
        """Draw an RGB triangle legend on the desired axis."""
        if loc not in range(1, 11):
            loc = 2

        from mpl_toolkits.axes_grid1.inset_locator import inset_axes

        inset_ax = inset_axes(ax, width=1, height=1, loc=loc)
        mesh = 35
        x = []
        y = []
        color = []
        for r in range(mesh):
            for g in range(mesh):
                for b in range(mesh):
                    if not (r == 0 and b == 0 and g == 0):
                        r1 = r / (r + g + b)
                        g1 = g / (r + g + b)
                        b1 = b / (r + g + b)
                        x.append(0.33 * (2.0 * g1 + r1) / (r1 + b1 + g1))
                        y.append(0.33 * np.sqrt(3) * r1 / (r1 + b1 + g1))
                        rc = math.sqrt(r**2 / (r**2 + g**2 + b**2))
                        gc = math.sqrt(g**2 / (r**2 + g**2 + b**2))
                        bc = math.sqrt(b**2 / (r**2 + g**2 + b**2))
                        color.append([rc, gc, bc])

        # x = [n + 0.25 for n in x]  # nudge x coordinates
        # y = [n + (max_y - 1) for n in y]  # shift y coordinates to top
        # plot the triangle
        inset_ax.scatter(x, y, s=7, marker=".", edgecolor=color)
        inset_ax.set_xlim([-0.35, 1.00])
        inset_ax.set_ylim([-0.35, 1.00])

        # add the labels
        inset_ax.text(
            0.70,
            -0.2,
            g_label,
            fontsize=13,
            family="Times New Roman",
            color=(0, 0, 0),
            horizontalalignment="left",
        )
        inset_ax.text(
            0.325,
            0.70,
            r_label,
            fontsize=13,
            family="Times New Roman",
            color=(0, 0, 0),
            horizontalalignment="center",
        )
        inset_ax.text(
            -0.05,
            -0.2,
            b_label,
            fontsize=13,
            family="Times New Roman",
            color=(0, 0, 0),
            horizontalalignment="right",
        )

        inset_ax.axis("off")

    @staticmethod
    def _rb_line(ax, r_label, b_label, loc) -> None:
        # Draw an rb bar legend on the desired axis

        if loc not in range(1, 11):
            loc = 2
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes

        inset_ax = inset_axes(ax, width=1.2, height=0.4, loc=loc)

        x, y, color = [], [], []
        for idx in range(1000):
            x.append(idx / 1800.0 + 0.55)
            y.append(0)
            color.append([math.sqrt(c) for c in [1 - (idx / 1000) ** 2, 0, (idx / 1000) ** 2]])

        # plot the bar
        inset_ax.scatter(x, y, s=250.0, marker="s", c=color)
        inset_ax.set_xlim([-0.1, 1.7])
        inset_ax.text(
            1.35,
            0,
            b_label,
            fontsize=13,
            family="Times New Roman",
            color=(0, 0, 0),
            horizontalalignment="left",
            verticalalignment="center",
        )
        inset_ax.text(
            0.30,
            0,
            r_label,
            fontsize=13,
            family="Times New Roman",
            color=(0, 0, 0),
            horizontalalignment="right",
            verticalalignment="center",
        )

        inset_ax.axis("off")


class BoltztrapPlotter:
    # TODO: We need a unittest for this. Come on folks.
    """Plot Boltztrap data."""

    def __init__(self, bz) -> None:
        """
        Args:
            bz: a BoltztrapAnalyzer object.
        """
        self._bz = bz

    def _plot_doping(self, plt, temp) -> None:
        if len(self._bz.doping) != 0:
            limit = 2.21e15
            plt.axvline(self._bz.mu_doping["n"][temp][0], linewidth=3.0, linestyle="--")
            plt.text(
                self._bz.mu_doping["n"][temp][0] + 0.01,
                limit,
                f"$n$=10^{{{math.log10(self._bz.doping['n'][0])}}}$",
                color="b",
            )
            plt.axvline(self._bz.mu_doping["n"][temp][-1], linewidth=3.0, linestyle="--")
            plt.text(
                self._bz.mu_doping["n"][temp][-1] + 0.01,
                limit,
                f"$n$=10^{{{math.log10(self._bz.doping['n'][-1])}}}$",
                color="b",
            )
            plt.axvline(self._bz.mu_doping["p"][temp][0], linewidth=3.0, linestyle="--")
            plt.text(
                self._bz.mu_doping["p"][temp][0] + 0.01,
                limit,
                f"$p$=10^{{{math.log10(self._bz.doping['p'][0])}}}$",
                color="b",
            )
            plt.axvline(self._bz.mu_doping["p"][temp][-1], linewidth=3.0, linestyle="--")
            plt.text(
                self._bz.mu_doping["p"][temp][-1] + 0.01,
                limit,
                f"$p$=10^{{{math.log10(self._bz.doping['p'][-1])}}}$",
                color="b",
            )

    def _plot_bg_limits(self, plt) -> None:
        plt.axvline(0.0, color="k", linewidth=3.0)
        plt.axvline(self._bz.gap, color="k", linewidth=3.0)

    def plot_seebeck_eff_mass_mu(self, temps=(300,), output="average", Lambda=0.5):
        """Plot respect to the chemical potential of the Seebeck effective mass
        calculated as explained in Ref.
        Gibbs, Z. M. et al., Effective mass and fermi surface complexity factor
        from ab initio band structure calculations.
        npj Computational Materials 3, 8 (2017).

        Args:
            output: 'average' returns the seebeck effective mass calculated
                using the average of the three diagonal components of the
                seebeck tensor. 'tensor' returns the seebeck effective mass
                respect to the three diagonal components of the seebeck tensor.
            temps: list of temperatures of calculated seebeck.
            Lambda: fitting parameter used to model the scattering (0.5 means
                constant relaxation time).

        Returns:
            a matplotlib object
        """
        ax = pretty_plot(9, 7)
        for temp in temps:
            sbk_mass = self._bz.get_seebeck_eff_mass(output=output, temp=temp, Lambda=0.5)
            # remove noise inside the gap
            start = self._bz.mu_doping["p"][temp][0]
            stop = self._bz.mu_doping["n"][temp][0]
            mu_steps_1 = []
            mu_steps_2 = []
            sbk_mass_1 = []
            sbk_mass_2 = []
            for i, mu in enumerate(self._bz.mu_steps):
                if mu <= start:
                    mu_steps_1.append(mu)
                    sbk_mass_1.append(sbk_mass[i])
                elif mu >= stop:
                    mu_steps_2.append(mu)
                    sbk_mass_2.append(sbk_mass[i])

            ax.plot(mu_steps_1, sbk_mass_1, label=f"{temp}K", linewidth=3)
            ax.plot(mu_steps_2, sbk_mass_2, linewidth=3.0)
            if output == "average":
                ax.get_lines()[1].set_c(ax.get_lines()[0].get_c())
            elif output == "tensor":
                ax.get_lines()[3].set_c(ax.get_lines()[0].get_c())
                ax.get_lines()[4].set_c(ax.get_lines()[1].get_c())
                ax.get_lines()[5].set_c(ax.get_lines()[2].get_c())

        ax.set_xlabel("E-E$_f$ (eV)", fontsize=30)
        ax.set_ylabel("Seebeck effective mass", fontsize=30)
        ax.set_xticks(fontsize=25)
        ax.set_yticks(fontsize=25)
        if output == "tensor":
            ax.legend(
                [f"{dim}_{T}K" for T in temps for dim in ("x", "y", "z")],
                fontsize=20,
            )
        elif output == "average":
            ax.legend(fontsize=20)
        plt.tight_layout()
        return ax

    def plot_complexity_factor_mu(self, temps=(300,), output="average", Lambda=0.5):
        """Plot respect to the chemical potential of the Fermi surface complexity
        factor calculated as explained in Ref.
        Gibbs, Z. M. et al., Effective mass and fermi surface complexity factor
        from ab initio band structure calculations.
        npj Computational Materials 3, 8 (2017).

        Args:
            output: 'average' returns the complexity factor calculated using the average
                of the three diagonal components of the seebeck and conductivity tensors.
                'tensor' returns the complexity factor respect to the three
                diagonal components of seebeck and conductivity tensors.
            temps: list of temperatures of calculated seebeck and conductivity.
            Lambda: fitting parameter used to model the scattering (0.5 means constant
                relaxation time).

        Returns:
            a matplotlib object
        """
        ax = pretty_plot(9, 7)
        for T in temps:
            cmplx_fact = self._bz.get_complexity_factor(output=output, temp=T, Lambda=Lambda)
            start = self._bz.mu_doping["p"][T][0]
            stop = self._bz.mu_doping["n"][T][0]
            mu_steps_1 = []
            mu_steps_2 = []
            cmplx_fact_1 = []
            cmplx_fact_2 = []
            for i, mu in enumerate(self._bz.mu_steps):
                if mu <= start:
                    mu_steps_1.append(mu)
                    cmplx_fact_1.append(cmplx_fact[i])
                elif mu >= stop:
                    mu_steps_2.append(mu)
                    cmplx_fact_2.append(cmplx_fact[i])

            ax.plot(mu_steps_1, cmplx_fact_1, label=f"{T}K", linewidth=3.0)
            ax.plot(mu_steps_2, cmplx_fact_2, linewidth=3.0)
            if output == "average":
                ax.gca().get_lines()[1].set_c(ax.gca().get_lines()[0].get_c())
            elif output == "tensor":
                ax.gca().get_lines()[3].set_c(ax.gca().get_lines()[0].get_c())
                ax.gca().get_lines()[4].set_c(ax.gca().get_lines()[1].get_c())
                ax.gca().get_lines()[5].set_c(ax.gca().get_lines()[2].get_c())

        ax.set_xlabel("E-E$_f$ (eV)", fontsize=30)
        ax.set_ylabel("Complexity Factor", fontsize=30)
        ax.set_xticks(fontsize=25)
        ax.set_yticks(fontsize=25)
        if output == "tensor":
            ax.legend(
                [f"{dim}_{T}K" for T in temps for dim in ("x", "y", "z")],
                fontsize=20,
            )
        elif output == "average":
            ax.legend(fontsize=20)
        plt.tight_layout()
        return ax

    def plot_seebeck_mu(
        self,
        temp: float = 600,
        output: str = "eig",
        xlim: Sequence[float] | None = None,
    ):
        """Plot the seebeck coefficient in function of Fermi level.

        Args:
            temp (float): the temperature
            output (str): "eig" or "average"
            xlim (tuple[float, float]): a 2-tuple of min and max fermi energy. Defaults to (0, band gap)

        Returns:
            a matplotlib object
        """
        ax = pretty_plot(9, 7)
        seebeck = self._bz.get_seebeck(output=output, doping_levels=False)[temp]
        ax.plot(self._bz.mu_steps, seebeck, linewidth=3.0)

        self._plot_bg_limits(ax)
        self._plot_doping(ax, temp)
        if output == "eig":
            ax.legend(["S$_1$", "S$_2$", "S$_3$"])
        if xlim is None:
            ax.set_xlim(-0.5, self._bz.gap + 0.5)
        else:
            ax.set_xlim(xlim[0], xlim[1])
        ax.set_ylabel("Seebeck \n coefficient  ($\\mu$V/K)", fontsize=30.0)
        ax.set_xlabel("E-E$_f$ (eV)", fontsize=30)
        ax.set_xticks(fontsize=25)
        ax.set_yticks(fontsize=25)
        plt.tight_layout()
        return ax

    def plot_conductivity_mu(
        self,
        temp: float = 600,
        output: str = "eig",
        relaxation_time: float = 1e-14,
        xlim: Sequence[float] | None = None,
    ):
        """Plot the conductivity in function of Fermi level. Semi-log plot.

        Args:
            temp (float): the temperature
            output (str): "eig" or "average"
            relaxation_time (float): A relaxation time in s. Defaults to 1e-14 and the plot is in
               units of relaxation time
            xlim (tuple[float, float]): a 2-tuple of min and max fermi energy. Defaults to (0, band gap)

        Returns:
            a matplotlib object
        """
        cond = self._bz.get_conductivity(relaxation_time=relaxation_time, output=output, doping_levels=False)[temp]
        ax = pretty_plot(9, 7)
        ax.semilogy(self._bz.mu_steps, cond, linewidth=3.0)
        self._plot_bg_limits(ax)
        self._plot_doping(ax, temp)
        if output == "eig":
            ax.legend(["$\\Sigma_1$", "$\\Sigma_2$", "$\\Sigma_3$"])
        if xlim is None:
            ax.set_xlim(-0.5, self._bz.gap + 0.5)
        else:
            ax.set_xlim(xlim)
        ax.set_ylim([1e13 * relaxation_time, 1e20 * relaxation_time])
        ax.set_ylabel("conductivity,\n $\\Sigma$ (1/($\\Omega$ m))", fontsize=30.0)
        ax.set_xlabel("E-E$_f$ (eV)", fontsize=30.0)
        ax.set_xticks(fontsize=25)
        ax.set_yticks(fontsize=25)
        plt.tight_layout()
        return ax

    def plot_power_factor_mu(
        self,
        temp: float = 600,
        output: str = "eig",
        relaxation_time: float = 1e-14,
        xlim: Sequence[float] | None = None,
    ):
        """Plot the power factor in function of Fermi level. Semi-log plot.

        Args:
            temp (float): the temperature
            output (str): "eig" or "average"
            relaxation_time (float): A relaxation time in s. Defaults to 1e-14 and the plot is in
               units of relaxation time
            xlim (tuple[float, float]): a 2-tuple of min and max fermi energy. Defaults to (0, band gap)

        Returns:
            a matplotlib object
        """
        ax = pretty_plot(9, 7)
        pow_factor = self._bz.get_power_factor(relaxation_time=relaxation_time, output=output, doping_levels=False)[
            temp
        ]
        ax.semilogy(self._bz.mu_steps, pow_factor, linewidth=3.0)
        self._plot_bg_limits(ax)
        self._plot_doping(ax, temp)
        if output == "eig":
            ax.legend(["PF$_1$", "PF$_2$", "PF$_3$"])
        if xlim is None:
            ax.set_xlim(-0.5, self._bz.gap + 0.5)
        else:
            ax.set_xlim(xlim)
        ax.set_ylabel("Power factor, ($\\mu$W/(mK$^2$))", fontsize=30.0)
        ax.set_xlabel("E-E$_f$ (eV)", fontsize=30.0)
        ax.set_xticks(fontsize=25)
        ax.set_yticks(fontsize=25)
        plt.tight_layout()
        return ax

    def plot_zt_mu(
        self,
        temp: float = 600,
        output: str = "eig",
        relaxation_time: float = 1e-14,
        xlim: Sequence[float] | None = None,
    ) -> plt.Axes:
        """Plot the ZT as function of Fermi level.

        Args:
            temp (float): the temperature
            output (str): "eig" or "average"
            relaxation_time (float): A relaxation time in s. Defaults to 1e-14 and the plot is in
               units of relaxation time
            xlim (tuple[float, float]): a 2-tuple of min and max fermi energy. Defaults to (0, band gap)

        Returns:
            plt.Axes: matplotlib axes object
        """
        ax = pretty_plot(9, 7)
        zt = self._bz.get_zt(relaxation_time=relaxation_time, output=output, doping_levels=False)[temp]
        ax.plot(self._bz.mu_steps, zt, linewidth=3.0)
        self._plot_bg_limits(ax)
        self._plot_doping(ax, temp)
        if output == "eig":
            ax.legend(["ZT$_1$", "ZT$_2$", "ZT$_3$"])
        if xlim is None:
            ax.set_xlim(-0.5, self._bz.gap + 0.5)
        else:
            ax.set_xlim(xlim)
        ax.set_ylabel("ZT", fontsize=30.0)
        ax.set_xlabel("E-E$_f$ (eV)", fontsize=30.0)
        ax.set_xticks(fontsize=25)
        ax.set_yticks(fontsize=25)
        plt.tight_layout()
        return ax

    def plot_seebeck_temp(self, doping="all", output: Literal["average", "eigs"] = "average"):
        """Plot the Seebeck coefficient in function of temperature for different
        doping levels.

        Args:
            doping (str): the default 'all' plots all the doping levels in the analyzer.
                Specify a list of doping levels if you want to plot only some.
            output: with 'average' you get an average of the three directions
                with 'eigs' you get all the three directions.

        Returns:
            a matplotlib object
        """
        sbk = self._bz.get_seebeck(output="average") if output == "average" else self._bz.get_seebeck(output="eigs")

        ax = pretty_plot(22, 14)
        temperatures = sorted(sbk["n"])
        doping = self._bz.doping["n"] if doping == "all" else doping
        for idx, dop_type in enumerate("np"):
            plt.subplot(121 + idx)
            for dop in doping:
                dop_idx = self._bz.doping[dop_type].index(dop)
                sbk_temp = []
                for temp in temperatures:
                    sbk_temp.append(sbk[dop_type][temp][dop_idx])
                if output == "average":
                    ax.plot(temperatures, sbk_temp, marker="s", label=f"{dop} $cm^{-3}$")
                elif output == "eigs":
                    for xyz in range(3):
                        ax.plot(
                            temperatures,
                            list(zip(*sbk_temp, strict=True))[xyz],
                            marker="s",
                            label=f"{xyz} {dop} $cm^{{-3}}$",
                        )
            ax.set_title(f"{dop_type}-type", fontsize=20)
            if idx == 0:
                ax.set_ylabel("Seebeck \n coefficient  ($\\mu$V/K)", fontsize=30.0)
            ax.set_xlabel("Temperature (K)", fontsize=30.0)

            ax.legend(loc="best", fontsize=15)
            ax.grid()
            ax.set_xticks(fontsize=25)
            ax.set_yticks(fontsize=25)

        plt.tight_layout()

        return ax

    def plot_conductivity_temp(
        self,
        doping="all",
        output: Literal["average", "eigs"] = "average",
        relaxation_time=1e-14,
    ):
        """Plot the conductivity in function of temperature for different doping levels.

        Args:
            doping (str): the default 'all' plots all the doping levels in the analyzer.
                Specify a list of doping levels if you want to plot only some.
            output: with 'average' you get an average of the three directions
                with 'eigs' you get all the three directions.
            relaxation_time: specify a constant relaxation time value

        Returns:
            a matplotlib object
        """
        if output == "average":
            cond = self._bz.get_conductivity(relaxation_time=relaxation_time, output="average")
        else:
            cond = self._bz.get_conductivity(relaxation_time=relaxation_time, output="eigs")

        ax = pretty_plot(22, 14)
        temperatures = sorted(cond["n"])
        doping = self._bz.doping["n"] if doping == "all" else doping
        for idx, dop_type in enumerate("np"):
            plt.subplot(121 + idx)
            for dop in doping:
                dop_idx = self._bz.doping[dop_type].index(dop)
                cond_temp = []
                for temp in temperatures:
                    cond_temp.append(cond[dop_type][temp][dop_idx])
                if output == "average":
                    ax.plot(temperatures, cond_temp, marker="s", label=f"{dop} $cm^{-3}$")
                elif output == "eigs":
                    for xyz in range(3):
                        ax.plot(
                            temperatures,
                            list(zip(*cond_temp, strict=True))[xyz],
                            marker="s",
                            label=f"{xyz} {dop} $cm^{{-3}}$",
                        )
            ax.set_title(f"{dop_type}-type", fontsize=20)
            if idx == 0:
                ax.set_ylabel("conductivity $\\sigma$ (1/($\\Omega$ m))", fontsize=30.0)
            ax.set_xlabel("Temperature (K)", fontsize=30.0)

            ax.legend(loc="best", fontsize=15)
            ax.grid()
            ax.set_xticks(fontsize=25)
            ax.set_yticks(fontsize=25)
            ax.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))

        plt.tight_layout()

        return ax

    def plot_power_factor_temp(
        self,
        doping="all",
        output: Literal["average", "eigs"] = "average",
        relaxation_time=1e-14,
    ):
        """Plot the Power Factor in function of temperature for different doping levels.

        Args:
            doping (str): the default 'all' plots all the doping levels in the analyzer.
                Specify a list of doping levels if you want to plot only some.
            output: with 'average' you get an average of the three directions
                with 'eigs' you get all the three directions.
            relaxation_time: specify a constant relaxation time value

        Returns:
            a matplotlib object
        """
        if output == "average":
            pow_factor = self._bz.get_power_factor(relaxation_time=relaxation_time, output="average")
        else:
            pow_factor = self._bz.get_power_factor(relaxation_time=relaxation_time, output="eigs")

        ax = pretty_plot(22, 14)
        temperatures = sorted(pow_factor["n"])
        doping = self._bz.doping["n"] if doping == "all" else doping
        for idx, dop_type in enumerate("np"):
            plt.subplot(121 + idx)
            for dop in doping:
                dop_idx = self._bz.doping[dop_type].index(dop)
                pf_temp = []
                for temp in temperatures:
                    pf_temp.append(pow_factor[dop_type][temp][dop_idx])
                if output == "average":
                    ax.plot(temperatures, pf_temp, marker="s", label=f"{dop} $cm^{-3}$")
                elif output == "eigs":
                    for xyz in range(3):
                        ax.plot(
                            temperatures,
                            list(zip(*pf_temp, strict=True))[xyz],
                            marker="s",
                            label=f"{xyz} {dop} $cm^{{-3}}$",
                        )
            ax.set_title(f"{dop_type}-type", fontsize=20)
            if idx == 0:
                ax.set_ylabel("Power Factor ($\\mu$W/(mK$^2$))", fontsize=30.0)
            ax.set_xlabel("Temperature (K)", fontsize=30.0)

            ax.legend(loc="best", fontsize=15)
            ax.grid()
            ax.set_xticks(fontsize=25)
            ax.set_yticks(fontsize=25)
            ax.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))

        plt.tight_layout()
        return ax

    def plot_zt_temp(
        self,
        doping="all",
        output: Literal["average", "eigs"] = "average",
        relaxation_time=1e-14,
    ):
        """Plot the figure of merit zT in function of temperature for different doping levels.

        Args:
            doping (str): the default 'all' plots all the doping levels in the analyzer.
                Specify a list of doping levels if you want to plot only some.
            output: with 'average' you get an average of the three directions
                with 'eigs' you get all the three directions.
            relaxation_time: specify a constant relaxation time value

        Raises:
            ValueError: if output is not 'average' or 'eigs'

        Returns:
            a matplotlib object
        """
        if output not in {"average", "eigs"}:
            raise ValueError(f"{output=} must be 'average' or 'eigs'")
        zt = self._bz.get_zt(relaxation_time=relaxation_time, output=output)

        ax = pretty_plot(22, 14)
        temperatures = sorted(zt["n"])
        doping = self._bz.doping["n"] if doping == "all" else doping
        for idx, dop_type in enumerate("np"):
            plt.subplot(121 + idx)
            for dop in doping:
                dop_idx = self._bz.doping[dop_type].index(dop)
                zt_temp = []
                for temp in temperatures:
                    zt_temp.append(zt[dop_type][temp][dop_idx])
                if output == "average":
                    ax.plot(temperatures, zt_temp, marker="s", label=f"{dop} $cm^{-3}$")
                elif output == "eigs":
                    for xyz in range(3):
                        ax.plot(
                            temperatures,
                            list(zip(*zt_temp, strict=True))[xyz],
                            marker="s",
                            label=f"{xyz} {dop} $cm^{{-3}}$",
                        )
            ax.set_title(f"{dop_type}-type", fontsize=20)
            if idx == 0:
                ax.set_ylabel("zT", fontsize=30.0)
            ax.set_xlabel("Temperature (K)", fontsize=30.0)

            ax.legend(loc="best", fontsize=15)
            ax.grid()
            ax.set_xticks(fontsize=25)
            ax.set_yticks(fontsize=25)

        plt.tight_layout()
        return ax

    def plot_eff_mass_temp(self, doping="all", output: Literal["average", "eigs"] = "average"):
        """Plot the average effective mass in function of temperature
        for different doping levels.

        Args:
            doping (str): the default 'all' plots all the doping levels in the analyzer.
                Specify a list of doping levels if you want to plot only some.
            output ('average' | 'eigs'): with 'average' you get an average of the three directions
                with 'eigs' you get all the three directions.

        Returns:
            a matplotlib Axes object
        """
        if output == "average":
            eff_mass = self._bz.get_average_eff_mass(output="average")
        else:
            eff_mass = self._bz.get_average_eff_mass(output="eigs")

        ax_main = pretty_plot(22, 14)
        temperatures = sorted(eff_mass["n"])
        doping = self._bz.doping["n"] if doping == "all" else doping
        for idx, dop_type in enumerate("np"):
            ax = plt.subplot(121 + idx)
            for dop in doping:
                dop_idx = self._bz.doping[dop_type].index(dop)
                em_temp = [eff_mass[dop_type][temp][dop_idx] for temp in temperatures]
                if output == "average":
                    ax.plot(temperatures, em_temp, marker="s", label=f"{dop} $cm^{{-3}}$")
                elif output == "eigs":
                    for xyz in range(3):
                        ax.plot(
                            temperatures,
                            list(zip(*em_temp, strict=True))[xyz],
                            marker="s",
                            label=f"{xyz} {dop} $cm^{{-3}}$",
                        )
            ax.set_title(f"{dop_type}-type", fontsize=20)
            if idx == 0:
                ax.set_ylabel("Effective mass (m$_e$)", fontsize=30.0)
            ax.set_xlabel("Temperature (K)", fontsize=30.0)

            ax.legend(loc="best", fontsize=15)
            ax.grid()
            ax.tick_params(labelsize=25)

        plt.tight_layout()
        return ax_main

    def plot_seebeck_dop(self, temps="all", output: Literal["average", "eigs"] = "average"):
        """Plot the Seebeck in function of doping levels for different temperatures.

        Args:
            temps: the default 'all' plots all the temperatures in the analyzer.
                Specify a list of temperatures if you want to plot only some.
            output: with 'average' you get an average of the three directions
                with 'eigs' you get all the three directions.

        Returns:
            a matplotlib object
        """
        sbk = self._bz.get_seebeck(output="average") if output == "average" else self._bz.get_seebeck(output="eigs")

        temperatures = sorted(sbk["n"]) if temps == "all" else temps
        ax = pretty_plot(22, 14)
        for idx, dop_type in enumerate("np"):
            plt.subplot(121 + idx)
            for temp in temperatures:
                if output == "eigs":
                    for xyz in range(3):
                        ax.semilogx(
                            self._bz.doping[dop_type],
                            list(zip(*sbk[dop_type][temp], strict=True))[xyz],
                            marker="s",
                            label=f"{xyz} {temp} K",
                        )
                elif output == "average":
                    ax.semilogx(
                        self._bz.doping[dop_type],
                        sbk[dop_type][temp],
                        marker="s",
                        label=f"{temp} K",
                    )
            ax.set_title(f"{dop_type}-type", fontsize=20)
            if idx == 0:
                ax.set_ylabel("Seebeck coefficient ($\\mu$V/K)", fontsize=30.0)
            ax.set_xlabel("Doping concentration (cm$^{-3}$)", fontsize=30.0)

            p = "lower right" if idx == 0 else "best"
            ax.legend(loc=p, fontsize=15)
            ax.grid()
            ax.set_xticks(fontsize=25)
            ax.set_yticks(fontsize=25)

        plt.tight_layout()

        return ax

    def plot_conductivity_dop(
        self,
        temps="all",
        output: Literal["average", "eigs"] = "average",
        relaxation_time=1e-14,
    ):
        """Plot the conductivity in function of doping levels for different
        temperatures.

        Args:
            temps: the default 'all' plots all the temperatures in the analyzer.
                Specify a list of temperatures if you want to plot only some.
            output: with 'average' you get an average of the three directions
                with 'eigs' you get all the three directions.
            relaxation_time: specify a constant relaxation time value

        Returns:
            a matplotlib object
        """
        if output == "average":
            cond = self._bz.get_conductivity(relaxation_time=relaxation_time, output="average")
        else:
            cond = self._bz.get_conductivity(relaxation_time=relaxation_time, output="eigs")

        temperatures = sorted(cond["n"]) if temps == "all" else temps
        ax = pretty_plot(22, 14)
        for idx, dop_type in enumerate("np"):
            plt.subplot(121 + idx)
            for temp in temperatures:
                if output == "eigs":
                    for xyz in range(3):
                        ax.semilogx(
                            self._bz.doping[dop_type],
                            list(zip(*cond[dop_type][temp], strict=True))[xyz],
                            marker="s",
                            label=f"{xyz} {temp} K",
                        )
                elif output == "average":
                    ax.semilogx(
                        self._bz.doping[dop_type],
                        cond[dop_type][temp],
                        marker="s",
                        label=f"{temp} K",
                    )
            ax.set_title(f"{dop_type}-type", fontsize=20)
            if idx == 0:
                ax.set_ylabel("conductivity $\\sigma$ (1/($\\Omega$ m))", fontsize=30.0)
            ax.set_xlabel("Doping concentration ($cm^{-3}$)", fontsize=30.0)
            ax.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
            ax.legend(fontsize=15)
            ax.grid()
            ax.set_xticks(fontsize=25)
            ax.set_yticks(fontsize=25)

        plt.tight_layout()

        return ax

    def plot_power_factor_dop(
        self,
        temps="all",
        output: Literal["average", "eigs"] = "average",
        relaxation_time=1e-14,
    ):
        """Plot the Power Factor in function of doping levels for different temperatures.

        Args:
            temps: the default 'all' plots all the temperatures in the analyzer.
                Specify a list of temperatures if you want to plot only some.
            output: with 'average' you get an average of the three directions
                with 'eigs' you get all the three directions.
            relaxation_time: specify a constant relaxation time value

        Returns:
            a matplotlib object
        """
        if output == "average":
            pow_factor = self._bz.get_power_factor(relaxation_time=relaxation_time, output="average")
        else:
            pow_factor = self._bz.get_power_factor(relaxation_time=relaxation_time, output="eigs")

        temperatures = sorted(pow_factor["n"]) if temps == "all" else temps
        ax = pretty_plot(22, 14)
        for idx, dop_type in enumerate("np"):
            plt.subplot(121 + idx)
            for temp in temperatures:
                if output == "eigs":
                    for xyz in range(3):
                        ax.semilogx(
                            self._bz.doping[dop_type],
                            list(zip(*pow_factor[dop_type][temp], strict=True))[xyz],
                            marker="s",
                            label=f"{xyz} {temp} K",
                        )
                elif output == "average":
                    ax.semilogx(
                        self._bz.doping[dop_type],
                        pow_factor[dop_type][temp],
                        marker="s",
                        label=f"{temp} K",
                    )
            ax.set_title(f"{dop_type}-type", fontsize=20)
            if idx == 0:
                ax.set_ylabel("Power Factor  ($\\mu$W/(mK$^2$))", fontsize=30.0)
            ax.set_xlabel("Doping concentration ($cm^{-3}$)", fontsize=30.0)
            ax.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
            p = "best"  # 'lower right' if i == 0 else ''
            ax.legend(loc=p, fontsize=15)
            ax.grid()
            ax.set_xticks(fontsize=25)
            ax.set_yticks(fontsize=25)

        plt.tight_layout()

        return ax

    def plot_zt_dop(
        self,
        temps="all",
        output: Literal["average", "eigs"] = "average",
        relaxation_time=1e-14,
    ):
        """Plot the figure of merit zT in function of doping levels for different
        temperatures.

        Args:
            temps: the default 'all' plots all the temperatures in the analyzer.
                Specify a list of temperatures if you want to plot only some.
            output: with 'average' you get an average of the three directions
                with 'eigs' you get all the three directions.
            relaxation_time: specify a constant relaxation time value

        Returns:
            a matplotlib object
        """
        if output == "average":
            zt = self._bz.get_zt(relaxation_time=relaxation_time, output="average")
        else:
            zt = self._bz.get_zt(relaxation_time=relaxation_time, output="eigs")

        temperatures = sorted(zt["n"]) if temps == "all" else temps
        ax = pretty_plot(22, 14)
        for idx, dop_type in enumerate("np"):
            plt.subplot(121 + idx)
            for temp in temperatures:
                if output == "eigs":
                    for xyz in range(3):
                        ax.semilogx(
                            self._bz.doping[dop_type],
                            list(zip(*zt[dop_type][temp], strict=True))[xyz],
                            marker="s",
                            label=f"{xyz} {temp} K",
                        )
                elif output == "average":
                    ax.semilogx(
                        self._bz.doping[dop_type],
                        zt[dop_type][temp],
                        marker="s",
                        label=f"{temp} K",
                    )
            ax.set_title(f"{dop_type}-type", fontsize=20)
            if idx == 0:
                ax.set_ylabel("zT", fontsize=30.0)
            ax.set_xlabel("Doping concentration ($cm^{-3}$)", fontsize=30.0)

            p = "lower right" if idx == 0 else "best"
            ax.legend(loc=p, fontsize=15)
            ax.grid()
            ax.set_xticks(fontsize=25)
            ax.set_yticks(fontsize=25)

        plt.tight_layout()

        return ax

    def plot_eff_mass_dop(self, temps="all", output: Literal["average", "eigs"] = "average"):
        """Plot the average effective mass in function of doping levels
        for different temperatures.

        Args:
            temps: the default 'all' plots all the temperatures in the analyzer.
                Specify a list of temperatures if you want to plot only some.
            output: with 'average' you get an average of the three directions
                with 'eigs' you get all the three directions.
            relaxation_time: specify a constant relaxation time value

        Returns:
            a matplotlib object
        """
        if output == "average":
            em = self._bz.get_average_eff_mass(output="average")
        else:
            em = self._bz.get_average_eff_mass(output="eigs")

        temperatures = sorted(em["n"]) if temps == "all" else temps
        ax = pretty_plot(22, 14)
        for idx, dop_type in enumerate("np"):
            plt.subplot(121 + idx)
            for temp in temperatures:
                if output == "eigs":
                    for xyz in range(3):
                        ax.semilogx(
                            self._bz.doping[dop_type],
                            list(zip(*em[dop_type][temp], strict=True))[xyz],
                            marker="s",
                            label=f"{xyz} {temp} K",
                        )
                elif output == "average":
                    ax.semilogx(
                        self._bz.doping[dop_type],
                        em[dop_type][temp],
                        marker="s",
                        label=f"{temp} K",
                    )
            ax.set_title(f"{dop_type}-type", fontsize=20)
            if idx == 0:
                ax.set_ylabel("Effective mass (m$_e$)", fontsize=30.0)
            ax.set_xlabel("Doping concentration ($cm^{-3}$)", fontsize=30.0)

            p = "lower right" if idx == 0 else "best"
            ax.legend(loc=p, fontsize=15)
            ax.grid()
            ax.set_xticks(fontsize=25)
            ax.set_yticks(fontsize=25)

        plt.tight_layout()

        return ax

    def plot_dos(self, sigma=0.05):
        """Plot dos.

        Args:
            sigma: a smearing

        Returns:
            a matplotlib object
        """
        plotter = DosPlotter(sigma=sigma)
        plotter.add_dos("t", self._bz.dos)
        return plotter.get_plot()

    def plot_carriers(self, temp=300):
        """Plot the carrier concentration in function of Fermi level.

        Args:
            temp: the temperature

        Returns:
            a matplotlib object
        """
        ax = pretty_plot(9, 7)
        carriers = [abs(c / (self._bz.vol * 1e-24)) for c in self._bz._carrier_conc[temp]]
        ax.semilogy(self._bz.mu_steps, carriers, linewidth=3.0, color="r")
        self._plot_bg_limits(ax)
        self._plot_doping(ax, temp)
        ax.set_xlim(-0.5, self._bz.gap + 0.5)
        ax.set_ylim(1e14, 1e22)
        ax.set_ylabel("carrier concentration (cm-3)", fontsize=30.0)
        ax.set_xlabel("E-E$_f$ (eV)", fontsize=30)
        ax.set_xticks(fontsize=25)
        ax.set_yticks(fontsize=25)
        plt.tight_layout()
        return ax

    def plot_hall_carriers(self, temp=300):
        """Plot the Hall carrier concentration in function of Fermi level.

        Args:
            temp: the temperature

        Returns:
            a matplotlib object
        """
        ax = pretty_plot(9, 7)
        hall_carriers = [abs(i) for i in self._bz.get_hall_carrier_concentration()[temp]]
        ax.semilogy(self._bz.mu_steps, hall_carriers, linewidth=3.0, color="r")
        self._plot_bg_limits(ax)
        self._plot_doping(ax, temp)
        ax.set_xlim(-0.5, self._bz.gap + 0.5)
        ax.set_ylim(1e14, 1e22)
        ax.set_ylabel("Hall carrier concentration (cm-3)", fontsize=30.0)
        ax.set_xlabel("E-E$_f$ (eV)", fontsize=30)
        ax.set_xticks(fontsize=25)
        ax.set_yticks(fontsize=25)
        plt.tight_layout()
        return ax


class CohpPlotter:
    """Plot crystal orbital Hamilton populations (COHPs) or crystal orbital overlap
    populations (COOPs). It is modeled after the DosPlotter object.
    """

    def __init__(self, zero_at_efermi=True, are_coops=False, are_cobis=False) -> None:
        """
        Args:
            zero_at_efermi: Whether to shift all populations to have zero
                energy at the Fermi level. Defaults to True.
            are_coops: Switch to indicate that these are COOPs, not COHPs.
                Defaults to False for COHPs.
            are_cobis: Switch to indicate that these are COBIs or multi-center COBIs, not COHPs/COOPs.
                Defaults to False for COHPs.
        """
        self.zero_at_efermi = zero_at_efermi
        self.are_coops = are_coops
        self.are_cobis = are_cobis
        self._cohps: dict[str, dict[str, np.ndarray | dict[Spin, np.ndarray] | float]] = {}

    def add_cohp(self, label, cohp) -> None:
        """Add a COHP for plotting.

        Args:
            label: Label for the COHP. Must be unique.

            cohp: COHP object.
        """
        energies = cohp.energies - cohp.efermi if self.zero_at_efermi else cohp.energies
        populations = cohp.get_cohp()
        int_populations = cohp.get_icohp()
        self._cohps[label] = {
            "energies": energies,
            "COHP": populations,
            "ICOHP": int_populations,
            "efermi": cohp.efermi,
        }

    def add_cohp_dict(self, cohp_dict, key_sort_func=None) -> None:
        """Add a dictionary of COHPs with an optional sorting function
        for the keys.

        Args:
            cohp_dict: dict of the form {label: Cohp}

            key_sort_func: function used to sort the cohp_dict keys.
        """
        keys = sorted(cohp_dict, key=key_sort_func) if key_sort_func else list(cohp_dict)
        for label in keys:
            self.add_cohp(label, cohp_dict[label])

    def get_cohp_dict(self):
        """Get the added COHPs as a json-serializable dict. Note that if you
        have specified smearing for the COHP plot, the populations returned
        will be the smeared and not the original populations.

        Returns:
            dict: Dict of COHP data of the form {label: {"efermi": efermi,
            "energies": ..., "COHP": {Spin.up: ...}, "ICOHP": ...}}.
        """
        return jsanitize(self._cohps)

    def get_plot(
        self,
        xlim=None,
        ylim=None,
        plot_negative=None,
        integrated=False,
        invert_axes=True,
    ):
        """Get a matplotlib plot showing the COHP.

        Args:
            xlim: Specifies the x-axis limits. Defaults to None for
                automatic determination.

            ylim: Specifies the y-axis limits. Defaults to None for
                automatic determination.

            plot_negative: It is common to plot -COHP(E) so that the
                sign means the same for COOPs and COHPs. Defaults to None
                for automatic determination: If are_coops is True, this
                will be set to False, else it will be set to True.

            integrated: Switch to plot ICOHPs. Defaults to False.

            invert_axes: Put the energies onto the y-axis, which is
                common in chemistry.

        Returns:
            A matplotlib object.
        """
        if self.are_coops:
            cohp_label = "COOP"
        elif self.are_cobis:
            cohp_label = "COBI"
        else:
            cohp_label = "COHP"

        if plot_negative is None:
            plot_negative = (not self.are_coops) and (not self.are_cobis)

        if integrated:
            cohp_label = f"I{cohp_label} (eV)"

        if plot_negative:
            cohp_label = f"-{cohp_label}"

        energy_label = "$E - E_f$ (eV)" if self.zero_at_efermi else "$E$ (eV)"

        ncolors = max(3, len(self._cohps))
        ncolors = min(9, ncolors)

        colors = palettable.colorbrewer.qualitative.Set1_9.mpl_colors

        ax = pretty_plot(12, 8)

        allpts = []
        keys = list(self._cohps)
        idx = key = None

        for idx, key in enumerate(keys):
            energies = self._cohps[key]["energies"]
            populations = self._cohps[key]["COHP"] if not integrated else self._cohps[key]["ICOHP"]
            for spin in [Spin.up, Spin.down]:
                if spin in populations:
                    if invert_axes:
                        x = -populations[spin] if plot_negative else populations[spin]
                        y = energies
                    else:
                        x = energies
                        y = -populations[spin] if plot_negative else populations[spin]
                    allpts.extend(list(zip(x, y, strict=True)))
                    if spin == Spin.up:
                        ax.plot(
                            x,
                            y,
                            color=colors[idx % ncolors],
                            linestyle="-",
                            label=str(key),
                            linewidth=3,
                        )
                    else:
                        ax.plot(
                            x,
                            y,
                            color=colors[idx % ncolors],
                            linestyle="--",
                            linewidth=3,
                        )

        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)
        elif not invert_axes:
            xlim = ax.get_xlim()
            relevant_y = [p[1] for p in allpts if xlim[0] < p[0] < xlim[1]]
            ax.set_ylim((min(relevant_y), max(relevant_y)))
        if not xlim and invert_axes:
            ylim = ax.get_ylim()
            relevant_y = [p[0] for p in allpts if ylim[0] < p[1] < ylim[1]]
            ax.set_xlim((min(relevant_y), max(relevant_y)))

        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        if not invert_axes:
            ax.axhline(y=0, color="k", linewidth=2)
            if self.zero_at_efermi:
                ax.plot([0, 0], ylim, "k--", linewidth=2)
            else:
                ax.plot(
                    [self._cohps[key]["efermi"], self._cohps[key]["efermi"]],
                    ylim,
                    color=colors[idx % ncolors],
                    linestyle="--",
                    linewidth=2,
                )
        else:
            ax.axvline(x=0, color="k", linewidth=2)
            if self.zero_at_efermi:
                ax.plot(xlim, [0, 0], "k--", linewidth=2)
            else:
                ax.plot(
                    xlim,
                    [self._cohps[key]["efermi"], self._cohps[key]["efermi"]],
                    color=colors[idx % ncolors],
                    linestyle="--",
                    linewidth=2,
                )

        if invert_axes:
            ax.set_xlabel(cohp_label)
            ax.set_ylabel(energy_label)
        else:
            ax.set_xlabel(energy_label)
            ax.set_ylabel(cohp_label)

        ax.legend()
        legend_text = ax.legend().get_texts()
        plt.setp(legend_text, fontsize=30)
        plt.tight_layout()
        return ax

    def save_plot(self, filename: str, xlim=None, ylim=None) -> None:
        """Save matplotlib plot to a file.

        Args:
            filename (str): File name to write to. Must include extension to specify image format.
            xlim: Specifies the x-axis limits. Defaults to None for
                automatic determination.
            ylim: Specifies the y-axis limits. Defaults to None for
                automatic determination.
        """
        self.get_plot(xlim, ylim)
        plt.savefig(filename)

    def show(self, xlim=None, ylim=None) -> None:
        """Show the plot using matplotlib.

        Args:
            xlim: Specifies the x-axis limits. Defaults to None for
                automatic determination.
            ylim: Specifies the y-axis limits. Defaults to None for
                automatic determination.
        """
        self.get_plot(xlim, ylim)
        plt.show()


@requires(mlab is not None, "MayAvi mlab not installed! Please install mayavi.")
def plot_fermi_surface(
    data,
    structure,
    cbm,
    energy_levels=None,
    multiple_figure=True,
    mlab_figure=None,
    kpoints_dict=None,
    colors=None,
    transparency_factor=None,
    labels_scale_factor=0.05,
    points_scale_factor=0.02,
    interactive=True,
):
    """Plot the Fermi surface at specific energy value using Boltztrap 1 FERMI
    mode.

    The easiest way to use this plotter is:

        1. Run boltztrap in 'FERMI' mode using BoltztrapRunner,
        2. Load BoltztrapAnalyzer using your method of choice (e.g., from_files)
        3. Pass in your BoltztrapAnalyzer's fermi_surface_data as this
            function's data argument.

    Args:
        data: energy values in a 3D grid from a CUBE file via read_cube_file
            function, or from a BoltztrapAnalyzer.fermi_surface_data
        structure: structure object of the material
        energy_levels ([float]): Energy values for plotting the fermi surface(s)
            By default 0 eV correspond to the VBM, as in the plot of band
            structure along symmetry line.
            Default: One surface, with max energy value + 0.01 eV
        cbm (bool): True if the considered band is a conduction band or not.
        multiple_figure (bool): If True a figure for each energy level will be
            shown. If False all the surfaces will be shown in the same figure.
            In this last case, tune the transparency factor.
        mlab_figure (mayavi.mlab.figure): A previous figure to plot a new
            surface on.
        kpoints_dict (dict): dictionary of kpoints to label in the plot.
            Example: {"K":[0.5,0.0,0.5]}, coords are fractional
        colors ([tuple]): Iterable of 3-tuples (r,g,b) of integers to define
            the colors of each surface (one per energy level).
            Should be the same length as the number of surfaces being plotted.
            Example (3 surfaces): colors=[(1,0,0), (0,1,0), (0,0,1)]
            Example (2 surfaces): colors=[(0, 0.5, 0.5)]
        transparency_factor (float): Values in the range [0,1] to tune the
            opacity of each surface. Should be one transparency_factor per
            surface.
        labels_scale_factor (float): factor to tune size of the kpoint labels
        points_scale_factor (float): factor to tune size of the kpoint points
        interactive (bool): if True an interactive figure will be shown.
            If False a non interactive figure will be shown, but it is possible
            to plot other surfaces on the same figure. To make it interactive,
            run mlab.show().

    Returns:
        tuple[mlab.figure, mlab]: The mlab plotter and an interactive
            figure to control the plot.

    Note: Experimental.
        Please, double check the surface shown by using some other software and report issues.
    """
    bz = structure.lattice.reciprocal_lattice.get_wigner_seitz_cell()
    cell = structure.lattice.reciprocal_lattice.matrix

    fact = 1 if not cbm else -1
    data_1d = data.ravel()
    en_min = np.min(fact * data_1d)
    en_max = np.max(fact * data_1d)

    if energy_levels is None:
        energy_levels = [en_min + 0.01] if cbm else [en_max - 0.01]
        print(f"Energy level set to: {energy_levels[0]} eV")

    else:
        for e in energy_levels:
            if e > en_max or e < en_min:
                raise BoltztrapError(f"energy level {e} not in the range of possible energies: [{en_min}, {en_max}]")

    n_surfaces = len(energy_levels)
    if colors is None:
        colors = [(0, 0, 1)] * n_surfaces

    if transparency_factor is None:
        transparency_factor = [1] * n_surfaces

    fig = mlab_figure or None

    if kpoints_dict is None:
        kpoints_dict = {}

    if mlab_figure is None and not multiple_figure:
        fig = mlab.figure(size=(1024, 768), bgcolor=(1, 1, 1))
        for iface in range(len(bz)):
            for line in itertools.combinations(bz[iface], 2):
                for jface in range(len(bz)):
                    if (
                        iface < jface
                        and any(np.all(line[0] == x) for x in bz[jface])
                        and any(np.all(line[1] == x) for x in bz[jface])
                    ):
                        mlab.plot3d(
                            *zip(line[0], line[1], strict=True),
                            color=(0, 0, 0),
                            tube_radius=None,
                            figure=fig,
                        )
        for key, coords in kpoints_dict.items():
            label_coords = structure.lattice.reciprocal_lattice.get_cartesian_coords(coords)
            mlab.points3d(
                *label_coords,
                scale_factor=points_scale_factor,
                color=(0, 0, 0),
                figure=fig,
            )
            mlab.text3d(
                *label_coords,
                text=key,
                scale=labels_scale_factor,
                color=(0, 0, 0),
                figure=fig,
            )

    for i, isolevel in enumerate(energy_levels):
        alpha = transparency_factor[i]
        color = colors[i]
        if multiple_figure:
            fig = mlab.figure(size=(1024, 768), bgcolor=(1, 1, 1))

            for iface in range(len(bz)):
                for line in itertools.combinations(bz[iface], 2):
                    for jface in range(len(bz)):
                        if (
                            iface < jface
                            and any(np.all(line[0] == x) for x in bz[jface])
                            and any(np.all(line[1] == x) for x in bz[jface])
                        ):
                            mlab.plot3d(
                                *zip(line[0], line[1], strict=True),
                                color=(0, 0, 0),
                                tube_radius=None,
                                figure=fig,
                            )

            for key, coords in kpoints_dict.items():
                label_coords = structure.lattice.reciprocal_lattice.get_cartesian_coords(coords)
                mlab.points3d(
                    *label_coords,
                    scale_factor=points_scale_factor,
                    color=(0, 0, 0),
                    figure=fig,
                )
                mlab.text3d(
                    *label_coords,
                    text=key,
                    scale=labels_scale_factor,
                    color=(0, 0, 0),
                    figure=fig,
                )

        cp = mlab.contour3d(
            fact * data,
            contours=[isolevel],
            transparent=True,
            colormap="hot",
            color=color,
            opacity=alpha,
            figure=fig,
        )

        polydata = cp.actor.actors[0].mapper.input
        pts = np.array(polydata.points)  # - 1
        polydata.points = np.dot(pts, cell / np.array(data.shape)[:, None])

        cx, cy, cz = (np.mean(np.array(polydata.points)[:, i]) for i in range(3))

        polydata.points = (np.array(polydata.points) - [cx, cy, cz]) * 2

        # mlab.view(distance='auto')
        fig.scene.isometric_view()

    if interactive:
        mlab.show()

    return fig, mlab


def plot_wigner_seitz(lattice, ax: plt.Axes = None, **kwargs):
    """Add the skeleton of the Wigner-Seitz cell of the lattice to a matplotlib Axes.

    Args:
        lattice: Lattice object
        ax: matplotlib Axes or None if a new figure should be created.
        kwargs: kwargs passed to the matplotlib function 'plot'. Color defaults to black
            and linewidth to 1.

    Returns:
        matplotlib figure and matplotlib ax
    """
    ax, fig = get_ax3d_fig(ax)

    kwargs.setdefault("color", "k")
    kwargs.setdefault("linewidth", 1)

    bz = lattice.get_wigner_seitz_cell()
    for iface in range(len(bz)):
        for line in itertools.combinations(bz[iface], 2):
            for jface in range(len(bz)):
                if (
                    iface < jface
                    and any(np.all(line[0] == x) for x in bz[jface])
                    and any(np.all(line[1] == x) for x in bz[jface])
                ):
                    ax.plot(*zip(line[0], line[1], strict=True), **kwargs)

    return fig, ax


def plot_lattice_vectors(lattice, ax: plt.Axes = None, **kwargs):
    """Add the basis vectors of the lattice provided to a matplotlib Axes.

    Args:
        lattice: Lattice object
        ax: matplotlib Axes or None if a new figure should be created.
        kwargs: kwargs passed to the matplotlib function 'plot'. Color defaults to green
            and linewidth to 3.

    Returns:
        matplotlib figure and matplotlib ax
    """
    ax, fig = get_ax3d_fig(ax)

    if "color" not in kwargs:
        kwargs["color"] = "g"
    if "linewidth" not in kwargs:
        kwargs["linewidth"] = 3

    vertex1 = lattice.get_cartesian_coords([0.0, 0.0, 0.0])
    vertex2 = lattice.get_cartesian_coords([1.0, 0.0, 0.0])
    ax.plot(*zip(vertex1, vertex2, strict=True), **kwargs)
    vertex2 = lattice.get_cartesian_coords([0.0, 1.0, 0.0])
    ax.plot(*zip(vertex1, vertex2, strict=True), **kwargs)
    vertex2 = lattice.get_cartesian_coords([0.0, 0.0, 1.0])
    ax.plot(*zip(vertex1, vertex2, strict=True), **kwargs)

    return fig, ax


def plot_path(line, lattice=None, coords_are_cartesian=False, ax: plt.Axes = None, **kwargs):
    """Add a line passing through the coordinates listed in 'line' to a matplotlib Axes.

    Args:
        line: list of coordinates.
        lattice: Lattice object used to convert from reciprocal to Cartesian coordinates
        coords_are_cartesian: Set to True if you are providing
            coordinates in Cartesian coordinates. Defaults to False.
            Requires lattice if False.
        ax: matplotlib Axes or None if a new figure should be created.
        kwargs: kwargs passed to the matplotlib function 'plot'. Color defaults to red
            and linewidth to 3.

    Returns:
        matplotlib figure and matplotlib ax
    """
    ax, fig = get_ax3d_fig(ax)

    if "color" not in kwargs:
        kwargs["color"] = "r"
    if "linewidth" not in kwargs:
        kwargs["linewidth"] = 3

    for k in range(1, len(line)):
        vertex1 = line[k - 1]
        vertex2 = line[k]
        if not coords_are_cartesian:
            if lattice is None:
                raise ValueError("coords_are_cartesian False requires the lattice")
            vertex1 = lattice.get_cartesian_coords(vertex1)
            vertex2 = lattice.get_cartesian_coords(vertex2)
        ax.plot(*zip(vertex1, vertex2, strict=True), **kwargs)

    return fig, ax


def plot_labels(labels, lattice=None, coords_are_cartesian=False, ax: plt.Axes = None, **kwargs):
    """Add labels to a matplotlib Axes.

    Args:
        labels: dict containing the label as a key and the coordinates as value.
        lattice: Lattice object used to convert from reciprocal to Cartesian coordinates
        coords_are_cartesian: Set to True if you are providing.
            coordinates in Cartesian coordinates. Defaults to False.
            Requires lattice if False.
        ax: matplotlib Axes or None if a new figure should be created.
        kwargs: kwargs passed to the matplotlib function 'text'. Color defaults to blue
            and size to 25.

    Returns:
        matplotlib figure and matplotlib ax
    """
    ax, fig = get_ax3d_fig(ax)

    if "color" not in kwargs:
        kwargs["color"] = "b"
    if "size" not in kwargs:
        kwargs["size"] = 25

    for k, coords in labels.items():
        label = k
        if k.startswith("\\") or k.find("_") != -1:
            label = f"${k}$"
        off = 0.01
        if coords_are_cartesian:
            coords = np.array(coords)
        else:
            if lattice is None:
                raise ValueError("coords_are_cartesian False requires the lattice")
            coords = lattice.get_cartesian_coords(coords)
        ax.text(*(coords + off), s=label, **kwargs)

    return fig, ax


def fold_point(p, lattice, coords_are_cartesian=False):
    """Folds a point with coordinates p inside the first Brillouin zone of the lattice.

    Args:
        p: coordinates of one point
        lattice: Lattice object used to convert from reciprocal to Cartesian coordinates
        coords_are_cartesian: Set to True if you are providing
            coordinates in Cartesian coordinates. Defaults to False.

    Returns:
        The Cartesian coordinates folded inside the first Brillouin zone
    """
    p = lattice.get_fractional_coords(p) if coords_are_cartesian else np.array(p)

    p = np.mod(p + 0.5 - 1e-10, 1) - 0.5 + 1e-10
    p = lattice.get_cartesian_coords(p)

    closest_lattice_point = None
    smallest_distance = 10000
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            for k in (-1, 0, 1):
                lattice_point = np.dot((i, j, k), lattice.matrix)
                dist = np.linalg.norm(p - lattice_point)
                if closest_lattice_point is None or dist < smallest_distance:
                    closest_lattice_point = lattice_point
                    smallest_distance = dist

    if not np.allclose(closest_lattice_point, (0, 0, 0)):
        p -= closest_lattice_point

    return p


def plot_points(
    points,
    lattice=None,
    coords_are_cartesian=False,
    fold=False,
    ax: plt.Axes = None,
    **kwargs,
):
    """Add Points to a matplotlib Axes.

    Args:
        points: list of coordinates
        lattice: Lattice object used to convert from reciprocal to Cartesian coordinates
        coords_are_cartesian: Set to True if you are providing
            coordinates in Cartesian coordinates. Defaults to False.
            Requires lattice if False.
        fold: whether the points should be folded inside the first Brillouin Zone.
            Defaults to False. Requires lattice if True.
        ax: matplotlib Axes or None if a new figure should be created.
        kwargs: kwargs passed to the matplotlib function 'scatter'. Color defaults to blue

    Returns:
        matplotlib figure and matplotlib ax
    """
    ax, fig = get_ax3d_fig(ax)

    if "color" not in kwargs:
        kwargs["color"] = "b"

    if (not coords_are_cartesian or fold) and lattice is None:
        raise ValueError("coords_are_cartesian False or fold True require the lattice")

    for p in points:
        if fold:
            p = fold_point(p, lattice, coords_are_cartesian=coords_are_cartesian)

        elif not coords_are_cartesian:
            p = lattice.get_cartesian_coords(p)

        ax.scatter(*p, **kwargs)

    return fig, ax


@add_fig_kwargs
def plot_brillouin_zone_from_kpath(kpath, ax: plt.Axes = None, **kwargs) -> plt.Axes:
    """Get the plot (as a matplotlib object) of the symmetry line path in
        the Brillouin Zone.

    Args:
        kpath (HighSymmKpath): a HighSymmKPath object
        ax: matplotlib Axes or None if a new figure should be created.
        **kwargs: provided by add_fig_kwargs decorator

    Returns:
        plt.Axes: matplotlib Axes
    """
    lines = [[kpath.kpath["kpoints"][k] for k in p] for p in kpath.kpath["path"]]
    return plot_brillouin_zone(
        bz_lattice=kpath.prim_rec,
        lines=lines,
        ax=ax,
        labels=kpath.kpath["kpoints"],
        **kwargs,
    )


@add_fig_kwargs
def plot_brillouin_zone(
    bz_lattice,
    lines=None,
    labels=None,
    kpoints=None,
    fold=False,
    coords_are_cartesian: bool = False,
    ax: plt.Axes = None,
    **kwargs,
):
    """Plots a 3D representation of the Brillouin zone of the structure.
    Can add to the plot paths, labels and kpoints.

    Args:
        bz_lattice: Lattice object of the Brillouin zone
        lines: list of lists of coordinates. Each list represent a different path
        labels: dict containing the label as a key and the coordinates as value.
        kpoints: list of coordinates
        fold: whether the points should be folded inside the first Brillouin Zone.
            Defaults to False. Requires lattice if True.
        coords_are_cartesian: Set to True if you are providing
            coordinates in Cartesian coordinates. Defaults to False.
        ax: matplotlib Axes or None if a new figure should be created.
        kwargs: provided by add_fig_kwargs decorator

    Returns:
        matplotlib figure
    """
    fig, ax = plot_lattice_vectors(bz_lattice, ax=ax)
    plot_wigner_seitz(bz_lattice, ax=ax)
    if lines is not None:
        for line in lines:
            plot_path(line, bz_lattice, coords_are_cartesian=coords_are_cartesian, ax=ax)

    if labels is not None:
        plot_labels(labels, bz_lattice, coords_are_cartesian=coords_are_cartesian, ax=ax)
        plot_points(
            labels.values(),
            bz_lattice,
            coords_are_cartesian=coords_are_cartesian,
            fold=False,
            ax=ax,
        )

    if kpoints is not None:
        plot_points(
            kpoints,
            bz_lattice,
            coords_are_cartesian=coords_are_cartesian,
            ax=ax,
            fold=fold,
        )

    ax.set_xlim3d(-1, 1)
    ax.set_ylim3d(-1, 1)
    ax.set_zlim3d(-1, 1)

    # ax.set_aspect('equal')
    ax.axis("off")

    return fig


def plot_ellipsoid(
    hessian,
    center,
    lattice=None,
    rescale=1.0,
    ax: plt.Axes = None,
    coords_are_cartesian=False,
    arrows=False,
    **kwargs,
):
    """Plots a 3D ellipsoid rappresenting the Hessian matrix in input.
    Useful to get a graphical visualization of the effective mass
    of a band in a single k-point.

    Args:
        hessian: the Hessian matrix
        center: the center of the ellipsoid in reciprocal coords (Default)
        lattice: Lattice object of the Brillouin zone
        rescale: factor for size scaling of the ellipsoid
        ax: matplotlib Axes or None if a new figure should be created.
        coords_are_cartesian: Set to True if you are providing a center in
            Cartesian coordinates. Defaults to False.
        arrows: whether to plot arrows for the principal axes of the ellipsoid. Defaults to False.
        **kwargs: passed to the matplotlib function 'plot_wireframe'.
            Color defaults to blue, rstride and cstride
            default to 4, alpha defaults to 0.2.

    Returns:
        matplotlib figure and matplotlib ax

    Example of use:
        fig,ax=plot_wigner_seitz(struct.reciprocal_lattice)
        plot_ellipsoid(hessian,[0.0,0.0,0.0], struct.reciprocal_lattice,ax=ax)
    """
    if (not coords_are_cartesian) and lattice is None:
        raise ValueError("coords_are_cartesian False or fold True require the lattice")

    if not coords_are_cartesian:
        center = lattice.get_cartesian_coords(center)

    if "color" not in kwargs:
        kwargs["color"] = "b"
    if "rstride" not in kwargs:
        kwargs["rstride"] = 4
    if "cstride" not in kwargs:
        kwargs["cstride"] = 4
    if "alpha" not in kwargs:
        kwargs["alpha"] = 0.2

    # calculate the ellipsoid
    # find the rotation matrix and radii of the axes
    _U, s, rotation = np.linalg.svd(hessian)
    radii = 1.0 / np.sqrt(s)

    # from polar coordinates
    u = np.linspace(0.0, 2.0 * np.pi, 100)
    v = np.linspace(0.0, np.pi, 100)
    x = radii[0] * np.outer(np.cos(u), np.sin(v))
    y = radii[1] * np.outer(np.sin(u), np.sin(v))
    z = radii[2] * np.outer(np.ones_like(u), np.cos(v))
    for i in range(len(x)):
        for j in range(len(x)):
            [x[i, j], y[i, j], z[i, j]] = np.dot([x[i, j], y[i, j], z[i, j]], rotation) * rescale + center

    # add the ellipsoid to the current axes
    ax, fig = get_ax3d_fig(ax)
    ax.plot_wireframe(x, y, z, **kwargs)

    if arrows:
        color = ("b", "g", "r")
        em = np.zeros((3, 3))
        for i in range(3):
            em[i, :] = rotation[i, :] / np.linalg.norm(rotation[i, :])
        for i in range(3):
            ax.quiver3D(
                center[0],
                center[1],
                center[2],
                em[i, 0],
                em[i, 1],
                em[i, 2],
                pivot="tail",
                arrow_length_ratio=0.2,
                length=radii[i] * rescale,
                color=color[i],
            )

    return fig, ax
