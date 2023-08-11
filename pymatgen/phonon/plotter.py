"""This module implements plotter for DOS and band structure."""

from __future__ import annotations

import logging
from collections import namedtuple

import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as const
from monty.json import jsanitize

from pymatgen.electronic_structure.plotter import plot_brillouin_zone
from pymatgen.phonon.bandstructure import PhononBandStructureSymmLine
from pymatgen.phonon.gruneisen import GruneisenPhononBandStructureSymmLine
from pymatgen.util.plotting import add_fig_kwargs, get_ax_fig_plt, pretty_plot

logger = logging.getLogger(__name__)

FreqUnits = namedtuple("FreqUnits", ["factor", "label"])


def freq_units(units):
    """

    Args:
        units: str, accepted values: thz, ev, mev, ha, cm-1, cm^-1.

    Returns:
        Returns conversion factor from THz to the required units and the label in the form of a namedtuple

    """
    d = {
        "thz": FreqUnits(1, "THz"),
        "ev": FreqUnits(const.value("hertz-electron volt relationship") * const.tera, "eV"),
        "mev": FreqUnits(
            const.value("hertz-electron volt relationship") * const.tera / const.milli,
            "meV",
        ),
        "ha": FreqUnits(const.value("hertz-hartree relationship") * const.tera, "Ha"),
        "cm-1": FreqUnits(
            const.value("hertz-inverse meter relationship") * const.tera * const.centi,
            "cm^{-1}",
        ),
        "cm^-1": FreqUnits(
            const.value("hertz-inverse meter relationship") * const.tera * const.centi,
            "cm^{-1}",
        ),
    }
    try:
        return d[units.lower().strip()]
    except KeyError:
        raise KeyError(f"Value for units `{units}` unknown\nPossible values are:\n {list(d)}")


class PhononDosPlotter:
    """
    Class for plotting phonon DOSs. The interface is extremely flexible given there are many
    different ways in which people want to view DOS.
    Typical usage is:
        # Initializes plotter with some optional args. Defaults are usually fine
        plotter = PhononDosPlotter().

        # Add DOS with a label
        plotter.add_dos("Total DOS", dos)

        # Alternatively, you can add a dict of DOSes. This is the typical form
        # returned by CompletePhononDos.get_element_dos().
    """

    def __init__(self, stack=False, sigma=None):
        """
        Args:
            stack: Whether to plot the DOS as a stacked area graph
            sigma: A float specifying a standard deviation for Gaussian smearing
                the DOS for nicer looking plots. Defaults to None for no smearing.
        """
        self.stack = stack
        self.sigma = sigma
        self._doses = {}

    def add_dos(self, label, dos):
        """
        Adds a dos for plotting.

        Args:
            label:
                label for the DOS. Must be unique.
            dos:
                PhononDos object
        """
        densities = dos.get_smeared_densities(self.sigma) if self.sigma else dos.densities
        self._doses[label] = {"frequencies": dos.frequencies, "densities": densities}

    def add_dos_dict(self, dos_dict, key_sort_func=None):
        """
        Add a dictionary of doses, with an optional sorting function for the
        keys.

        Args:
            dos_dict: dict of {label: Dos}
            key_sort_func: function used to sort the dos_dict keys.
        """
        keys = sorted(dos_dict, key=key_sort_func) if key_sort_func else list(dos_dict)
        for label in keys:
            self.add_dos(label, dos_dict[label])

    def get_dos_dict(self):
        """
        Returns the added doses as a json-serializable dict. Note that if you
        have specified smearing for the DOS plot, the densities returned will
        be the smeared densities, not the original densities.

        Returns:
            Dict of dos data. Generally of the form, {label: {'frequencies':..,
            'densities': ...}}
        """
        return jsanitize(self._doses)

    def get_plot(self, xlim=None, ylim=None, units="thz"):
        """
        Get a matplotlib plot showing the DOS.

        Args:
            xlim: Specifies the x-axis limits. Set to None for automatic
                determination.
            ylim: Specifies the y-axis limits.
            units: units for the frequencies. Accepted values thz, ev, mev, ha, cm-1, cm^-1.
        """
        u = freq_units(units)

        ncolors = max(3, len(self._doses))
        ncolors = min(9, ncolors)

        import palettable

        colors = palettable.colorbrewer.qualitative.Set1_9.mpl_colors  # pylint: disable=E1101

        y = None
        alldensities = []
        allfrequencies = []
        plt = pretty_plot(12, 8)

        # Note that this complicated processing of frequencies is to allow for
        # stacked plots in matplotlib.
        for dos in self._doses.values():
            frequencies = dos["frequencies"] * u.factor
            densities = dos["densities"]
            if y is None:
                y = np.zeros(frequencies.shape)
            if self.stack:
                y += densities
                newdens = y.copy()
            else:
                newdens = densities
            allfrequencies.append(frequencies)
            alldensities.append(newdens)

        keys = list(self._doses)
        keys.reverse()
        alldensities.reverse()
        allfrequencies.reverse()
        allpts = []
        for i, (key, frequencies, densities) in enumerate(zip(keys, allfrequencies, alldensities)):
            allpts.extend(list(zip(frequencies, densities)))
            if self.stack:
                plt.fill(frequencies, densities, color=colors[i % ncolors], label=str(key))
            else:
                plt.plot(
                    frequencies,
                    densities,
                    color=colors[i % ncolors],
                    label=str(key),
                    linewidth=3,
                )

        if xlim:
            plt.xlim(xlim)
        if ylim:
            plt.ylim(ylim)
        else:
            xlim = plt.xlim()
            relevanty = [p[1] for p in allpts if xlim[0] < p[0] < xlim[1]]
            plt.ylim((min(relevanty), max(relevanty)))

        ylim = plt.ylim()
        plt.plot([0, 0], ylim, "k--", linewidth=2)

        plt.xlabel(rf"$\mathrm{{Frequencies\ ({u.label})}}$")
        plt.ylabel(r"$\mathrm{Density\ of\ states}$")

        plt.legend()
        leg = plt.gca().get_legend()
        ltext = leg.get_texts()  # all the text.Text instance in the legend
        plt.setp(ltext, fontsize=30)
        plt.tight_layout()
        return plt

    def save_plot(self, filename, img_format="eps", xlim=None, ylim=None, units="thz"):
        """
        Save matplotlib plot to a file.

        Args:
            filename: Filename to write to.
            img_format: Image format to use. Defaults to EPS.
            xlim: Specifies the x-axis limits. Set to None for automatic
                determination.
            ylim: Specifies the y-axis limits.
            units: units for the frequencies. Accepted values thz, ev, mev, ha, cm-1, cm^-1
        """
        plt = self.get_plot(xlim, ylim, units=units)
        plt.savefig(filename, format=img_format)
        plt.close()

    def show(self, xlim=None, ylim=None, units="thz"):
        """
        Show the plot using matplotlib.

        Args:
            xlim: Specifies the x-axis limits. Set to None for automatic
                determination.
            ylim: Specifies the y-axis limits.
            units: units for the frequencies. Accepted values thz, ev, mev, ha, cm-1, cm^-1.
        """
        plt = self.get_plot(xlim, ylim, units=units)
        plt.show()


class PhononBSPlotter:
    """Class to plot or get data to facilitate the plot of band structure objects."""

    def __init__(self, bs):
        """
        Args:
            bs: A PhononBandStructureSymmLine object.
        """
        if not isinstance(bs, PhononBandStructureSymmLine):
            raise ValueError(
                "PhononBSPlotter only works with PhononBandStructureSymmLine objects. "
                "A PhononBandStructure object (on a uniform grid for instance and "
                "not along symmetry lines won't work)"
            )
        self._bs = bs
        self._nb_bands = self._bs.nb_bands

    def _maketicks(self, plt):
        """Utility private method to add ticks to a band structure."""
        ticks = self.get_ticks()
        # Sanitize only plot the uniq values
        uniq_d = []
        uniq_l = []
        temp_ticks = list(zip(ticks["distance"], ticks["label"]))
        for i, tt in enumerate(temp_ticks):
            if i == 0:
                uniq_d.append(tt[0])
                uniq_l.append(tt[1])
                logger.debug(f"Adding label {tt[0]} at {tt[1]}")
            elif tt[1] == temp_ticks[i - 1][1]:
                logger.debug(f"Skipping label {tt[1]}")
            else:
                logger.debug(f"Adding label {tt[0]} at {tt[1]}")
                uniq_d.append(tt[0])
                uniq_l.append(tt[1])

        logger.debug(f"Unique labels are {list(zip(uniq_d, uniq_l))}")
        plt.gca().set_xticks(uniq_d)
        plt.gca().set_xticklabels(uniq_l)

        for i in range(len(ticks["label"])):
            if ticks["label"][i] is not None:
                # don't print the same label twice
                if i != 0:
                    if ticks["label"][i] == ticks["label"][i - 1]:
                        logger.debug(f"already print label... skipping label {ticks['label'][i]}")
                    else:
                        logger.debug(f"Adding a line at {ticks['distance'][i]} for label {ticks['label'][i]}")
                        plt.axvline(ticks["distance"][i], color="k")
                else:
                    logger.debug(f"Adding a line at {ticks['distance'][i]} for label {ticks['label'][i]}")
                    plt.axvline(ticks["distance"][i], color="k")
        return plt

    def bs_plot_data(self):
        """
        Get the data nicely formatted for a plot.

        Returns:
            A dict of the following format:
            ticks: A dict with the 'distances' at which there is a qpoint (the
            x axis) and the labels (None if no label)
            frequencies: A list (one element for each branch) of frequencies for
            each qpoint: [branch][qpoint][mode]. The data is
            stored by branch to facilitate the plotting
            lattice: The reciprocal lattice.
        """
        distance = []
        frequency = []

        ticks = self.get_ticks()

        for b in self._bs.branches:
            frequency.append([])
            distance.append([self._bs.distance[j] for j in range(b["start_index"], b["end_index"] + 1)])

            for i in range(self._nb_bands):
                frequency[-1].append([self._bs.bands[i][j] for j in range(b["start_index"], b["end_index"] + 1)])

        return {
            "ticks": ticks,
            "distances": distance,
            "frequency": frequency,
            "lattice": self._bs.lattice_rec.as_dict(),
        }

    def get_plot(self, ylim=None, units="thz"):
        """
        Get a matplotlib object for the bandstructure plot.

        Args:
            ylim: Specify the y-axis (frequency) limits; by default None let
                the code choose.
            units: units for the frequencies. Accepted values thz, ev, mev, ha, cm-1, cm^-1.
        """
        u = freq_units(units)

        plt = pretty_plot(12, 8)

        band_linewidth = 1

        data = self.bs_plot_data()
        for d in range(len(data["distances"])):
            for i in range(self._nb_bands):
                plt.plot(
                    data["distances"][d],
                    [data["frequency"][d][i][j] * u.factor for j in range(len(data["distances"][d]))],
                    "b-",
                    linewidth=band_linewidth,
                )

        self._maketicks(plt)

        # plot y=0 line
        plt.axhline(0, linewidth=1, color="k")

        # Main X and Y Labels
        plt.xlabel(r"$\mathrm{Wave\ Vector}$", fontsize=30)
        ylabel = rf"$\mathrm{{Frequencies\ ({u.label})}}$"
        plt.ylabel(ylabel, fontsize=30)

        # X range (K)
        # last distance point
        x_max = data["distances"][-1][-1]
        plt.xlim(0, x_max)

        if ylim is not None:
            plt.ylim(ylim)

        plt.tight_layout()

        return plt

    def _get_weight(self, vec: np.ndarray, indices: list[list[int]]) -> np.ndarray:
        """
        compute the weight for each combination of sites according to the
        eigenvector.
        """
        num_atom = int(self._nb_bands / 3)
        new_vec = np.zeros(num_atom)
        for i in range(num_atom):
            new_vec[i] = np.linalg.norm(vec[i * 3 : i * 3 + 3])
        # get the projectors for each group
        gw = []
        norm_f = 0
        for comb in indices:
            projector = np.zeros(len(new_vec))
            for idx in range(len(projector)):
                if idx in comb:
                    projector[idx] = 1
            group_weight = np.dot(projector, new_vec)
            gw.append(group_weight)
            norm_f += group_weight
        return np.array(gw, dtype=float) / norm_f

    @staticmethod
    def _make_color(colors: list[int]) -> list[int]:
        """Convert the eigendisplacements to rgb colors."""
        # if there are two groups, use red and blue
        if len(colors) == 2:
            return [colors[0], 0, colors[1]]
        if len(colors) == 3:
            return colors
        # if there are four groups, use cyan, magenta, yellow and black
        if len(colors) == 4:
            r = (1 - colors[0]) * (1 - colors[3])
            g = (1 - colors[1]) * (1 - colors[3])
            b = (1 - colors[2]) * (1 - colors[3])
            return [r, g, b]
        raise ValueError(f"Expected 2, 3 or 4 colors, got {len(colors)}")

    def get_proj_plot(
        self,
        site_comb: str | list[list[int]] = "element",
        ylim: tuple[None | float, None | float] | None = None,
        units: str = "thz",
        rgb_labels: tuple[None | str] | None = None,
    ) -> plt.Axes:
        """
        Get a matplotlib object for the bandstructure plot projected along atomic
        sites.

        Args:
            site_comb: a list of list, for example, [[0],[1],[2,3,4]];
                the numbers in each sublist represents the indices of atoms;
                the atoms in a same sublist will be plotted in a same color;
                if not specified, unique elements are automatically grouped.
            ylim: Specify the y-axis (frequency) limits; by default None let
                the code choose.
            units: units for the frequencies. Accepted values thz, ev, mev, ha, cm-1, cm^-1.
                Defaults to "thz".
            rgb_labels: a list of rgb colors for the labels; if not specified,
                the colors will be automatically generated.
        """
        from matplotlib.collections import LineCollection

        from pymatgen.electronic_structure.plotter import BSDOSPlotter

        elements = [e.symbol for e in self._bs.structure.elements]
        if site_comb == "element":
            assert 2 <= len(elements) <= 4, "the compound must have 2, 3 or 4 unique elements"
            indices: list[list[int]] = [[] for _ in range(len(elements))]
            for i, ele in enumerate(self._bs.structure.species):
                for j, unique_species in enumerate(self._bs.structure.elements):
                    if ele == unique_species:
                        indices[j].append(i)
        else:
            assert isinstance(site_comb, list)
            assert 2 <= len(site_comb) <= 4, "the length of site_comb must be 2, 3 or 4"
            all_sites = self._bs.structure.sites
            all_indices = {*range(len(all_sites))}
            for comb in site_comb:
                for idx in comb:
                    assert 0 <= idx < len(all_sites), "one or more indices in site_comb does not exist"
                    all_indices.remove(idx)
            if len(all_indices) != 0:
                raise Exception(f"not all {len(all_sites)} indices are included in site_comb")
            indices = site_comb  # type: ignore[assignment]
        assert rgb_labels is None or len(rgb_labels) == len(indices), "wrong number of rgb_labels"

        u = freq_units(units)
        fig, ax = plt.subplots(figsize=(12, 8), dpi=300)
        self._maketicks(plt)

        data = self.bs_plot_data()
        k_dist = np.array(data["distances"]).flatten()
        for d in range(1, len(k_dist)):
            # consider 2 k points each time so they connect
            colors = []
            for idx in range(self._nb_bands):
                eigenvec_1 = self._bs.eigendisplacements[idx][d - 1].flatten()
                eigenvec_2 = self._bs.eigendisplacements[idx][d].flatten()
                colors1 = self._get_weight(eigenvec_1, indices)
                colors2 = self._get_weight(eigenvec_2, indices)
                colors.append(self._make_color((colors1 + colors2) / 2))
            seg = np.zeros((self._nb_bands, 2, 2))
            seg[:, :, 1] = self._bs.bands[:, d - 1 : d + 1] * u.factor
            seg[:, 0, 0] = k_dist[d - 1]
            seg[:, 1, 0] = k_dist[d]
            ls = LineCollection(seg, colors=colors, linestyles="-", linewidths=2.5)
            ax.add_collection(ls)
        if ylim is None:
            y_max: float = max(max(b) for b in self._bs.bands) * u.factor
            y_min: float = min(min(b) for b in self._bs.bands) * u.factor
            y_margin = (y_max - y_min) * 0.05
            ylim = (y_min - y_margin, y_max + y_margin)
        ax.set_ylim(ylim)
        xlim = [min(k_dist), max(k_dist)]
        ax.set_xlim(xlim)
        ax.set_xlabel(r"$\mathrm{Wave\ Vector}$", fontsize=28)
        ylabel = rf"$\mathrm{{Frequencies\ ({u.label})}}$"
        ax.set_ylabel(ylabel, fontsize=28)
        ax.tick_params(labelsize=28)
        # make color legend
        labels: list[str]
        if rgb_labels is not None:
            labels = rgb_labels  # type: ignore[assignment]
        elif site_comb == "element":
            labels = [e.symbol for e in self._bs.structure.elements]
        else:
            labels = [f"{i}" for i in range(len(site_comb))]
        if len(indices) == 2:
            BSDOSPlotter._rb_line(ax, labels[0], labels[1], "best")
        elif len(indices) == 3:
            BSDOSPlotter._rgb_triangle(ax, labels[0], labels[1], labels[2], "best")
        else:
            # for 4 combinations, build a color square?
            pass
        return ax

    def show(self, ylim=None, units="thz"):
        """
        Show the plot using matplotlib.

        Args:
            ylim: Specify the y-axis (frequency) limits; by default None let
                the code choose.
            units: units for the frequencies. Accepted values thz, ev, mev, ha, cm-1, cm^-1.
        """
        plt = self.get_plot(ylim, units=units)
        plt.show()

    def save_plot(self, filename, img_format="eps", ylim=None, units="thz"):
        """
        Save matplotlib plot to a file.

        Args:
            filename: Filename to write to.
            img_format: Image format to use. Defaults to EPS.
            ylim: Specifies the y-axis limits.
            units: units for the frequencies. Accepted values thz, ev, mev, ha, cm-1, cm^-1.
        """
        plt = self.get_plot(ylim=ylim, units=units)
        plt.savefig(filename, format=img_format)
        plt.close()

    def show_proj(
        self,
        site_comb: str | list[list[int]] = "element",
        ylim: tuple[None | float, None | float] | None = None,
        units: str = "thz",
        rgb_labels: tuple[str] | None = None,
    ):
        """
        Show the projected plot using matplotlib.

        Args:
            site_comb: A list of list of indices of sites to combine. For example,
                [[0, 1], [2, 3]] will combine the projections of sites 0 and 1,
                and sites 2 and 3. Defaults to "element", which will combine
                sites by element.
            ylim: Specify the y-axis (frequency) limits; by default None let
                the code choose.
            units: units for the frequencies. Accepted values thz, ev, mev, ha, cm-1, cm^-1.
                Defaults to "thz".
            rgb_labels: A list of labels for the rgb triangle. Defaults to None,
                which will use the element symbols.
        """
        self.get_proj_plot(site_comb=site_comb, ylim=ylim, units=units, rgb_labels=rgb_labels)
        plt.show()

    def get_ticks(self):
        """
        Get all ticks and labels for a band structure plot.

        Returns:
            A dict with 'distance': a list of distance at which ticks should
            be set and 'label': a list of label for each of those ticks.
        """
        tick_distance = []
        tick_labels = []
        previous_label = self._bs.qpoints[0].label
        previous_branch = self._bs.branches[0]["name"]
        for i, c in enumerate(self._bs.qpoints):
            if c.label is not None:
                tick_distance.append(self._bs.distance[i])
                this_branch = None
                for b in self._bs.branches:
                    if b["start_index"] <= i <= b["end_index"]:
                        this_branch = b["name"]
                        break
                if c.label != previous_label and previous_branch != this_branch:
                    label1 = c.label
                    if label1.startswith("\\") or label1.find("_") != -1:
                        label1 = "$" + label1 + "$"
                    label0 = previous_label
                    if label0.startswith("\\") or label0.find("_") != -1:
                        label0 = "$" + label0 + "$"
                    tick_labels.pop()
                    tick_distance.pop()
                    tick_labels.append(label0 + "$\\mid$" + label1)
                elif c.label.startswith("\\") or c.label.find("_") != -1:
                    tick_labels.append("$" + c.label + "$")
                else:
                    tick_labels.append(c.label)
                previous_label = c.label
                previous_branch = this_branch
        return {"distance": tick_distance, "label": tick_labels}

    def plot_compare(self, other_plotter, units="thz"):
        """
        plot two band structure for comparison. One is in red the other in blue.
        The two band structures need to be defined on the same symmetry lines!
        and the distance between symmetry lines is the one of the band structure
        used to build the PhononBSPlotter.

        Args:
            other_plotter: another PhononBSPlotter object defined along the same symmetry lines
            units: units for the frequencies. Accepted values thz, ev, mev, ha, cm-1, cm^-1.
                Defaults to 'thz'.

        Returns:
            a matplotlib object with both band structures
        """
        u = freq_units(units)

        data_orig = self.bs_plot_data()
        data = other_plotter.bs_plot_data()

        if len(data_orig["distances"]) != len(data["distances"]):
            raise ValueError("The two objects are not compatible.")

        plt = self.get_plot(units=units)
        band_linewidth = 1
        for i in range(other_plotter._nb_bands):
            for d in range(len(data_orig["distances"])):
                plt.plot(
                    data_orig["distances"][d],
                    [data["frequency"][d][i][j] * u.factor for j in range(len(data_orig["distances"][d]))],
                    "r-",
                    linewidth=band_linewidth,
                )

        return plt

    def plot_brillouin(self):
        """Plot the Brillouin zone."""
        # get labels and lines
        labels = {}
        for q in self._bs.qpoints:
            if q.label:
                labels[q.label] = q.frac_coords

        lines = []
        for b in self._bs.branches:
            lines.append(
                [
                    self._bs.qpoints[b["start_index"]].frac_coords,
                    self._bs.qpoints[b["end_index"]].frac_coords,
                ]
            )

        plot_brillouin_zone(self._bs.lattice_rec, lines=lines, labels=labels)


class ThermoPlotter:
    """
    Plotter for thermodynamic properties obtained from phonon DOS.
    If the structure corresponding to the DOS, it will be used to extract the formula unit and provide
    the plots in units of mol instead of mole-cell.
    """

    def __init__(self, dos, structure=None):
        """
        Args:
            dos: A PhononDos object.
            structure: A Structure object corresponding to the structure used for the calculation.
        """
        self.dos = dos
        self.structure = structure

    def _plot_thermo(self, func, temperatures, factor=1, ax=None, ylabel=None, label=None, ylim=None, **kwargs):
        """
        Plots a thermodynamic property for a generic function from a PhononDos instance.

        Args:
            func: the thermodynamic function to be used to calculate the property
            temperatures: a list of temperatures
            factor: a multiplicative factor applied to the thermodynamic property calculated. Used to change
                the units.
            ax: matplotlib :class:`Axes` or None if a new figure should be created.
            ylabel: label for the y axis
            label: label of the plot
            ylim: tuple specifying the y-axis limits.
            kwargs: kwargs passed to the matplotlib function 'plot'.

        Returns:
            matplotlib figure
        """
        ax, fig, plt = get_ax_fig_plt(ax)

        values = []

        for t in temperatures:
            values.append(func(t, structure=self.structure) * factor)

        ax.plot(temperatures, values, label=label, **kwargs)

        if ylim:
            ax.set_ylim(ylim)

        ax.set_xlim((np.min(temperatures), np.max(temperatures)))
        ylim = plt.ylim()
        if ylim[0] < 0 < ylim[1]:
            plt.plot(plt.xlim(), [0, 0], "k-", linewidth=1)

        ax.set_xlabel(r"$T$ (K)")
        if ylabel:
            ax.set_ylabel(ylabel)

        return fig

    @add_fig_kwargs
    def plot_cv(self, tmin, tmax, ntemp, ylim=None, **kwargs):
        """
        Plots the constant volume specific heat C_v in a temperature range.

        Args:
            tmin: minimum temperature
            tmax: maximum temperature
            ntemp: number of steps
            ylim: tuple specifying the y-axis limits.
            kwargs: kwargs passed to the matplotlib function 'plot'.

        Returns:
            matplotlib figure
        """
        temperatures = np.linspace(tmin, tmax, ntemp)

        ylabel = "$C_v$ (J/K/mol)" if self.structure else "$C_v$ (J/K/mol-c)"

        return self._plot_thermo(self.dos.cv, temperatures, ylabel=ylabel, ylim=ylim, **kwargs)

    @add_fig_kwargs
    def plot_entropy(self, tmin, tmax, ntemp, ylim=None, **kwargs):
        """
        Plots the vibrational entrpy in a temperature range.

        Args:
            tmin: minimum temperature
            tmax: maximum temperature
            ntemp: number of steps
            ylim: tuple specifying the y-axis limits.
            kwargs: kwargs passed to the matplotlib function 'plot'.

        Returns:
            matplotlib figure
        """
        temperatures = np.linspace(tmin, tmax, ntemp)

        ylabel = "$S$ (J/K/mol)" if self.structure else "$S$ (J/K/mol-c)"

        return self._plot_thermo(self.dos.entropy, temperatures, ylabel=ylabel, ylim=ylim, **kwargs)

    @add_fig_kwargs
    def plot_internal_energy(self, tmin, tmax, ntemp, ylim=None, **kwargs):
        """
        Plots the vibrational internal energy in a temperature range.

        Args:
            tmin: minimum temperature
            tmax: maximum temperature
            ntemp: number of steps
            ylim: tuple specifying the y-axis limits.
            kwargs: kwargs passed to the matplotlib function 'plot'.

        Returns:
            matplotlib figure
        """
        temperatures = np.linspace(tmin, tmax, ntemp)

        ylabel = "$\\Delta E$ (kJ/mol)" if self.structure else "$\\Delta E$ (kJ/mol-c)"

        return self._plot_thermo(
            self.dos.internal_energy, temperatures, ylabel=ylabel, ylim=ylim, factor=1e-3, **kwargs
        )

    @add_fig_kwargs
    def plot_helmholtz_free_energy(self, tmin, tmax, ntemp, ylim=None, **kwargs):
        """
        Plots the vibrational contribution to the Helmoltz free energy in a temperature range.

        Args:
            tmin: minimum temperature
            tmax: maximum temperature
            ntemp: number of steps
            ylim: tuple specifying the y-axis limits.
            kwargs: kwargs passed to the matplotlib function 'plot'.

        Returns:
            matplotlib figure
        """
        temperatures = np.linspace(tmin, tmax, ntemp)

        ylabel = "$\\Delta F$ (kJ/mol)" if self.structure else "$\\Delta F$ (kJ/mol-c)"

        return self._plot_thermo(
            self.dos.helmholtz_free_energy, temperatures, ylabel=ylabel, ylim=ylim, factor=1e-3, **kwargs
        )

    @add_fig_kwargs
    def plot_thermodynamic_properties(self, tmin, tmax, ntemp, ylim=None, **kwargs):
        """
        Plots all the thermodynamic properties in a temperature range.

        Args:
            tmin: minimum temperature
            tmax: maximum temperature
            ntemp: number of steps
            ylim: tuple specifying the y-axis limits.
            kwargs: kwargs passed to the matplotlib function 'plot'.

        Returns:
            matplotlib figure
        """
        temperatures = np.linspace(tmin, tmax, ntemp)

        mol = "" if self.structure else "-c"

        fig = self._plot_thermo(
            self.dos.cv,
            temperatures,
            ylabel="Thermodynamic properties",
            ylim=ylim,
            label=rf"$C_v$ (J/K/mol{mol})",
            **kwargs,
        )
        self._plot_thermo(
            self.dos.entropy, temperatures, ylim=ylim, ax=fig.axes[0], label=rf"$S$ (J/K/mol{mol})", **kwargs
        )
        self._plot_thermo(
            self.dos.internal_energy,
            temperatures,
            ylim=ylim,
            ax=fig.axes[0],
            factor=1e-3,
            label=rf"$\Delta E$ (kJ/mol{mol})",
            **kwargs,
        )
        self._plot_thermo(
            self.dos.helmholtz_free_energy,
            temperatures,
            ylim=ylim,
            ax=fig.axes[0],
            factor=1e-3,
            label=rf"$\Delta F$ (kJ/mol{mol})",
            **kwargs,
        )

        fig.axes[0].legend(loc="best")

        return fig


class GruneisenPlotter:
    """Class to plot Gruneisenparameter Object."""

    def __init__(self, gruneisen):
        """
        Class to plot information from Gruneisenparameter Object
        Args:
            gruneisen: GruneisenParameter Object.
        """
        self._gruneisen = gruneisen

    def get_plot(self, marker="o", markersize=None, units="thz"):
        """
        will produce a plot
        Args:
            marker: marker for the depiction
            markersize: size of the marker
            units: unit for the plots, accepted units: thz, ev, mev, ha, cm-1, cm^-1.

        Returns: plot
        """
        u = freq_units(units)

        xs = self._gruneisen.frequencies.flatten() * u.factor
        ys = self._gruneisen.gruneisen.flatten()

        plt = pretty_plot(12, 8)

        plt.xlabel(rf"$\mathrm{{Frequency\ ({u.label})}}$")
        plt.ylabel(r"$\mathrm{Grüneisen\ parameter}$")

        n = len(ys) - 1
        for i, (x, y) in enumerate(zip(xs, ys)):
            color = (1.0 / n * i, 0, 1.0 / n * (n - i))

            if markersize:
                plt.plot(x, y, marker, color=color, markersize=markersize)
            else:
                plt.plot(x, y, marker, color=color)

        plt.tight_layout()

        return plt

    def show(self, units="thz"):
        """
        will show the plot
        Args:
            units: units for the plot, accepted units: thz, ev, mev, ha, cm-1, cm^-1.

        Returns: plot
        """
        plt = self.get_plot(units=units)
        plt.show()

    def save_plot(self, filename, img_format="pdf", units="thz"):
        """
        Will save the plot to a file
        Args:
            filename: name of the filename
            img_format: format of the saved plot
            units: accepted units: thz, ev, mev, ha, cm-1, cm^-1.
        """
        plt = self.get_plot(units=units)
        plt.savefig(filename, format=img_format)
        plt.close()


class GruneisenPhononBSPlotter(PhononBSPlotter):
    """Class to plot or get data to facilitate the plot of band structure objects."""

    def __init__(self, bs):
        """
        Args:
            bs: A GruneisenPhononBandStructureSymmLine object.
        """
        if not isinstance(bs, GruneisenPhononBandStructureSymmLine):
            raise ValueError(
                "GruneisenPhononBSPlotter only works with GruneisenPhononBandStructureSymmLine objects. "
                "A GruneisenPhononBandStructure object (on a uniform grid for instance and "
                "not along symmetry lines won't work)"
            )
        super().__init__(bs)

    def bs_plot_data(self):
        """
        Get the data nicely formatted for a plot.

        Returns:
            A dict of the following format:
            ticks: A dict with the 'distances' at which there is a qpoint (the
            x axis) and the labels (None if no label)
            frequencies: A list (one element for each branch) of frequencies for
            each qpoint: [branch][qpoint][mode]. The data is
            stored by branch to facilitate the plotting
            gruneisen: GruneisenPhononBandStructureSymmLine
            lattice: The reciprocal lattice.
        """
        distance, frequency, gruneisen = ([] for _ in range(3))

        ticks = self.get_ticks()

        for b in self._bs.branches:
            frequency.append([])
            gruneisen.append([])
            distance.append([self._bs.distance[j] for j in range(b["start_index"], b["end_index"] + 1)])

            for i in range(self._nb_bands):
                frequency[-1].append([self._bs.bands[i][j] for j in range(b["start_index"], b["end_index"] + 1)])
                gruneisen[-1].append([self._bs.gruneisen[i][j] for j in range(b["start_index"], b["end_index"] + 1)])

        return {
            "ticks": ticks,
            "distances": distance,
            "frequency": frequency,
            "gruneisen": gruneisen,
            "lattice": self._bs.lattice_rec.as_dict(),
        }

    def get_plot_gs(self, ylim=None):
        """
        Get a matplotlib object for the gruneisen bandstructure plot.

        Args:
            ylim: Specify the y-axis (gruneisen) limits; by default None let
                the code choose.
        """
        plt = pretty_plot(12, 8)

        # band_linewidth = 1

        data = self.bs_plot_data()
        for d in range(len(data["distances"])):
            for i in range(self._nb_bands):
                plt.plot(
                    data["distances"][d],
                    [data["gruneisen"][d][i][j] for j in range(len(data["distances"][d]))],
                    "b-",
                    # linewidth=band_linewidth)
                    marker="o",
                    markersize=2,
                    linewidth=2,
                )

        self._maketicks(plt)

        # plot y=0 line
        plt.axhline(0, linewidth=1, color="k")

        # Main X and Y Labels
        plt.xlabel(r"$\mathrm{Wave\ Vector}$", fontsize=30)
        plt.ylabel(r"$\mathrm{Grüneisen\ Parameter}$", fontsize=30)

        # X range (K)
        # last distance point
        x_max = data["distances"][-1][-1]
        plt.xlim(0, x_max)

        if ylim is not None:
            plt.ylim(ylim)

        plt.tight_layout()

        return plt

    def show_gs(self, ylim=None):
        """
        Show the plot using matplotlib.

        Args:
            ylim: Specifies the y-axis limits.
        """
        plt = self.get_plot_gs(ylim)
        plt.show()

    def save_plot_gs(self, filename, img_format="eps", ylim=None):
        """
        Save matplotlib plot to a file.

        Args:
            filename: Filename to write to.
            img_format: Image format to use. Defaults to EPS.
            ylim: Specifies the y-axis limits.
        """
        plt = self.get_plot_gs(ylim=ylim)
        plt.savefig(filename, format=img_format)
        plt.close()

    def plot_compare_gs(self, other_plotter: GruneisenPhononBSPlotter) -> plt:
        """
        plot two band structure for comparison. One is in red the other in blue.
        The two band structures need to be defined on the same symmetry lines!
        and the distance between symmetry lines is
        the one of the band structure used to build the PhononBSPlotter.

        Args:
            other_plotter (GruneisenPhononBSPlotter): another phonon DOS plotter defined along
                the same symmetry lines.

        Returns:
            a matplotlib object with both band structures
        """
        data_orig = self.bs_plot_data()
        data = other_plotter.bs_plot_data()

        if len(data_orig["distances"]) != len(data["distances"]):
            raise ValueError("The two objects are not compatible.")

        plt = self.get_plot()
        band_linewidth = 1
        for i in range(other_plotter._nb_bands):
            for d in range(len(data_orig["distances"])):
                plt.plot(
                    data_orig["distances"][d],
                    [data["gruneisen"][d][i][j] for j in range(len(data_orig["distances"][d]))],
                    "r-",
                    linewidth=band_linewidth,
                )

        return plt
