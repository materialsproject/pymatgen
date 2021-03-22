# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


"""
This module provides plotting capabilities for battery related applications.
"""


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Jul 12, 2012"


from collections import OrderedDict

import plotly.graph_objects as go
from pymatgen.util.plotting import pretty_plot


class VoltageProfilePlotter:
    """
    A plotter to make voltage profile plots for batteries.
    """

    def __init__(self, xaxis="capacity", hide_negative=False):
        """
        Args:
            xaxis: The quantity to use as the xaxis. Can be either
            - capacity_grav: the graviometric capcity
            - capacity_vol: the volumetric capacity
            - x_form: the number of working ions per formula unit of the host
            - frac_x: the atomic fraction of the working ion
            hide_negative: If True only plot the voltage steps above zero
        """
        self._electrodes = OrderedDict()
        self.xaxis = xaxis
        self.hide_negative = hide_negative

    def add_electrode(self, electrode, label=None):
        """
        Add an electrode to the plot.

        Args:
            electrode: An electrode. All electrodes satisfying the
                AbstractElectrode interface should work.
            label: A label for the electrode. If None, defaults to a counting
                system, i.e. 'Electrode 1', 'Electrode 2', ...
        """
        if not label:
            label = "Electrode {}".format(len(self._electrodes) + 1)
        self._electrodes[label] = electrode

    def get_plot_data(self, electrode, term_zero=True):
        """
        Args:
            electrode: Electrode object
            term_zero: If True append zero voltage point at the end

        Returns:
            Plot data in x, y.
        """
        x = []
        y = []
        cap = 0

        for sub_electrode in electrode.get_sub_electrodes(adjacent_only=True):
            if self.hide_negative and sub_electrode.get_average_voltage() < 0:
                continue
            if self.xaxis in {"capacity_grav", "capacity"}:
                x.append(cap)
                cap += sub_electrode.get_capacity_grav()
                x.append(cap)
            elif self.xaxis == "capacity_vol":
                x.append(cap)
                cap += sub_electrode.get_capacity_vol()
                x.append(cap)
            elif self.xaxis == "x_form":
                x.append(sub_electrode.x_charge)
                x.append(sub_electrode.x_discharge)
            elif self.xaxis == "frac_x":
                x.append(sub_electrode.voltage_pairs[0].frac_charge)
                x.append(sub_electrode.voltage_pairs[0].frac_discharge)
            else:
                raise NotImplementedError("x_axis must be capacity_grav/capacity_vol/x_form/frac_x")
            y.extend([sub_electrode.get_average_voltage()] * 2)

        if term_zero:
            x.append(x[-1])
            y.append(0)
        return x, y

    def get_plot(self, width=8, height=8, term_zero=True):
        """
        Returns a plot object.

        Args:
            width: Width of the plot. Defaults to 8 in.
            height: Height of the plot. Defaults to 6 in.
            term_zero: If True append zero voltage point at the end

        Returns:
            A matplotlib plot object.
        """
        plt = pretty_plot(width, height)
        wion_symbol = set()
        formula = set()

        for label, electrode in self._electrodes.items():
            (x, y) = self.get_plot_data(electrode, term_zero=term_zero)
            wion_symbol.add(electrode.working_ion.symbol)
            formula.add(electrode._framework_formula)
            plt.plot(x, y, "-", linewidth=2, label=label)

        plt.legend()
        plt.xlabel(self._choose_best_x_lable(formula=formula, wion_symbol=wion_symbol))
        plt.ylabel("Voltage (V)")
        plt.tight_layout()
        return plt

    def get_plotly_figure(
        self,
        width=800,
        height=600,
        font_dict=None,
        term_zero=True,
        **kwargs,
    ):
        """
        Return plotly Figure object
        Args:
            width: Width of the plot. Defaults to 800 px.
            height: Height of the plot. Defaults to 600 px.
            font: dictionary that defines the font
            term_zero: If True append zero voltage point at the end
            **kwargs:

        Returns:

        """
        font_dict = dict(family="Arial", size=24, color="#000000") if font_dict is None else font_dict
        hover_temp = "Voltage : %{y:.2f} V"

        data = []
        wion_symbol = set()
        formula = set()
        for label, electrode in self._electrodes.items():
            (x, y) = self.get_plot_data(electrode, term_zero=term_zero)
            wion_symbol.add(electrode.working_ion.symbol)
            formula.add(electrode._framework_formula)
            data.append(go.Scatter(x=x, y=y, name=label, hovertemplate=hover_temp))

        fig = go.Figure(
            data=data,
            layout=go.Layout(
                title="Voltage vs. Capacity",
                width=width,
                height=height,
                font=font_dict,
                xaxis=dict(title=self._choose_best_x_lable(formula=formula, wion_symbol=wion_symbol)),
                yaxis=dict(title="Voltage (V)"),
                **kwargs,
            ),
        )

        fig.update_layout(template="plotly_white", title_x=0.5)
        return fig

    def _choose_best_x_lable(self, formula, wion_symbol):
        if self.xaxis in {"capacity", "capacity_grav"}:
            return "Capacity (mAh/g)"
        if self.xaxis == "capacity_vol":
            return "Capacity (Ah/l)"

        if len(formula) == 1:
            formula = formula.pop()
        else:
            formula = None

        if len(wion_symbol) == 1:
            wion_symbol = wion_symbol.pop()
        else:
            wion_symbol = None

        if self.xaxis == "x_form":
            if formula and wion_symbol:
                return f"x in {wion_symbol}<sub>x</sub>{formula}"
            return "x Workion Ion per Host F.U."

        if self.xaxis == "frac_x":
            if wion_symbol:
                return f"Atomic Fraction of {wion_symbol}"
            return "Atomic Fraction of Working Ion"
        raise RuntimeError("No xaxis label can be determined")

    def show(self, width=8, height=6):
        """
        Show the voltage profile plot.

        Args:
            width: Width of the plot. Defaults to 8 in.
            height: Height of the plot. Defaults to 6 in.
        """
        self.get_plot(width, height).show()

    def save(self, filename, image_format="eps", width=8, height=6):
        """
        Save the plot to an image file.

        Args:
            filename: Filename to save to.
            image_format: Format to save to. Defaults to eps.
        """
        self.get_plot(width, height).savefig(filename, format=image_format)
