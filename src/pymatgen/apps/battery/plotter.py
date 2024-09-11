"""This module provides plotting capabilities for battery related applications."""

from __future__ import annotations

from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import plotly.graph_objects as go

from pymatgen.util.plotting import pretty_plot

if TYPE_CHECKING:
    from pymatgen.apps.battery.battery_abc import AbstractElectrode

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Jul 12, 2012"


class VoltageProfilePlotter:
    """A plotter to make voltage profile plots for batteries."""

    def __init__(self, xaxis: str = "capacity", hide_negative: bool = False) -> None:
        """
        Args:
            xaxis: The quantity to use as the xaxis. Can be either
            - capacity_grav: the gravimetric capacity
            - capacity_vol: the volumetric capacity
            - x_form: the number of working ions per formula unit of the host
            - frac_x: the atomic fraction of the working ion
            hide_negative: If True only plot the voltage steps above zero.
        """
        self._electrodes: dict[str, AbstractElectrode] = {}
        self.xaxis = xaxis
        self.hide_negative = hide_negative

    def add_electrode(self, electrode: AbstractElectrode, label: str | None = None) -> None:
        """Add an electrode to the plot.

        Args:
            electrode: An electrode. All electrodes satisfying the
                AbstractElectrode interface should work.
            label: A label for the electrode. If None, defaults to a counting
                system, i.e. 'Electrode 1', 'Electrode 2', ...
        """
        if label is None:
            label = f"Electrode {len(self._electrodes) + 1}"
        self._electrodes[label] = electrode

    def get_plot_data(self, electrode: AbstractElectrode, term_zero: bool = True) -> tuple[list, list]:
        """
        Args:
            electrode: Electrode object
            term_zero: If True append zero voltage point at the end.

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
                x.extend((sub_electrode.x_charge, sub_electrode.x_discharge))
            elif self.xaxis == "frac_x":
                x.extend((sub_electrode.voltage_pairs[0].frac_charge, sub_electrode.voltage_pairs[0].frac_discharge))
            else:
                raise NotImplementedError("x_axis must be capacity_grav/capacity_vol/x_form/frac_x")
            y.extend([sub_electrode.get_average_voltage()] * 2)

        if term_zero:
            x.append(x[-1])
            y.append(0)
        return x, y

    def get_plot(self, width: float = 8, height: float = 8, term_zero: bool = True, ax: plt.Axes = None) -> plt.Axes:
        """Get a plot object.

        Args:
            width: Width of the plot. Defaults to 8 in.
            height: Height of the plot. Defaults to 6 in.
            term_zero: If True append zero voltage point at the end
            ax (plt.Axes): matplotlib axes object. Defaults to None.

        Returns:
            plt.Axes: matplotlib axes object.
        """
        ax = ax or pretty_plot(width, height)
        working_ion_symbols = set()
        formula = set()

        for key, electrode in self._electrodes.items():
            x, y = self.get_plot_data(electrode, term_zero=term_zero)
            working_ion_symbols.add(electrode.working_ion.symbol)
            formula.add(electrode.framework_formula)
            ax.plot(x, y, "-", linewidth=2, label=key)

        ax.legend()
        ax.set_xlabel(self._choose_best_x_label(formula=formula, work_ion_symbol=working_ion_symbols))
        ax.set_ylabel("Voltage (V)")
        plt.tight_layout()
        return ax

    def get_plotly_figure(
        self,
        width: float = 800,
        height: float = 600,
        font_dict: dict | None = None,
        term_zero: bool = True,
        **kwargs,
    ) -> plt.Figure:
        """Return plotly Figure object.

        Args:
            width: Width of the plot. Defaults to 800 px.
            height: Height of the plot. Defaults to 600 px.
            font_dict: define the font. Defaults to {"family": "Arial", "size": 24, "color": "#000000"}
            term_zero: If True append zero voltage point at the end
            **kwargs: passed to plotly.graph_objects.Layout
        """
        font_dict = font_dict or {"family": "Arial", "size": 24, "color": "#000000"}
        hover_temp = "Voltage (V): %{y:.2f}<br>x: %{x:.3f}"

        data = []
        working_ion_symbols = set()
        formula = set()
        for key, electrode in self._electrodes.items():
            x, y = self.get_plot_data(electrode, term_zero=term_zero)
            working_ion_symbols.add(electrode.working_ion.symbol)
            formula.add(electrode.framework_formula)
            # add Nones to x and y so vertical connecting lines are not plotted
            plot_x, plot_y = [x[0]], [y[0]]
            for i in range(1, len(x)):
                if x[i - 1] == x[i]:
                    plot_x.append(None)
                    plot_y.append(None)
                plot_x.append(x[i])
                plot_y.append(y[i])
            data.append(go.Scatter(x=plot_x, y=plot_y, name=key, hovertemplate=hover_temp))

        fig = go.Figure(
            data=data,
            layout=dict(
                title="Voltage vs. Capacity",
                width=width,
                height=height,
                font=font_dict,
                xaxis={"title": self._choose_best_x_label(formula=formula, work_ion_symbol=working_ion_symbols)},
                yaxis={"title": "Voltage (V)"},
                **kwargs,
            ),
        )

        fig.update_layout(template="plotly_white", title_x=0.5)
        return fig

    def _choose_best_x_label(self, formula: set[str], work_ion_symbol: set[str]) -> str:
        if self.xaxis in {"capacity", "capacity_grav"}:
            return "Capacity (mAh/g)"
        if self.xaxis == "capacity_vol":
            return "Capacity (Ah/l)"

        _formula: str | None = formula.pop() if len(formula) == 1 else None

        _work_ion_symbol: str | None = work_ion_symbol.pop() if len(work_ion_symbol) == 1 else None

        if self.xaxis == "x_form":
            if _formula and _work_ion_symbol:
                return f"x in {_work_ion_symbol}<sub>x</sub>{_formula}"
            return "x Work Ion per Host F.U."

        if self.xaxis == "frac_x":
            if _work_ion_symbol:
                return f"Atomic Fraction of {_work_ion_symbol}"
            return "Atomic Fraction of Working Ion"
        raise RuntimeError("No xaxis label can be determined")

    def show(self, width: float = 8, height: float = 6) -> None:
        """Show the voltage profile plot.

        Args:
            width: Width of the plot. Defaults to 8 in.
            height: Height of the plot. Defaults to 6 in.
        """
        self.get_plot(width, height).show()

    def save(self, filename: str, width: float = 8, height: float = 6) -> None:
        """Save the plot to an image file.

        Args:
            filename (str): Filename to save to. Must include extension to specify image format.
            width: Width of the plot. Defaults to 8 in.
            height: Height of the plot. Defaults to 6 in.
        """
        self.get_plot(width, height).savefig(filename)
