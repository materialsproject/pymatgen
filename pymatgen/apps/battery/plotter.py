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
from pymatgen.util.plotting import pretty_plot


class VoltageProfilePlotter:
    """
    A plotter to make voltage profile plots for batteries.

    Args:
        xaxis: The quantity to use as the xaxis. Can be either capacity (the
            default), or the frac_x.
    """

    def __init__(self, xaxis="capacity"):
        self._electrodes = OrderedDict()
        self.xaxis = xaxis

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

    def get_plot_data(self, electrode):
        x = []
        y = []
        cap = 0
        most_discharged = electrode[-1].frac_discharge
        norm = most_discharged / (1 - most_discharged)
        for vpair in electrode:
            if self.xaxis == "capacity":
                x.append(cap)
                cap += vpair.mAh / electrode.normalization_mass
                x.append(cap)
            else:
                x.append(vpair.frac_charge / (1 - vpair.frac_charge) / norm)
                x.append(vpair.frac_discharge / (1 - vpair.frac_discharge)
                         / norm)
            y.extend([vpair.voltage] * 2)

        x.append(x[-1])
        y.append(0)
        return x, y

    def get_plot(self, width=8, height=8):
        """
        Returns a plot object.

        Args:
            width: Width of the plot. Defaults to 8 in.
            height: Height of the plot. Defaults to 6 in.

        Returns:
            A matplotlib plot object.
        """
        plt = pretty_plot(width, height)
        for label, electrode in self._electrodes.items():
            (x, y) = self.get_plot_data(electrode)
            plt.plot(x, y, '-', linewidth=2, label=label)

        plt.legend()
        if self.xaxis == "capacity":
            plt.xlabel('Capacity (mAh/g)')
        else:
            plt.xlabel('Fraction')
        plt.ylabel('Voltage (V)')
        plt.tight_layout()
        return plt

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
