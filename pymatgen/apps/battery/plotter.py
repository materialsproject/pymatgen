#!/usr/bin/env python

'''
This module provides plotting capabilities for battery related applications.
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Jul 12, 2012"


from collections import OrderedDict
from pymatgen.util.plotting_utils import get_publication_quality_plot

class VoltageProfilePlotter(object):
    """
    A plotter to make voltage profile plots for batteries.
    """

    def __init__(self):
        self._electrodes = OrderedDict()

    def add_electrode(self, electrode, label=None):
        """
        Add an electrode to the plot.
        
        Args:
            electrode:
                An electrode. All electrodes satisfying the AbstractElectrode
                interface should work.
            label:
                A label for the electrode. If None, defaults to a counting
                system, i.e. 'Electrode 1', 'Electrode 2', ...
        """
        if not label:
            label = "Electrode {}".format(len(self._electrodes) + 1)
        self._electrodes[label] = electrode

    def get_plot(self, width, height):
        """
        Returns a plot object.
        
        Args:
            width:
                Width of the plot.
            height:
                Height of the plot.
                
        Returns:
            A matplotlib plot object.
        """
        plt = get_publication_quality_plot(width, height)
        for label, electrode in self._electrodes.items():
            xcoords = []
            ycoords = []
            cap = 0
            for vpair in electrode:
                xcoords.append(cap)

                cap += vpair.mAh / electrode.normalization_mass
                xcoords.append(cap)
                ycoords.extend([vpair.voltage] * 2)

            xcoords.append(cap)
            ycoords.append(0)

            plt.plot(xcoords, ycoords, '-', linewidth=2, label=label)

        plt.legend(self._electrodes.keys())
        plt.xlabel('Capacity (mAh/g)')
        plt.ylabel('Voltage (V)')
        plt.tight_layout()
        return plt

    def show(self, width=8, height=6):
        """
        Show the voltage profile plot.
        
        Args:
            width:
                Width of the plot. Defaults to 8 in.
            height:
                Height of the plot. Defaults to 6 in.
        """
        self.get_plot(width, height).show()

    def save(self, filename, image_format="eps", width=8, height=6):
        """
        Save the plot to an image file.
        
        Args:
            filename:
                Filename to save to.
            image_format:
                Format to save to. Defaults to eps.
        """
        self.get_plot(width, height).savefig(filename, format=image_format)
