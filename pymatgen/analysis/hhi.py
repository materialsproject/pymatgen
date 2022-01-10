# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module is used to estimate the Herfindahl-Hirschman Index, or HHI, of
chemical compounds. The HHI is a measure of how geographically confined or
dispersed the elements comprising a compound are. A low HHI is desirable
because it means the component elements are geographically dispersed.

Data/strategy from "Data-Driven Review of Thermoelectric Materials:
Performance and Resource Considerations" by Gaultois et al., published
in Chemistry of Materials (2013).
"""

import os

from monty.design_patterns import singleton

from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element

__author__ = "Anubhav Jain"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Anubhav Jain"
__email__ = "ajain@lbl.gov"
__date__ = "Oct 27, 2014"


csv_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "hhi_data.csv")


@singleton
class HHIModel:
    """
    HHI calculator.
    """

    def __init__(self):
        """
        Init for HHIModel.
        """
        self.symbol_hhip_hhir = {}  # symbol->(HHI_production, HHI reserve)

        with open(csv_path) as f:
            for line in f:
                if line[0] != "#":
                    symbol, hhi_production, hhi_reserve = line.split(",")
                    self.symbol_hhip_hhir[symbol] = (
                        float(hhi_production),
                        float(hhi_reserve),
                    )

    def _get_hhi_el(self, el_or_symbol):
        """
        Returns the tuple of HHI_production, HHI reserve for a single element only
        """
        if isinstance(el_or_symbol, Element):
            el_or_symbol = el_or_symbol.symbol

        return (
            self.symbol_hhip_hhir[el_or_symbol][0],
            self.symbol_hhip_hhir[el_or_symbol][1],
        )

    def get_hhi(self, comp_or_form):
        """
        Gets the reserve and production HHI for a compound.

        Args:
            comp_or_form (Composition or String): A Composition or String formula

        Returns:
            A tuple representing the (HHI_production, HHI_reserve)
        """

        try:
            if not isinstance(comp_or_form, Composition):
                comp_or_form = Composition(comp_or_form)

            hhi_p = 0
            hhi_r = 0

            for e in comp_or_form.elements:
                percent = comp_or_form.get_wt_fraction(e)
                dp, dr = self._get_hhi_el(e)
                hhi_p += dp * percent
                hhi_r += dr * percent
            return hhi_p, hhi_r

        except Exception:
            return None, None

    def get_hhi_production(self, comp_or_form):
        """
        Gets the production HHI for a compound.

        Args:
            comp_or_form (Composition or String): A Composition or String formula

        Returns:
            The HHI production value
        """
        return self.get_hhi(comp_or_form)[0]

    def get_hhi_reserve(self, comp_or_form):
        """
        Gets the reserve HHI for a compound.

        Args:
            comp_or_form (Composition or String): A Composition or String formula

        Returns:
            The HHI reserve value
        """
        return self.get_hhi(comp_or_form)[1]

    @staticmethod
    def get_hhi_designation(hhi):
        """
        Gets a designation for low, medium, high HHI, as specified in "U.S.
        Department of Justice and the Federal Trade Commission, Horizontal
        merger guidelines; 2010."

        Args:
            hhi (float): HHI value

        Returns:
            The designation as String
        """

        if hhi is None:
            return None

        if 0 <= hhi < 1500:
            return "low"

        if 1500 <= hhi <= 2500:
            return "medium"

        return "high"
