"""This module defines Entry classes for containing experimental data."""

from __future__ import annotations

from monty.json import MSONable

from pymatgen.analysis.phase_diagram import PDEntry
from pymatgen.analysis.thermochemistry import ThermoData
from pymatgen.core.composition import Composition

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Jun 27, 2012"


class ExpEntry(PDEntry, MSONable):
    """
    An lightweight ExpEntry object containing experimental data for a
    composition for many purposes. Extends a PDEntry so that it can be used for
    phase diagram generation and reaction calculation.

    Current version works only with solid phases and at 298K. Further
    extensions for temperature dependence are planned.
    """

    def __init__(self, composition, thermodata, temperature=298):
        """
        Args:
            composition: Composition of the entry. For flexibility, this can take
                the form of all the typical input taken by a Composition, including
                a {symbol: amt} dict, a string formula, and others.
            thermodata: A sequence of ThermoData associated with the entry.
            temperature: A temperature for the entry in Kelvin. Defaults to 298K.
        """
        comp = Composition(composition)
        self._thermodata = thermodata
        found = False
        enthalpy = float("inf")
        for data in self._thermodata:
            if data.type == "fH" and data.value < enthalpy and (data.phaseinfo not in ("gas", "liquid")):
                enthalpy = data.value
                found = True
        if not found:
            raise ValueError("List of Thermodata does not contain enthalpy values.")
        self.temperature = temperature
        super().__init__(comp, enthalpy)

    def __repr__(self):
        return f"ExpEntry {self.composition.formula}, Energy = {self.energy:.4f}"

    @classmethod
    def from_dict(cls, d):
        """
        :param d: Dict representation.
        :return: ExpEntry
        """
        thermodata = [ThermoData.from_dict(td) for td in d["thermodata"]]
        return cls(d["composition"], thermodata, d["temperature"])

    def as_dict(self):
        """MSONable dict."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "thermodata": [td.as_dict() for td in self._thermodata],
            "composition": self.composition.as_dict(),
            "temperature": self.temperature,
        }
