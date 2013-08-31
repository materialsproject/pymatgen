#!/usr/bin/env python

"""
A module to perform experimental thermochemical data analysis.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Jun 10, 2012"


from pymatgen.core.composition import Composition


STANDARD_TEMP = 298.0


class ThermoData(object):
    """
    A object container for an experimental Thermochemical Data.
    """
    def __init__(self, data_type, cpdname, phaseinfo, formula, value, ref="",
                 method="", temp_range=(298, 298), uncertainty=None):
        """
        Args:
            data_type:
                The thermochemical data type. Should be one of the following:
                    fH - Formation enthalpy
                    S - Entropy
                    A, B, C, D, E, F, G, H - variables for use in the various
                    equations for generating formation enthaplies or Cp at
                    various temperatures.
            cpdname:
                A name for the compound. For example, hematite for Fe2O3.
            phaseinfo:
                String denoting the phase. For example, "solid", "liquid",
                "gas" or "tetragonal".
            formula:
                A proper string formula, e.g., Fe2O3
            value:
                The value of the data.
            ref:
                A reference, if any, for the data.
            method:
                The method by which the data was determined, if available.
            temp_range:
                Temperature range of validity for the data in Kelvin. Defaults
                to 298 K only.
            uncertainty:
                An uncertainty for the data, if available.
        """
        self.type = data_type
        self.formula = formula
        self.composition = Composition(self.formula)
        self.reduced_formula = self.composition.reduced_formula
        self.compound_name = cpdname
        self.phaseinfo = phaseinfo
        self.value = value
        self.temp_range = temp_range
        self.method = method
        self.ref = ref
        self.uncertainty = uncertainty

    @classmethod
    def from_dict(cls, d):
        return ThermoData(d["type"], d["compound_name"],
                          d["phaseinfo"], d["formula"],
                          d["value"], d["ref"], d["method"],
                          d["temp_range"], d.get("uncertainty", None))

    @property
    def to_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__, "type": self.type,
                "formula": self.formula, "compound_name": self.compound_name,
                "phaseinfo": self.phaseinfo, "value": self.value,
                "temp_range": self.temp_range, "method": self.method,
                "ref": self.ref, "uncertainty": self.uncertainty}

    def __repr__(self):
        props = ["formula", "compound_name", "phaseinfo", "type", "temp_range",
                 "value", "method", "ref", "uncertainty"]
        output = ["{} : {}".format(k, getattr(self, k)) for k in props]
        return "\n".join(output)

    def __str__(self):
        return "{}_{}_{} = {}, Valid T : {}, Ref = {}".format(
            self.type, self.formula, self.phaseinfo, self.value,
            self.temp_range, self.ref
        )
