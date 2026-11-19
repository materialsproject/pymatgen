from __future__ import annotations

import pytest

from pymatgen.analysis.thermochemistry import ThermoData
from pymatgen.core import Composition


class TestThermoData:
    def setup_method(self):
        self.data = ThermoData(
            data_type="fH",
            cpdname="hematite",
            phaseinfo="solid",
            formula="Fe2O3",
            value=-825.5,
            ref="Kubaschewski",
            method="calorimetry",
            temp_range=(298, 1000),
            uncertainty=1.3,
        )

    def test_init(self):
        assert self.data.type == "fH"
        assert self.data.compound_name == "hematite"
        assert self.data.phaseinfo == "solid"
        assert self.data.formula == "Fe2O3"
        assert self.data.composition == Composition("Fe2O3")
        assert self.data.reduced_formula == "Fe2O3"
        assert self.data.value == pytest.approx(-825.5)
        assert self.data.temp_range == (298, 1000)
        assert self.data.uncertainty == pytest.approx(1.3)

    def test_defaults(self):
        data = ThermoData("S", "water", "liquid", "H4O2", 70.0)
        assert data.ref == ""
        assert data.method == ""
        assert data.temp_range == (298, 298)
        assert data.uncertainty is None
        # reduced_formula should reduce a non-reduced input formula
        assert data.reduced_formula == "H2O"

    def test_as_from_dict(self):
        dct = self.data.as_dict()
        assert dct["@module"] == "pymatgen.analysis.thermochemistry"
        assert dct["@class"] == "ThermoData"

        data = ThermoData.from_dict(dct)
        assert data.as_dict() == dct

    def test_from_dict_without_uncertainty(self):
        dct = self.data.as_dict()
        del dct["uncertainty"]
        assert ThermoData.from_dict(dct).uncertainty is None

    def test_str(self):
        assert str(self.data) == "fH_Fe2O3_solid = -825.5, Valid T : (298, 1000), Ref = Kubaschewski"

    def test_repr(self):
        lines = repr(self.data).splitlines()
        assert lines == [
            "formula : Fe2O3",
            "compound_name : hematite",
            "phaseinfo : solid",
            "type : fH",
            "temp_range : (298, 1000)",
            "value : -825.5",
            "method : calorimetry",
            "ref : Kubaschewski",
            "uncertainty : 1.3",
        ]
