from __future__ import annotations

import dataclasses
import json

import matplotlib.pyplot as plt
import pytest
from monty.json import MontyDecoder
from pytest import approx

from pymatgen.apps.battery.conversion_battery import ConversionElectrode
from pymatgen.apps.battery.insertion_battery import InsertionElectrode
from pymatgen.apps.battery.plotter import VoltageProfilePlotter
from pymatgen.core.composition import Composition
from pymatgen.entries.computed_entries import ComputedEntry
from tests.testing import TEST_FILES_DIR

TEST_DIR = f"{TEST_FILES_DIR}/apps/battery"


class TestVoltageProfilePlotter:
    def setup_method(self):
        entry_Li = ComputedEntry("Li", -1.90753119)

        with open(f"{TEST_DIR}/LiTiO2_batt.json", encoding="utf-8") as file:
            entries_LTO = json.load(file, cls=MontyDecoder)
        self.ie_LTO = InsertionElectrode.from_entries(entries_LTO, entry_Li)

        with open(f"{TEST_DIR}/FeF3_batt.json", encoding="utf-8") as file:
            entries = json.load(file, cls=MontyDecoder)
        self.ce_FF = ConversionElectrode.from_composition_and_entries(Composition("FeF3"), entries)

    def test_name(self):
        plotter = VoltageProfilePlotter(xaxis="frac_x")
        plotter.add_electrode(self.ie_LTO, "LTO insertion")
        plotter.add_electrode(self.ce_FF, "FeF3 conversion")
        assert plotter.get_plot_data(self.ie_LTO) is not None
        assert plotter.get_plot_data(self.ce_FF) is not None

    def test_plotly(self):
        plotter = VoltageProfilePlotter(xaxis="frac_x")
        plotter.add_electrode(self.ie_LTO, "LTO insertion")
        plotter.add_electrode(self.ce_FF, "FeF3 conversion")
        fig = plotter.get_plotly_figure()
        assert fig.layout.xaxis.title.text == "Atomic Fraction of Li"

        plotter = VoltageProfilePlotter(xaxis="x_form")
        plotter.add_electrode(self.ce_FF, "FeF3 conversion")
        fig = plotter.get_plotly_figure()
        assert fig.layout.xaxis.title.text == "x in Li<sub>x</sub>FeF3"

        plotter.add_electrode(self.ie_LTO, "LTO insertion")
        fig = plotter.get_plotly_figure()
        assert fig.layout.xaxis.title.text == "x Work Ion per Host F.U."


# Expected stair-step profiles: x values for LiTiO2 (insertion) and FeF3 (conversion),
# keyed by x-axis option. The voltages do not depend on the x-axis.
LTO_VOLTAGES = [2.78583901, 2.78583901, 0.89702381, 0.89702381]
FEF3_VOLTAGES = [3.07574726, 3.07574726, 3.05482526, 3.05482526]
EXPECTED_X = {
    "capacity": ([0, 160.803162, 160.803162, 315.177481], [0, 223.753636, 223.753636, 624.783676]),
    "capacity_grav": ([0, 160.803162, 160.803162, 315.177481], [0, 223.753636, 223.753636, 624.783676]),
    "capacity_vol": ([0, 586.861586, 586.861586, 1189.858567], [0, 811.966315, 811.966315, 2233.437651]),
    "x_form": ([0, 0.5, 0.5, 1], [0, 1, 1, 3]),
    "frac_x": ([0, 1 / 7, 1 / 7, 0.25], [0, 0.2, 0.2, 3 / 7]),
}
EXPECTED_X_LABEL = {
    "capacity": "Capacity (mAh/g)",
    "capacity_grav": "Capacity (mAh/g)",
    "capacity_vol": "Capacity (Ah/l)",
    "x_form": "x Work Ion per Host F.U.",
    "frac_x": "Atomic Fraction of Li",
}


@pytest.fixture(scope="module")
def ie_lto() -> InsertionElectrode:
    with open(f"{TEST_DIR}/LiTiO2_batt.json", encoding="utf-8") as file:
        entries = json.load(file, cls=MontyDecoder)
    return InsertionElectrode.from_entries(entries, ComputedEntry("Li", -1.90753119))


@pytest.fixture(scope="module")
def ce_fef3() -> ConversionElectrode:
    with open(f"{TEST_DIR}/FeF3_batt.json", encoding="utf-8") as file:
        entries = json.load(file, cls=MontyDecoder)
    return ConversionElectrode.from_composition_and_entries(Composition("FeF3"), entries)


@pytest.fixture
def close_figures():
    yield
    plt.close("all")


@pytest.mark.parametrize("xaxis", list(EXPECTED_X))
def test_get_plot_data(xaxis, ie_lto, ce_fef3):
    plotter = VoltageProfilePlotter(xaxis=xaxis)
    lto_x, fef3_x = EXPECTED_X[xaxis]

    x, y = plotter.get_plot_data(ie_lto)
    assert x == approx([*lto_x, lto_x[-1]], abs=1e-5)
    assert y == approx([*LTO_VOLTAGES, 0])

    x, y = plotter.get_plot_data(ce_fef3, term_zero=False)
    assert x == approx(fef3_x, abs=1e-5)
    assert y == approx(FEF3_VOLTAGES)


def test_get_plot_data_invalid_xaxis(ie_lto):
    plotter = VoltageProfilePlotter(xaxis="voltage")
    with pytest.raises(NotImplementedError, match="x_axis must be"):
        plotter.get_plot_data(ie_lto)
    # With no electrodes there is no data to fail on, but no label can be chosen either.
    with pytest.raises(RuntimeError, match="No xaxis label can be determined"):
        plotter.get_plotly_figure()


def test_hide_negative(ce_fef3):
    negative_pair = dataclasses.replace(ce_fef3[1], voltage=-0.5)
    electrode = dataclasses.replace(ce_fef3, voltage_pairs=[ce_fef3[0], negative_pair])

    x, y = VoltageProfilePlotter(xaxis="x_form").get_plot_data(electrode)
    assert y == approx([3.07574726, 3.07574726, -0.5, -0.5, 0])

    x, y = VoltageProfilePlotter(xaxis="x_form", hide_negative=True).get_plot_data(electrode)
    assert x == approx([0, 1, 1])
    assert y == approx([3.07574726, 3.07574726, 0])


def test_frac_x_label_mixed_working_ions(ce_fef3):
    with open(f"{TEST_DIR}/MnO2_batt.json", encoding="utf-8") as file:
        entries = json.load(file, cls=MontyDecoder)
    ce_mno2 = ConversionElectrode.from_composition_and_entries(Composition("MnO2"), entries, working_ion_symbol="Mg")
    plotter = VoltageProfilePlotter(xaxis="frac_x")
    plotter.add_electrode(ce_fef3, "FeF3 (Li)")
    plotter.add_electrode(ce_mno2, "MnO2 (Mg)")
    fig = plotter.get_plotly_figure()
    assert fig.layout.xaxis.title.text == "Atomic Fraction of Working Ion"


def test_add_electrode_default_labels(ie_lto, ce_fef3):
    plotter = VoltageProfilePlotter()
    plotter.add_electrode(ie_lto)
    plotter.add_electrode(ce_fef3)
    plotter.add_electrode(ie_lto, "custom")
    fig = plotter.get_plotly_figure()
    assert [trace.name for trace in fig.data] == ["Electrode 1", "Electrode 2", "custom"]


@pytest.mark.parametrize("xaxis", list(EXPECTED_X))
def test_get_plot(xaxis, ie_lto, ce_fef3, close_figures):
    plotter = VoltageProfilePlotter(xaxis=xaxis)
    plotter.add_electrode(ie_lto, "LTO")
    plotter.add_electrode(ce_fef3, "FeF3")
    ax = plotter.get_plot()
    assert isinstance(ax, plt.Axes)
    assert ax.get_xlabel() == EXPECTED_X_LABEL[xaxis]
    assert ax.get_ylabel() == "Voltage (V)"
    assert [text.get_text() for text in ax.get_legend().get_texts()] == ["LTO", "FeF3"]

    lto_line, fef3_line = ax.get_lines()
    assert lto_line.get_xdata() == approx([*EXPECTED_X[xaxis][0], EXPECTED_X[xaxis][0][-1]], abs=1e-5)
    assert lto_line.get_ydata() == approx([*LTO_VOLTAGES, 0])
    assert fef3_line.get_ydata() == approx([*FEF3_VOLTAGES, 0])


def test_get_plot_existing_ax_no_term_zero(ce_fef3, close_figures):
    _fig, ax = plt.subplots()
    plotter = VoltageProfilePlotter(xaxis="x_form")
    plotter.add_electrode(ce_fef3, "FeF3")
    assert plotter.get_plot(term_zero=False, ax=ax) is ax
    assert ax.get_xlabel() == "x in Li<sub>x</sub>FeF3"
    (line,) = ax.get_lines()
    assert line.get_ydata() == approx(FEF3_VOLTAGES)


@pytest.mark.parametrize("xaxis", list(EXPECTED_X))
def test_get_plotly_figure(xaxis, ie_lto, ce_fef3):
    plotter = VoltageProfilePlotter(xaxis=xaxis)
    plotter.add_electrode(ie_lto, "LTO")
    plotter.add_electrode(ce_fef3, "FeF3")
    fig = plotter.get_plotly_figure(width=500, height=400, term_zero=False, paper_bgcolor="red")
    assert fig.layout.width == 500
    assert fig.layout.height == 400
    assert fig.layout.paper_bgcolor == "red"
    assert fig.layout.font.family == "Arial"
    assert fig.layout.title.text == "Voltage vs. Capacity"
    assert fig.layout.xaxis.title.text == EXPECTED_X_LABEL[xaxis]
    assert fig.layout.yaxis.title.text == "Voltage (V)"
    lto_trace, fef3_trace = fig.data
    assert list(lto_trace.x) == approx(EXPECTED_X[xaxis][0], abs=1e-5)
    assert list(lto_trace.y) == approx(LTO_VOLTAGES)
    assert list(fef3_trace.x) == approx(EXPECTED_X[xaxis][1], abs=1e-5)
    assert list(fef3_trace.y) == approx(FEF3_VOLTAGES)


def test_get_plotly_figure_font(ce_fef3):
    plotter = VoltageProfilePlotter(xaxis="frac_x")
    plotter.add_electrode(ce_fef3)
    fig = plotter.get_plotly_figure(font_dict={"family": "Helvetica", "size": 12})
    assert fig.layout.font.family == "Helvetica"
    assert fig.layout.font.size == 12
    # term_zero appends a final zero-voltage point at the last x value
    assert list(fig.data[0].y)[-1] == 0
    assert list(fig.data[0].x)[-1] == approx(3 / 7)


def test_save(ie_lto, tmp_path, close_figures):
    plotter = VoltageProfilePlotter()
    plotter.add_electrode(ie_lto, "LTO")
    filename = tmp_path / "voltage_profile.png"
    plotter.save(str(filename))
    assert filename.is_file()
    assert filename.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"


def test_show(ie_lto, monkeypatch, close_figures):
    calls = []
    monkeypatch.setattr(plt, "show", lambda *args, **kwargs: calls.append(plt.gca()))
    plotter = VoltageProfilePlotter()
    plotter.add_electrode(ie_lto, "LTO")
    plotter.show()
    assert len(calls) == 1
    assert calls[0].get_xlabel() == "Capacity (mAh/g)"
