from __future__ import annotations

import matplotlib.pyplot as plt

from pymatgen.util.plotting import periodic_table_heatmap, van_arkel_triangle
from pymatgen.util.testing import PymatgenTest

try:
    import pymatviz
    from plotly.graph_objects import Figure
except ImportError:
    pymatviz = None


class TestFunc(PymatgenTest):
    def test_plot_periodic_heatmap(self):
        random_data = {"Te": 0.11083, "Au": 0.75756, "Th": 1.24758, "Ni": -2.0354}
        fig = periodic_table_heatmap(random_data)
        if pymatviz:
            assert isinstance(fig, Figure)
        else:
            assert isinstance(fig, plt.Axes)

        # Test all keywords
        periodic_table_heatmap(
            random_data,
            cmap="plasma",
            max_row=10,
            cbar_label_size=18,
            cmap_range=[0, 1],
            cbar_label="Hello World",
            blank_color="white",
            value_format=".4f",
            edge_color="black",
            value_fontsize=12,
            symbol_fontsize=18,
            readable_fontcolor=True,
        )

    def test_van_arkel_triangle(self):
        random_list = [("Fe", "C"), ("Ni", "F")]
        ax = van_arkel_triangle(random_list)
        assert isinstance(ax, plt.Axes)
        assert ax.get_title() == ""
        assert ax.get_xlabel() == r"$\frac{\chi_{A}+\chi_{B}}{2}$"
        assert ax.get_ylabel() == r"$|\chi_{A}-\chi_{B}|$"
        ax = van_arkel_triangle(random_list, annotate=True)
        assert isinstance(ax, plt.Axes)
