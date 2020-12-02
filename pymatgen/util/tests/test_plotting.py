import unittest

import matplotlib

from pymatgen.util.plotting import periodic_table_heatmap, van_arkel_triangle
from pymatgen.util.testing import PymatgenTest


class FuncTestCase(PymatgenTest):
    def test_plot_periodic_heatmap(self):
        random_data = {
            "Te": 0.11083818874391202,
            "Au": 0.7575629917425387,
            "Th": 1.2475885304040335,
            "Ni": -2.0354391922547705,
        }
        plt = periodic_table_heatmap(random_data, cmap="plasma")
        plt = periodic_table_heatmap(random_data)
        plt = periodic_table_heatmap(random_data, max_row=7)
        plt = periodic_table_heatmap(random_data, max_row=10)
        plt = periodic_table_heatmap(random_data, cbar_label_size=18)
        plt = periodic_table_heatmap(random_data, cmap_range=[0, 1])
        plt = periodic_table_heatmap(random_data, cbar_label="Hello World")
        plt = periodic_table_heatmap(random_data, blank_color="white")
        plt = periodic_table_heatmap(random_data, value_format="%.4f")

    def test_van_arkel_triangle(self):
        random_list = [("Fe", "C"), ("Ni", "F")]
        plt = van_arkel_triangle(random_list)
        plt = van_arkel_triangle(random_list, annotate=True)


if __name__ == "__main__":
    unittest.main()
