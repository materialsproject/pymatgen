import unittest
from pymatgen.util.plotting import periodic_table_heatmap
from pymatgen.util.testing import PymatgenTest
import matplotlib


class FuncTestCase(PymatgenTest):

    def test_plot_periodic_heatmap(self):
        random_data = {'Te': 0.11083818874391202, 'Au': 0.7575629917425387,
                       'Th': 1.2475885304040335, 'Ni': -2.0354391922547705}
        plt = periodic_table_heatmap(random_data, cmap="plasma")
        plt = periodic_table_heatmap(random_data)
        plt = periodic_table_heatmap(random_data, max_row=7)
        plt = periodic_table_heatmap(random_data, max_row=10)


if __name__ == "__main__":
    unittest.main()