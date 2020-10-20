""" This module tests an f* diagram generator."""

import os
import unittest
import numpy as np
import pandas as pd
import plotly.express as px

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.Fstar.fstar import FStarDiagram
from pymatgen.util.testing import PymatgenTest
from pymatgen.io.cif import CifParser

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'tests')


class FStarDiagram_test(PymatgenTest):

    def setUp(self):
        from pymatgen.core.periodic_table import Element

        def cust_scatter(element, occupancy, ind, ind2):
            scat = np.log10(Element(element).Z) * occupancy  # no physical meaning, just an example of use.
            return scat

        c = cust_scatter
        cif_list = [i for i in os.listdir(os.getcwd()) if i.endswith('.cif')]
        struct_list = [CifParser(file).get_structures(primitive=False)[0] for file in cif_list]
        self.fstar_default = FStarDiagram(struct_list)
        self.fstar_neutron = FStarDiagram(struct_list, scattering_type='Neutron')
        self.fstar_xray = FStarDiagram(struct_list, scattering_type='X-ray')
        self.fstar_custom = FStarDiagram(struct_list, scattering_type='Custom', custom_scatter=c)

    def test_edit_fstar_diagram(self):
        self.assertEqual(len(self.fstar_default.site_labels), 3)
        new = FStarDiagram(self.struct_list)
        self.assertEqual(self.fstar_default.plot, new.plot)
        new.edit_fstar_diagram(combine_list=[['1Co', '0Li']])
        self.assertEqual(len(new.site_labels), 4)
        self.assertEqual(self.fstar_default.plot, new.plot)
        new.edit_fstar_diagram(plot_list=['1Co', '2O', '0Li'])
        self.assertTrue(self.fstar_default.plot != new.plot)

    def test_get_site_labels(self):

        self.assertEqual(len(self.fstar_default.get_site_labels()), len(self.fstar_default._equiv_inds[0]))
        self.assertTrue(
            self.fstar_default.get_site_labels()[0] != self.fstar_default.get_site_labels()[1] !=
            self.fstar_default.get_site_labels()[2])
        self.assertTrue(len(self.fstar_default.get_site_labels()) == len(self.fstar_neutron.get_site_labels()) ==
                        len(self.fstar_xray.get_site_labels()) == len(self.fstar_custom.get_site_labels()))

    def test_get_fstar_coords(self):
        self.assertEqual(len(self.fstar_default.get_fstar_coords()[0]), len(self.fstar_default.get_site_labels()))
        coord_sum_default = [sum(coord) for coord in self.fstar_default.get_fstar_coords()]
        self.assertEqual(sum(coord_sum_default), len(self.fstar_default.get_fstar_coords()))
        self.assertEqual(len(self.fstar_neutron.get_fstar_coords()[0]), len(self.fstar_neutron.get_site_labels()))
        coord_sum_neutron = [sum(coord) for coord in self.fstar_neutron.get_fstar_coords()]
        self.assertEqual(sum(coord_sum_neutron), len(self.fstar_neutron.get_fstar_coords()))
        self.assertEqual(len(self.fstar_xray.get_fstar_coords()[0]), len(self.fstar_xray.get_site_labels()))
        coord_sum_xray = [sum(coord) for coord in self.fstar_xray.get_fstar_coords()]
        self.assertEqual(sum(coord_sum_xray), len(self.fstar_xray.get_fstar_coords()))
        self.assertEqual(len(self.fstar_custom.get_fstar_coords()[0]), len(self.fstar_custom.get_site_labels()))
        coord_sum_custom = [sum(coord) for coord in self.fstar_custom.get_fstar_coords()]
        self.assertEqual(sum(coord_sum_custom), len(self.fstar_custom.get_fstar_coords()))


if __name__ == '__main__':
    unittest.main()
