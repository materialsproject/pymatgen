# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest2 as unittest
import os

from pymatgen.analysis.pourbaix.maker import PourbaixDiagram
from pymatgen.analysis.pourbaix.entry import PourbaixEntryIO

try:
    from pymatgen.analysis.pourbaix.analyzer import PourbaixAnalyzer
except ImportError:
    PourbaixAnalyzer = None


@unittest.skipIf(PourbaixAnalyzer is None, "ImportError while importing PourbaixAnalyzer")
class TestPourbaixAnalyzer(unittest.TestCase):

    def setUp(self):
        module_dir = os.path.dirname(os.path.abspath(__file__))
        (elements, entries) = PourbaixEntryIO.from_csv(os.path.join(module_dir,
                                                    "test_entries.csv"))
        self.num_simplices = {"Zn(s)": 7, "ZnO2(s)": 7, "Zn[2+]": 4, "ZnO2[2-]": 4, "ZnHO2[-]": 4}
        self.e_above_hull_test = {"ZnHO[+]": 0.0693, "ZnO(aq)": 0.0624}
        self.decomp_test = {"ZnHO[+]": {"ZnO(s)": 0.5, "Zn[2+]": 0.5}, "ZnO(aq)": {"ZnO(s)": 1.0}}
        self.pd = PourbaixDiagram(entries)
        self.analyzer = PourbaixAnalyzer(self.pd)

    def test_get_facet_chempots(self):
        range_map = self.analyzer.get_chempot_range_map()
        range_map_dict = {}
        for PourEntry in range_map.keys():
            range_map_dict[PourEntry.name] = range_map[PourEntry]
        for entry in self.num_simplices.keys():
            self.assertEqual(len(range_map_dict[entry]), self.num_simplices[entry])

    def test_get_decomp(self):
        for entry in [entry for entry in self.pd.all_entries if entry not in self.pd.stable_entries]:
            decomp_entries = self.analyzer.get_decomposition(entry)
            for entr in decomp_entries:
                self.assertEqual(decomp_entries[entr], self.decomp_test[entry.name][entr.name])
            e_above_hull = self.analyzer.get_e_above_hull(entry)
            self.assertAlmostEqual(e_above_hull, self.e_above_hull_test[entry.name], 3)

if __name__ == '__main__':
    unittest.main()
