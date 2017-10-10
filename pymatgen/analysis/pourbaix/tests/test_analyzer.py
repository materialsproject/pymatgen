# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest
import os

from pymatgen.analysis.pourbaix.maker import PourbaixDiagram
from pymatgen.analysis.pourbaix.entry import PourbaixEntryIO, PourbaixEntry
from monty.serialization import loadfn

try:
    from pymatgen.analysis.pourbaix.analyzer import PourbaixAnalyzer
except ImportError:
    PourbaixAnalyzer = None

test_dir = os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', 'test_files')

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
        self.multi_data = loadfn(os.path.join(test_dir, 'multicomp_pbx.json'))
        
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
    
    def test_binary(self):
        # Test get_all_decomp_and_e_above_hull
        pd_binary = PourbaixDiagram(self.multi_data['binary'], 
                                    comp_dict = {"Ag": 0.5, "Te": 0.5})
        analyzer_binary = PourbaixAnalyzer(pd_binary)

        te_entry = pd_binary._unprocessed_entries[4]
        de, hull_e, entries = analyzer_binary.get_all_decomp_and_e_above_hull(te_entry)
        self.assertEqual(len(de), 10)
        self.assertEqual(len(hull_e), 10)
        self.assertAlmostEqual(hull_e[0], 5.4419792326439893)

    def test_ternary(self):
        # Ternary
        te_entry = self.multi_data['ternary'][20]
        pd_ternary = PourbaixDiagram(self.multi_data['ternary'],
                                     comp_dict = {"Ag": 0.33333, 
                                                  "Te": 0.33333,
                                                  "N": 0.33333})
        analyzer_ternary = PourbaixAnalyzer(pd_ternary)
        de, hull_e, entries = analyzer_ternary.get_all_decomp_and_e_above_hull(te_entry)
        self.assertEqual(len(de), 116)
        self.assertEqual(len(hull_e), 116)
        self.assertAlmostEqual(hull_e[0], 29.2520325229)

    def test_v_fe(self):
        v_fe_entries = self.multi_data['v_fe']
        pd = PourbaixDiagram(v_fe_entries, comp_dict = {"Fe": 0.5, "V": 0.5})
        analyzer_vfe = PourbaixAnalyzer(pd)
        entries = sorted(pd._all_entries, key=lambda x: x.energy)
        self.assertAlmostEqual(entries[1].weights[1], 0.6666666666)
        self.assertAlmostEqual(entries[1].energy, -110.77582628499995)
        self.assertAlmostEqual(entries[100].energy, -20.685496454465955)


if __name__ == '__main__':
    unittest.main()
