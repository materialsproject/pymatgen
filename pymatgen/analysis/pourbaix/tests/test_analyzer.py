# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import os
import warnings
from pymatgen.analysis.pourbaix.maker import PourbaixDiagram
from pymatgen.analysis.pourbaix.entry import PourbaixEntryIO, PourbaixEntry
from pymatgen.analysis.pourbaix.analyzer import PourbaixAnalyzer
from pymatgen.util.testing import PymatgenTest

from monty.serialization import loadfn

test_dir = os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', 'test_files')

# TODO: refactor with more sensible binary/ternary test data
class TestPourbaixAnalyzer(PymatgenTest):

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
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.resetwarnings()

    def test_chempot_range_map(self):
        range_map = self.analyzer.get_chempot_range_map()
        range_map_dict = {}
        for PourEntry in range_map.keys():
            range_map_dict[PourEntry.name] = range_map[PourEntry]
        for entry in self.num_simplices.keys():
            self.assertEqual(len(range_map_dict[entry]), self.num_simplices[entry])
        ZnO2_entry = [e for e in self.pd.all_entries if e.name == "ZnO2(s)"][0]
        test_vertices = self.analyzer.pourbaix_domain_vertices[ZnO2_entry]
        self.assertArrayAlmostEqual(test_vertices[0], [16, 4])

    def test_get_decomp(self):
        for entry in [entry for entry in self.pd.all_entries
                      if entry not in self.pd.stable_entries]:
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

        te_entry = [e for e in pd_binary._unprocessed_entries
                    if e.composition.formula == "Te3"][0]
        de, hull_e, entries = analyzer_binary.get_all_decomp_and_e_above_hull(te_entry)
        self.assertEqual(len(de), 10)
        self.assertEqual(len(hull_e), 10)
        # Find a specific multientry to test
        tuples = zip(de, hull_e, entries)
        test_tuple = [t for t in tuples if t[2].name=='Te(s) + Ag[2+]'][0]
        self.assertAlmostEqual(test_tuple[1], 5.1396968548627315)

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
        tuples = zip(de, hull_e, entries)
        test_tuple = [t for t in tuples if t[2].name=='N2(s) + TeO4[2-] + Ag[2+]'][0]
        self.assertAlmostEqual(test_tuple[1], 50.337069095866745)

    def test_get_entry_stability(self):
        stab = self.analyzer.get_entry_stability(self.pd.all_entries[0], pH=0, V=1)
        self.assertAlmostEqual(stab, 3.88159999)

        # binary
        pd_binary = PourbaixDiagram(self.multi_data['binary'],
                                    comp_dict = {"Ag": 0.5, "Te": 0.5})
        analyzer_binary = PourbaixAnalyzer(pd_binary)
        ghull = analyzer_binary.get_entry_stability(self.multi_data['binary'][16], 
                                                    0, -5)

if __name__ == '__main__':
    unittest.main()
