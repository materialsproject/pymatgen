# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals


import unittest
import os

from pymatgen.analysis.pourbaix.maker import PourbaixDiagram
from pymatgen.analysis.pourbaix.entry import PourbaixEntryIO


class TestPourbaixDiagram(unittest.TestCase):
    def setUp(self):
        module_dir = os.path.dirname(os.path.abspath(__file__))
        (elements, entries) = PourbaixEntryIO.from_csv(os.path.join(module_dir,
                                                    "test_entries.csv"))
        self.entries = entries

    def test_pourbaix_diagram(self):
        pbx = PourbaixDiagram(self.entries)
        self.assertEqual(len(pbx.facets), 3, "Incorrect number of facets")
        self.assertEqual(set([e.name for e in pbx.stable_entries]),
                         {"ZnO(s)", "Zn[2+]", "ZnHO2[-]", "ZnO2[2-]", "Zn(s)"},
                         "List of stable entries does not match")

        pbx_nofilter = PourbaixDiagram(self.entries, filter_solids=False)
        self.assertEqual(len(pbx_nofilter.facets), 6,
                         "Incorrect number of facets in solid-unfiltered Pbx")
        self.assertEqual(set([e.name for e in pbx_nofilter.stable_entries]),
                         {"ZnO(s)", "Zn[2+]", "ZnHO2[-]", "ZnO2[2-]", "Zn(s)", "ZnO2(s)"},
                         "List of stable entries for unfiltered pbx does not match")

        pbx_lowconc = PourbaixDiagram(self.entries, conc_dict={"Zn": 1e-8})
        self.assertEqual(set([e.name for e in pbx_lowconc.stable_entries]),
                         {"ZnO(aq)", "Zn[2+]", "ZnHO2[-]", "ZnO2[2-]", "Zn(s)"})

    def test_get_pourbaix_domains(self):
        pbx = PourbaixDiagram(self.entries)
        domains = pbx.get_pourbaix_domains()
        import nose; nose.tools.set_trace()

if __name__ == '__main__':
    unittest.main()
