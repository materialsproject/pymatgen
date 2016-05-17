# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest2 as unittest
import os

from pymatgen.analysis.pourbaix.maker import PourbaixDiagram
from pymatgen.analysis.pourbaix.entry import PourbaixEntryIO


class TestPourbaixDiagram(unittest.TestCase):

    def setUp(self):
        module_dir = os.path.dirname(os.path.abspath(__file__))
        (elements, entries) = PourbaixEntryIO.from_csv(os.path.join(module_dir,
                                                    "test_entries.csv"))
        self._entries = entries
        self._pd = PourbaixDiagram(entries)
        self.list_of_stable_entries = ["ZnO(s)", "Zn[2+]", "ZnO2(s)", "ZnHO2[-]", "ZnO2[2-]", "Zn(s)"]

        
    def test_pourbaix_diagram(self):
        self.assertEqual(len(self._pd.facets), 6, "Incorrect number of facets")
        for entry in self._pd.stable_entries:
            self.assertIn(entry.name, self.list_of_stable_entries, "List of stable entries does not match")

if __name__ == '__main__':
    unittest.main()
