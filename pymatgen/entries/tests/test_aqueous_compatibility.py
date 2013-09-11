#!/usr/bin/env python
"""
"""

__author__ = "Sai Jayaraman" 
__copyright__ = "Copyright 2013, The Materials Project" 
__version__ = "1.0"
__maintainer__ = "Sai Jayaraman"
__email__ = "sjayaram@mit.edu"
__date__ = "Mar 7, 2013"

import unittest

from pymatgen.entries.aqueouscompatibility import  MITAqueousCompatibility
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.core import Composition


class AqueousCompatibilityTest(unittest.TestCase):

    def setUp(self):
        self.compat = MITAqueousCompatibility()

    def test_compound_energy(self):
        entry = ComputedEntry(Composition("H2O"), -16)
        entry = self.compat.process_entry(entry)
        self.assertAlmostEqual(entry.energy, -15.1005, 4)

        entry = ComputedEntry(Composition("H2O"), -24)
        entry = self.compat.process_entry(entry)
        self.assertAlmostEqual(entry.energy, -15.1005, 4)

        entry = ComputedEntry(Composition("Cl"), -24)
        entry = self.compat.process_entry(entry)
        self.assertAlmostEqual(entry.energy, -1.6092, 4)


if __name__ == "__main__":
    unittest.main()
