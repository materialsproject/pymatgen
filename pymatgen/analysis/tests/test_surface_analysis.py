import unittest
import os
import random
import glob
<<<<<<< HEAD
import json
=======
>>>>>>> 8d1bdb462de8d36d908b8f58eef65379482abab6

import numpy as np

from pymatgen import SETTINGS
from pymatgen.io.vasp.outputs import Vasprun
<<<<<<< HEAD
from pymatgen.analysis.surface_analysis import SurfaceEnergyCalculator, \
    SurfaceEnergyPlotter, Composition
from pymatgen.util.testing import PymatgenTest
from pymatgen.entries.computed_entries import ComputedStructureEntry
=======
from pymatgen.analysis.surface_analysis import SurfaceEnergyAnalyzer
from pymatgen.util.testing import PymatgenTest
>>>>>>> 8d1bdb462de8d36d908b8f58eef65379482abab6

__author__ = "Richard Tran"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Richard Tran"
__email__ = "rit001@eng.ucsd.edu"
__date__ = "Aug 24, 2017"


def get_path(path_str):
    cwd = os.path.abspath(os.path.dirname(__file__))
<<<<<<< HEAD
    path = os.path.join(cwd, "..", "..", "..", "test_files",
                        "surface_tests", path_str)
    return path

@unittest.skipIf(not SETTINGS.get("PMG_MAPI_KEY"),
                 "PMG_MAPI_KEY environment variable not set.")

class SurfaceEnergyCalculatorTest(PymatgenTest):

    def setUp(self):

        self.entry_dict = get_entry_dict(os.path.join(get_path(""),
                                                      "Cu_entries.txt"))
        c = Composition("Cu")
        self.Cu_analyzer = SurfaceEnergyCalculator("mp-30", c, c)
=======
    path = os.path.join(cwd, "..", "..", "..", "test_files", "surface_tests",
                        path_str)
    return path

@unittest.skipIf(not SETTINGS.get("PMG_MAPI_KEY"), "PMG_MAPI_KEY environment variable not set.")
class SurfaceEnergyAnalyzerTest(PymatgenTest):

    def setUp(self):



        vasprun_dict = {}
        for v in glob.glob(os.path.join(get_path(""), "*")):
            if ".xml.Cu" in v:
                vasprun_dict[tuple([int(i) for i in v[-6:].strip(".gz")])] = [Vasprun(v)]
        self.vasprun_dict = vasprun_dict
        self.Cu_analyzer = SurfaceEnergyAnalyzer("mp-30", self.vasprun_dict, "Cu")
>>>>>>> 8d1bdb462de8d36d908b8f58eef65379482abab6

    def test_gamma_calculator(self):

        # make sure we've loaded all our files correctly
<<<<<<< HEAD
        self.assertEqual(len(self.entry_dict.keys()), 13)
        all_se = []
        for hkl, entries in self.entry_dict.items():
            u1, u2 = self.Cu_analyzer.chempot_range()
            se1 = self.Cu_analyzer.calculate_gamma_at_u(entries[0], u1)
            se2 = self.Cu_analyzer.calculate_gamma_at_u(entries[0], u2)
            # For a stoichiometric system, we expect surface
            # energy to be independent of chemical potential
            self.assertEqual(se1, se2)
            all_se.append(se1)

        # The (111) facet should be the most stable
        self.assertEqual(min(all_se),
                         self.Cu_analyzer.calculate_gamma_at_u(\
                             self.entry_dict[(1,1,1)][0], 0))

    def test_chempot_range(self):

        # just check we always get two values in a list
        chempot = self.Cu_analyzer.chempot_range()
        self.assertEqual(len(chempot), 2)

class SurfaceEnergyPlotterTest(PymatgenTest):

    def setUp(self):

        entry_dict = get_entry_dict(os.path.join(get_path(""),
                                                 "Cu_entries.txt"))
        c = Composition("Cu")
        calculator = SurfaceEnergyCalculator("mp-30", c, c)
        self.Cu_analyzer = SurfaceEnergyPlotter(entry_dict, calculator)
=======
        self.assertEqual(len(self.vasprun_dict.keys()), 13)
        for hkl, vaspruns in self.vasprun_dict.items():
            se_range = self.Cu_analyzer.calculate_gamma(vaspruns[0])
            # For a stoichiometric system, we expect surface
            # energy to be independent of chemical potential
            self.assertEqual(se_range[0], se_range[1])

    def test_get_intersections(self):

        # Just test if its working, for Cu, there are no different
        #  terminations so everything should be a nonetype
        for hkl in self.vasprun_dict.keys():
            self.assertFalse(self.Cu_analyzer.get_intersections(hkl))
>>>>>>> 8d1bdb462de8d36d908b8f58eef65379482abab6

    def test_get_wulff_shape_dict(self):

        # for pure Cu, all facets are independent of chemical potential,
        # so we assume Wulff does not change wrt chemical potential
        wulff_dict = self.Cu_analyzer.wulff_shape_dict(at_intersections=True)
        self.assertEqual(len(wulff_dict.keys()), 1)
        wulffshape = list(wulff_dict.values())[0]
        # The Wulff shape of Cu should have at least 70% (100) and (111) facets
        area_fraction_dict = wulffshape.area_fraction_dict
        self.assertGreater(area_fraction_dict[(1,0,0)]+\
                           area_fraction_dict[(1,1,1)], 0.7)
        # test out self.wulff_shape_from_chempot(), all Wulff
        # shapes should the same regardless of chemical potential
        wulff1 = self.Cu_analyzer.wulff_shape_from_chempot(\
            min(self.Cu_analyzer.chempot_range))
        wulff2 = self.Cu_analyzer.wulff_shape_from_chempot(\
            max(self.Cu_analyzer.chempot_range))
<<<<<<< HEAD
        for hkl in self.Cu_analyzer.entry_dict.keys():
            self.assertEqual(wulff1.area_fraction_dict[hkl],
                             wulff2.area_fraction_dict[hkl])

    def test_get_intersections(self):

        # Just test if its working, for Cu, there are no different
        #  terminations so everything should be a nonetype
        for hkl in self.Cu_analyzer.entry_dict.keys():
            self.assertFalse(self.Cu_analyzer.get_intersections(hkl))


if __name__ == "__main__":
    unittest.main()

def get_entry_dict(filename):
    # helper to generate an entry_dict

    entry_dict = {}
    with open(filename) as entries:
        entries = json.loads(entries.read())
    for k in entries.keys():
        n = k[25:]
        miller_index = []
        for i, s in enumerate(n):
            if i > 2:
                break
            miller_index.append(int(s))
        hkl = tuple(miller_index)
        if hkl not in entry_dict.keys():
            entry_dict[hkl] = []
        entry_dict[hkl].append(ComputedStructureEntry.from_dict(entries[k]))

    return entry_dict
=======
        for hkl in self.Cu_analyzer.vasprun_dict.keys():
            self.assertEqual(wulff1.area_fraction_dict[hkl],
                             wulff2.area_fraction_dict[hkl])

if __name__ == "__main__":
    unittest.main()
>>>>>>> 8d1bdb462de8d36d908b8f58eef65379482abab6
