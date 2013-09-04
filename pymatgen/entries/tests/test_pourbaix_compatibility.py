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

from pymatgen.entries.pourbaixcompatibility import  MITPourbaixCompatibility, oxide_type
from pymatgen.core.periodic_table import Element
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.entries.computed_entries import ComputedStructureEntry, ComputedEntry
from pymatgen.core import Composition


class PourbaixCompatibilityTest(unittest.TestCase):

    def setUp(self):
        self.compat = MITPourbaixCompatibility()

    def test_no_struct_compat(self):
        lio2_entry_nostruct = ComputedEntry(Composition("Li2O4"), -29.1943757, data={"oxide_type": "superoxide"})
        lio2_entry_corrected = self.compat.process_entry(lio2_entry_nostruct)
        self.assertAlmostEqual(lio2_entry_corrected.energy, -27.0367157, 4)

        lio2_entry_nostruct = ComputedEntry(Composition("Li2(OH)2"), -29.85285134, data={"oxide_type": "hydroxide"})
        lio2_entry_corrected = self.compat.process_entry(lio2_entry_nostruct)
        self.assertAlmostEqual(lio2_entry_corrected.energy, -28.9047, 4)

    def test_process_entry_superoxide(self):
        el_li = Element("Li")
        el_o = Element("O")
        latt = Lattice([[3.985034, 0.0, 0.0],
                        [0.0, 4.881506, 0.0],
                        [0.0, 0.0, 2.959824]])
        elts = [el_li, el_li, el_o, el_o, el_o, el_o]
        coords = list()
        coords.append([0.500000, 0.500000, 0.500000])
        coords.append([0.0, 0.0, 0.0])
        coords.append([0.632568, 0.085090, 0.500000])
        coords.append([0.367432, 0.914910, 0.500000])
        coords.append([0.132568, 0.414910, 0.000000])
        coords.append([0.867432, 0.585090, 0.000000])
        struct = Structure(latt, elts, coords)
        lio2_entry = ComputedStructureEntry(struct, -29.1943757)
        lio2_entry_corrected = self.compat.process_entry(lio2_entry)
        self.assertAlmostEqual(lio2_entry_corrected.energy, -27.0367157, 4)

    
    def test_process_entry_peroxide(self):
        latt = Lattice.from_parameters(3.159597, 3.159572, 7.685205, 89.999884, 89.999674, 60.000510)
        el_li = Element("Li")
        el_o = Element("O")
        elts = [el_li, el_li, el_li, el_li, el_o, el_o, el_o, el_o]
        coords = [[0.666656, 0.666705, 0.750001], 
                  [0.333342, 0.333378, 0.250001], 
                  [0.000001, 0.000041, 0.500001],
                  [0.000001, 0.000021, 0.000001], 
                  [0.333347, 0.333332, 0.649191], 
                  [0.333322, 0.333353, 0.850803], 
                  [0.666666, 0.666686, 0.350813], 
                  [0.666665, 0.666684, 0.149189]]
        struct = Structure(latt, elts, coords)
        li2o2_entry = ComputedStructureEntry(struct, -38.80421304)
        li2o2_entry_corrected = self.compat.process_entry(li2o2_entry)
        self.assertAlmostEqual(li2o2_entry_corrected.energy, -37.80087, 4)
        
    def test_process_entry_ozonide(self):
        el_li = Element("Li")
        el_o = Element("O")
        elts = [el_li, el_o, el_o, el_o]
        latt = Lattice.from_parameters(3.999911, 3.999911, 3.999911, 133.847504, 102.228244, 95.477342)
        coords = [[0.513004, 0.513004, 1.000000],
                  [0.017616, 0.017616, 0.000000], 
                  [0.649993, 0.874790, 0.775203], 
                  [0.099587, 0.874790, 0.224797]]
        struct = Structure(latt, elts, coords)
        lio3_entry = ComputedStructureEntry(struct, -19.20762604)
        lio3_entry_corrected = self.compat.process_entry(lio3_entry)
        self.assertEqual(lio3_entry_corrected, None)

    def test_process_entry_hydroxide(self):
        el_li = Element("Li")
        el_o = Element("O")
        el_h = Element("H")
        latt = Lattice.from_parameters(3.565276, 3.565276, 4.384277, 90.000000, 90.000000, 90.000000)
        elts = [el_h, el_h, el_li, el_li, el_o, el_o]
        coords = [[0.000000, 0.500000, 0.413969],
                  [0.500000, 0.000000, 0.586031], 
                  [0.000000, 0.000000, 0.000000],
                  [0.500000, 0.500000, 0.000000],
                  [0.000000, 0.500000, 0.192672],
                  [0.500000, 0.000000, 0.807328]]
        struct = Structure(latt, elts, coords)
        lioh_entry = ComputedStructureEntry(struct, -29.85285134)
        lioh_entry_corrected = self.compat.process_entry(lioh_entry)
        self.assertAlmostEqual(lioh_entry_corrected.energy, -28.9047, 4)
        
    def test_oxide_type(self):
        el_li = Element("Li")
        el_o = Element("O")
        latt = Lattice([[3.985034, 0.0, 0.0],
                        [0.0, 4.881506, 0.0],
                        [0.0, 0.0, 2.959824]])
        elts = [el_li, el_li, el_o, el_o, el_o, el_o]
        coords = list()
        coords.append([0.500000, 0.500000, 0.500000])
        coords.append([0.0, 0.0, 0.0])
        coords.append([0.632568, 0.085090, 0.500000])
        coords.append([0.367432, 0.914910, 0.500000])
        coords.append([0.132568, 0.414910, 0.000000])
        coords.append([0.867432, 0.585090, 0.000000])
        struct = Structure(latt, elts, coords)
        self.assertEqual(oxide_type(struct, 1.1), "superoxide")
        
        el_li = Element("Li")
        el_o = Element("O")
        elts = [el_li, el_o, el_o, el_o]
        latt = Lattice.from_parameters(3.999911, 3.999911, 3.999911, 133.847504, 102.228244, 95.477342)
        coords = [[0.513004, 0.513004, 1.000000],
                  [0.017616, 0.017616, 0.000000], 
                  [0.649993, 0.874790, 0.775203], 
                  [0.099587, 0.874790, 0.224797]]
        struct = Structure(latt, elts, coords)
        self.assertEqual(oxide_type(struct, 1.1), "ozonide")

        latt = Lattice.from_parameters(3.159597, 3.159572, 7.685205, 89.999884, 89.999674, 60.000510)
        el_li = Element("Li")
        el_o = Element("O")
        elts = [el_li, el_li, el_li, el_li, el_o, el_o, el_o, el_o]
        coords = [[0.666656, 0.666705, 0.750001], 
                  [0.333342, 0.333378, 0.250001], 
                  [0.000001, 0.000041, 0.500001],
                  [0.000001, 0.000021, 0.000001], 
                  [0.333347, 0.333332, 0.649191], 
                  [0.333322, 0.333353, 0.850803], 
                  [0.666666, 0.666686, 0.350813], 
                  [0.666665, 0.666684, 0.149189]]
        struct = Structure(latt, elts, coords)
        self.assertEqual(oxide_type(struct, 1.1), "peroxide")
    
        el_li = Element("Li")
        el_o = Element("O")
        el_h = Element("H")
        latt = Lattice.from_parameters(3.565276, 3.565276, 4.384277, 90.000000, 90.000000, 90.000000)
        elts = [el_h, el_h, el_li, el_li, el_o, el_o]
        coords = [[0.000000, 0.500000, 0.413969],
                  [0.500000, 0.000000, 0.586031], 
                  [0.000000, 0.000000, 0.000000],
                  [0.500000, 0.500000, 0.000000],
                  [0.000000, 0.500000, 0.192672],
                  [0.500000, 0.000000, 0.807328]]
        struct = Structure(latt, elts, coords)
        self.assertEqual(oxide_type(struct, 1.1), "hydroxide")
        
        el_li = Element("Li")
        el_n = Element("N")
        el_h = Element("H")
        latt = Lattice.from_parameters(3.565276, 3.565276, 4.384277, 90.000000, 90.000000, 90.000000)
        elts = [el_h, el_h, el_li, el_li, el_n, el_n]
        coords = [[0.000000, 0.500000, 0.413969],
                  [0.500000, 0.000000, 0.586031], 
                  [0.000000, 0.000000, 0.000000],
                  [0.500000, 0.500000, 0.000000],
                  [0.000000, 0.500000, 0.192672],
                  [0.500000, 0.000000, 0.807328]]
        struct = Structure(latt, elts, coords)
        self.assertEqual(oxide_type(struct, 1.1), "None")
        


if __name__ == "__main__":
    unittest.main()
