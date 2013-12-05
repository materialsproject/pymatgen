#!/usr/bin/env python

'''
'''

from __future__ import division

__author__ = "Bharat Medasani"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "bkmedasani@lbl.gov"
__date__ = "Aug 2, 2013"

import unittest
import os
import re

from pymatgen.io.cifio import CifParser
from pymatgen.io.zeoio import *
from pymatgen.io.vaspio.vasp_input import Poscar
from pymatgen.core.structure import Structure, Molecule
from pymatgen.analysis.bond_valence import BVAnalyzer
from pymatgen.core.periodic_table import Specie

try:
    import zeo
except ImportError:
    zeo = None

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


@unittest.skipIf(not zeo, "zeo not present.")
class ZeoCssrTest(unittest.TestCase):
    def setUp(self):
        filepath = os.path.join(test_dir, 'POSCAR')
        p = Poscar.from_file(filepath)
        self.zeocssr = ZeoCssr(p.structure)

    def test_str(self):
        expected_string = """4.7595 10.4118 6.0672
90.00 90.00 90.00 SPGR =  1 P 1    OPT = 1
24 0
0 Fe4 P4 O16
1 Fe 0.4749 0.2187 0.7500 0 0 0 0 0 0 0 0 0.0000
2 Fe 0.9749 0.2813 0.2500 0 0 0 0 0 0 0 0 0.0000
3 Fe 0.0251 0.7187 0.7500 0 0 0 0 0 0 0 0 0.0000
4 Fe 0.5251 0.7813 0.2500 0 0 0 0 0 0 0 0 0.0000
5 P 0.4182 0.0946 0.2500 0 0 0 0 0 0 0 0 0.0000
6 P 0.9182 0.4054 0.7500 0 0 0 0 0 0 0 0 0.0000
7 P 0.0818 0.5946 0.2500 0 0 0 0 0 0 0 0 0.0000
8 P 0.5818 0.9054 0.7500 0 0 0 0 0 0 0 0 0.0000
9 O 0.7071 0.0434 0.7500 0 0 0 0 0 0 0 0 0.0000
10 O 0.7413 0.0966 0.2500 0 0 0 0 0 0 0 0 0.0000
11 O 0.2854 0.1657 0.0461 0 0 0 0 0 0 0 0 0.0000
12 O 0.2854 0.1657 0.4539 0 0 0 0 0 0 0 0 0.0000
13 O 0.7854 0.3343 0.5461 0 0 0 0 0 0 0 0 0.0000
14 O 0.7854 0.3343 0.9539 0 0 0 0 0 0 0 0 0.0000
15 O 0.2413 0.4034 0.7500 0 0 0 0 0 0 0 0 0.0000
16 O 0.2071 0.4566 0.2500 0 0 0 0 0 0 0 0 0.0000
17 O 0.7929 0.5434 0.7500 0 0 0 0 0 0 0 0 0.0000
18 O 0.7587 0.5966 0.2500 0 0 0 0 0 0 0 0 0.0000
19 O 0.2146 0.6657 0.0461 0 0 0 0 0 0 0 0 0.0000
20 O 0.2146 0.6657 0.4539 0 0 0 0 0 0 0 0 0.0000
21 O 0.7146 0.8343 0.5461 0 0 0 0 0 0 0 0 0.0000
22 O 0.7146 0.8343 0.9539 0 0 0 0 0 0 0 0 0.0000
23 O 0.2587 0.9034 0.7500 0 0 0 0 0 0 0 0 0.0000
24 O 0.2929 0.9566 0.2500 0 0 0 0 0 0 0 0 0.0000"""
        self.assertEqual(str(self.zeocssr), expected_string)

    def test_from_file(self):
        filename = os.path.join(test_dir, "EDI.cssr")
        zeocssr = ZeoCssr.from_file(filename)
        self.assertIsInstance(zeocssr.structure, Structure)

class ZeoVoronoiXYZTest(unittest.TestCase):
    def setUp(self):
        coords = [
                [0.000000, 0.000000, 0.000000],
                [0.000000, 0.000000, 1.089000],
                [1.026719, 0.000000, -0.363000],
                [-0.513360, -0.889165, -0.363000],
                [-0.513360, 0.889165, -0.363000]]
        prop = [0.4, 0.2, 0.2, 0.2, 0.2]
        self.mol = Molecule(
                ["C", "H", "H", "H", "H"], coords, 
                site_properties={"voronoi_radius":prop})
        self.xyz = ZeoVoronoiXYZ(self.mol)

    def test_str(self):
        ans = """5
H4 C1
C 0.000000 0.000000 0.000000 0.400000
H 1.089000 0.000000 0.000000 0.200000
H -0.363000 1.026719 0.000000 0.200000
H -0.363000 -0.513360 -0.889165 0.200000
H -0.363000 -0.513360 0.889165 0.200000"""
        self.assertEqual(str(self.xyz), ans)
        self.assertEqual(str(self.xyz), ans)

    def test_from_file(self):
        filename = os.path.join(test_dir, "EDI_voro.xyz")
        vorXYZ = ZeoVoronoiXYZ.from_file(filename)
        self.assertIsInstance(vorXYZ.molecule, Molecule)


@unittest.skipIf(not zeo, "zeo not present.")
class GetVoronoiNodesTest(unittest.TestCase):
    def setUp(self):
        filepath = os.path.join(test_dir, 'POSCAR')
        p = Poscar.from_file(filepath)
        self.structure = p.structure
        bv = BVAnalyzer()
        valences = bv.get_valences(self.structure)
        el = [site.species_string for site in self.structure.sites]
        valence_dict = dict(zip(el, valences))
        self.rad_dict = {}
        for k, v in valence_dict.items():
            self.rad_dict[k] = float(Specie(k,v).ionic_radius)

        assert len(self.rad_dict) == len(self.structure.composition)

    def test_get_voronoi_nodes(self):
        vor_struct = get_voronoi_nodes(self.structure, self.rad_dict)
        self.assertIsInstance(vor_struct, Structure)


@unittest.skip("The function is deprecated")
class GetVoidVolumeSurfaceTest(unittest.TestCase):
    def setUp(self):
        filepath1 = os.path.join(test_dir, 'Li2O.cif')
        p = CifParser(filepath1).get_structures(False)[0]
        bv = BVAnalyzer()
        valences = bv.get_valences(p)
        el = [site.species_string for site in p.sites]
        val_dict = dict(zip(el, valences))
        self._radii = {}
        for k,v in val_dict.items():
            k1 = re.sub('[1-9,+,\-]', '', k)
            self._radii[k1] = float(Specie(k1, v).ionic_radius)
        p.remove(0)
        self._vac_struct = p
    
    def test_void_volume_surface_area(self):
        pass
        vol, sa = get_void_volume_surfarea(
                self._vac_struct, self._radii
                )
        #print "vol:  ", vol, "sa:  ", sa
        self.assertIsInstance(vol, float)
        self.assertIsInstance(sa, float)

if __name__ == "__main__":
    unittest.main()
