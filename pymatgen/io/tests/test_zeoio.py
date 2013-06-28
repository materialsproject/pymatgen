#!/usr/bin/env python

'''
'''

from __future__ import division


import unittest
import os

from pymatgen.io.cifio import CifParser
from pymatgen.io.zeoio import ZeoCssr, ZeoVoronoiXYZ, get_voronoi_nodes
from pymatgen.io.zeoio import get_void_volume_surfarea
from pymatgen.io.vaspio.vasp_input import Poscar
from pymatgen.core.structure import Structure, Molecule
from pymatgen.defects.point_defects import Vacancy
from pymatgen.analysis.bond_valence import BVAnalyzer
from pymatgen.core.periodic_table import Specie

test_dir = os.path.join('/Users/mbkumar/Research/Defects/pymatgen',
                        'test_files')


class ZeoCssrTest(unittest.TestCase):
    def setUp(self):
        filepath = os.path.join(test_dir, 'POSCAR')
        p = Poscar.from_file(filepath)
        self.zeocssr = ZeoCssr(p.structure)

    def test_str(self):
        expected_string = """10.4118 6.0672 4.7595
90.00 90.00 90.00 SPGR =  1 P 1    OPT = 1
24 0
0 Fe4 P4 O16
1 Fe 0.2187 0.7500 0.4749 0 0 0 0 0 0 0 0 0.0000
2 Fe 0.2813 0.2500 0.9749 0 0 0 0 0 0 0 0 0.0000
3 Fe 0.7187 0.7500 0.0251 0 0 0 0 0 0 0 0 0.0000
4 Fe 0.7813 0.2500 0.5251 0 0 0 0 0 0 0 0 0.0000
5 P 0.0946 0.2500 0.4182 0 0 0 0 0 0 0 0 0.0000
6 P 0.4054 0.7500 0.9182 0 0 0 0 0 0 0 0 0.0000
7 P 0.5946 0.2500 0.0818 0 0 0 0 0 0 0 0 0.0000
8 P 0.9054 0.7500 0.5818 0 0 0 0 0 0 0 0 0.0000
9 O 0.0434 0.7500 0.7071 0 0 0 0 0 0 0 0 0.0000
10 O 0.0966 0.2500 0.7413 0 0 0 0 0 0 0 0 0.0000
11 O 0.1657 0.0461 0.2854 0 0 0 0 0 0 0 0 0.0000
12 O 0.1657 0.4539 0.2854 0 0 0 0 0 0 0 0 0.0000
13 O 0.3343 0.5461 0.7854 0 0 0 0 0 0 0 0 0.0000
14 O 0.3343 0.9539 0.7854 0 0 0 0 0 0 0 0 0.0000
15 O 0.4034 0.7500 0.2413 0 0 0 0 0 0 0 0 0.0000
16 O 0.4566 0.2500 0.2071 0 0 0 0 0 0 0 0 0.0000
17 O 0.5434 0.7500 0.7929 0 0 0 0 0 0 0 0 0.0000
18 O 0.5966 0.2500 0.7587 0 0 0 0 0 0 0 0 0.0000
19 O 0.6657 0.0461 0.2146 0 0 0 0 0 0 0 0 0.0000
20 O 0.6657 0.4539 0.2146 0 0 0 0 0 0 0 0 0.0000
21 O 0.8343 0.5461 0.7146 0 0 0 0 0 0 0 0 0.0000
22 O 0.8343 0.9539 0.7146 0 0 0 0 0 0 0 0 0.0000
23 O 0.9034 0.7500 0.2587 0 0 0 0 0 0 0 0 0.0000
24 O 0.9566 0.2500 0.2929 0 0 0 0 0 0 0 0 0.0000"""
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
H 0.000000 0.000000 1.089000 0.200000
H 1.026719 0.000000 -0.363000 0.200000
H -0.513360 -0.889165 -0.363000 0.200000
H -0.513360 0.889165 -0.363000 0.200000"""
        self.assertEqual(str(self.xyz), ans)

    def test_from_file(self):
        filename = os.path.join(test_dir, "EDI_voro.xyz")
        vorXYZ = ZeoVoronoiXYZ.from_file(filename)
        self.assertIsInstance(vorXYZ.molecule, Molecule)

class GetVoronoiNodesTest(unittest.TestCase):
    def setUp(self):
        filepath = os.path.join(test_dir, 'POSCAR')
        p = Poscar.from_file(filepath)
        self.structure = p.structure
        bv = BVAnalyzer()
        valences = bv.get_valences(self.structure)
        self.rad_dict = {}
        valence_dict = {}
        sites = self.structure.sites
        for i in range(len(sites)):
            if sites[i].species_string in valence_dict.keys():
                continue
            else:
                valence_dict[sites[i].species_string] = valences[i]
                if len(valence_dict) == len(self.structure.composition):
                    break

        self.rad_dict = {}
        for el in valence_dict.keys():
            val = valence_dict[el]
            self.rad_dict[el] = Specie(el, val).ionic_radius
            if not self.rad_dict[el]: #get covalent radii
                pass
        #print self._rad_dict
        assert len(self.rad_dict) == len(self.structure.composition)

    def test_get_voronoi_nodes(self):
        vor_struct = get_voronoi_nodes(self.structure, self.rad_dict)
        self.assertIsInstance(vor_struct, Structure)

class GetVoidVolumeSurfaceTest(unittest.TestCase):
    def setUp(self):
        filepath = os.path.join(test_dir, 'POSCAR')
        p = Poscar.from_file(filepath)
        filepath1 = os.path.join(test_dir, 'LiMn2O4.cif')
        p1 = CifParser(filepath1).get_structures()[0]
        #vacancy = Vacancy(p.structure)
        self._vacancy = Vacancy(p1)

        um = [[1,0,0],[0,1,0],[0,0,1]]
        self._vac_struct = self._vacancy.make_supercells_with_defects(um)[1]
        #print self._vac_struct
    
    def test_void_volume_surface_area(self):
        #for struct in self._vac_struct:
        #   vol, sa = get_void_volume_surfarea(struct)
        #   self.assertIsInstance(vol, float)
        #   self.assertIsInstance(sa, float)
        vol, sa = get_void_volume_surfarea(
                self._vac_struct, self._vacancy.ionic_radii
                )
        print "vol:  ", vol, "sa:  ", sa
        self.assertIsInstance(vol, float)
        self.assertIsInstance(sa, float)

if __name__ == "__main__":
    #unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(GetVoronoiNodesTest)
    unittest.TextTestRunner().run(suite)
