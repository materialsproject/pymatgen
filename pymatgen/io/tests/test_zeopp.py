# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


__author__ = "Bharat Medasani"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "bkmedasani@lbl.gov"
__date__ = "Aug 2, 2013"

import unittest
import os
import re

from pymatgen.core.periodic_table import Species
from pymatgen.core.structure import Structure, Molecule
from pymatgen.io.cif import CifParser
from pymatgen.io.zeopp import ZeoCssr, ZeoVoronoiXYZ, get_voronoi_nodes, \
    get_high_accuracy_voronoi_nodes, get_void_volume_surfarea, \
    get_free_sphere_params
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.analysis.bond_valence import BVAnalyzer

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


# @unittest.skipIf(not zeo, "zeo not present.")
class ZeoCssrOxiTest(unittest.TestCase):
    def setUp(self):
        filepath = os.path.join(test_dir, 'POSCAR')
        p = Poscar.from_file(filepath)
        structure = BVAnalyzer().get_oxi_state_decorated_structure(p.structure)
        self.zeocssr = ZeoCssr(structure)

    def test_str(self):
        expected_string = """4.7595 10.4118 6.0672
90.00 90.00 90.00 SPGR =  1 P 1    OPT = 1
24 0
0 Fe4 P4 O16
1 Fe3+ 0.4749 0.2187 0.7500 0 0 0 0 0 0 0 0 0.0000
2 Fe3+ 0.9749 0.2813 0.2500 0 0 0 0 0 0 0 0 0.0000
3 Fe3+ 0.0251 0.7187 0.7500 0 0 0 0 0 0 0 0 0.0000
4 Fe3+ 0.5251 0.7813 0.2500 0 0 0 0 0 0 0 0 0.0000
5 P5+ 0.4182 0.0946 0.2500 0 0 0 0 0 0 0 0 0.0000
6 P5+ 0.9182 0.4054 0.7500 0 0 0 0 0 0 0 0 0.0000
7 P5+ 0.0818 0.5946 0.2500 0 0 0 0 0 0 0 0 0.0000
8 P5+ 0.5818 0.9054 0.7500 0 0 0 0 0 0 0 0 0.0000
9 O2- 0.7071 0.0434 0.7500 0 0 0 0 0 0 0 0 0.0000
10 O2- 0.7413 0.0966 0.2500 0 0 0 0 0 0 0 0 0.0000
11 O2- 0.2854 0.1657 0.0461 0 0 0 0 0 0 0 0 0.0000
12 O2- 0.2854 0.1657 0.4539 0 0 0 0 0 0 0 0 0.0000
13 O2- 0.7854 0.3343 0.5461 0 0 0 0 0 0 0 0 0.0000
14 O2- 0.7854 0.3343 0.9539 0 0 0 0 0 0 0 0 0.0000
15 O2- 0.2413 0.4034 0.7500 0 0 0 0 0 0 0 0 0.0000
16 O2- 0.2071 0.4566 0.2500 0 0 0 0 0 0 0 0 0.0000
17 O2- 0.7929 0.5434 0.7500 0 0 0 0 0 0 0 0 0.0000
18 O2- 0.7587 0.5966 0.2500 0 0 0 0 0 0 0 0 0.0000
19 O2- 0.2146 0.6657 0.0461 0 0 0 0 0 0 0 0 0.0000
20 O2- 0.2146 0.6657 0.4539 0 0 0 0 0 0 0 0 0.0000
21 O2- 0.7146 0.8343 0.5461 0 0 0 0 0 0 0 0 0.0000
22 O2- 0.7146 0.8343 0.9539 0 0 0 0 0 0 0 0 0.0000
23 O2- 0.2587 0.9034 0.7500 0 0 0 0 0 0 0 0 0.0000
24 O2- 0.2929 0.9566 0.2500 0 0 0 0 0 0 0 0 0.0000"""
        self.assertEqual(str(self.zeocssr), expected_string)

    def test_from_file(self):
        filename = os.path.join(test_dir, "EDI_oxistate_decorated.cssr")
        zeocssr = ZeoCssr.from_file(filename)
        self.assertIsInstance(zeocssr.structure, Structure)


@unittest.skipIf(not zeo, "zeo not present.")
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
            site_properties={"voronoi_radius": prop})
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
        vor = ZeoVoronoiXYZ.from_file(filename)
        self.assertIsInstance(vor.molecule, Molecule)


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
            self.rad_dict[k] = float(Species(k, v).ionic_radius)

        assert len(self.rad_dict) == len(self.structure.composition)

    def test_get_voronoi_nodes(self):
        vor_node_struct, vor_edge_center_struct, vor_face_center_struct = \
            get_voronoi_nodes(self.structure, self.rad_dict)
        self.assertIsInstance(vor_node_struct, Structure)
        self.assertIsInstance(vor_edge_center_struct, Structure)
        self.assertIsInstance(vor_face_center_struct, Structure)
        print(len(vor_node_struct.sites))
        print(len(vor_face_center_struct.sites))


@unittest.skip("file free_sph.cif not present")
class GetFreeSphereParamsTest(unittest.TestCase):
    def setUp(self):
        filepath = os.path.join(test_dir, 'free_sph.cif')
        self.structure = Structure.from_file(filepath)
        self.rad_dict = {'Ge': 0.67, 'P': 0.52, 'S': 1.7,
                         'La': 1.17, 'Zr': 0.86, 'O': 1.26}

    def test_get_free_sphere_params(self):
        free_sph_params = get_free_sphere_params(self.structure,
                                                 rad_dict=self.rad_dict)
        # Zeo results can change in future. Hence loose comparison
        self.assertAlmostEqual(
            free_sph_params['inc_sph_max_dia'], 2.58251, places=1)
        self.assertAlmostEqual(
            free_sph_params['free_sph_max_dia'], 1.29452, places=1)
        self.assertAlmostEqual(
            free_sph_params['inc_sph_along_free_sph_path_max_dia'],
            2.58251, places=1)


@unittest.skipIf(not zeo, "zeo not present.")
class GetHighAccuracyVoronoiNodesTest(unittest.TestCase):
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
            self.rad_dict[k] = float(Species(k, v).ionic_radius)

        assert len(self.rad_dict) == len(self.structure.composition)

    def test_get_voronoi_nodes(self):
        # vor_node_struct, vor_ec_struct, vor_fc_struct = \
        #    get_high_accuracy_voronoi_nodes(self.structure, self.rad_dict)
        vor_node_struct = \
            get_high_accuracy_voronoi_nodes(self.structure, self.rad_dict)
        self.assertIsInstance(vor_node_struct, Structure)
        # self.assertIsInstance(vor_ec_struct, Structure)
        # self.assertIsInstance(vor_fc_struct, Structure)
        print(len(vor_node_struct.sites))
        # print(len(vor_fc_struct.sites))


@unittest.skipIf(not zeo, "zeo not present.")
class GetVoronoiNodesMultiOxiTest(unittest.TestCase):
    def setUp(self):
        filepath = os.path.join(test_dir, 'POSCAR')
        p = Poscar.from_file(filepath)
        self.structure = p.structure
        bv = BVAnalyzer()
        self.structure = bv.get_oxi_state_decorated_structure(self.structure)
        valences = bv.get_valences(self.structure)
        radii = []
        for i in range(len(valences)):
            el = self.structure.sites[i].specie.symbol
            radius = Species(el, valences[i]).ionic_radius
            radii.append(radius)
        el = [site.species_string for site in self.structure.sites]
        self.rad_dict = dict(zip(el, radii))
        for el in self.rad_dict.keys():
            print((el, self.rad_dict[el].real))

    def test_get_voronoi_nodes(self):
        vor_node_struct, vor_edge_center_struct, vor_face_center_struct = \
            get_voronoi_nodes(self.structure, self.rad_dict)
        self.assertIsInstance(vor_node_struct, Structure)
        self.assertIsInstance(vor_edge_center_struct, Structure)
        self.assertIsInstance(vor_face_center_struct, Structure)


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
        for k, v in val_dict.items():
            k1 = re.sub(r'[1-9,+,\-]', '', k)
            self._radii[k1] = float(Species(k1, v).ionic_radius)
        p.remove(0)
        self._vac_struct = p

    def test_void_volume_surface_area(self):
        pass
        vol, sa = get_void_volume_surfarea(self._vac_struct, self._radii)
        # print "vol:  ", vol, "sa:  ", sa
        self.assertIsInstance(vol, float)
        self.assertIsInstance(sa, float)


if __name__ == "__main__":
    unittest.main()
