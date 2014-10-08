#!/usr/bin/python


import unittest

from pymatgen.core.lattice import Lattice
from pymatgen.io.smartio import CifParser
from pymatgen import write_structure
from pymatgen.core.surface import Slab, SurfaceGenerator
import os
from pymatgen.core import Structure
import itertools



def get_path(path_str):
    file_name = "pymatgen/pymatgen/core/tests/surface_tests/" + path_str
    path = os.path.join(os.path.expanduser("~"), file_name)
    return path

class SlabTest(unittest.TestCase):

    def setUp(self):
        self.cu = Structure(Lattice.cubic(3), ["Cu", "Cu", "Cu", "Cu"],
                            [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5],
                             [0, 0.5, 0.5]])

        self.libcc = Structure(Lattice.cubic(3.51004), ["Li", "Li"],
                               [[0, 0, 0], [0.5, 0.5, 0.5]])

        Li_Fe_P_O4 = CifParser(get_path("LiFePO4.cif"))
        self.lifepo4 = (Li_Fe_P_O4.get_structures(primitive = False)[0])

        Li_Fe_P_O4_compare = CifParser(get_path("LiFePO4_010_original_3.005.cif"))
        self.lifepo4_compare = (Li_Fe_P_O4_compare.get_structures(primitive = False)[0])


    def test_init(self):
        for hkl in itertools.product(xrange(4), xrange(4), xrange(4)):
            if any(hkl):
                ssize = 6
                vsize = 10
                gen = SurfaceGenerator(self.cu, hkl, ssize, vsize, standardize=False)
                s = gen.get_slab([0])[0]
                if hkl == [0, 1, 1]:
                    self.assertEqual(len(s), 13)
                    self.assertAlmostEqual(s.surface_area, 12.727922061357855)

                manual = self.cu.copy()
                manual.make_supercell(s.scale_factor)
                self.assertEqual(manual.lattice.lengths_and_angles,
                                 s.lattice.lengths_and_angles)

        # # For visual debugging
        # from pymatgen import write_structure
        # write_structure(s.parent, "cu.cif")
        # write_structure(s, "cu_slab_%s_%.3f_%.3f.cif" %
        #                  (str(hkl), ssize, vsize))

    def test_adsorb_atom(self):
        gen = SurfaceGenerator(self.cu,[0, 0, 1], 5, 5, standardize=False)
        s001 = gen.get_slab([0])[0]
        # print s001
        # O adsorb on 4Cu[0.5, 0.5, 0.25], abc = [3, 3, 12]
        # 1. test site_a = abc input
        s001_ad1 = Slab.adsorb_atom(structure_a=s001, site_a=[0.5, 0.5, 0.25], atom= ['O'],
                                    distance=2)
        self.assertEqual(len(s001_ad1), 9)
        for i in xrange(len(s001_ad1)):
            if str(s001_ad1[i].specie) == 'O':
                print s001_ad1[i].frac_coords
                self.assertAlmostEqual(s001_ad1[i].a, 0.5)
                self.assertAlmostEqual(s001_ad1[i].b, 0.5)
                self.assertAlmostEqual(s001_ad1[i].c, 0.4166667)
        self.assertEqual(s001_ad1.lattice.lengths_and_angles,
                                 s001.lattice.lengths_and_angles)
        # 2. test site_a = xyz input
        s001_ad2 = Slab.adsorb_atom(structure_a=s001, site_a=[1.5, 1.5, 3], atom= ['O'],
                                    distance=2, xyz=1)
        self.assertEqual(len(s001_ad2), 9)
        for i in xrange(len(s001_ad2)):
            if str(s001_ad2[i].specie) == 'O':
                print s001_ad2[i].frac_coords
                self.assertAlmostEqual(s001_ad2[i].a, 0.5)
                self.assertAlmostEqual(s001_ad2[i].b, 0.5)
                self.assertAlmostEqual(s001_ad2[i].c, 0.4166667)

    def test_make_terminations(self):
        gen = SurfaceGenerator(self.lifepo4, [0,1,0], 10, 10, lll_reduce=False)
        manual_lifepo4 = gen.get_slab([3.005])[0]
        dist_lifepo4 = gen.get_all_slab()[0]
        check_bonds = gen.get_non_bond_breaking_slab('P5+', 'O2-')[0]
        a = manual_lifepo4.organize_along_c()
        b = dist_lifepo4.organize_along_c()
        lifepo4_compare = Slab([], [], self.lifepo4_compare, [], self.lifepo4_compare, [], [])
        c = lifepo4_compare.organize_along_c()
        d = check_bonds.organize_along_c()

        # Compares the sites of a (010) slab for LiFePO4 generated with the original
        # surface.py algorithm to the same slab generated using manual, check_bonds,
        #  and distance with a shift of 3.005A.

        for i in range(0, len(self.lifepo4_compare)):
            self.assertAlmostEqual(b[i].species_string, c[i].species_string)
            self.assertAlmostEqual(b[i].species_string, c[i].species_string)
            for ii in range(0, 2):
                self.assertAlmostEqual(b.frac_coords[i][ii], c.frac_coords[i][ii])
                self.assertAlmostEqual(a.frac_coords[i][ii], c.frac_coords[i][ii])


    # def test_make_interface(self):
    #     gen = Surface_Generator(self.lifepo4, [0, 0, 1], 10, 10)
    #     slab = gen.distance()[0]
    #     interface = slab.make_interface(self.libcc)
    #
    #     # For visual debugging
    #     # write_structure(interface, "Li_LiFePO4_interface.cif")
    #
    #     for i in range(0, len(slab)):
    #         self.assertEqual(slab[i].coords[0], interface[i].coords[0])
    #         self.assertEqual(slab[i].coords[1], interface[i].coords[1])
    #         self.assertEqual(slab[i].coords[2], interface[i].coords[2])

if __name__ == "__main__":
    unittest.main()
