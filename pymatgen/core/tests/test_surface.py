#!/usr/bin/python


import unittest


from pymatgen.core import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.io.smartio import CifParser
from pymatgen.core.surface import Slab, SurfaceGenerator
import os
import numpy as np

from pymatgen.util.testing import PymatgenTest


def get_path(path_str):
    cwd = os.path.abspath(os.path.dirname(__file__))
    path = os.path.join(cwd, "..", "..", "..", "test_files", "surface_tests",
                        path_str)
    return path


class SlabTest(PymatgenTest):

    def setUp(self):
        C1 = CifParser(get_path("ZnO-wz.cif"))
        zno = C1.get_structures(primitive=False)
        zno1 = zno[0]
        zno55 = SurfaceGenerator(zno1, [1, 0, 0], 5, 5, lll_reduce=False,
                                 standardize=False).get_slab()
        #print zno55[0]
        #write_structure(zno55[0], "surface_tests/ZnO55.cif")
        self.zno1 = zno1
        self.zno55 = zno55
        self.h = Structure(Lattice.cubic(3), ["H"],
                            [[0, 0, 0]])
        self.libcc = Structure(Lattice.cubic(3.51004), ["Li", "Li"],
                               [[0, 0, 0], [0.5, 0.5, 0.5]])

    def test_init(self):
        zno_slab = Slab(self.zno55, self.zno55.miller_index, 0,
                        self.zno55.normal, self.zno55.scale_factor,
                        self.zno1, self.zno55.min_slab_size,
                        self.zno55.min_vac_size)
        m =self.zno55.lattice.matrix
        area = np.linalg.norm(np.cross(m[0], m[1]))
        self.assertAlmostEqual(zno_slab.surface_area, area)
        self.assertEqual(zno_slab.lattice.lengths_and_angles,
                         self.zno55.lattice.lengths_and_angles)
        self.assertEqual(zno_slab.normal.all, self.zno55.normal.all)
        self.assertEqual(zno_slab.parent.composition, self.zno1.composition)
        self.assertEqual(len(zno_slab), 8)

    def test_add_adsorbate_atom(self):
        zno_slab = Slab(self.zno55, self.zno55.miller_index,
                        0, self.zno55.normal, self.zno55.scale_factor,
                        self.zno1, self.zno55.min_slab_size, self.zno55.min_vac_size)
        zno_slab.add_adsorbate_atom([1], 'H', 1)

        self.assertEqual(len(zno_slab), 9)
        self.assertEqual(str(zno_slab[8].specie),'H')
        self.assertAlmostEqual(zno_slab.get_distance(1, 8), 1.0)
        self.assertTrue(zno_slab[8].c > zno_slab[0].c)
        m =self.zno55.lattice.matrix
        area = np.linalg.norm(np.cross(m[0], m[1]))
        self.assertAlmostEqual(zno_slab.surface_area, area)
        self.assertEqual(zno_slab.lattice.lengths_and_angles,
                         self.zno55.lattice.lengths_and_angles)
        self.assertEqual(zno_slab.normal.all, self.zno55.normal.all)


class SurfaceGeneratorTest(PymatgenTest):

    def setup(self):
        Li_Fe_P_O4 = CifParser(get_path("LiFePO4.cif"))
        self.lifepo4 = (Li_Fe_P_O4.get_structures(primitive = False)[0])
        Li_Fe_P_O4_compare = CifParser(get_path("LiFePO4_010_original_3.005.cif"))
        self.lifepo4_compare = (Li_Fe_P_O4_compare.get_structures(primitive = False)[0])

    def test_get_slab(self):
        s = self.get_structure("LiFePO4")
        gen = SurfaceGenerator(s, [0, 0, 1], 10, 10)
        s = gen.get_slab(0.25)
        self.assertAlmostEqual(s.lattice.abc[2], 20.820740000000001)

    def test_get_all_slabs(self):
        gen = SurfaceGenerator(self.get_structure("CsCl"), [0, 0, 1], 10, 10)
        self.assertEqual(len(gen.get_slabs()), 2)
        s = self.get_structure("LiFePO4")
        gen = SurfaceGenerator(s, [0, 0, 1], 10, 10)
        self.assertEqual(len(gen.get_slabs()), 18)

        self.assertEqual(len(gen.get_slabs(bonds={("P", "O"): 3})), 6)

        # At this threshold, only the origin and center Li results in clustering.
        # All other sites are non-clustered. So the # of slabs is # of sites
        # in LiFePO4 unit cell - 2.
        self.assertEqual(len(gen.get_slabs(tol=1e-4)), 26)

    def test_CsCl(self):
        # Checks the species and coordinates in every site in three
        # equivalent orientations of CsCl to see if they are the same
        miller_indices = [[0,0,1], [0,1,0], [1,0,0]]
        config1, config2 = [], []
        for index in miller_indices:
            cscl = SurfaceGenerator(self.get_structure("CsCl"), index, 10, 10)
            all_cscl = cscl.get_slabs()
            self.assertEqual(len(all_cscl), 2)
            config1.append(all_cscl[0])
            config2.append(all_cscl[1])
        for i in range(0,1):
            config1[2].organize_along_c
            config1[i].organize_along_c
            config2[2].organize_along_c
            config2[i].organize_along_c
            for ii in xrange(len(config1[2])):
                self.assertEqual(config1[i][ii].species_string,
                                 config1[2][ii].species_string)
                for iii in range(0, 3):
                    self.assertEqual(config1[i].frac_coords[ii][iii],
                                     config1[2].frac_coords[ii][iii])
            for ii in xrange(len(config2[2])):

                self.assertEqual(config1[i][ii].species_string,
                                 config1[2][ii].species_string)
                for iii in range(0, 3):
                    self.assertEqual(config2[i].frac_coords[ii][iii],
                                     config2[2].frac_coords[ii][iii])







if __name__ == "__main__":
    unittest.main()
