#!/usr/bin/python


import unittest


from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.core.surface import Slab, SlabGenerator, generate_all_slabs, \
    get_symmetrically_distinct_miller_indices
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
        zno1 = Structure.from_file(get_path("ZnO-wz.cif"), primitive=False)
        zno55 = SlabGenerator(zno1, [1, 0, 0], 5, 5, lll_reduce=False,
                              center_slab=False).get_slab()
        #print zno55[0]
        #write_structure(zno55[0], "surface_tests/ZnO55.cif")
        self.zno1 = zno1
        self.zno55 = zno55
        self.h = Structure(Lattice.cubic(3), ["H"],
                            [[0, 0, 0]])
        self.libcc = Structure(Lattice.cubic(3.51004), ["Li", "Li"],
                               [[0, 0, 0], [0.5, 0.5, 0.5]])

    def test_init(self):
        zno_slab = Slab(self.zno55.lattice, self.zno55.species, self.zno55.frac_coords,
                        self.zno55.miller_index, self.zno55.oriented_unit_cell, 0,
                        self.zno55.scale_factor)
        m =self.zno55.lattice.matrix
        area = np.linalg.norm(np.cross(m[0], m[1]))
        self.assertAlmostEqual(zno_slab.surface_area, area)
        self.assertEqual(zno_slab.lattice.lengths_and_angles,
                         self.zno55.lattice.lengths_and_angles)
        self.assertEqual(zno_slab.oriented_unit_cell.composition,
                         self.zno1.composition)
        self.assertEqual(len(zno_slab), 8)

    def test_add_adsorbate_atom(self):
        zno_slab = Slab(self.zno55.lattice, self.zno55.species, self.zno55.frac_coords,
                        self.zno55.miller_index, self.zno55.oriented_unit_cell, 0,
                        self.zno55.scale_factor)
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


class SlabGeneratorTest(PymatgenTest):
    def test_get_slab(self):
        s = self.get_structure("LiFePO4")
        gen = SlabGenerator(s, [0, 0, 1], 10, 10)
        s = gen.get_slab(0.25)
        self.assertAlmostEqual(s.lattice.abc[2], 20.820740000000001)

    def test_get_slabs(self):
        gen = SlabGenerator(self.get_structure("CsCl"), [0, 0, 1], 10, 10)
        self.assertEqual(len(gen.get_slabs()), 2)
        s = self.get_structure("LiFePO4")
        gen = SlabGenerator(s, [0, 0, 1], 10, 10)
        self.assertEqual(len(gen.get_slabs()), 18)

        self.assertEqual(len(gen.get_slabs(bonds={("P", "O"): 3})), 6)

        # There are no slabs in LFP that does not break either P-O or Fe-O
        # bonds for a miller index of [0, 0, 1].
        self.assertEqual(len(gen.get_slabs(bonds={("P", "O"): 3,
                                                  ("Fe", "O"): 3})), 0)

        # At this threshold, only the origin and center Li results in clustering.
        # All other sites are non-clustered. So the # of slabs is # of sites
        # in LiFePO4 unit cell - 2.
        self.assertEqual(len(gen.get_slabs(tol=1e-4)), 26)

        self.LiCoO2 = Structure.from_file(get_path("icsd_LiCoO2.cif"), primitive=False)
        licoo2 = SlabGenerator(self.LiCoO2, [0,0,1], 10, 10)
        lico = licoo2.get_slabs(bonds={("Co", "O"): 3})
        self.assertEqual(len(lico), 3)

    def test_triclinic_TeI(self):
        # Test case for a triclinic structure of TeI. Only these three
        # Miller indices are used because it is easier to identify which
        # atoms should be in a surface together. The closeness of the sites
        # in other Miller indices can cause some ambiguity when choosing a
        # higher tolerance.
        mill = [[0, 0, 1],[0, 1, 0],[1, 0, 0]]
        numb_slabs = {'[0, 0, 1]':6, '[0, 1, 0]':3, '[1, 0, 0]':8}
        TeI = Structure.from_file(get_path("icsd_TeI.cif"), primitive=False)
        for i in mill:
            trclnc_TeI = SlabGenerator(TeI, i, 10, 10)
            TeI_slabs = trclnc_TeI.get_slabs(tol=0.05)
            self.assertEqual(numb_slabs[str(i)], len(TeI_slabs))


class FuncTest(PymatgenTest):

    def setUp(self):
        self.cscl = self.get_structure("CsCl")
        self.lifepo4 = self.get_structure("LiFePO4")
        self.tei = Structure.from_file(get_path("icsd_TeI.cif"),
                                       primitive=False)
        self.LiCoO2 = Structure.from_file(get_path("icsd_LiCoO2.cif"), primitive=False)

        self.p1 = Structure(Lattice.from_parameters(3, 4, 5, 31, 43, 50),
                            ["H", "He"], [[0, 0, 0], [0.1, 0.2, 0.3]])

    def test_get_symmetrically_distinct_miller_indices(self):
        indices = get_symmetrically_distinct_miller_indices(self.cscl, 1)
        self.assertEqual(len(indices), 3)
        indices = get_symmetrically_distinct_miller_indices(self.cscl, 2)
        self.assertEqual(len(indices), 6)

        self.assertEqual(len(get_symmetrically_distinct_miller_indices(
            self.lifepo4, 1)), 7)

        # The TeI P-1 structure should have 13 unique millers (only inversion
        # symmetry eliminates pairs)
        indices = get_symmetrically_distinct_miller_indices(self.tei, 1)
        self.assertEqual(len(indices), 13)

        # P1 and P-1 should have the same # of miller indices since surfaces
        # always have inversion symmetry.
        indices = get_symmetrically_distinct_miller_indices(self.p1, 1)
        self.assertEqual(len(indices), 13)

    def test_generate_all_slabs(self):
        slabs = generate_all_slabs(self.cscl, 1, 10, 10)
        #Only three possible slabs, one each in (100), (110) and (111).
        self.assertEqual(len(slabs), 3)

        slabs1 = generate_all_slabs(self.lifepo4, 1, 10, 10,
                                   bonds={("P", "O"): 3})
        self.assertEqual(len(slabs1), 5)

        slabs2 = generate_all_slabs(self.lifepo4, 1, 10, 10,
                                   bonds={("P", "O"): 3,
                                          ("Fe", "O"): 3})
        self.assertEqual(len(slabs2), 0)

        # There should be only three possible stable surfaces, all of which are
        # in the (001) oriented unit cell
        slabs3 = generate_all_slabs(self.LiCoO2, 1, 10, 10,
                           bonds={("Co", "O"): 3})
        self.assertEqual(len(slabs3), 3)
        mill = [0, 0, 1]
        for lico in slabs3:
            for i in xrange(3):
                self.assertEqual(lico.miller_index[i], mill[i])


if __name__ == "__main__":
    unittest.main()
