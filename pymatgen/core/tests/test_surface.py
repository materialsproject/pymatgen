#!/usr/bin/python


import unittest2 as unittest
import os
import random

import numpy as np

from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.core.surface import Slab, SlabGenerator, generate_all_slabs, \
    get_symmetrically_distinct_miller_indices
from pymatgen.symmetry.groups import SpaceGroup
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
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

        Ti = Structure(Lattice.hexagonal(4.6, 2.82), ["Ti", "Ti", "Ti"],
                       [[0.000000, 0.000000, 0.000000],
                       [0.333333, 0.666667, 0.500000],
                       [0.666667, 0.333333, 0.500000]])

        Ag_fcc = Structure(Lattice.cubic(4.06), ["Ag", "Ag", "Ag", "Ag"],
                           [[0.000000, 0.000000, 0.000000],
                           [0.000000, 0.500000, 0.500000],
                           [0.500000, 0.000000, 0.500000],
                           [0.500000, 0.500000, 0.000000]])

        laue_groups = ["-1", "2/m", "mmm", "4/m",
                       "4/mmm", "-3", "-3m", "6/m",
                       "6/mmm", "m-3", "m-3m"]

        self.ti = Ti
        self.agfcc = Ag_fcc
        self.zno1 = zno1
        self.zno55 = zno55
        self.h = Structure(Lattice.cubic(3), ["H"],
                            [[0, 0, 0]])
        self.libcc = Structure(Lattice.cubic(3.51004), ["Li", "Li"],
                               [[0, 0, 0], [0.5, 0.5, 0.5]])
        self.laue_groups = laue_groups

    def test_init(self):
        zno_slab = Slab(self.zno55.lattice, self.zno55.species,
                        self.zno55.frac_coords,
                        self.zno55.miller_index,
                        self.zno55.oriented_unit_cell,
                        0, self.zno55.scale_factor)
        m =self.zno55.lattice.matrix
        area = np.linalg.norm(np.cross(m[0], m[1]))
        self.assertAlmostEqual(zno_slab.surface_area, area)
        self.assertEqual(zno_slab.lattice.lengths_and_angles,
                         self.zno55.lattice.lengths_and_angles)
        self.assertEqual(zno_slab.oriented_unit_cell.composition,
                         self.zno1.composition)
        self.assertEqual(len(zno_slab), 8)

    def test_add_adsorbate_atom(self):
        zno_slab = Slab(self.zno55.lattice, self.zno55.species,
                        self.zno55.frac_coords,
                        self.zno55.miller_index,
                        self.zno55.oriented_unit_cell,
                        0, self.zno55.scale_factor)
        zno_slab.add_adsorbate_atom([1], 'H', 1)

        self.assertEqual(len(zno_slab), 9)
        self.assertEqual(str(zno_slab[8].specie), 'H')
        self.assertAlmostEqual(zno_slab.get_distance(1, 8), 1.0)
        self.assertTrue(zno_slab[8].c > zno_slab[0].c)
        m = self.zno55.lattice.matrix
        area = np.linalg.norm(np.cross(m[0], m[1]))
        self.assertAlmostEqual(zno_slab.surface_area, area)
        self.assertEqual(zno_slab.lattice.lengths_and_angles,
                         self.zno55.lattice.lengths_and_angles)

    def test_get_sorted_structure(self):
        species = [str(site.specie) for site in
                   self.zno55.get_sorted_structure()]
        self.assertEqual(species, ["Zn2+"] * 4 + ["O2-"] * 4)

    def test_methods(self):
        #Test various structure methods
        self.zno55.get_primitive_structure()

    def test_as_from_dict(self):
        d = self.zno55.as_dict()
        obj = Slab.from_dict(d)
        self.assertEqual(obj.miller_index, (1, 0, 0))

    def test_dipole_and_is_polar(self):
        self.assertArrayAlmostEqual(self.zno55.dipole, [0, 0, 0])
        self.assertFalse(self.zno55.is_polar())
        cscl = self.get_structure("CsCl")
        cscl.add_oxidation_state_by_element({"Cs": 1, "Cl": -1})
        slab = SlabGenerator(cscl, [1, 0, 0], 5, 5,
                             lll_reduce=False, center_slab=False).get_slab()
        self.assertArrayAlmostEqual(slab.dipole, [-4.209, 0, 0])
        self.assertTrue(slab.is_polar())

    def test_symmetrization(self):

        # Restricted to elemental materials due to the risk of
        # broken stoichiometry. For compound materials, use is_polar()

        # Get all slabs for P6/mmm Ti and Fm-3m Ag up to index of 2

        all_Ti_slabs = generate_all_slabs(self.ti, 2, 10, 10, bonds=None,
                                          tol=1e-3, max_broken_bonds=0,
                                          lll_reduce=False, center_slab=False,
                                          primitive=True, max_normal_search=2,
                                          symmetrize=True)

        all_Ag_fcc_slabs = generate_all_slabs(self.agfcc, 2, 10, 10, bonds=None,
                                              tol=1e-3, max_broken_bonds=0,
                                              lll_reduce=False, center_slab=False,
                                              primitive=True, max_normal_search=2,
                                              symmetrize=True)

        all_slabs = [all_Ti_slabs, all_Ag_fcc_slabs]

        for i, slabs in enumerate(all_slabs):

            assymetric_count = 0
            symmetric_count = 0

            for i, slab in enumerate(slabs):
                sg = SpacegroupAnalyzer(slab)
                pg = sg.get_point_group_symbol()

                # Check if a slab is symmetric
                if str(pg) not in self.laue_groups:
                    assymetric_count += 1
                else:
                    symmetric_count += 1

            # Check if slabs are all symmetric
            self.assertEqual(assymetric_count, 0)
            self.assertEqual(symmetric_count, len(slabs))

class SlabGeneratorTest(PymatgenTest):

    def test_get_slab(self):
        s = self.get_structure("LiFePO4")
        gen = SlabGenerator(s, [0, 0, 1], 10, 10)
        s = gen.get_slab(0.25)
        self.assertAlmostEqual(s.lattice.abc[2], 20.820740000000001)

        fcc = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3), ["Fe"],
                                        [[0, 0, 0]])
        gen = SlabGenerator(fcc, [1, 1, 1], 10, 10)
        slab = gen.get_slab()
        gen = SlabGenerator(fcc, [1, 1, 1], 10, 10, primitive=False)
        slab_non_prim = gen.get_slab()
        self.assertEqual(len(slab), 6)
        self.assertEqual(len(slab_non_prim), len(slab) * 4)

        #Some randomized testing of cell vectors
        for i in range(1, 231):
            i = random.randint(1, 230)
            sg = SpaceGroup.from_int_number(i)
            if sg.crystal_system == "hexagonal" or (sg.crystal_system == \
                    "trigonal" and sg.symbol.endswith("H")):
                latt = Lattice.hexagonal(5, 10)
            else:
                #Cubic lattice is compatible with all other space groups.
                latt = Lattice.cubic(5)
            s = Structure.from_spacegroup(i, latt, ["H"], [[0, 0, 0]])
            miller = (0, 0, 0)
            while miller == (0, 0, 0):
                miller = (random.randint(0, 6), random.randint(0, 6),
                          random.randint(0, 6))
            gen = SlabGenerator(s, miller, 10, 10)
            a, b, c = gen.oriented_unit_cell.lattice.matrix
            self.assertAlmostEqual(np.dot(a, gen._normal), 0)
            self.assertAlmostEqual(np.dot(b, gen._normal), 0)

    def test_normal_search(self):
        fcc = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3), ["Fe"],
                                        [[0, 0, 0]])
        for miller in [(1, 0, 0), (1, 1, 0), (1, 1, 1), (2, 1, 1)]:
            gen = SlabGenerator(fcc, miller, 10, 10)
            gen_normal = SlabGenerator(fcc, miller, 10, 10,
                                       max_normal_search=max(miller))
            slab = gen_normal.get_slab()
            self.assertAlmostEqual(slab.lattice.alpha, 90)
            self.assertAlmostEqual(slab.lattice.beta, 90)
            self.assertGreaterEqual(len(gen_normal.oriented_unit_cell),
                                    len(gen.oriented_unit_cell))

        graphite = self.get_structure("Graphite")
        for miller in [(1, 0, 0), (1, 1, 0), (0, 0, 1), (2, 1, 1)]:
            gen = SlabGenerator(graphite, miller, 10, 10)
            gen_normal = SlabGenerator(graphite, miller, 10, 10,
                                       max_normal_search=max(miller))
            self.assertGreaterEqual(len(gen_normal.oriented_unit_cell),
                                    len(gen.oriented_unit_cell))

        sc = Structure(Lattice.hexagonal(3.32, 5.15), ["Sc", "Sc"],
                       [[1/3, 2/3, 0.25], [2/3, 1/3, 0.75]])
        gen = SlabGenerator(sc, (1, 1, 1), 10, 10, max_normal_search=1)
        self.assertAlmostEqual(gen.oriented_unit_cell.lattice.angles[1], 90)

    def test_get_slabs(self):
        gen = SlabGenerator(self.get_structure("CsCl"), [0, 0, 1], 10, 10)

        #Test orthogonality of some internal variables.
        a, b, c = gen.oriented_unit_cell.lattice.matrix
        self.assertAlmostEqual(np.dot(a, gen._normal), 0)
        self.assertAlmostEqual(np.dot(b, gen._normal), 0)

        self.assertEqual(len(gen.get_slabs()), 1)

        s = self.get_structure("LiFePO4")
        gen = SlabGenerator(s, [0, 0, 1], 10, 10)
        self.assertEqual(len(gen.get_slabs()), 5)

        self.assertEqual(len(gen.get_slabs(bonds={("P", "O"): 3})), 2)

        # There are no slabs in LFP that does not break either P-O or Fe-O
        # bonds for a miller index of [0, 0, 1].
        self.assertEqual(len(gen.get_slabs(
            bonds={("P", "O"): 3, ("Fe", "O"): 3})), 0)

        #If we allow some broken bonds, there are a few slabs.
        self.assertEqual(len(gen.get_slabs(
            bonds={("P", "O"): 3, ("Fe", "O"): 3},
            max_broken_bonds=2)), 2)

        # At this threshold, only the origin and center Li results in
        # clustering. All other sites are non-clustered. So the of
        # slabs is of sites in LiFePO4 unit cell - 2 + 1.
        self.assertEqual(len(gen.get_slabs(tol=1e-4)), 15)

        LiCoO2=Structure.from_file(get_path("icsd_LiCoO2.cif"),
                                          primitive=False)
        gen = SlabGenerator(LiCoO2, [0, 0, 1], 10, 10)
        lco = gen.get_slabs(bonds={("Co", "O"): 3})
        self.assertEqual(len(lco), 1)
        a, b, c = gen.oriented_unit_cell.lattice.matrix
        self.assertAlmostEqual(np.dot(a, gen._normal), 0)
        self.assertAlmostEqual(np.dot(b, gen._normal), 0)

        scc = Structure.from_spacegroup("Pm-3m", Lattice.cubic(3), ["Fe"],
                                        [[0, 0, 0]])
        gen = SlabGenerator(scc, [0, 0, 1], 10, 10)
        slabs = gen.get_slabs()
        self.assertEqual(len(slabs), 1)
        gen = SlabGenerator(scc, [1, 1, 1], 10, 10, max_normal_search=1)
        slabs = gen.get_slabs()
        self.assertEqual(len(slabs), 1)

    def test_triclinic_TeI(self):
        # Test case for a triclinic structure of TeI. Only these three
        # Miller indices are used because it is easier to identify which
        # atoms should be in a surface together. The closeness of the sites
        # in other Miller indices can cause some ambiguity when choosing a
        # higher tolerance.
        numb_slabs = {(0, 0, 1): 5, (0, 1, 0): 3, (1, 0, 0): 7}
        TeI = Structure.from_file(get_path("icsd_TeI.cif"),
                                  primitive=False)
        for k, v in numb_slabs.items():
            trclnc_TeI = SlabGenerator(TeI, k, 10, 10)
            TeI_slabs = trclnc_TeI.get_slabs()
            self.assertEqual(v, len(TeI_slabs))

    def test_get_orthogonal_c_slab(self):
        TeI = Structure.from_file(get_path("icsd_TeI.cif"),
                                  primitive=False)
        trclnc_TeI = SlabGenerator(TeI, (0, 0, 1), 10, 10)
        TeI_slabs = trclnc_TeI.get_slabs()
        slab = TeI_slabs[0]
        norm_slab = slab.get_orthogonal_c_slab()
        self.assertAlmostEqual(norm_slab.lattice.angles[0], 90)
        self.assertAlmostEqual(norm_slab.lattice.angles[1], 90)


class MillerIndexFinderTests(PymatgenTest):

    def setUp(self):
        self.cscl = Structure.from_spacegroup(
            "Pm-3m", Lattice.cubic(4.2), ["Cs", "Cl"],
            [[0, 0, 0], [0.5, 0.5, 0.5]])

        self.lifepo4 = self.get_structure("LiFePO4")
        self.tei = Structure.from_file(get_path("icsd_TeI.cif"),
                                       primitive=False)
        self.LiCoO2 = Structure.from_file(get_path("icsd_LiCoO2.cif"),
                                          primitive=False)

        self.p1 = Structure(Lattice.from_parameters(3, 4, 5, 31, 43, 50),
                            ["H", "He"], [[0, 0, 0], [0.1, 0.2, 0.3]])
        self.graphite = self.get_structure("Graphite")

    def test_get_symmetrically_distinct_miller_indices(self):

        # Tests to see if the function obtains the known number of unique slabs

        indices = get_symmetrically_distinct_miller_indices(self.cscl, 1)
        self.assertEqual(len(indices), 3)
        indices = get_symmetrically_distinct_miller_indices(self.cscl, 2)
        self.assertEqual(len(indices), 6)

        self.assertEqual(len(get_symmetrically_distinct_miller_indices(self.lifepo4, 1)), 7)

        # The TeI P-1 structure should have 13 unique millers (only inversion
        # symmetry eliminates pairs)
        indices = get_symmetrically_distinct_miller_indices(self.tei, 1)
        self.assertEqual(len(indices), 13)

        # P1 and P-1 should have the same # of miller indices since surfaces
        # always have inversion symmetry.
        indices = get_symmetrically_distinct_miller_indices(self.p1, 1)
        self.assertEqual(len(indices), 13)

        indices = get_symmetrically_distinct_miller_indices(self.graphite, 2)
        self.assertEqual(len(indices), 12)

    def test_generate_all_slabs(self):

        slabs = generate_all_slabs(self.cscl, 1, 10, 10)

        # Only three possible slabs, one each in (100), (110) and (111).
        self.assertEqual(len(slabs), 3)

        slabs = generate_all_slabs(self.cscl, 1, 10, 10,
                                   bonds={("Cs", "Cl"): 4})
        # No slabs if we don't allow broken Cs-Cl
        self.assertEqual(len(slabs), 0)

        slabs = generate_all_slabs(self.cscl, 1, 10, 10,
                                   bonds={("Cs", "Cl"): 4},
                                   max_broken_bonds=100)
        self.assertEqual(len(slabs), 3)

        slabs1 = generate_all_slabs(self.lifepo4, 1, 10, 10, tol=0.1,
                                    bonds={("P", "O"): 3})
        self.assertEqual(len(slabs1), 4)

        slabs2 = generate_all_slabs(self.lifepo4, 1, 10, 10,
                                    bonds={("P", "O"): 3, ("Fe", "O"): 3})
        self.assertEqual(len(slabs2), 0)

        # There should be only one possible stable surfaces, all of which are
        # in the (001) oriented unit cell
        slabs3 = generate_all_slabs(self.LiCoO2, 1, 10, 10,
                                    bonds={("Co", "O"): 3})
        self.assertEqual(len(slabs3), 1)
        mill = (0, 0, 1)
        for s in slabs3:
            self.assertEqual(s.miller_index, mill)

if __name__ == "__main__":
    unittest.main()
