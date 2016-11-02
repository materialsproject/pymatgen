# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import unittest2 as unittest
import os
import json
import numpy as np

from pymatgen.analysis.structure_matcher import StructureMatcher, \
    ElementComparator, FrameworkComparator, OrderDisorderElementComparator
from monty.json import MontyDecoder
from pymatgen.core.operations import SymmOp
from pymatgen import Structure, Element, Lattice
from pymatgen.util.coord_utils import find_in_coord_list_pbc
from pymatgen.util.testing import PymatgenTest

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class StructureMatcherTest(PymatgenTest):

    def setUp(self):
        with open(os.path.join(test_dir, "TiO2_entries.json"), 'r') as fp:
            entries = json.load(fp, cls=MontyDecoder)
        self.struct_list = [e.structure for e in entries]
        self.oxi_structs = [self.get_structure("Li2O"),
                            Structure.from_file(os.path.join(
                                test_dir, "POSCAR.Li2O"))]

    def test_ignore_species(self):
        s1 = Structure.from_file(os.path.join(test_dir, "LiFePO4.cif"))
        s2 = Structure.from_file(os.path.join(test_dir, "POSCAR"))
        m = StructureMatcher(ignored_species=["Li"], primitive_cell=False,
                             attempt_supercell=True)
        self.assertTrue(m.fit(s1, s2))
        self.assertTrue(m.fit_anonymous(s1, s2))
        groups = m.group_structures([s1, s2])
        self.assertEqual(len(groups), 1)
        s2.make_supercell((2, 1, 1))
        ss1 = m.get_s2_like_s1(s2, s1, include_ignored_species=True)
        self.assertAlmostEqual(ss1.lattice.a, 20.820740000000001)
        self.assertEqual(ss1.composition.reduced_formula, "LiFePO4")

        self.assertEqual({
            k.symbol: v.symbol for k, v in
            m.get_best_electronegativity_anonymous_mapping(s1, s2).items()},
                         {"Fe": "Fe", "P": "P", "O": "O"})



    def test_get_supercell_size(self):
        l = Lattice.cubic(1)
        l2 = Lattice.cubic(0.9)
        s1 = Structure(l, ['Mg', 'Cu', 'Ag', 'Cu', 'Ag'], [[0]*3]*5)
        s2 = Structure(l2, ['Cu', 'Cu', 'Ag'], [[0]*3]*3)

        sm = StructureMatcher(supercell_size='volume')
        self.assertEqual(sm._get_supercell_size(s1, s2),
                         (1, True))
        self.assertEqual(sm._get_supercell_size(s2, s1),
                         (1, True))

        sm = StructureMatcher(supercell_size='num_sites')
        self.assertEqual(sm._get_supercell_size(s1, s2),
                         (2, False))
        self.assertEqual(sm._get_supercell_size(s2, s1),
                         (2, True))

        sm = StructureMatcher(supercell_size='Ag')
        self.assertEqual(sm._get_supercell_size(s1, s2),
                         (2, False))
        self.assertEqual(sm._get_supercell_size(s2, s1),
                         (2, True))

        sm = StructureMatcher(supercell_size='wfieoh')
        self.assertRaises(ValueError, sm._get_supercell_size, s1, s2)

    def test_cmp_fstruct(self):
        sm = StructureMatcher()

        s1 = np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]])
        s2 = np.array([[0.11, 0.22, 0.33]])
        frac_tol = np.array([0.02, 0.03, 0.04])
        mask = np.array([[False, False]])
        mask2 = np.array([[True, False]])

        self.assertRaises(ValueError, sm._cmp_fstruct, s2, s1, frac_tol, mask.T)
        self.assertRaises(ValueError, sm._cmp_fstruct, s1, s2, frac_tol, mask.T)

        self.assertTrue(sm._cmp_fstruct(s1, s2, frac_tol, mask))
        self.assertFalse(sm._cmp_fstruct(s1, s2, frac_tol/2, mask))
        self.assertFalse(sm._cmp_fstruct(s1, s2, frac_tol, mask2))

    def test_cart_dists(self):
        sm = StructureMatcher()
        l = Lattice.orthorhombic(1, 2, 3)

        s1 = np.array([[0.13, 0.25, 0.37], [0.1, 0.2, 0.3]])
        s2 = np.array([[0.11, 0.22, 0.33]])
        s3 = np.array([[0.1, 0.2, 0.3], [0.11, 0.2, 0.3]])
        s4 = np.array([[0.1, 0.2, 0.3], [0.1, 0.6, 0.7]])
        mask = np.array([[False, False]])
        mask2 = np.array([[False, True]])
        mask3 = np.array([[False, False], [False, False]])
        mask4 = np.array([[False, True], [False, True]])

        n1 = (len(s1) / l.volume) ** (1/3)
        n2 = (len(s2) / l.volume) ** (1/3)

        self.assertRaises(ValueError, sm._cart_dists, s2, s1, l, mask.T, n2)
        self.assertRaises(ValueError, sm._cart_dists, s1, s2, l, mask.T, n1)

        d, ft, s = sm._cart_dists(s1, s2, l, mask, n1)
        self.assertTrue(np.allclose(d, [0]))
        self.assertTrue(np.allclose(ft, [-0.01, -0.02, -0.03]))
        self.assertTrue(np.allclose(s, [1]))

        #check that masking best value works
        d, ft, s = sm._cart_dists(s1, s2, l, mask2, n1)
        self.assertTrue(np.allclose(d, [0]))
        self.assertTrue(np.allclose(ft, [0.02, 0.03, 0.04]))
        self.assertTrue(np.allclose(s, [0]))

        #check that averaging of translation is done properly
        d, ft, s = sm._cart_dists(s1, s3, l, mask3, n1)
        self.assertTrue(np.allclose(d, [0.08093341]*2))
        self.assertTrue(np.allclose(ft, [0.01, 0.025, 0.035]))
        self.assertTrue(np.allclose(s, [1, 0]))

        #check distances are large when mask allows no 'real' mapping
        d, ft, s = sm._cart_dists(s1, s4, l, mask4, n1)
        self.assertTrue(np.min(d) > 1e8)
        self.assertTrue(np.min(ft) > 1e8)

    def test_get_mask(self):
        sm = StructureMatcher(comparator=ElementComparator())
        l = Lattice.cubic(1)
        s1 = Structure(l, ['Mg', 'Cu', 'Ag', 'Cu'], [[0]*3]*4)
        s2 = Structure(l, ['Cu', 'Cu', 'Ag'], [[0]*3]*3)

        result = [[True, False, True, False],
                  [True, False, True, False],
                  [True, True, False, True]]
        m, inds, i = sm._get_mask(s1, s2, 1, True)
        self.assertTrue(np.all(m == result))
        self.assertTrue(i == 2)
        self.assertEqual(inds, [2])

        #test supercell with match
        result = [[1, 1, 0, 0, 1, 1, 0, 0],
                  [1, 1, 0, 0, 1, 1, 0, 0],
                  [1, 1, 1, 1, 0, 0, 1, 1]]
        m, inds, i = sm._get_mask(s1, s2, 2, True)
        self.assertTrue(np.all(m == result))
        self.assertTrue(i == 2)
        self.assertTrue(np.allclose(inds, np.array([4])))

        #test supercell without match
        result = [[1, 1, 1, 1, 1, 1],
                  [0, 0, 0, 0, 1, 1],
                  [1, 1, 1, 1, 0, 0],
                  [0, 0, 0, 0, 1, 1]]
        m, inds, i = sm._get_mask(s2, s1, 2, True)
        self.assertTrue(np.all(m == result))
        self.assertTrue(i == 0)
        self.assertTrue(np.allclose(inds, np.array([])))

        #test s2_supercell
        result = [[1, 1, 1], [1, 1, 1],
                  [0, 0, 1], [0, 0, 1],
                  [1, 1, 0], [1, 1, 0],
                  [0, 0, 1], [0, 0, 1]]
        m, inds, i = sm._get_mask(s2, s1, 2, False)
        self.assertTrue(np.all(m == result))
        self.assertTrue(i == 0)
        self.assertTrue(np.allclose(inds, np.array([])))

        #test for multiple translation indices
        s1 = Structure(l, ['Cu', 'Ag', 'Cu', 'Ag', 'Ag'], [[0]*3]*5)
        s2 = Structure(l, ['Ag', 'Cu', 'Ag'], [[0]*3]*3)
        result = [[1, 0, 1, 0, 0],
                  [0, 1, 0, 1, 1],
                  [1, 0, 1, 0, 0]]
        m, inds, i = sm._get_mask(s1, s2, 1, True)

        self.assertTrue(np.all(m == result))
        self.assertTrue(i == 1)
        self.assertTrue(np.allclose(inds, [0, 2]))

    def test_get_supercells(self):
        sm = StructureMatcher(comparator=ElementComparator())
        l = Lattice.cubic(1)
        l2 = Lattice.cubic(0.5)
        s1 = Structure(l, ['Mg', 'Cu', 'Ag', 'Cu'], [[0]*3]*4)
        s2 = Structure(l2, ['Cu', 'Cu', 'Ag'], [[0]*3]*3)
        scs = list(sm._get_supercells(s1, s2, 8, False))
        for x in scs:
            self.assertAlmostEqual(abs(np.linalg.det(x[3])), 8)
            self.assertEqual(len(x[0]), 4)
            self.assertEqual(len(x[1]), 24)
        self.assertEqual(len(scs), 48)

        scs = list(sm._get_supercells(s2, s1, 8, True))
        for x in scs:
            self.assertAlmostEqual(abs(np.linalg.det(x[3])), 8)
            self.assertEqual(len(x[0]), 24)
            self.assertEqual(len(x[1]), 4)
        self.assertEqual(len(scs), 48)

    def test_fit(self):
        """
        Take two known matched structures
            1) Ensure match
            2) Ensure match after translation and rotations
            3) Ensure no-match after large site translation
            4) Ensure match after site shuffling
            """
        sm = StructureMatcher()

        self.assertTrue(sm.fit(self.struct_list[0], self.struct_list[1]))

        # Test rotational/translational invariance
        op = SymmOp.from_axis_angle_and_translation([0, 0, 1], 30, False,
                                                    np.array([0.4, 0.7, 0.9]))
        self.struct_list[1].apply_operation(op)
        self.assertTrue(sm.fit(self.struct_list[0], self.struct_list[1]))

        #Test failure under large atomic translation
        self.struct_list[1].translate_sites([0], [.4, .4, .2],
                                            frac_coords=True)
        self.assertFalse(sm.fit(self.struct_list[0], self.struct_list[1]))

        self.struct_list[1].translate_sites([0], [-.4, -.4, -.2],
                                            frac_coords=True)
        # random.shuffle(editor._sites)
        self.assertTrue(sm.fit(self.struct_list[0], self.struct_list[1]))
        #Test FrameworkComporator
        sm2 = StructureMatcher(comparator=FrameworkComparator())
        lfp = self.get_structure("LiFePO4")
        nfp = self.get_structure("NaFePO4")
        self.assertTrue(sm2.fit(lfp, nfp))
        self.assertFalse(sm.fit(lfp, nfp))

        #Test anonymous fit.
        self.assertEqual(sm.fit_anonymous(lfp, nfp), True)
        self.assertAlmostEqual(sm.get_rms_anonymous(lfp, nfp)[0],
                               0.060895871160262717)

        #Test partial occupancies.
        s1 = Structure(Lattice.cubic(3),
                       [{"Fe": 0.5}, {"Fe": 0.5}, {"Fe": 0.5}, {"Fe": 0.5}],
                       [[0, 0, 0], [0.25, 0.25, 0.25],
                        [0.5, 0.5, 0.5], [0.75, 0.75, 0.75]])
        s2 = Structure(Lattice.cubic(3),
                       [{"Fe": 0.25}, {"Fe": 0.5}, {"Fe": 0.5}, {"Fe": 0.75}],
                       [[0, 0, 0], [0.25, 0.25, 0.25],
                        [0.5, 0.5, 0.5], [0.75, 0.75, 0.75]])
        self.assertFalse(sm.fit(s1, s2))
        self.assertFalse(sm.fit(s2, s1))
        s2 = Structure(Lattice.cubic(3),
                       [{"Mn": 0.5}, {"Mn": 0.5}, {"Mn": 0.5},
                        {"Mn": 0.5}],
                       [[0, 0, 0], [0.25, 0.25, 0.25],
                        [0.5, 0.5, 0.5], [0.75, 0.75, 0.75]])
        self.assertEqual(sm.fit_anonymous(s1, s2), True)

        self.assertAlmostEqual(sm.get_rms_anonymous(s1, s2)[0], 0)

    def test_oxi(self):
        """Test oxidation state removal matching"""
        sm = StructureMatcher()
        self.assertFalse(sm.fit(self.oxi_structs[0], self.oxi_structs[1]))
        sm = StructureMatcher(comparator=ElementComparator())
        self.assertTrue(sm.fit(self.oxi_structs[0], self.oxi_structs[1]))

    def test_primitive(self):
        """Test primitive cell reduction"""
        sm = StructureMatcher(primitive_cell=True)
        self.struct_list[1].make_supercell([[2, 0, 0], [0, 3, 0], [0, 0, 1]])
        self.assertTrue(sm.fit(self.struct_list[0], self.struct_list[1]))

    def test_class(self):
        # Tests entire class as single working unit
        sm = StructureMatcher()
        # Test group_structures and find_indices
        out = sm.group_structures(self.struct_list)
        self.assertEqual(list(map(len, out)), [4, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1])
        self.assertEqual(sum(map(len, out)), len(self.struct_list))
        for s in self.struct_list[::2]:
            s.replace_species({'Ti': 'Zr', 'O':'Ti'})
        out = sm.group_structures(self.struct_list, anonymous=True)
        self.assertEqual(list(map(len, out)), [4, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1])

    def test_mix(self):
        structures = [self.get_structure("Li2O"),
                      self.get_structure("Li2O2"),
                      self.get_structure("LiFePO4")]
        for fname in ["POSCAR.Li2O", "POSCAR.LiFePO4"]:
            structures.append(Structure.from_file(os.path.join(test_dir, fname)))
        sm = StructureMatcher(comparator=ElementComparator())
        groups = sm.group_structures(structures)
        for g in groups:
            formula = g[0].composition.reduced_formula
            if formula in ["Li2O", "LiFePO4"]:
                self.assertEqual(len(g), 2)
            else:
                self.assertEqual(len(g), 1)

    def test_left_handed_lattice(self):
        """Ensure Left handed lattices are accepted"""
        sm = StructureMatcher()
        s = Structure.from_file(os.path.join(test_dir, "Li3GaPCO7.json"))
        self.assertTrue(sm.fit(s, s))

    def test_as_dict_and_from_dict(self):
        sm = StructureMatcher(ltol=0.1, stol=0.2, angle_tol=2,
                              primitive_cell=False, scale=False,
                              comparator=FrameworkComparator())
        d = sm.as_dict()
        sm2 = StructureMatcher.from_dict(d)
        self.assertEqual(sm2.as_dict(), d)

    def test_no_scaling(self):
        sm = StructureMatcher(ltol=0.1, stol=0.1, angle_tol=2,
                              scale=False, comparator=ElementComparator())
        self.assertTrue(sm.fit(self.struct_list[0], self.struct_list[1]))

        self.assertTrue(sm.get_rms_dist(self.struct_list[0],
                                        self.struct_list[1])[0] < 0.0008)

    def test_supercell_fit(self):
        sm = StructureMatcher(attempt_supercell=False)
        s1 = Structure.from_file(os.path.join(test_dir, "Al3F9.json"))
        s2 = Structure.from_file(os.path.join(test_dir, "Al3F9_distorted.json"))

        self.assertFalse(sm.fit(s1, s2))

        sm = StructureMatcher(attempt_supercell=True)

        self.assertTrue(sm.fit(s1, s2))
        self.assertTrue(sm.fit(s2, s1))

    def test_get_lattices(self):
        sm = StructureMatcher(ltol=0.2, stol=0.3, angle_tol=5,
                              primitive_cell=True, scale=True,
                              attempt_supercell=False)
        l1 = Lattice.from_lengths_and_angles([1, 2.1, 1.9] , [90, 89, 91])
        l2 = Lattice.from_lengths_and_angles([1.1, 2, 2] , [89, 91, 90])
        s1 = Structure(l1, [], [])
        s2 = Structure(l2, [], [])

        lattices = list(sm._get_lattices(s=s1, target_lattice=s2.lattice))
        self.assertEqual(len(lattices), 16)

        l3 = Lattice.from_lengths_and_angles([1.1, 2, 20] , [89, 91, 90])
        s3 = Structure(l3, [], [])

        lattices = list(sm._get_lattices(s=s1, target_lattice=s3.lattice))
        self.assertEqual(len(lattices), 0)

    def test_find_match1(self):
        sm = StructureMatcher(ltol=0.2, stol=0.3, angle_tol=5,
                              primitive_cell=True, scale=True,
                              attempt_supercell=False)
        l = Lattice.orthorhombic(1, 2, 3)
        s1 = Structure(l, ['Si', 'Si', 'Ag'],
                       [[0,0,0.1],[0,0,0.2],[.7,.4,.5]])
        s2 = Structure(l, ['Si', 'Si', 'Ag'],
                       [[0,0.1,0],[0,0.1,-0.95],[.7,.5,.375]])

        s1, s2, fu, s1_supercell = sm._preprocess(s1, s2, False)
        match = sm._strict_match(s1, s2, fu, s1_supercell = True, use_rms = True, break_on_match = False)
        scale_matrix = match[2]
        s2.make_supercell(scale_matrix)
        fc = s2.frac_coords + match[3]
        fc -= np.round(fc)
        self.assertAlmostEqual(np.sum(fc), 0.9)
        self.assertAlmostEqual(np.sum(fc[:,:2]), 0.1)
        cart_dist = np.sum(match[1] * (l.volume/3) ** (1/3))
        self.assertAlmostEqual(cart_dist, 0.15)

    def test_find_match2(self):
        sm = StructureMatcher(ltol=0.2, stol=0.3, angle_tol=5,
                              primitive_cell=True, scale=True,
                              attempt_supercell=False)
        l = Lattice.orthorhombic(1, 2, 3)
        s1 = Structure(l, ['Si', 'Si'], [[0,0,0.1],[0,0,0.2]])
        s2 = Structure(l, ['Si', 'Si'], [[0,0.1,0],[0,0.1,-0.95]])

        s1, s2, fu, s1_supercell = sm._preprocess(s1, s2, False)

        match = sm._strict_match(s1, s2, fu, s1_supercell = False,
                          use_rms = True, break_on_match = False)
        scale_matrix = match[2]
        s2.make_supercell(scale_matrix)
        s2.translate_sites(range(len(s2)), match[3])

        self.assertAlmostEqual(np.sum(s2.frac_coords), 0.3)
        self.assertAlmostEqual(np.sum(s2.frac_coords[:,:2]), 0)

    def test_supercell_subsets(self):
        sm = StructureMatcher(ltol=0.2, stol=0.3, angle_tol=5,
                              primitive_cell=False, scale=True,
                              attempt_supercell=True, allow_subset=True,
                              supercell_size='volume')
        sm_no_s = StructureMatcher(ltol=0.2, stol=0.3, angle_tol=5,
                              primitive_cell=False, scale=True,
                              attempt_supercell=True, allow_subset=False,
                              supercell_size='volume')
        l = Lattice.orthorhombic(1, 2, 3)
        s1 = Structure(l, ['Ag', 'Si', 'Si'],
                       [[.7,.4,.5],[0,0,0.1],[0,0,0.2]])
        s1.make_supercell([2,1,1])
        s2 = Structure(l, ['Si', 'Si', 'Ag'],
                       [[0,0.1,-0.95],[0,0.1,0],[-.7,.5,.375]])

        shuffle = [0,2,1,3,4,5]
        s1 = Structure.from_sites([s1[i] for i in shuffle])

        #test when s1 is exact supercell of s2
        result = sm.get_s2_like_s1(s1, s2)
        for a, b in zip(s1, result):
            self.assertTrue(a.distance(b) < 0.08)
            self.assertEqual(a.species_and_occu, b.species_and_occu)

        self.assertTrue(sm.fit(s1, s2))
        self.assertTrue(sm.fit(s2, s1))
        self.assertTrue(sm_no_s.fit(s1, s2))
        self.assertTrue(sm_no_s.fit(s2, s1))

        rms = (0.048604032430991401, 0.059527539448807391)
        self.assertTrue(np.allclose(sm.get_rms_dist(s1, s2), rms))
        self.assertTrue(np.allclose(sm.get_rms_dist(s2, s1), rms))

        #test when the supercell is a subset of s2
        subset_supercell = s1.copy()
        del subset_supercell[0]
        result = sm.get_s2_like_s1(subset_supercell, s2)
        self.assertEqual(len(result), 6)
        for a, b in zip(subset_supercell, result):
            self.assertTrue(a.distance(b) < 0.08)
            self.assertEqual(a.species_and_occu, b.species_and_occu)

        self.assertTrue(sm.fit(subset_supercell, s2))
        self.assertTrue(sm.fit(s2, subset_supercell))
        self.assertFalse(sm_no_s.fit(subset_supercell, s2))
        self.assertFalse(sm_no_s.fit(s2, subset_supercell))

        rms = (0.053243049896333279, 0.059527539448807336)
        self.assertTrue(np.allclose(sm.get_rms_dist(subset_supercell, s2), rms))
        self.assertTrue(np.allclose(sm.get_rms_dist(s2, subset_supercell), rms))

        #test when s2 (once made a supercell) is a subset of s1
        s2_missing_site = s2.copy()
        del s2_missing_site[1]
        result = sm.get_s2_like_s1(s1, s2_missing_site)
        for a, b in zip((s1[i] for i in (0, 2, 4, 5)), result):
            self.assertTrue(a.distance(b) < 0.08)
            self.assertEqual(a.species_and_occu, b.species_and_occu)

        self.assertTrue(sm.fit(s1, s2_missing_site))
        self.assertTrue(sm.fit(s2_missing_site, s1))
        self.assertFalse(sm_no_s.fit(s1, s2_missing_site))
        self.assertFalse(sm_no_s.fit(s2_missing_site, s1))

        rms = (0.029763769724403633, 0.029763769724403987)
        self.assertTrue(np.allclose(sm.get_rms_dist(s1, s2_missing_site), rms))
        self.assertTrue(np.allclose(sm.get_rms_dist(s2_missing_site, s1), rms))

    def test_get_s2_large_s2(self):
        sm = StructureMatcher(ltol=0.2, stol=0.3, angle_tol=5,
                              primitive_cell=False, scale=False,
                              attempt_supercell=True, allow_subset=False,
                              supercell_size='volume')

        l = Lattice.orthorhombic(1, 2, 3)
        s1 = Structure(l, ['Ag', 'Si', 'Si'],
                       [[.7,.4,.5],[0,0,0.1],[0,0,0.2]])

        l2 = Lattice.orthorhombic(1.01, 2.01, 3.01)
        s2 = Structure(l2, ['Si', 'Si', 'Ag'],
                       [[0,0.1,-0.95],[0,0.1,0],[-.7,.5,.375]])
        s2.make_supercell([[0,-1,0],[1,0,0],[0,0,1]])

        result = sm.get_s2_like_s1(s1, s2)

        for x,y in zip(s1, result):
            self.assertLess(x.distance(y), 0.08)

    def test_get_mapping(self):
        sm = StructureMatcher(ltol=0.2, stol=0.3, angle_tol=5,
                              primitive_cell=False, scale=True,
                              attempt_supercell=False,
                              allow_subset = True)
        l = Lattice.orthorhombic(1, 2, 3)
        s1 = Structure(l, ['Ag', 'Si', 'Si'],
                       [[.7,.4,.5],[0,0,0.1],[0,0,0.2]])
        s1.make_supercell([2,1,1])
        s2 = Structure(l, ['Si', 'Si', 'Ag'],
                       [[0,0.1,-0.95],[0,0.1,0],[-.7,.5,.375]])

        shuffle = [2,0,1,3,5,4]
        s1 = Structure.from_sites([s1[i] for i in shuffle])
        #test the mapping
        s2.make_supercell([2,1,1])
        #equal sizes
        for i, x in enumerate(sm.get_mapping(s1, s2)):
            self.assertEqual(s1[x].species_and_occu,
                             s2[i].species_and_occu)

        del s1[0]
        #s1 is subset of s2
        for i, x in enumerate(sm.get_mapping(s2, s1)):
            self.assertEqual(s1[i].species_and_occu,
                             s2[x].species_and_occu)
        #s2 is smaller than s1
        del s2[0]
        del s2[1]
        self.assertRaises(ValueError, sm.get_mapping, s2, s1)

    def test_get_supercell_matrix(self):
        sm = StructureMatcher(ltol=0.1, stol=0.3, angle_tol=2,
                              primitive_cell=False, scale=True,
                              attempt_supercell=True)

        l = Lattice.orthorhombic(1, 2, 3)

        s1 = Structure(l, ['Si', 'Si', 'Ag'],
                       [[0,0,0.1],[0,0,0.2],[.7,.4,.5]])
        s1.make_supercell([2,1,1])
        s2 = Structure(l, ['Si', 'Si', 'Ag'],
                       [[0,0.1,0],[0,0.1,-0.95],[-.7,.5,.375]])
        result = sm.get_supercell_matrix(s1, s2)
        self.assertTrue((result == [[-2,0,0],[0,1,0],[0,0,1]]).all())

        s1 = Structure(l, ['Si', 'Si', 'Ag'],
                       [[0,0,0.1],[0,0,0.2],[.7,.4,.5]])
        s1.make_supercell([[1, -1, 0],[0, 0, -1],[0, 1, 0]])

        s2 = Structure(l, ['Si', 'Si', 'Ag'],
                       [[0,0.1,0],[0,0.1,-0.95],[-.7,.5,.375]])
        result = sm.get_supercell_matrix(s1, s2)
        self.assertTrue((result == [[-1,-1,0],[0,0,-1],[0,1,0]]).all())

        #test when the supercell is a subset
        sm = StructureMatcher(ltol=0.1, stol=0.3, angle_tol=2,
                              primitive_cell=False, scale=True,
                              attempt_supercell=True, allow_subset=True)
        del s1[0]
        result = sm.get_supercell_matrix(s1, s2)
        self.assertTrue((result == [[-1,-1,0],[0,0,-1],[0,1,0]]).all())

    def test_subset(self):
        sm = StructureMatcher(ltol=0.2, stol=0.3, angle_tol=5,
                              primitive_cell=False, scale=True,
                              attempt_supercell=False,
                              allow_subset=True)
        l = Lattice.orthorhombic(10, 20, 30)
        s1 = Structure(l, ['Si', 'Si', 'Ag'],
                       [[0,0,0.1],[0,0,0.2],[.7,.4,.5]])
        s2 = Structure(l, ['Si', 'Ag'],
                       [[0,0.1,0],[-.7,.5,.4]])
        result = sm.get_s2_like_s1(s1, s2)

        self.assertEqual(len(find_in_coord_list_pbc(result.frac_coords,
                                                    [0,0,0.1])), 1)
        self.assertEqual(len(find_in_coord_list_pbc(result.frac_coords,
                                                    [0.7,0.4,0.5])), 1)

        #test with fewer species in s2
        s1 = Structure(l, ['Si', 'Ag', 'Si'],
                       [[0,0,0.1],[0,0,0.2],[.7,.4,.5]])
        s2 = Structure(l, ['Si', 'Si'],
                       [[0,0.1,0],[-.7,.5,.4]])
        result = sm.get_s2_like_s1(s1, s2)
        mindists = np.min(s1.lattice.get_all_distances(
            s1.frac_coords, result.frac_coords), axis=0)
        self.assertLess(np.max(mindists), 1e-6)

        self.assertEqual(len(find_in_coord_list_pbc(result.frac_coords,
                                                    [0,0,0.1])), 1)
        self.assertEqual(len(find_in_coord_list_pbc(result.frac_coords,
                                                    [0.7,0.4,0.5])), 1)

        #test with not enough sites in s1
        #test with fewer species in s2
        s1 = Structure(l, ['Si', 'Ag', 'Cl'],
                       [[0,0,0.1],[0,0,0.2],[.7,.4,.5]])
        s2 = Structure(l, ['Si', 'Si'],
                       [[0,0.1,0],[-.7,.5,.4]])
        self.assertEqual(sm.get_s2_like_s1(s1, s2), None)

    def test_out_of_cell_s2_like_s1(self):
        l = Lattice.cubic(5)
        s1 = Structure(l, ['Si', 'Ag', 'Si'],
                       [[0,0,-0.02],[0,0,0.001],[.7,.4,.5]])
        s2 = Structure(l, ['Si', 'Ag', 'Si'],
                       [[0,0,0.98],[0,0,0.99],[.7,.4,.5]])
        new_s2 = StructureMatcher(primitive_cell=False).get_s2_like_s1(s1, s2)
        dists = np.sum((s1.cart_coords - new_s2.cart_coords) ** 2, axis=-1) ** 0.5
        self.assertLess(np.max(dists), 0.1)

    def test_disordered_primitive_to_ordered_supercell(self):
        sm_atoms = StructureMatcher(ltol=0.2, stol=0.3, angle_tol=5,
                                    primitive_cell=False, scale=True,
                                    attempt_supercell=True,
                                    allow_subset=True,
                                    supercell_size = 'num_atoms',
                                    comparator=OrderDisorderElementComparator())
        sm_sites = StructureMatcher(ltol=0.2, stol=0.3, angle_tol=5,
                                    primitive_cell=False, scale=True,
                                    attempt_supercell=True,
                                    allow_subset=True,
                                    supercell_size = 'num_sites',
                                    comparator=OrderDisorderElementComparator())
        lp = Lattice.orthorhombic(10, 20, 30)
        pcoords = [[0,   0,   0],
                   [0.5, 0.5, 0.5]]
        ls = Lattice.orthorhombic(20,20,30)
        scoords = [[0,    0,   0],
                   [0.75, 0.5, 0.5]]
        prim = Structure(lp, [{'Na':0.5}, {'Cl':0.5}], pcoords)
        supercell = Structure(ls, ['Na', 'Cl'], scoords)
        supercell.make_supercell([[-1,1,0],[0,1,1],[1,0,0]])

        self.assertFalse(sm_sites.fit(prim, supercell))
        self.assertTrue(sm_atoms.fit(prim, supercell))

        self.assertRaises(ValueError, sm_atoms.get_s2_like_s1, prim, supercell)
        self.assertEqual(len(sm_atoms.get_s2_like_s1(supercell, prim)), 4)

    def test_ordered_primitive_to_disordered_supercell(self):
        sm_atoms = StructureMatcher(ltol=0.2, stol=0.3, angle_tol=5,
                                    primitive_cell=False, scale=True,
                                    attempt_supercell=True,
                                    allow_subset=True,
                                    supercell_size = 'num_atoms',
                                    comparator=OrderDisorderElementComparator())
        sm_sites = StructureMatcher(ltol=0.2, stol=0.3, angle_tol=5,
                                    primitive_cell=False, scale=True,
                                    attempt_supercell=True,
                                    allow_subset=True,
                                    supercell_size = 'num_sites',
                                    comparator=OrderDisorderElementComparator())
        lp = Lattice.orthorhombic(10, 20, 30)
        pcoords = [[0,   0,   0],
                   [0.5, 0.5, 0.5]]
        ls = Lattice.orthorhombic(20,20,30)
        scoords = [[0,    0,   0],
                   [0.5,  0,   0],
                   [0.25, 0.5, 0.5],
                   [0.75, 0.5, 0.5]]
        s1 = Structure(lp, ['Na', 'Cl'], pcoords)
        s2 = Structure(ls, [{'Na':0.5}, {'Na':0.5}, {'Cl':0.5}, {'Cl':0.5}], scoords)

        self.assertTrue(sm_sites.fit(s1, s2))
        self.assertFalse(sm_atoms.fit(s1, s2))

    def test_electronegativity(self):
        sm = StructureMatcher(ltol=0.2, stol=0.3, angle_tol=5)

        s1 = Structure.from_file(os.path.join(test_dir, "Na2Fe2PAsO4S4.json"))
        s2 = Structure.from_file(os.path.join(test_dir, "Na2Fe2PNO4Se4.json"))
        self.assertEqual(sm.get_best_electronegativity_anonymous_mapping(s1, s2),
                    {Element('S'): Element('Se'),
                     Element('As'): Element('N'),
                     Element('Fe'): Element('Fe'),
                     Element('Na'): Element('Na'),
                     Element('P'): Element('P'),
                     Element('O'): Element('O'),})
        self.assertEqual(len(sm.get_all_anonymous_mappings(s1, s2)), 2)

        # test include_dist
        dists = {Element('N'): 0, Element('P'): 0.0010725064}
        for mapping, d in sm.get_all_anonymous_mappings(s1, s2, include_dist=True):
            self.assertAlmostEqual(dists[mapping[Element('As')]], d)

    def test_rms_vs_minimax(self):
        # This tests that structures with adjusted RMS less than stol, but minimax
        # greater than stol are treated properly
        # stol=0.3 gives exactly an ftol of 0.1 on the c axis
        sm = StructureMatcher(ltol=0.2, stol=0.301, angle_tol=1, primitive_cell=False)
        l = Lattice.orthorhombic(1, 2, 12)

        sp = ["Si", "Si", "Al"]
        s1 = Structure(l, sp, [[0.5, 0, 0], [0, 0, 0], [0, 0, 0.5]])
        s2 = Structure(l, sp, [[0.5, 0, 0], [0, 0, 0], [0, 0, 0.6]])
        self.assertArrayAlmostEqual(sm.get_rms_dist(s1, s2),
                                    (0.32 ** 0.5 / 2, 0.4))

        self.assertEqual(sm.fit(s1, s2), False)
        self.assertEqual(sm.fit_anonymous(s1, s2), False)
        self.assertEqual(sm.get_mapping(s1, s2), None)


if __name__ == '__main__':
    unittest.main()
