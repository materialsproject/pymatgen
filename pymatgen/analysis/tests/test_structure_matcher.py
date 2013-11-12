from __future__ import division
import unittest
import os
import json
import numpy as np

from pymatgen.analysis.structure_matcher import StructureMatcher, \
    ElementComparator, FrameworkComparator, OrderDisorderElementComparator
from pymatgen.serializers.json_coders import PMGJSONDecoder
from pymatgen.core.operations import SymmOp
from pymatgen.io.smartio import read_structure
from pymatgen.core import Structure, Composition, Lattice
from pymatgen.util.coord_utils import find_in_coord_list_pbc

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class StructureMatcherTest(unittest.TestCase):

    def setUp(self):
        with open(os.path.join(test_dir, "TiO2_entries.json"), 'r') as fp:
            entries = json.load(fp, cls=PMGJSONDecoder)
        self.struct_list = [e.structure for e in entries]
        self.oxi_structs = [read_structure(os.path.join(test_dir, fname))
                            for fname in ["Li2O.cif", "POSCAR.Li2O"]]

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
        lfp = read_structure(os.path.join(test_dir, "LiFePO4.cif"))
        nfp = read_structure(os.path.join(test_dir, "NaFePO4.cif"))
        self.assertTrue(sm2.fit(lfp, nfp))
        self.assertFalse(sm.fit(lfp, nfp))

        #Test anonymous fit.
        self.assertEqual(sm.fit_anonymous(lfp, nfp),
                         {Composition("Li"): Composition("Na")})
        self.assertAlmostEqual(sm.get_minimax_rms_anonymous(lfp, nfp)[0],
                               0.096084154118549828)

        #Test partial occupancies.
        s1 = Structure([[3, 0, 0], [0, 3, 0], [0, 0, 3]],
                       [{"Fe": 0.5}, {"Fe": 0.5}, {"Fe": 0.5}, {"Fe": 0.5}],
                       [[0, 0, 0], [0.25, 0.25, 0.25],
                        [0.5, 0.5, 0.5], [0.75, 0.75, 0.75]])
        s2 = Structure([[3, 0, 0], [0, 3, 0], [0, 0, 3]],
                       [{"Fe": 0.25}, {"Fe": 0.5}, {"Fe": 0.5}, {"Fe": 0.75}],
                       [[0, 0, 0], [0.25, 0.25, 0.25],
                        [0.5, 0.5, 0.5], [0.75, 0.75, 0.75]])
        self.assertFalse(sm.fit(s1, s2))
        self.assertFalse(sm.fit(s2, s1))
        s2 = Structure([[3, 0, 0], [0, 3, 0], [0, 0, 3]],
                       [{"Fe": 0.25}, {"Fe": 0.25}, {"Fe": 0.25},
                        {"Fe": 0.25}],
                       [[0, 0, 0], [0.25, 0.25, 0.25],
                        [0.5, 0.5, 0.5], [0.75, 0.75, 0.75]])
        self.assertEqual(sm.fit_anonymous(s1, s2),
                         {Composition("Fe0.5"): Composition("Fe0.25")})

        self.assertAlmostEqual(sm.get_minimax_rms_anonymous(s1, s2)[0], 0)

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
        self.assertEqual(map(len, out), [4, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1])
        self.assertEqual(sum(map(len, out)), len(self.struct_list))

    def test_mix(self):
        structures = []
        for fname in ["POSCAR.Li2O", "Li2O.cif", "Li2O2.cif", "LiFePO4.cif",
                      "POSCAR.LiFePO4"]:
            structures.append(read_structure(os.path.join(test_dir, fname)))
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
        s = read_structure(os.path.join(test_dir, "Li3GaPCO7.cif"))
        self.assertTrue(sm.fit(s, s))

    def test_to_dict_and_from_dict(self):
        sm = StructureMatcher(ltol=0.1, stol=0.2, angle_tol=2,
                              primitive_cell=False, scale=False,
                              comparator=FrameworkComparator())
        d = sm.to_dict
        sm2 = StructureMatcher.from_dict(d)
        self.assertEqual(sm2.to_dict, d)

    def test_no_scaling(self):
        sm = StructureMatcher(ltol=0.1, stol=0.1, angle_tol=2,
                              scale=False, comparator=ElementComparator())
        self.assertTrue(sm.fit(self.struct_list[0], self.struct_list[1]))

        self.assertTrue(sm.get_rms_dist(self.struct_list[0],
                                        self.struct_list[1])[0] < 0.0008)

    def test_supercell_fit(self):
        sm = StructureMatcher(attempt_supercell=False)
        s1 = read_structure(os.path.join(test_dir, "Al3F9.cif"))
        s2 = read_structure(os.path.join(test_dir, "Al3F9_distorted.cif"))

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

        lattices = list(sm._get_lattices(s = s1, target_s = s2))
        self.assertEqual(len(lattices), 16)

        l3 = Lattice.from_lengths_and_angles([1.1, 2, 20] , [89, 91, 90])
        s3 = Structure(l3, [], [])

        lattices = list(sm._get_lattices(s = s1, target_s = s3))
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
        match = sm._find_match(s1, s2, break_on_match = False,
                               use_rms = True, niggli = False)
        scale_matrix = np.round(np.dot(match[2].matrix,
                            s2.lattice.inv_matrix)).astype('int')
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
        match = sm._find_match(s1, s2, break_on_match = False,
                               use_rms = True, niggli = False)
        scale_matrix = np.round(np.dot(match[2].matrix,
                                       s2.lattice.inv_matrix)).astype('int')
        s2.make_supercell(scale_matrix)
        fc = s2.frac_coords + match[3]
        fc -= np.round(fc)

        self.assertAlmostEqual(np.sum(fc), 0.3)
        self.assertAlmostEqual(np.sum(fc[:,:2]), 0)

    def test_get_s2_like_s1(self):
        sm = StructureMatcher(ltol=0.2, stol=0.3, angle_tol=5,
                              primitive_cell=False, scale=True,
                              attempt_supercell=True)
        l = Lattice.orthorhombic(1, 2, 3)
        s1 = Structure(l, ['Si', 'Si', 'Ag'],
                       [[0,0,0.1],[0,0,0.2],[.7,.4,.5]])
        s1.make_supercell([2,1,1])
        s2 = Structure(l, ['Si', 'Si', 'Ag'],
                       [[0,0.1,0],[0,0.1,-0.95],[-.7,.5,.375]])
        result = sm.get_s2_like_s1(s1, s2)

        self.assertEqual(len(find_in_coord_list_pbc(result.frac_coords,
                                                    [0.35,0.4,0.5])), 1)
        self.assertEqual(len(find_in_coord_list_pbc(result.frac_coords,
                                                    [0,0,0.125])), 1)
        self.assertEqual(len(find_in_coord_list_pbc(result.frac_coords,
                                                    [0,0,0.175])), 1)

    def test_get_supercell_matrix(self):
        sm = StructureMatcher(ltol=0.2, stol=0.3, angle_tol=5,
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
        s1 = Structure(lp, [{'Na':0.5}, {'Cl':0.5}], pcoords)
        s2 = Structure(ls, ['Na', 'Cl'], scoords)

        self.assertFalse(sm_sites.fit(s1, s2))
        self.assertTrue(sm_atoms.fit(s1, s2))
        self.assertRaises(ValueError, sm_atoms.get_s2_like_s1, s1, s2)
        self.assertEqual(len(sm_atoms.get_s2_like_s1(s2, s1)), 4)

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

        s1 = read_structure(os.path.join(test_dir, "Na2Fe2PAsO4S4.cif"))
        s2 = read_structure(os.path.join(test_dir, "Na2Fe2PNO4Se4.cif"))
        self.assertAlmostEqual(sm.fit_with_electronegativity(s1, s2),
                               {Composition('S'): Composition('Se'), Composition('As'): Composition('N')})


if __name__ == '__main__':
    unittest.main()
