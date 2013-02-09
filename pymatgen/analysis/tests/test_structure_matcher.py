import unittest
import os
import json
import numpy as np
import random
from pymatgen.analysis.structure_matcher import StructureMatcher, \
    ElementComparator, FrameworkComparator
from pymatgen.serializers.json_coders import PMGJSONDecoder
from pymatgen.core.operations import SymmOp
from pymatgen.core.structure_modifier import StructureEditor
from pymatgen.core.structure_modifier import SupercellMaker
from pymatgen.io.smartio import read_structure
from pymatgen.core.structure import Structure

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class StructureMatcherTest(unittest.TestCase):
    def setUp(self):
        with open(os.path.join(test_dir, "TiO2_entries.json"), 'rb') as fp:
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

        """Test rotational/translational invariance"""
        op = SymmOp.from_axis_angle_and_translation([0, 0, 1], 30, False,
                                                    np.array([0.4, 0.7, 0.9]))
        editor = StructureEditor(self.struct_list[1])
        editor.apply_operation(op)
        self.assertTrue(sm.fit(self.struct_list[0], editor.modified_structure))

        """Test failure under large atomic translation"""
        editor.translate_sites([0], [.4, .4, .2], frac_coords=True)
        self.assertFalse(sm.fit(self.struct_list[0],
                                editor.modified_structure))

        editor.translate_sites([0], [-.4, -.4, -.2], frac_coords=True)
        """Test match under shuffling of sites"""
        random.shuffle(editor._sites)
        self.assertTrue(sm.fit(self.struct_list[0], editor.modified_structure))
        """Test FrameworkComporator"""
        sm2 = StructureMatcher(comparator=FrameworkComparator())
        lfp = read_structure(os.path.join(test_dir, "LiFePO4.cif"))
        nfp = read_structure(os.path.join(test_dir, "NaFePO4.cif"))
        self.assertTrue(sm2.fit(lfp, nfp))
        self.assertFalse(sm.fit(lfp, nfp))

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

    def test_oxi(self):
        """Test oxidation state removal matching"""
        sm = StructureMatcher()
        self.assertFalse(sm.fit(self.oxi_structs[0], self.oxi_structs[1]))
        sm = StructureMatcher(comparator=ElementComparator())
        self.assertTrue(sm.fit(self.oxi_structs[0], self.oxi_structs[1]))

    def test_primitive(self):
        """Test primitive cell reduction"""
        sm = StructureMatcher(primitive_cell=True)
        mod = SupercellMaker(self.struct_list[1],
                             scaling_matrix=[[2, 0, 0], [0, 3, 0], [0, 0, 1]])
        super_cell = mod.modified_structure
        self.assertTrue(sm.fit(self.struct_list[0], super_cell))

    def test_class(self):
        """Tests entire class as single working unit"""
        sm = StructureMatcher()
        """ Test group_structures and find_indicies"""
        out = sm.group_structures(self.struct_list)

        self.assertEqual(sm.find_indexes(self.struct_list, out),
                         [0, 0, 0, 1, 2, 3, 4, 0, 5, 6, 7, 8, 8, 9, 9, 10])

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


if __name__ == '__main__':
    unittest.main()