# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import os
import unittest

import numpy as np

from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Molecule, Structure
from pymatgen.io.cif import CifParser
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.symmetry.analyzer import (
    PointGroupAnalyzer,
    SpacegroupAnalyzer,
    cluster_sites,
    iterative_symmetrize,
)
from pymatgen.util.testing import PymatgenTest

test_dir_mol = os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules")


class SpacegroupAnalyzerTest(PymatgenTest):
    def setUp(self):
        p = Poscar.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR"))
        self.structure = p.structure
        self.sg = SpacegroupAnalyzer(self.structure, 0.001)
        self.disordered_structure = self.get_structure("Li10GeP2S12")
        self.disordered_sg = SpacegroupAnalyzer(self.disordered_structure, 0.001)
        s = p.structure.copy()
        site = s[0]
        del s[0]
        s.append(site.species, site.frac_coords)
        self.sg3 = SpacegroupAnalyzer(s, 0.001)
        graphite = self.get_structure("Graphite")
        graphite.add_site_property("magmom", [0.1] * len(graphite))
        self.sg4 = SpacegroupAnalyzer(graphite, 0.001)
        self.structure4 = graphite

    def test_primitive(self):
        s = Structure.from_spacegroup("Fm-3m", np.eye(3) * 3, ["Cu"], [[0, 0, 0]])
        a = SpacegroupAnalyzer(s)
        self.assertEqual(len(s), 4)
        self.assertEqual(len(a.find_primitive()), 1)

    def test_is_laue(self):
        s = Structure.from_spacegroup("Fm-3m", np.eye(3) * 3, ["Cu"], [[0, 0, 0]])
        a = SpacegroupAnalyzer(s)
        self.assertTrue(a.is_laue())

    def test_magnetic(self):
        lfp = PymatgenTest.get_structure("LiFePO4")
        sg = SpacegroupAnalyzer(lfp, 0.1)
        self.assertEqual(sg.get_space_group_symbol(), "Pnma")
        magmoms = [0] * len(lfp)
        magmoms[4] = 1
        magmoms[5] = -1
        magmoms[6] = 1
        magmoms[7] = -1
        lfp.add_site_property("magmom", magmoms)
        sg = SpacegroupAnalyzer(lfp, 0.1)
        self.assertEqual(sg.get_space_group_symbol(), "Pnma")

    def test_get_space_symbol(self):
        self.assertEqual(self.sg.get_space_group_symbol(), "Pnma")
        self.assertEqual(self.disordered_sg.get_space_group_symbol(), "P4_2/nmc")
        self.assertEqual(self.sg3.get_space_group_symbol(), "Pnma")
        self.assertEqual(self.sg4.get_space_group_symbol(), "P6_3/mmc")

    def test_get_space_number(self):
        self.assertEqual(self.sg.get_space_group_number(), 62)
        self.assertEqual(self.disordered_sg.get_space_group_number(), 137)
        self.assertEqual(self.sg4.get_space_group_number(), 194)

    def test_get_hall(self):
        self.assertEqual(self.sg.get_hall(), "-P 2ac 2n")
        self.assertEqual(self.disordered_sg.get_hall(), "P 4n 2n -1n")

    def test_get_pointgroup(self):
        self.assertEqual(self.sg.get_point_group_symbol(), "mmm")
        self.assertEqual(self.disordered_sg.get_point_group_symbol(), "4/mmm")

    def test_get_symmetry_operations(self):

        for sg, structure in [(self.sg, self.structure), (self.sg4, self.structure4)]:

            pgops = sg.get_point_group_operations()
            fracsymmops = sg.get_symmetry_operations()
            symmops = sg.get_symmetry_operations(True)
            latt = structure.lattice
            for fop, op, pgop in zip(fracsymmops, symmops, pgops):
                # translation vector values should all be 0 or 0.5
                t = fop.translation_vector * 2
                self.assertArrayAlmostEqual(t - np.round(t), 0)

                self.assertArrayAlmostEqual(fop.rotation_matrix, pgop.rotation_matrix)
                for site in structure:
                    newfrac = fop.operate(site.frac_coords)
                    newcart = op.operate(site.coords)
                    self.assertTrue(np.allclose(latt.get_fractional_coords(newcart), newfrac))
                    found = False
                    newsite = PeriodicSite(site.species, newcart, latt, coords_are_cartesian=True)
                    for testsite in structure:
                        if newsite.is_periodic_image(testsite, 1e-3):
                            found = True
                            break
                    self.assertTrue(found)

                # Make sure this works for any position, not just the atomic
                # ones.
                random_fcoord = np.random.uniform(size=(3))
                random_ccoord = latt.get_cartesian_coords(random_fcoord)
                newfrac = fop.operate(random_fcoord)
                newcart = op.operate(random_ccoord)
                self.assertTrue(np.allclose(latt.get_fractional_coords(newcart), newfrac))

    def test_get_symmetry_dataset(self):
        ds = self.sg.get_symmetry_dataset()
        self.assertEqual(ds["international"], "Pnma")

    def test_get_crystal_system(self):
        crystal_system = self.sg.get_crystal_system()
        self.assertEqual("orthorhombic", crystal_system)
        self.assertEqual("tetragonal", self.disordered_sg.get_crystal_system())

        orig_spg = self.sg._space_group_data["number"]
        self.sg._space_group_data["number"] = 0
        try:
            crystal_system = self.sg.get_crystal_system()
        except ValueError as exc:
            self.assertEqual(str(exc), "Received invalid space group 0")
        finally:
            self.sg._space_group_data["number"] = orig_spg

    def test_get_refined_structure(self):
        for a in self.sg.get_refined_structure().lattice.angles:
            self.assertEqual(a, 90)
        refined = self.disordered_sg.get_refined_structure()
        for a in refined.lattice.angles:
            self.assertEqual(a, 90)
        self.assertEqual(refined.lattice.a, refined.lattice.b)

        structure = self.get_structure("Li2O")
        structure.add_site_property("magmom", [1.0] * len(structure))
        sg = SpacegroupAnalyzer(structure, 0.01)
        refined_struct = sg.get_refined_structure(keep_site_properties=True)
        self.assertEqual(refined_struct.site_properties["magmom"], [1.0] * len(refined_struct))

        structure = self.get_structure("Li2O")
        structure.add_site_property("magmom", [1.0] * len(structure))
        sg = SpacegroupAnalyzer(structure, 0.01)
        refined_struct = sg.get_refined_structure(keep_site_properties=False)
        self.assertEqual(refined_struct.site_properties.get("magmom", None), None)

    def test_get_symmetrized_structure(self):
        symm_struct = self.sg.get_symmetrized_structure()
        for a in symm_struct.lattice.angles:
            self.assertEqual(a, 90)
        self.assertEqual(len(symm_struct.equivalent_sites), 5)

        symm_struct = self.disordered_sg.get_symmetrized_structure()
        self.assertEqual(len(symm_struct.equivalent_sites), 8)
        self.assertEqual([len(i) for i in symm_struct.equivalent_sites], [16, 4, 8, 4, 2, 8, 8, 8])
        s1 = symm_struct.equivalent_sites[1][1]
        s2 = symm_struct[symm_struct.equivalent_indices[1][1]]
        self.assertEqual(s1, s2)
        self.assertEqual(self.sg4.get_symmetrized_structure()[0].magmom, 0.1)
        self.assertEqual(symm_struct.wyckoff_symbols[0], "16h")

        # Check copying
        self.assertEqual(symm_struct.copy(), symm_struct)
        d = symm_struct.as_dict()
        from pymatgen.symmetry.structure import SymmetrizedStructure

        ss = SymmetrizedStructure.from_dict(d)
        self.assertEqual(ss.wyckoff_symbols[0], "16h")
        self.assertIn("SymmetrizedStructure", str(ss))

    def test_find_primitive(self):
        """
        F m -3 m Li2O testing of converting to primitive cell
        """
        parser = CifParser(os.path.join(PymatgenTest.TEST_FILES_DIR, "Li2O.cif"))
        structure = parser.get_structures(False)[0]
        s = SpacegroupAnalyzer(structure)
        primitive_structure = s.find_primitive()
        self.assertEqual(primitive_structure.formula, "Li2 O1")
        self.assertTrue(primitive_structure.site_properties.get("magmom", None) is None)
        # This isn't what is expected. All the angles should be 60
        self.assertAlmostEqual(primitive_structure.lattice.alpha, 60)
        self.assertAlmostEqual(primitive_structure.lattice.beta, 60)
        self.assertAlmostEqual(primitive_structure.lattice.gamma, 60)
        self.assertAlmostEqual(primitive_structure.lattice.volume, structure.lattice.volume / 4.0)

        structure = parser.get_structures(False)[0]
        structure.add_site_property("magmom", [1.0] * len(structure))
        s = SpacegroupAnalyzer(structure)
        primitive_structure = s.find_primitive(keep_site_properties=True)
        self.assertEqual(primitive_structure.site_properties["magmom"], [1.0] * len(primitive_structure))

        structure = parser.get_structures(False)[0]
        structure.add_site_property("magmom", [1.0] * len(structure))
        s = SpacegroupAnalyzer(structure)
        primitive_structure = s.find_primitive(keep_site_properties=False)
        self.assertEqual(primitive_structure.site_properties.get("magmom", None), None)

    def test_get_ir_reciprocal_mesh(self):
        grid = self.sg.get_ir_reciprocal_mesh()
        self.assertEqual(len(grid), 216)
        self.assertAlmostEqual(grid[1][0][0], 0.1)
        self.assertAlmostEqual(grid[1][0][1], 0.0)
        self.assertAlmostEqual(grid[1][0][2], 0.0)
        self.assertAlmostEqual(grid[1][1], 2)

    def test_get_ir_reciprocal_mesh_map(self):
        mesh = (6, 6, 6)
        grid = self.sg.get_ir_reciprocal_mesh(mesh=mesh)
        full_grid, mapping = self.sg.get_ir_reciprocal_mesh_map(mesh=mesh)
        self.assertEqual(len(np.unique(mapping)), len(grid))
        for _, i in enumerate(np.unique(mapping)):
            self.assertAlmostEqual(full_grid[i][0], grid[_][0][0])
            self.assertAlmostEqual(full_grid[i][1], grid[_][0][1])
            self.assertAlmostEqual(full_grid[i][2], grid[_][0][2])

    def test_get_conventional_standard_structure(self):
        parser = CifParser(os.path.join(PymatgenTest.TEST_FILES_DIR, "bcc_1927.cif"))
        structure = parser.get_structures(False)[0]
        s = SpacegroupAnalyzer(structure, symprec=1e-2)
        conv = s.get_conventional_standard_structure()
        self.assertAlmostEqual(conv.lattice.alpha, 90)
        self.assertAlmostEqual(conv.lattice.beta, 90)
        self.assertAlmostEqual(conv.lattice.gamma, 90)
        self.assertAlmostEqual(conv.lattice.a, 9.1980270633769461)
        self.assertAlmostEqual(conv.lattice.b, 9.1980270633769461)
        self.assertAlmostEqual(conv.lattice.c, 9.1980270633769461)

        parser = CifParser(os.path.join(PymatgenTest.TEST_FILES_DIR, "btet_1915.cif"))
        structure = parser.get_structures(False)[0]
        s = SpacegroupAnalyzer(structure, symprec=1e-2)
        conv = s.get_conventional_standard_structure()
        self.assertAlmostEqual(conv.lattice.alpha, 90)
        self.assertAlmostEqual(conv.lattice.beta, 90)
        self.assertAlmostEqual(conv.lattice.gamma, 90)
        self.assertAlmostEqual(conv.lattice.a, 5.0615106678044235)
        self.assertAlmostEqual(conv.lattice.b, 5.0615106678044235)
        self.assertAlmostEqual(conv.lattice.c, 4.2327080177761687)

        parser = CifParser(os.path.join(PymatgenTest.TEST_FILES_DIR, "orci_1010.cif"))
        structure = parser.get_structures(False)[0]
        s = SpacegroupAnalyzer(structure, symprec=1e-2)
        conv = s.get_conventional_standard_structure()
        self.assertAlmostEqual(conv.lattice.alpha, 90)
        self.assertAlmostEqual(conv.lattice.beta, 90)
        self.assertAlmostEqual(conv.lattice.gamma, 90)
        self.assertAlmostEqual(conv.lattice.a, 2.9542233922299999)
        self.assertAlmostEqual(conv.lattice.b, 4.6330325651443296)
        self.assertAlmostEqual(conv.lattice.c, 5.373703587040775)

        parser = CifParser(os.path.join(PymatgenTest.TEST_FILES_DIR, "orcc_1003.cif"))
        structure = parser.get_structures(False)[0]
        s = SpacegroupAnalyzer(structure, symprec=1e-2)
        conv = s.get_conventional_standard_structure()
        self.assertAlmostEqual(conv.lattice.alpha, 90)
        self.assertAlmostEqual(conv.lattice.beta, 90)
        self.assertAlmostEqual(conv.lattice.gamma, 90)
        self.assertAlmostEqual(conv.lattice.a, 4.1430033493799998)
        self.assertAlmostEqual(conv.lattice.b, 31.437979757624728)
        self.assertAlmostEqual(conv.lattice.c, 3.99648651)

        parser = CifParser(os.path.join(PymatgenTest.TEST_FILES_DIR, "orac_632475.cif"))
        structure = parser.get_structures(False)[0]
        s = SpacegroupAnalyzer(structure, symprec=1e-2)
        conv = s.get_conventional_standard_structure()
        self.assertAlmostEqual(conv.lattice.alpha, 90)
        self.assertAlmostEqual(conv.lattice.beta, 90)
        self.assertAlmostEqual(conv.lattice.gamma, 90)
        self.assertAlmostEqual(conv.lattice.a, 3.1790663399999999)
        self.assertAlmostEqual(conv.lattice.b, 9.9032878699999998)
        self.assertAlmostEqual(conv.lattice.c, 3.5372412099999999)

        parser = CifParser(os.path.join(PymatgenTest.TEST_FILES_DIR, "monoc_1028.cif"))
        structure = parser.get_structures(False)[0]
        s = SpacegroupAnalyzer(structure, symprec=1e-2)
        conv = s.get_conventional_standard_structure()
        self.assertAlmostEqual(conv.lattice.alpha, 90)
        self.assertAlmostEqual(conv.lattice.beta, 117.53832420192903)
        self.assertAlmostEqual(conv.lattice.gamma, 90)
        self.assertAlmostEqual(conv.lattice.a, 14.033435583000625)
        self.assertAlmostEqual(conv.lattice.b, 3.96052850731)
        self.assertAlmostEqual(conv.lattice.c, 6.8743926325200002)
        parser = CifParser(os.path.join(PymatgenTest.TEST_FILES_DIR, "hex_1170.cif"))
        structure = parser.get_structures(False)[0]
        s = SpacegroupAnalyzer(structure, symprec=1e-2)
        conv = s.get_conventional_standard_structure()
        self.assertAlmostEqual(conv.lattice.alpha, 90)
        self.assertAlmostEqual(conv.lattice.beta, 90)
        self.assertAlmostEqual(conv.lattice.gamma, 120)
        self.assertAlmostEqual(conv.lattice.a, 3.699919902005897)
        self.assertAlmostEqual(conv.lattice.b, 3.699919902005897)
        self.assertAlmostEqual(conv.lattice.c, 6.9779585500000003)

        structure = Structure.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "tric_684654.json"))
        s = SpacegroupAnalyzer(structure, symprec=1e-2)
        conv = s.get_conventional_standard_structure()
        self.assertAlmostEqual(conv.lattice.alpha, 74.09581916308757)
        self.assertAlmostEqual(conv.lattice.beta, 75.72817279281173)
        self.assertAlmostEqual(conv.lattice.gamma, 63.63234318667333)
        self.assertAlmostEqual(conv.lattice.a, 3.741372924048738)
        self.assertAlmostEqual(conv.lattice.b, 3.9883228679270686)
        self.assertAlmostEqual(conv.lattice.c, 7.288495840048958)

        structure = Structure.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "tric_684654.json"))
        structure.add_site_property("magmom", [1.0] * len(structure))
        s = SpacegroupAnalyzer(structure, symprec=1e-2)
        conv = s.get_conventional_standard_structure(keep_site_properties=True)
        self.assertEqual(conv.site_properties["magmom"], [1.0] * len(conv))

        structure = Structure.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "tric_684654.json"))
        structure.add_site_property("magmom", [1.0] * len(structure))
        s = SpacegroupAnalyzer(structure, symprec=1e-2)
        conv = s.get_conventional_standard_structure(keep_site_properties=False)
        self.assertEqual(conv.site_properties.get("magmom", None), None)

    def test_get_primitive_standard_structure(self):
        parser = CifParser(os.path.join(PymatgenTest.TEST_FILES_DIR, "bcc_1927.cif"))
        structure = parser.get_structures(False)[0]
        s = SpacegroupAnalyzer(structure, symprec=1e-2)
        prim = s.get_primitive_standard_structure()
        self.assertAlmostEqual(prim.lattice.alpha, 109.47122063400001)
        self.assertAlmostEqual(prim.lattice.beta, 109.47122063400001)
        self.assertAlmostEqual(prim.lattice.gamma, 109.47122063400001)
        self.assertAlmostEqual(prim.lattice.a, 7.9657251015812145)
        self.assertAlmostEqual(prim.lattice.b, 7.9657251015812145)
        self.assertAlmostEqual(prim.lattice.c, 7.9657251015812145)

        parser = CifParser(os.path.join(PymatgenTest.TEST_FILES_DIR, "btet_1915.cif"))
        structure = parser.get_structures(False)[0]
        s = SpacegroupAnalyzer(structure, symprec=1e-2)
        prim = s.get_primitive_standard_structure()
        self.assertAlmostEqual(prim.lattice.alpha, 105.015053349)
        self.assertAlmostEqual(prim.lattice.beta, 105.015053349)
        self.assertAlmostEqual(prim.lattice.gamma, 118.80658411899999)
        self.assertAlmostEqual(prim.lattice.a, 4.1579321075608791)
        self.assertAlmostEqual(prim.lattice.b, 4.1579321075608791)
        self.assertAlmostEqual(prim.lattice.c, 4.1579321075608791)

        parser = CifParser(os.path.join(PymatgenTest.TEST_FILES_DIR, "orci_1010.cif"))
        structure = parser.get_structures(False)[0]
        s = SpacegroupAnalyzer(structure, symprec=1e-2)
        prim = s.get_primitive_standard_structure()
        self.assertAlmostEqual(prim.lattice.alpha, 134.78923546600001)
        self.assertAlmostEqual(prim.lattice.beta, 105.856239333)
        self.assertAlmostEqual(prim.lattice.gamma, 91.276341676000001)
        self.assertAlmostEqual(prim.lattice.a, 3.8428217771014852)
        self.assertAlmostEqual(prim.lattice.b, 3.8428217771014852)
        self.assertAlmostEqual(prim.lattice.c, 3.8428217771014852)

        parser = CifParser(os.path.join(PymatgenTest.TEST_FILES_DIR, "orcc_1003.cif"))
        structure = parser.get_structures(False)[0]
        s = SpacegroupAnalyzer(structure, symprec=1e-2)
        prim = s.get_primitive_standard_structure()
        self.assertAlmostEqual(prim.lattice.alpha, 90)
        self.assertAlmostEqual(prim.lattice.beta, 90)
        self.assertAlmostEqual(prim.lattice.gamma, 164.985257335)
        self.assertAlmostEqual(prim.lattice.a, 15.854897098324196)
        self.assertAlmostEqual(prim.lattice.b, 15.854897098324196)
        self.assertAlmostEqual(prim.lattice.c, 3.99648651)

        parser = CifParser(os.path.join(PymatgenTest.TEST_FILES_DIR, "orac_632475.cif"))
        structure = parser.get_structures(False)[0]
        s = SpacegroupAnalyzer(structure, symprec=1e-2)
        prim = s.get_primitive_standard_structure()
        self.assertAlmostEqual(prim.lattice.alpha, 90)
        self.assertAlmostEqual(prim.lattice.beta, 90)
        self.assertAlmostEqual(prim.lattice.gamma, 144.40557588533386)
        self.assertAlmostEqual(prim.lattice.a, 5.2005185662155391)
        self.assertAlmostEqual(prim.lattice.b, 5.2005185662155391)
        self.assertAlmostEqual(prim.lattice.c, 3.5372412099999999)

        parser = CifParser(os.path.join(PymatgenTest.TEST_FILES_DIR, "monoc_1028.cif"))
        structure = parser.get_structures(False)[0]
        s = SpacegroupAnalyzer(structure, symprec=1e-2)
        prim = s.get_primitive_standard_structure()
        self.assertAlmostEqual(prim.lattice.alpha, 63.579155761999999)
        self.assertAlmostEqual(prim.lattice.beta, 116.42084423747779)
        self.assertAlmostEqual(prim.lattice.gamma, 148.47965136208569)
        self.assertAlmostEqual(prim.lattice.a, 7.2908007159612325)
        self.assertAlmostEqual(prim.lattice.b, 7.2908007159612325)
        self.assertAlmostEqual(prim.lattice.c, 6.8743926325200002)

        parser = CifParser(os.path.join(PymatgenTest.TEST_FILES_DIR, "hex_1170.cif"))
        structure = parser.get_structures(False)[0]
        s = SpacegroupAnalyzer(structure, symprec=1e-2)
        prim = s.get_primitive_standard_structure()
        self.assertAlmostEqual(prim.lattice.alpha, 90)
        self.assertAlmostEqual(prim.lattice.beta, 90)
        self.assertAlmostEqual(prim.lattice.gamma, 120)
        self.assertAlmostEqual(prim.lattice.a, 3.699919902005897)
        self.assertAlmostEqual(prim.lattice.b, 3.699919902005897)
        self.assertAlmostEqual(prim.lattice.c, 6.9779585500000003)

        parser = CifParser(os.path.join(PymatgenTest.TEST_FILES_DIR, "rhomb_3478_conv.cif"))
        structure = parser.get_structures(False)[0]
        s = SpacegroupAnalyzer(structure, symprec=1e-2)
        prim = s.get_primitive_standard_structure()
        self.assertAlmostEqual(prim.lattice.alpha, 28.049186140546812)
        self.assertAlmostEqual(prim.lattice.beta, 28.049186140546812)
        self.assertAlmostEqual(prim.lattice.gamma, 28.049186140546812)
        self.assertAlmostEqual(prim.lattice.a, 5.9352627428399982)
        self.assertAlmostEqual(prim.lattice.b, 5.9352627428399982)
        self.assertAlmostEqual(prim.lattice.c, 5.9352627428399982)

        parser = CifParser(os.path.join(PymatgenTest.TEST_FILES_DIR, "rhomb_3478_conv.cif"))
        structure = parser.get_structures(False)[0]
        structure.add_site_property("magmom", [1.0] * len(structure))
        s = SpacegroupAnalyzer(structure, symprec=1e-2)
        prim = s.get_primitive_standard_structure(keep_site_properties=True)
        self.assertEqual(prim.site_properties["magmom"], [1.0] * len(prim))

        parser = CifParser(os.path.join(PymatgenTest.TEST_FILES_DIR, "rhomb_3478_conv.cif"))
        structure = parser.get_structures(False)[0]
        structure.add_site_property("magmom", [1.0] * len(structure))
        s = SpacegroupAnalyzer(structure, symprec=1e-2)
        prim = s.get_primitive_standard_structure(keep_site_properties=False)
        self.assertEqual(prim.site_properties.get("magmom", None), None)

    def test_tricky_structure(self):
        # for some reason this structure kills spglib1.9
        # 1.7 can't find symmetry either, but at least doesn't kill python
        s = Structure.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR.tricky_symmetry"))
        sa = SpacegroupAnalyzer(s, 0.1)
        sa.get_space_group_symbol()
        sa.get_space_group_number()
        sa.get_point_group_symbol()
        sa.get_crystal_system()
        sa.get_hall()


class SpacegroupTest(unittest.TestCase):
    def setUp(self):
        p = Poscar.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR"))
        self.structure = p.structure
        self.sg1 = SpacegroupAnalyzer(self.structure, 0.001).get_space_group_operations()

    def test_are_symmetrically_equivalent(self):
        sites1 = [self.structure[i] for i in [0, 1]]
        sites2 = [self.structure[i] for i in [2, 3]]
        self.assertTrue(self.sg1.are_symmetrically_equivalent(sites1, sites2, 1e-3))

        sites1 = [self.structure[i] for i in [0, 1]]
        sites2 = [self.structure[i] for i in [0, 2]]
        self.assertFalse(self.sg1.are_symmetrically_equivalent(sites1, sites2, 1e-3))


H2O2 = Molecule(
    ["O", "O", "H", "H"],
    [
        [0, 0.727403, -0.050147],
        [0, -0.727403, -0.050147],
        [0.83459, 0.897642, 0.401175],
        [-0.83459, -0.897642, 0.401175],
    ],
)

C2H2F2Br2 = Molecule(
    ["C", "C", "F", "Br", "H", "F", "H", "Br"],
    [
        [-0.752000, 0.001000, -0.141000],
        [0.752000, -0.001000, 0.141000],
        [-1.158000, 0.991000, 0.070000],
        [-1.240000, -0.737000, 0.496000],
        [-0.924000, -0.249000, -1.188000],
        [1.158000, -0.991000, -0.070000],
        [0.924000, 0.249000, 1.188000],
        [1.240000, 0.737000, -0.496000],
    ],
)

H2O = Molecule(
    ["H", "O", "H"],
    [[0, 0.780362, -0.456316], [0, 0, 0.114079], [0, -0.780362, -0.456316]],
)

C2H4 = Molecule(
    ["C", "C", "H", "H", "H", "H"],
    [
        [0.0000, 0.0000, 0.6695],
        [0.0000, 0.0000, -0.6695],
        [0.0000, 0.9289, 1.2321],
        [0.0000, -0.9289, 1.2321],
        [0.0000, 0.9289, -1.2321],
        [0.0000, -0.9289, -1.2321],
    ],
)

NH3 = Molecule(
    ["N", "H", "H", "H"],
    [
        [0.0000, 0.0000, 0.0000],
        [0.0000, -0.9377, -0.3816],
        [0.8121, 0.4689, -0.3816],
        [-0.8121, 0.4689, -0.3816],
    ],
)

BF3 = Molecule(
    ["B", "F", "F", "F"],
    [
        [0.0000, 0.0000, 0.0000],
        [0.0000, -0.9377, 0.00],
        [0.8121, 0.4689, 0],
        [-0.8121, 0.4689, 0],
    ],
)

CH4 = Molecule(
    ["C", "H", "H", "H", "H"],
    [
        [0.000000, 0.000000, 0.000000],
        [0.000000, 0.000000, 1.08],
        [1.026719, 0.000000, -0.363000],
        [-0.513360, -0.889165, -0.363000],
        [-0.513360, 0.889165, -0.363000],
    ],
)

PF6 = Molecule(
    ["P", "F", "F", "F", "F", "F", "F"],
    [[0, 0, 0], [0, 0, 1], [0, 0, -1], [0, 1, 0], [0, -1, 0], [1, 0, 0], [-1, 0, 0]],
)


class PointGroupAnalyzerTest(PymatgenTest):
    def test_spherical(self):
        a = PointGroupAnalyzer(CH4)
        self.assertEqual(a.sch_symbol, "Td")
        self.assertEqual(len(a.get_pointgroup()), 24)
        self.assertEqual(a.get_rotational_symmetry_number(), 12)
        a = PointGroupAnalyzer(H2O)
        self.assertEqual(a.get_rotational_symmetry_number(), 2)
        a = PointGroupAnalyzer(PF6)
        self.assertEqual(a.sch_symbol, "Oh")
        self.assertEqual(len(a.get_pointgroup()), 48)
        m = Molecule.from_file(os.path.join(test_dir_mol, "c60.xyz"))
        a = PointGroupAnalyzer(m)
        self.assertEqual(a.sch_symbol, "Ih")

        cube_species = ["C", "C", "C", "C", "C", "C", "C", "C"]
        cube_coords = [
            [0, 0, 0],
            [1, 0, 0],
            [0, 1, 0],
            [1, 1, 0],
            [0, 0, 1],
            [0, 1, 1],
            [1, 0, 1],
            [1, 1, 1],
        ]

        m = Molecule(cube_species, cube_coords)
        a = PointGroupAnalyzer(m, 0.1)
        self.assertEqual(a.sch_symbol, "Oh")

    def test_tricky(self):
        m = Molecule.from_file(os.path.join(test_dir_mol, "dh.xyz"))
        a = PointGroupAnalyzer(m, 0.1)
        self.assertEqual(a.sch_symbol, "D*h")

    def test_linear(self):
        coords = [
            [0.000000, 0.000000, 0.000000],
            [0.000000, 0.000000, 1.08],
            [0, 0.000000, -1.08],
        ]
        mol = Molecule(["C", "H", "H"], coords)
        a = PointGroupAnalyzer(mol)
        self.assertEqual(a.sch_symbol, "D*h")
        mol = Molecule(["C", "H", "N"], coords)
        a = PointGroupAnalyzer(mol)
        self.assertEqual(a.sch_symbol, "C*v")

    def test_asym_top(self):
        coords = [
            [0.000000, 0.000000, 0.000000],
            [0.000000, 0.000000, 1.08],
            [1.026719, 0.000000, -0.363000],
            [-0.513360, -0.889165, -0.363000],
            [-0.513360, 0.889165, -0.363000],
        ]
        mol = Molecule(["C", "H", "F", "Br", "Cl"], coords)
        a = PointGroupAnalyzer(mol)

        self.assertEqual(a.sch_symbol, "C1")
        self.assertEqual(len(a.get_pointgroup()), 1)
        coords = [
            [0.000000, 0.000000, 1.08],
            [1.026719, 0.000000, -0.363000],
            [-0.513360, -0.889165, -0.363000],
            [-0.513360, 0.889165, -0.363000],
        ]
        cs_mol = Molecule(["H", "F", "Cl", "Cl"], coords)
        a = PointGroupAnalyzer(cs_mol)
        self.assertEqual(a.sch_symbol, "Cs")
        self.assertEqual(len(a.get_pointgroup()), 2)
        a = PointGroupAnalyzer(C2H2F2Br2)
        self.assertEqual(a.sch_symbol, "Ci")
        self.assertEqual(len(a.get_pointgroup()), 2)

    def test_cyclic(self):
        a = PointGroupAnalyzer(H2O2)
        self.assertEqual(a.sch_symbol, "C2")
        self.assertEqual(len(a.get_pointgroup()), 2)
        a = PointGroupAnalyzer(H2O)
        self.assertEqual(a.sch_symbol, "C2v")
        self.assertEqual(len(a.get_pointgroup()), 4)
        a = PointGroupAnalyzer(NH3)
        self.assertEqual(a.sch_symbol, "C3v")
        self.assertEqual(len(a.get_pointgroup()), 6)
        cs2 = Molecule.from_file(os.path.join(test_dir_mol, "Carbon_Disulfide.xyz"))
        a = PointGroupAnalyzer(cs2, eigen_tolerance=0.001)
        self.assertEqual(a.sch_symbol, "C2v")

    def test_dihedral(self):
        a = PointGroupAnalyzer(C2H4)
        self.assertEqual(a.sch_symbol, "D2h")
        self.assertEqual(len(a.get_pointgroup()), 8)
        a = PointGroupAnalyzer(BF3)
        self.assertEqual(a.sch_symbol, "D3h")
        self.assertEqual(len(a.get_pointgroup()), 12)
        m = Molecule.from_file(os.path.join(test_dir_mol, "b12h12.xyz"))
        a = PointGroupAnalyzer(m)
        self.assertEqual(a.sch_symbol, "Ih")

    def test_symmetrize_molecule1(self):
        np.random.seed(77)
        distortion = np.random.randn(len(C2H4), 3) / 10
        dist_mol = Molecule(C2H4.species, C2H4.cart_coords + distortion)

        eq = iterative_symmetrize(dist_mol, max_n=100, epsilon=1e-7)
        sym_mol, eq_sets, ops = eq["sym_mol"], eq["eq_sets"], eq["sym_ops"]

        self.assertTrue({0, 1} in eq_sets.values())
        self.assertTrue({2, 3, 4, 5} in eq_sets.values())

        coords = sym_mol.cart_coords
        for i, eq_set in eq_sets.items():
            for j in eq_set:
                _ = np.dot(ops[i][j], coords[i])
                self.assertTrue(np.allclose(np.dot(ops[i][j], coords[i]), coords[j]))

    def test_symmetrize_molecule2(self):
        np.random.seed(77)
        distortion = np.random.randn(len(C2H2F2Br2), 3) / 20
        dist_mol = Molecule(C2H2F2Br2.species, C2H2F2Br2.cart_coords + distortion)
        PA1 = PointGroupAnalyzer(C2H2F2Br2, tolerance=0.1)
        self.assertTrue(PA1.get_pointgroup().sch_symbol == "Ci")
        PA2 = PointGroupAnalyzer(dist_mol, tolerance=0.1)
        self.assertTrue(PA2.get_pointgroup().sch_symbol == "C1")
        eq = iterative_symmetrize(dist_mol, tolerance=0.3)
        PA3 = PointGroupAnalyzer(eq["sym_mol"], tolerance=0.1)
        self.assertTrue(PA3.get_pointgroup().sch_symbol == "Ci")

    def test_get_kpoint_weights(self):
        for name in ["SrTiO3", "LiFePO4", "Graphite"]:
            s = PymatgenTest.get_structure(name)
            a = SpacegroupAnalyzer(s)
            ir_mesh = a.get_ir_reciprocal_mesh((4, 4, 4))
            weights = [i[1] for i in ir_mesh]
            weights = np.array(weights) / sum(weights)
            for i, w in zip(weights, a.get_kpoint_weights([i[0] for i in ir_mesh])):
                self.assertAlmostEqual(i, w)

        for name in ["SrTiO3", "LiFePO4", "Graphite"]:
            s = PymatgenTest.get_structure(name)
            a = SpacegroupAnalyzer(s)
            ir_mesh = a.get_ir_reciprocal_mesh((1, 2, 3))
            weights = [i[1] for i in ir_mesh]
            weights = np.array(weights) / sum(weights)
            for i, w in zip(weights, a.get_kpoint_weights([i[0] for i in ir_mesh])):
                self.assertAlmostEqual(i, w)

        v = Vasprun(os.path.join(PymatgenTest.TEST_FILES_DIR, "vasprun.xml"))
        a = SpacegroupAnalyzer(v.final_structure)
        wts = a.get_kpoint_weights(v.actual_kpoints)

        for w1, w2 in zip(v.actual_kpoints_weights, wts):
            self.assertAlmostEqual(w1, w2)

        kpts = [[0, 0, 0], [0.15, 0.15, 0.15], [0.2, 0.2, 0.2]]
        self.assertRaises(ValueError, a.get_kpoint_weights, kpts)


class FuncTest(unittest.TestCase):
    def test_cluster_sites(self):
        o, c = cluster_sites(CH4, 0.1)
        self.assertEqual(o.specie.symbol, "C")
        self.assertEqual(len(c), 1)
        o, c = cluster_sites(C2H2F2Br2.get_centered_molecule(), 0.1)
        self.assertIsNone(o)
        self.assertEqual(len(c), 4)


if __name__ == "__main__":
    unittest.main()
