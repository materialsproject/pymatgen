#!/usr/bin/env python

"""
Created on Mar 9, 2012
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Mar 9, 2012"

import unittest
import os

import numpy as np

from pymatgen.core.sites import PeriodicSite
from pymatgen.io.vaspio.vasp_input import Poscar
from pymatgen.symmetry.finder import SymmetryFinder
from pymatgen.io.cifio import CifParser
from pymatgen.io.vaspio.vasp_output import Vasprun

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class SymmetryFinderTest(unittest.TestCase):

    def setUp(self):
        p = Poscar.from_file(os.path.join(test_dir, 'POSCAR'))
        self.structure = p.structure
        self.sg = SymmetryFinder(self.structure, 0.001)
        parser = CifParser(os.path.join(test_dir, 'Li10GeP2S12.cif'))
        self.disordered_structure = parser.get_structures()[0]
        self.disordered_sg = SymmetryFinder(self.disordered_structure, 0.001)
        s = p.structure.copy()
        site = s[0]
        s.remove(0)
        s.append(site.species_and_occu, site.frac_coords)
        self.sg3 = SymmetryFinder(s, 0.001)
        parser = CifParser(os.path.join(test_dir, 'Graphite.cif'))
        graphite = parser.get_structures()[0]
        graphite.add_site_property("magmom", [0.1] * len(graphite))
        self.sg4 = SymmetryFinder(graphite, 0.001)

    def test_get_space_symbol(self):
        self.assertEqual(self.sg.get_spacegroup_symbol(), "Pnma")
        self.assertEqual(self.disordered_sg.get_spacegroup_symbol(),
                         "P4_2/nmc")
        self.assertEqual(self.sg3.get_spacegroup_symbol(), "Pnma")
        self.assertEqual(self.sg4.get_spacegroup_symbol(), "R-3m")

    def test_get_space_number(self):
        self.assertEqual(self.sg.get_spacegroup_number(), 62)
        self.assertEqual(self.disordered_sg.get_spacegroup_number(), 137)
        self.assertEqual(self.sg4.get_spacegroup_number(), 166)

    def test_get_hall(self):
        self.assertEqual(self.sg.get_hall(), '-P 2ac 2n')
        self.assertEqual(self.disordered_sg.get_hall(), 'P 4n 2n -1n')

    def test_get_pointgroup(self):
        self.assertEqual(self.sg.get_point_group(), 'mmm')
        self.assertEqual(self.disordered_sg.get_point_group(), '4/mmm')

    def test_get_symmetry_dataset(self):
        ds = self.sg.get_symmetry_dataset()
        self.assertEqual(ds['international'], 'Pnma')

    def test_get_crystal_system(self):
        crystal_system = self.sg.get_crystal_system()
        self.assertEqual('orthorhombic', crystal_system)
        self.assertEqual('tetragonal', self.disordered_sg.get_crystal_system())

    def test_get_symmetry_operations(self):
        fracsymmops = self.sg.get_symmetry_operations()
        symmops = self.sg.get_symmetry_operations(True)
        self.assertEqual(len(symmops), 8)
        latt = self.structure.lattice
        for fop, op in zip(fracsymmops, symmops):
            for site in self.structure:
                newfrac = fop.operate(site.frac_coords)
                newcart = op.operate(site.coords)
                self.assertTrue(np.allclose(latt.get_fractional_coords(newcart),
                                            newfrac))
                found = False
                newsite = PeriodicSite(site.species_and_occu, newcart, latt,
                                       coords_are_cartesian=True)
                for testsite in self.structure:
                    if newsite.is_periodic_image(testsite, 1e-3):
                        found = True
                        break
                self.assertTrue(found)

    def test_get_refined_structure(self):
        for a in self.sg.get_refined_structure().lattice.angles:
            self.assertEqual(a, 90)
        refined = self.disordered_sg.get_refined_structure()
        for a in refined.lattice.angles:
            self.assertEqual(a, 90)
        self.assertEqual(refined.lattice.a, refined.lattice.b)
        parser = CifParser(os.path.join(test_dir, 'Li2O.cif'))
        s = parser.get_structures()[0]
        sg = SymmetryFinder(s, 0.001)
        self.assertEqual(sg.get_refined_structure().num_sites, 4 * s.num_sites)

    def test_get_symmetrized_structure(self):
        symm_struct = self.sg.get_symmetrized_structure()
        for a in symm_struct.lattice.angles:
            self.assertEqual(a, 90)
        self.assertEqual(len(symm_struct.equivalent_sites), 5)

        symm_struct = self.disordered_sg.get_symmetrized_structure()
        self.assertEqual(len(symm_struct.equivalent_sites), 8)
        self.assertEqual(map(len, symm_struct.equivalent_sites), [16,4,8,4,2,8,8,8])
        s1 = symm_struct.equivalent_sites[1][1]
        s2 = symm_struct[symm_struct.equivalent_indices[1][1]]
        self.assertEqual(s1, s2)
        self.assertEqual(self.sg4.get_symmetrized_structure()[0].magmom, 0.1)

    def test_find_primitive(self):
        """
        F m -3 m Li2O testing of converting to primitive cell
        """
        parser = CifParser(os.path.join(test_dir, 'Li2O.cif'))
        structure = parser.get_structures(False)[0]
        s = SymmetryFinder(structure)
        primitive_structure = s.find_primitive()
        self.assertEqual(primitive_structure.formula, "Li2 O1")
        # This isn't what is expected. All the angles should be 60
        self.assertAlmostEqual(primitive_structure.lattice.alpha, 60)
        self.assertAlmostEqual(primitive_structure.lattice.beta, 60)
        self.assertAlmostEqual(primitive_structure.lattice.gamma, 60)
        self.assertAlmostEqual(primitive_structure.lattice.volume,
                               structure.lattice.volume / 4.0)

    def test_get_ir_kpoints(self):
        v = Vasprun(os.path.join(test_dir, 'vasprun_Si_bands.xml'))
        finder = SymmetryFinder(v.final_structure)
        actual_kpts = v.actual_kpoints
        mapping = finder.get_ir_kpoints_mapping(actual_kpts)
        irr_kpts = finder.get_ir_kpoints(actual_kpts)
        self.assertLess(len(irr_kpts), len(actual_kpts))
        self.assertEqual(len(irr_kpts), len(set(mapping)))

    def test_get_ir_reciprocal_mesh(self):
        grid=self.sg.get_ir_reciprocal_mesh()
        self.assertEquals(len(grid), 216)
        self.assertAlmostEquals(grid[1][0][0], 0.1)
        self.assertAlmostEquals(grid[1][0][1], 0.0)
        self.assertAlmostEquals(grid[1][0][2], 0.0)
        self.assertEquals(grid[1][1], 2)

    def test_get_conventional_standard_structure(self):
        parser = CifParser(os.path.join(test_dir, 'bcc_1927.cif'))
        structure = parser.get_structures(False)[0]
        s = SymmetryFinder(structure, symprec=1e-2)
        conv = s.get_conventional_standard_structure()
        self.assertAlmostEqual(conv.lattice.alpha, 90)
        self.assertAlmostEqual(conv.lattice.beta, 90)
        self.assertAlmostEqual(conv.lattice.gamma, 90)
        self.assertAlmostEqual(conv.lattice.a, 9.1980270633769461)
        self.assertAlmostEqual(conv.lattice.b, 9.1980270633769461)
        self.assertAlmostEqual(conv.lattice.c, 9.1980270633769461)

        parser = CifParser(os.path.join(test_dir, 'btet_1915.cif'))
        structure = parser.get_structures(False)[0]
        s = SymmetryFinder(structure, symprec=1e-2)
        conv = s.get_conventional_standard_structure()
        self.assertAlmostEqual(conv.lattice.alpha, 90)
        self.assertAlmostEqual(conv.lattice.beta, 90)
        self.assertAlmostEqual(conv.lattice.gamma, 90)
        self.assertAlmostEqual(conv.lattice.a, 5.0615106678044235)
        self.assertAlmostEqual(conv.lattice.b, 5.0615106678044235)
        self.assertAlmostEqual(conv.lattice.c, 4.2327080177761687)

        parser = CifParser(os.path.join(test_dir, 'orci_1010.cif'))
        structure = parser.get_structures(False)[0]
        s = SymmetryFinder(structure, symprec=1e-2)
        conv = s.get_conventional_standard_structure()
        self.assertAlmostEqual(conv.lattice.alpha, 90)
        self.assertAlmostEqual(conv.lattice.beta, 90)
        self.assertAlmostEqual(conv.lattice.gamma, 90)
        self.assertAlmostEqual(conv.lattice.a, 2.9542233922299999)
        self.assertAlmostEqual(conv.lattice.b, 4.6330325651443296)
        self.assertAlmostEqual(conv.lattice.c, 5.373703587040775)

        parser = CifParser(os.path.join(test_dir, 'orcc_1003.cif'))
        structure = parser.get_structures(False)[0]
        s = SymmetryFinder(structure, symprec=1e-2)
        conv = s.get_conventional_standard_structure()
        self.assertAlmostEqual(conv.lattice.alpha, 90)
        self.assertAlmostEqual(conv.lattice.beta, 90)
        self.assertAlmostEqual(conv.lattice.gamma, 90)
        self.assertAlmostEqual(conv.lattice.a, 4.1430033493799998)
        self.assertAlmostEqual(conv.lattice.b, 31.437979757624728)
        self.assertAlmostEqual(conv.lattice.c, 3.99648651)

        parser = CifParser(os.path.join(test_dir, 'monoc_1028.cif'))
        structure = parser.get_structures(False)[0]
        s = SymmetryFinder(structure, symprec=1e-2)
        conv = s.get_conventional_standard_structure()
        self.assertAlmostEqual(conv.lattice.alpha, 62.461675798000002)
        self.assertAlmostEqual(conv.lattice.beta, 90)
        self.assertAlmostEqual(conv.lattice.gamma, 90)
        self.assertAlmostEqual(conv.lattice.a, 3.9605285073099998)
        self.assertAlmostEqual(conv.lattice.b, 14.033435583000625)
        self.assertAlmostEqual(conv.lattice.c, 6.8743926325200002)

        parser = CifParser(os.path.join(test_dir, 'rhomb_1170.cif'))
        structure = parser.get_structures(False)[0]
        s = SymmetryFinder(structure, symprec=1e-2)
        conv = s.get_conventional_standard_structure()
        self.assertAlmostEqual(conv.lattice.alpha, 90)
        self.assertAlmostEqual(conv.lattice.beta, 90)
        self.assertAlmostEqual(conv.lattice.gamma, 120)
        self.assertAlmostEqual(conv.lattice.a, 3.699919902005897)
        self.assertAlmostEqual(conv.lattice.b, 3.699919902005897)
        self.assertAlmostEqual(conv.lattice.c, 6.9779585500000003)

    def test_get_primitive_standard_structure(self):
        parser = CifParser(os.path.join(test_dir, 'bcc_1927.cif'))
        structure = parser.get_structures(False)[0]
        s = SymmetryFinder(structure, symprec=1e-2)
        prim = s.get_primitive_standard_structure()
        self.assertAlmostEqual(prim.lattice.alpha, 109.47122063400001)
        self.assertAlmostEqual(prim.lattice.beta, 109.47122063400001)
        self.assertAlmostEqual(prim.lattice.gamma, 109.47122063400001)
        self.assertAlmostEqual(prim.lattice.a, 7.9657251015812145)
        self.assertAlmostEqual(prim.lattice.b, 7.9657251015812145)
        self.assertAlmostEqual(prim.lattice.c, 7.9657251015812145)

        parser = CifParser(os.path.join(test_dir, 'btet_1915.cif'))
        structure = parser.get_structures(False)[0]
        s = SymmetryFinder(structure, symprec=1e-2)
        prim = s.get_primitive_standard_structure()
        self.assertAlmostEqual(prim.lattice.alpha, 105.015053349)
        self.assertAlmostEqual(prim.lattice.beta, 105.015053349)
        self.assertAlmostEqual(prim.lattice.gamma, 118.80658411899999)
        self.assertAlmostEqual(prim.lattice.a, 4.1579321075608791)
        self.assertAlmostEqual(prim.lattice.b, 4.1579321075608791)
        self.assertAlmostEqual(prim.lattice.c, 4.1579321075608791)

        parser = CifParser(os.path.join(test_dir, 'orci_1010.cif'))
        structure = parser.get_structures(False)[0]
        s = SymmetryFinder(structure, symprec=1e-2)
        prim = s.get_primitive_standard_structure()
        self.assertAlmostEqual(prim.lattice.alpha, 134.78923546600001)
        self.assertAlmostEqual(prim.lattice.beta, 105.856239333)
        self.assertAlmostEqual(prim.lattice.gamma, 91.276341676000001)
        self.assertAlmostEqual(prim.lattice.a, 3.8428217771014852)
        self.assertAlmostEqual(prim.lattice.b, 3.8428217771014852)
        self.assertAlmostEqual(prim.lattice.c, 3.8428217771014852)

        parser = CifParser(os.path.join(test_dir, 'orcc_1003.cif'))
        structure = parser.get_structures(False)[0]
        s = SymmetryFinder(structure, symprec=1e-2)
        prim = s.get_primitive_standard_structure()
        self.assertAlmostEqual(prim.lattice.alpha, 90)
        self.assertAlmostEqual(prim.lattice.beta, 90)
        self.assertAlmostEqual(prim.lattice.gamma, 164.985257335)
        self.assertAlmostEqual(prim.lattice.a, 15.854897098324196)
        self.assertAlmostEqual(prim.lattice.b, 15.854897098324196)
        self.assertAlmostEqual(prim.lattice.c, 3.99648651)

        parser = CifParser(os.path.join(test_dir, 'monoc_1028.cif'))
        structure = parser.get_structures(False)[0]
        s = SymmetryFinder(structure, symprec=1e-2)
        prim = s.get_primitive_standard_structure()
        self.assertAlmostEqual(prim.lattice.alpha, 63.579155761999999)
        self.assertAlmostEqual(prim.lattice.beta, 63.579155761999999)
        self.assertAlmostEqual(prim.lattice.gamma, 31.520348638000002)
        self.assertAlmostEqual(prim.lattice.a, 7.2908007159612325)
        self.assertAlmostEqual(prim.lattice.b, 7.2908007159612325)
        self.assertAlmostEqual(prim.lattice.c, 6.8743926325200002)

        parser = CifParser(os.path.join(test_dir, 'rhomb_1170.cif'))
        structure = parser.get_structures(False)[0]
        s = SymmetryFinder(structure, symprec=1e-2)
        prim = s.get_primitive_standard_structure()
        self.assertAlmostEqual(prim.lattice.alpha, 90)
        self.assertAlmostEqual(prim.lattice.beta, 90)
        self.assertAlmostEqual(prim.lattice.gamma, 120)
        self.assertAlmostEqual(prim.lattice.a, 3.699919902005897)
        self.assertAlmostEqual(prim.lattice.b, 3.699919902005897)
        self.assertAlmostEqual(prim.lattice.c, 6.9779585500000003)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
