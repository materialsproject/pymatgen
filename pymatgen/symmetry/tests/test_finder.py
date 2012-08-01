#!/usr/bin/env python

'''
Created on Mar 9, 2012
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Mar 9, 2012"

import unittest
import os

import numpy as np

from pymatgen.core.sites import PeriodicSite
from pymatgen.io.vaspio.vasp_input import Poscar
from pymatgen.symmetry.finder import SymmetryFinder, get_pointgroup
from pymatgen.io.cifio import CifParser
from pymatgen.core.structure_modifier import StructureEditor

import pymatgen

test_dir = os.path.join(os.path.dirname(os.path.abspath(pymatgen.__file__)),
                        '..', 'test_files')


class SymmetryFinderTest(unittest.TestCase):

    def setUp(self):
        p = Poscar.from_file(os.path.join(test_dir, 'POSCAR'))
        self.structure = p.struct
        self.sg = SymmetryFinder(self.structure, 0.001)
        parser = CifParser(os.path.join(test_dir, 'Li10GeP2S12.cif'))
        self.disordered_structure = parser.get_structures()[0]
        self.disordered_sg = SymmetryFinder(self.disordered_structure, 0.001)
        s = p.struct
        editor = StructureEditor(p.struct)
        site = s[0]
        editor.delete_site(0)
        editor.append_site(site.species_and_occu, site.frac_coords)
        self.sg3 = SymmetryFinder(editor.modified_structure, 0.001)

        parser = CifParser(os.path.join(test_dir, 'Graphite.cif'))
        graphite = parser.get_structures()[0]
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
        self.assertEqual(self.sg.get_pointgroup(), 'mmm')
        self.assertEqual(self.disordered_sg.get_pointgroup(), '4/mmm')

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

    def test_get_symmetrized_structure(self):
        symm_struct = self.sg.get_symmetrized_structure()
        for a in symm_struct.lattice.angles:
            self.assertEqual(a, 90)
        self.assertEqual(len(symm_struct.equivalent_sites), 5)

        symm_struct = self.disordered_sg.get_symmetrized_structure()
        self.assertEqual(len(symm_struct.equivalent_sites), 8)

    def test_get_primitive(self):
        """
        F m -3 m Li2O testing of converting to primitive cell
        """
        self.assertIsNone(self.sg.find_primitive())
        parser = CifParser(os.path.join(test_dir, 'Li2O.cif'))
        structure = parser.get_structures(False)[0]
        s = SymmetryFinder(structure)
        primitive_structure = s.find_primitive()
        self.assertEqual(primitive_structure.formula, "Li2 O1")
        # This isn't what is expected. All the angles should be 60
        self.assertAlmostEqual(primitive_structure.lattice.alpha, 120)
        self.assertAlmostEqual(primitive_structure.lattice.beta, 60)
        self.assertAlmostEqual(primitive_structure.lattice.gamma, 120)
        self.assertAlmostEqual(primitive_structure.lattice.volume,
                               structure.lattice.volume / 4.0)


class HelperFunctionsTest(unittest.TestCase):

    def setUp(self):
        p = Poscar.from_file(os.path.join(test_dir, 'POSCAR'))
        self.sg = SymmetryFinder(p.struct, 0.1)

    def test_get_pointgroup(self):
        (rots, trans) = self.sg.get_symmetry()
        pg = get_pointgroup(rots)
        self.assertEqual(pg[0].strip(), "mmm")

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
