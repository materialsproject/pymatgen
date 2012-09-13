#!/usr/bin/python

import unittest

from pymatgen.core.periodic_table import Element, Specie
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure, Molecule
from pymatgen.core.structure_modifier import StructureEditor, SupercellMaker, \
    MoleculeEditor
import numpy as np


class StructureEditorTest(unittest.TestCase):

    def setUp(self):
        self.si = Element("Si")
        self.fe = Element("Fe")
        self.ge = Element("Ge")
        coords = list()
        coords.append(np.array([0, 0, 0]))
        coords.append(np.array([0.75, 0.5, 0.75]))
        lattice = Lattice.cubic(10)
        s = Structure(lattice, ["Si", "Fe"], coords)
        self.modifier = StructureEditor(s)
        coords = [[0.000000, 0.000000, 0.000000],
                  [0.000000, 0.000000, 1.089000],
                  [1.026719, 0.000000, -0.363000],
                  [-0.513360, -0.889165, -0.363000],
                  [-0.513360, 0.889165, -0.363000]]

    def test_translate_sites(self):
        self.modifier.translate_sites([0, 1], [0.5, 0.5, 0.5],
                                      frac_coords=True)
        self.assertTrue(np.array_equal(self.modifier.modified_structure
                                       .frac_coords[0],
                                       np.array([0.5, 0.5, 0.5])))

        self.modifier.translate_sites([0], [0.5, 0.5, 0.5], frac_coords=False)
        self.assertTrue(np.array_equal(self.modifier.modified_structure
                                       .cart_coords[0],
                                       np.array([5.5, 5.5, 5.5])))

    def test_append_site(self):
        self.modifier.append_site(self.si, [0, 0.5, 0])
        self.assertEqual(self.modifier.modified_structure.formula, "Fe1 Si2",
                         "Wrong formula!")
        self.assertRaises(ValueError, self.modifier.append_site, self.si,
                          np.array([0, 0.5, 0]))

    def test_modified_structure(self):
        self.modifier.insert_site(1, self.si, [0, 0.25, 0])
        self.assertEqual(self.modifier.modified_structure.formula, "Fe1 Si2",
                         "Wrong formula!")

        self.modifier.delete_site(0)
        self.assertEqual(self.modifier.modified_structure.formula, "Fe1 Si1",
                         "Wrong formula!")

        self.modifier.replace_site(0, self.ge)
        self.assertEqual(self.modifier.modified_structure.formula, "Fe1 Ge1",
                         "Wrong formula!")

        self.modifier.append_site(self.si, [0, 0.75, 0])
        self.modifier.replace_species({self.si: self.ge})
        self.assertEqual(self.modifier.modified_structure.formula, "Fe1 Ge2",
                         "Wrong formula!")

        self.modifier.replace_species({self.ge: {self.ge: 0.5, self.si: 0.5}})
        self.assertEqual(self.modifier.modified_structure.formula,
                         "Fe1 Si1 Ge1", "Wrong formula!")

        #this should change the .5Si .5Ge sites to .75Si .25Ge
        self.modifier.replace_species({self.ge: {self.ge: 0.5, self.si: 0.5}})
        self.assertEqual(self.modifier.modified_structure.formula,
                         "Fe1 Si1.5 Ge0.5", "Wrong formula!")

        d = 0.1
        pre_perturbation_sites = self.modifier.modified_structure.sites
        self.modifier.perturb_structure(distance=d)
        post_perturbation_sites = self.modifier.modified_structure.sites

        for i, x in enumerate(pre_perturbation_sites):
            self.assertAlmostEqual(x.distance(post_perturbation_sites[i]), d,
                                   3, "Bad perturbation distance")

    def test_add_site_property(self):
        self.modifier.add_site_property("charge", [4.1, 5])
        s = self.modifier.modified_structure
        self.assertEqual(s[0].charge, 4.1)
        self.assertEqual(s[1].charge, 5)

        #test adding multiple properties.
        mod2 = StructureEditor(s)
        mod2.add_site_property("magmom", [3, 2])
        s = mod2.modified_structure
        self.assertEqual(s[0].charge, 4.1)
        self.assertEqual(s[0].magmom, 3)

    def test_add_oxidation_states(self):
        si = Element("Si")
        fe = Element("Fe")
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        lattice = Lattice.cubic(10)
        s = Structure(lattice, [si, fe], coords)
        oxidation_states = {"Fe": 2, "Si":-4}
        mod = StructureEditor(s)
        mod.add_oxidation_state_by_element(oxidation_states)
        mod_s = mod.modified_structure
        for site in mod_s:
            for k in site.species_and_occu.keys():
                self.assertEqual(k.oxi_state, oxidation_states[k.symbol],
                                 "Wrong oxidation state assigned!")
        oxidation_states = {"Fe": 2}
        self.assertRaises(ValueError, mod.add_oxidation_state_by_element,
                          oxidation_states)
        mod.add_oxidation_state_by_site([2, -4])
        mod_s = mod.modified_structure
        self.assertRaises(ValueError, mod.add_oxidation_state_by_site,
                          [1])

    def test_remove_oxidation_states(self):
        co_elem = Element("Co")
        o_elem = Element("O")
        co_specie = Specie("Co", 2)
        o_specie = Specie("O", -2)
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        lattice = Lattice.cubic(10)
        s_elem = Structure(lattice, [co_elem, o_elem], coords)
        s_specie = Structure(lattice, [co_specie, o_specie], coords)
        mod = StructureEditor(s_specie)
        mod.remove_oxidation_states()
        mod_s = mod.modified_structure
        self.assertEqual(s_elem, mod_s, "Oxidation state remover failed")


class SupercellMakerTest(unittest.TestCase):

    def setUp(self):
        si = Element("Si")
        fe = Element("Fe")
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        lattice = Lattice.cubic(10)
        s = Structure(lattice, [si, fe], coords)
        self.mod = SupercellMaker(s, [[1, 1, 0], [-1, 1, 0], [0, 0, 2]])

    def test_modified_structure(self):
        self.assertEquals(self.mod.modified_structure.formula, "Fe4 Si4",
                          "Wrong formula!")


class MoleculeEditorTest(unittest.TestCase):

    def setUp(self):
        coords = [[0.000000, 0.000000, 0.000000],
                  [0.000000, 0.000000, 1.089000],
                  [1.026719, 0.000000, -0.363000],
                  [-0.513360, -0.889165, -0.363000],
                  [-0.513360, 0.889165, -0.363000]]
        mol = Molecule(["C", "H", "H", "H", "H"], coords)
        self.modifier = MoleculeEditor(mol)

    def test_translate_sites(self):
        self.modifier.translate_sites([0, 1], [0.5, 0.5, 0.5])
        self.assertTrue(np.array_equal(self.modifier.modified_structure
                                       .cart_coords[0],
                                       np.array([0.5, 0.5, 0.5])))

    def test_append_site(self):
        self.modifier.append_site("Si", [0, 0.5, 0])
        self.assertEqual(self.modifier.modified_structure.formula, "Si1 H4 C1",
                         "Wrong formula!")
        self.assertRaises(ValueError, self.modifier.append_site, Element("Si"),
                          np.array([0, 0.5, 0]))

    def test_modified_structure(self):
        self.modifier.insert_site(1, "Si", [0, 0.25, 0])
        self.assertEqual(self.modifier.modified_structure.formula, "Si1 H4 C1", "Wrong formula!")

        self.modifier.delete_site(0)
        self.assertEqual(self.modifier.modified_structure.formula, "Si1 H4", "Wrong formula!")

        self.modifier.replace_site(0, "Ge")
        self.assertEqual(self.modifier.modified_structure.formula, "Ge1 H4", "Wrong formula!")

        self.modifier.append_site("Si", [0, 0.75, 0])
        self.modifier.replace_species({Element("Si"): Element("Ge")})

        self.assertEqual(self.modifier.modified_structure.formula, "Ge2 H4", "Wrong formula!")

        self.modifier.replace_species({Element("Ge"): {Element("Ge"):0.5, Element("Si"):0.5}})
        self.assertEqual(self.modifier.modified_structure.formula, "Si1 Ge1 H4", "Wrong formula!")

        #this should change the .5Si .5Ge sites to .75Si .25Ge
        self.modifier.replace_species({Element("Ge"): {Element("Ge"):0.5, Element("Si"):0.5}})
        self.assertEqual(self.modifier.modified_structure.formula, "Si1.5 Ge0.5 H4", "Wrong formula!")

        d = 0.1
        pre_perturbation_sites = self.modifier.modified_structure.sites
        self.modifier.perturb_structure(distance=d)
        post_perturbation_sites = self.modifier.modified_structure.sites

        for i, x in enumerate(pre_perturbation_sites):
            self.assertAlmostEqual(x.distance(post_perturbation_sites[i]), d, 3, "Bad perturbation distance")

    def test_add_site_property(self):
        self.modifier.add_site_property("charge", [4.1, -2, -2, -2, -2])
        s = self.modifier.modified_structure
        self.assertEqual(s[0].charge, 4.1)
        self.assertEqual(s[1].charge, -2)

        #test adding multiple properties.
        mod2 = MoleculeEditor(s)
        mod2.add_site_property("magmom", [3, 2, 2, 2, 2])
        s = mod2.modified_structure
        self.assertEqual(s[0].charge, 4.1)
        self.assertEqual(s[0].magmom, 3)


if __name__ == '__main__':
    unittest.main()

