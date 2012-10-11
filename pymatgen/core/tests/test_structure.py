#!/usr/bin/python

import unittest

from pymatgen.core.periodic_table import Element, Specie
from pymatgen.core.structure import Structure, Molecule, Composition, \
    StructureError
from pymatgen.core.lattice import Lattice
import numpy as np
import random


class StructureTest(unittest.TestCase):

    def setUp(self):
        self.si = Element("Si")
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        self.lattice = Lattice([[3.8401979337, 0.00, 0.00],
                                [1.9200989668, 3.3257101909, 0.00],
                                [0.00, -2.2171384943, 3.1355090603]])
        self.struct = Structure(self.lattice, [self.si, self.si], coords)
        self.assertEqual(len(self.struct), 2,
                         "Wrong number of sites in structure!")
        self.assertTrue(self.struct.is_ordered)
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0., 0, 0.0000001])
        self.assertRaises(StructureError, Structure, self.lattice,
                          [self.si, self.si], coords, True)
        self.propertied_structure = Structure(self.lattice, [self.si, self.si],
                                              coords,
                                              site_properties={'magmom':
                                                               [5, -5]})

    def test_volume_and_density(self):
        self.assertAlmostEqual(self.struct.volume, 40.04, 2, "Volume wrong!")
        self.assertAlmostEqual(self.struct.density, 2.33, 2,
                               "Incorrect density")

    def test_specie_init(self):
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        s = Structure(self.lattice, [{Specie('O', -2):1.0},
                                     {Specie('Mg', 2):0.8}], coords)
        self.assertEqual(str(s.composition), 'Mg0.8 O1')

    def test_get_sorted_structure(self):
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        s = Structure(self.lattice, ["O", "Li"], coords,
                      site_properties={'charge': [-2, 1]})
        sorted_s = s.get_sorted_structure()
        self.assertEqual(sorted_s[0].species_and_occu, {Element("Li"): 1})
        self.assertEqual(sorted_s[1].species_and_occu, {Element("O"): 1})
        self.assertEqual(sorted_s[0].charge, 1)
        self.assertEqual(sorted_s[1].charge, -2)

    def test_fractional_occupations(self):
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        s = Structure(self.lattice, [{Element('O'):1.0}, {Element('Mg'):0.8}],
                      coords)
        self.assertEqual(str(s.composition), 'Mg0.8 O1')
        self.assertFalse(s.is_ordered)

    def test_get_distance(self):
        self.assertAlmostEqual(self.struct.get_distance(0, 1), 2.35, 2,
                               "Distance calculated wrongly!")
        pt = [0.9, 0.9, 0.8]
        self.assertAlmostEqual(self.struct[0].distance_from_point(pt),
                               1.50332963784, 2,
                               "Distance calculated wrongly!")

    def test_to_dict(self):
        si = Specie("Si", 4)
        mn = Element("Mn")
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        struct = Structure(self.lattice, [{si:0.5, mn:0.5}, {si:0.5}], coords)
        self.assertIn("lattice", struct.to_dict)
        self.assertIn("sites", struct.to_dict)
        d = self.propertied_structure.to_dict
        self.assertEqual(d['sites'][0]['properties']['magmom'], 5)
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        s = Structure(self.lattice, [{Specie('O', -2,
                                             properties={"spin":3}):1.0},
                                     {Specie('Mg', 2,
                                             properties={"spin":2}):0.8}],
                      coords, site_properties={'magmom': [5, -5]})
        d = s.to_dict
        self.assertEqual(d['sites'][0]['properties']['magmom'], 5)
        self.assertEqual(d['sites'][0]['species'][0]['properties']['spin'], 3)

    def test_from_dict(self):
        d = self.propertied_structure.to_dict
        s = Structure.from_dict(d)
        self.assertEqual(s[0].magmom, 5)
        d = {'lattice': {'a': 3.8401979337, 'volume': 40.044794644251596,
                         'c': 3.8401979337177736, 'b': 3.840198994344244,
                         'matrix': [[3.8401979337, 0.0, 0.0],
                                    [1.9200989668, 3.3257101909, 0.0],
                                    [0.0, -2.2171384943, 3.1355090603]],
                         'alpha': 119.9999908639842, 'beta': 90.0,
                         'gamma': 60.000009137322195},
             'sites': [{'properties': {'magmom': 5}, 'abc': [0.0, 0.0, 0.0],
                        'occu': 1.0, 'species': [{'occu': 1.0,
                                                  'oxidation_state':-2,
                                                  'properties': {'spin': 3},
                                                  'element': 'O'}],
                        'label': 'O2-', 'xyz': [0.0, 0.0, 0.0]},
                       {'properties': {'magmom':-5}, 'abc': [0.75, 0.5, 0.75],
                        'occu': 0.8, 'species': [{'occu': 0.8,
                                                  'oxidation_state': 2,
                                                  'properties': {'spin': 2},
                                                  'element': 'Mg'}],
                        'label': 'Mg2+:0.800',
                        'xyz': [3.8401979336749994, 1.2247250003039056e-06,
                                2.351631795225]}]}
        s = Structure.from_dict(d)
        self.assertEqual(s[0].magmom, 5)
        self.assertEqual(s[0].specie.spin, 3)

    def test_site_properties(self):
        site_props = self.propertied_structure.site_properties
        self.assertEqual(site_props['magmom'], [5, -5])
        self.assertEqual(self.propertied_structure[0].magmom, 5)
        self.assertEqual(self.propertied_structure[1].magmom, -5)

    def test_copy(self):
        new_struct = self.propertied_structure.copy(site_properties={'charge':
                                                                     [2, 3]})
        self.assertEqual(new_struct[0].magmom, 5)
        self.assertEqual(new_struct[1].magmom, -5)
        self.assertEqual(new_struct[0].charge, 2)
        self.assertEqual(new_struct[1].charge, 3)

        coords = list()
        coords.append([0, 0, 0])
        coords.append([0., 0, 0.0000001])

        structure = Structure(self.lattice, ["O", "Si"], coords,
                              site_properties={'magmom': [5, -5]})

        new_struct = structure.copy(site_properties={'charge': [2, 3]},
                                    sanitize=True)
        self.assertEqual(new_struct[0].magmom, -5)
        self.assertEqual(new_struct[1].magmom, 5)
        self.assertEqual(new_struct[0].charge, 3)
        self.assertEqual(new_struct[1].charge, 2)
        self.assertAlmostEqual(new_struct.volume, structure.volume)

    def test_interpolate(self):
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        struct = Structure(self.lattice, [self.si, self.si], coords)
        coords2 = list()
        coords2.append([0, 0, 0])
        coords2.append([0.5, 0.5, 0.5])
        struct2 = Structure(self.struct.lattice, [self.si, self.si], coords2)
        int_s = struct.interpolate(struct2, 10)
        for s in int_s:
            self.assertIsNotNone(s, "Interpolation Failed!")
        self.assertTrue((int_s[1][1].frac_coords == [0.725, 0.5, 0.725]).all())

        badlattice = [[1, 0.00, 0.00], [0, 1, 0.00], [0.00, 0, 1]]
        struct2 = Structure(badlattice, [self.si, self.si], coords2)
        self.assertRaises(ValueError, struct.interpolate, struct2)

        coords2 = list()
        coords2.append([0, 0, 0])
        coords2.append([0.5, 0.5, 0.5])
        struct2 = Structure(self.struct.lattice, [self.si, Element("Fe")],
                            coords2)
        self.assertRaises(ValueError, struct.interpolate, struct2)

    def test_get_all_neighbors_and_get_neighbors(self):
        s = self.struct
        r = random.uniform(3, 6)
        all_nn = s.get_all_neighbors(r)
        for i in range(len(s)):
            self.assertEqual(len(all_nn[i]), len(s.get_neighbors(s[i], r)))

    def test_get_dist_matrix(self):
        ans = [[0., 2.3516318],
               [2.3516318, 0.]]
        self.assertTrue(np.allclose(self.struct.distance_matrix, ans))


class MoleculeTest(unittest.TestCase):

    def setUp(self):
        coords = [[0.000000, 0.000000, 0.000000],
                  [0.000000, 0.000000, 1.089000],
                  [1.026719, 0.000000, -0.363000],
                  [-0.513360, -0.889165, -0.363000],
                  [-0.513360, 0.889165, -0.363000]]
        self.coords = coords
        self.mol = Molecule(["C", "H", "H", "H", "H"], coords)

    def test_get_angle_dihedral(self):
        self.assertAlmostEqual(self.mol.get_angle(1, 0, 2), 109.47122144618737)
        self.assertAlmostEqual(self.mol.get_angle(3, 1, 2), 60.00001388659683)
        self.assertAlmostEqual(self.mol.get_dihedral(0, 1, 2, 3),
                               - 35.26438851071765)

        coords = list()
        coords.append([0, 0, 0])
        coords.append([0, 0, 1])
        coords.append([0, 1, 1])
        coords.append([1, 1, 1])
        self.mol2 = Molecule(["C", "O", "N", "S"], coords)
        self.assertAlmostEqual(self.mol2.get_dihedral(0, 1, 2, 3), -90)

    def test_get_covalent_bonds(self):
        self.assertEqual(len(self.mol.get_covalent_bonds()), 4)

    def test_properties(self):
        self.assertEqual(len(self.mol), 5)
        self.assertTrue(self.mol.is_ordered)
        self.assertEqual(self.mol.formula, "H4 C1")

    def test_repr_str(self):
        ans = """Molecule Summary (H4 C1)
Reduced Formula: H4C
Sites (5)
1 C     0.000000     0.000000     0.000000
2 H     0.000000     0.000000     1.089000
3 H     1.026719     0.000000    -0.363000
4 H    -0.513360    -0.889165    -0.363000
5 H    -0.513360     0.889165    -0.363000"""
        self.assertEqual(str(self.mol), ans)
        ans = """Molecule Summary
Non-periodic Site
xyz        : (0.0000, 0.0000, 0.0000)
element    : C
occupation : 1.00
Non-periodic Site
xyz        : (0.0000, 0.0000, 1.0890)
element    : H
occupation : 1.00
Non-periodic Site
xyz        : (1.0267, 0.0000, -0.3630)
element    : H
occupation : 1.00
Non-periodic Site
xyz        : (-0.5134, -0.8892, -0.3630)
element    : H
occupation : 1.00
Non-periodic Site
xyz        : (-0.5134, 0.8892, -0.3630)
element    : H
occupation : 1.00"""
        self.assertEqual(repr(self.mol), ans)

    def test_site_properties(self):
        propertied_mol = Molecule(["C", "H", "H", "H", "H"], self.coords,
                                  site_properties={'magmom':
                                                   [0.5, -0.5, 1, 2, 3]})
        self.assertEqual(propertied_mol[0].magmom, 0.5)
        self.assertEqual(propertied_mol[1].magmom, -0.5)

    def test_to_from_dict(self):
        propertied_mol = Molecule(["C", "H", "H", "H", "H"], self.coords,
                                  site_properties={'magmom':
                                                   [0.5, -0.5, 1, 2, 3]})
        d = propertied_mol.to_dict
        self.assertEqual(d['sites'][0]['properties']['magmom'], 0.5)
        mol = Molecule.from_dict(d)
        self.assertEqual(mol[0].magmom, 0.5)
        self.assertEqual(mol.formula, "H4 C1")

    def test_get_boxed_structure(self):
        s = self.mol.get_boxed_structure(9, 9, 9)
        self.assertTrue(np.allclose(s[1].frac_coords, [0.000000, 0.000000,
                                                       0.121000]))
        self.assertRaises(ValueError, self.mol.get_boxed_structure, 1, 1, 1)

    def test_get_distance(self):
        self.assertAlmostEqual(self.mol.get_distance(0, 1), 1.089)

    def test_get_neighbors(self):
        nn = self.mol.get_neighbors(self.mol[0], 1)
        self.assertEqual(len(nn), 0)
        nn = self.mol.get_neighbors(self.mol[0], 2)
        self.assertEqual(len(nn), 4)

    def test_get_neighbors_in_shell(self):
        nn = self.mol.get_neighbors_in_shell([0, 0, 0], 0, 1)
        self.assertEqual(len(nn), 1)
        nn = self.mol.get_neighbors_in_shell([0, 0, 0], 1, 0.9)
        self.assertEqual(len(nn), 4)
        nn = self.mol.get_neighbors_in_shell([0, 0, 0], 2, 0.1)
        self.assertEqual(len(nn), 0)

    def test_get_dist_matrix(self):
        ans = [[0.0, 1.089, 1.08899995636, 1.08900040717, 1.08900040717],
               [1.089, 0.0, 1.77832952654, 1.7783298026, 1.7783298026],
               [1.08899995636, 1.77832952654, 0.0, 1.77833003783,
                1.77833003783],
               [1.08900040717, 1.7783298026, 1.77833003783, 0.0, 1.77833],
               [1.08900040717, 1.7783298026, 1.77833003783, 1.77833, 0.0]]
        self.assertTrue(np.allclose(self.mol.distance_matrix, ans))

    def test_break_bond(self):
        (mol1, mol2) = self.mol.break_bond(0, 1)
        self.assertEqual(mol1.formula, "H3 C1")
        self.assertEqual(mol2.formula, "H1")


class CompositionTest(unittest.TestCase):

    def setUp(self):
        self.comp = list()
        self.comp.append(Composition.from_formula("Li3Fe2(PO4)3"))
        self.comp.append(Composition.from_formula("Li3Fe(PO4)O"))
        self.comp.append(Composition.from_formula("LiMn2O4"))
        self.comp.append(Composition.from_formula("Li4O4"))
        self.comp.append(Composition.from_formula("Li3Fe2Mo3O12"))
        self.comp.append(Composition.from_formula("Li3Fe2((PO4)3(CO3)5)2"))
        self.comp.append(Composition.from_formula("Li1.5Si0.5"))

        self.indeterminate_comp = []
        self.indeterminate_comp.append(
            Composition.ranked_compositions_from_indeterminate_formula("Co1",
                                                                       True)
        )
        self.indeterminate_comp.append(
            Composition.ranked_compositions_from_indeterminate_formula("Co1",
                                                                       False)
        )
        self.indeterminate_comp.append(
            Composition.ranked_compositions_from_indeterminate_formula("co2o3")
        )
        self.indeterminate_comp.append(
            Composition.ranked_compositions_from_indeterminate_formula("ncalu")
        )
        self.indeterminate_comp.append(
            Composition.ranked_compositions_from_indeterminate_formula("calun")
        )
        self.indeterminate_comp.append(
            Composition.ranked_compositions_from_indeterminate_formula(
                "liCoo2n (pO4)2")
        )
        self.indeterminate_comp.append(
            Composition.ranked_compositions_from_indeterminate_formula(
                "(co)2 (PO)4")
        )
        self.indeterminate_comp.append(
            Composition.ranked_compositions_from_indeterminate_formula("Fee3"))

    def test_init_(self):
        self.assertRaises(ValueError, Composition, {Element("H"):-0.1})
        f = {'Fe': 4, 'Li': 4, 'O': 16, 'P': 4}
        self.assertEqual("Li4 Fe4 P4 O16", Composition(f).formula)
        f = {None: 4, 'Li': 4, 'O': 16, 'P': 4}
        self.assertRaises(ValueError, Composition, f)
        f = {1: 2, 8: 1}
        self.assertEqual("H2 O1", Composition(f).formula)
        self.assertEqual("Na2 O1", Composition(Na=2, O=1).formula)

    def test_formula(self):
        correct_formulas = ['Li3 Fe2 P3 O12', 'Li3 Fe1 P1 O5', 'Li1 Mn2 O4',
                            'Li4 O4', 'Li3 Fe2 Mo3 O12', 'Li3 Fe2 P6 C10 O54',
                            'Li1.5 Si0.5']
        all_formulas = [c.formula for c in self.comp]
        self.assertEqual(all_formulas, correct_formulas)
        self.assertRaises(ValueError, Composition.from_formula, "(co2)(po4)2")

    def test_mixed_valence(self):
        comp = Composition({"Fe2+": 2, "Fe3+": 4, "Li+": 8})
        self.assertEqual(comp.reduced_formula, "Li4Fe3")
        self.assertEqual(comp.alphabetical_formula, "Fe6 Li8")
        self.assertEqual(comp.formula, "Li8 Fe6")

    def test_indeterminate_formula(self):
        correct_formulas = []
        correct_formulas.append(["Co1"])
        correct_formulas.append(["Co1", "C1 O1"])
        correct_formulas.append(["Co2 O3", "C1 O5"])
        correct_formulas.append(["N1 Ca1 Lu1", "U1 Al1 C1 N1"])
        correct_formulas.append(["N1 Ca1 Lu1", "U1 Al1 C1 N1"])
        correct_formulas.append(["Li1 Co1 P2 N1 O10", "Li1 P2 C1 N1 O11",
                                 "Li1 Co1 Po8 N1 O2", "Li1 Po8 C1 N1 O3"])
        correct_formulas.append(["Co2 P4 O4", "P4 C2 O6", "Co2 Po4",
                                 "Po4 C2 O2"])
        correct_formulas.append([])
        for i, c in enumerate(correct_formulas):
            self.assertEqual([Composition.from_formula(comp) for comp in c],
                             self.indeterminate_comp[i])

    def test_alphabetical_formula(self):
        correct_formulas = ['Fe2 Li3 O12 P3', 'Fe1 Li3 O5 P1', 'Li1 Mn2 O4',
                            'Li4 O4', 'Fe2 Li3 Mo3 O12', 'C10 Fe2 Li3 O54 P6',
                            'Li1.5 Si0.5']
        all_formulas = [c.alphabetical_formula for c in self.comp]
        self.assertEqual(all_formulas, correct_formulas)

    def test_reduced_composition(self):
        correct_reduced_formulas = ['Li3Fe2(PO4)3', 'Li3FePO5', 'LiMn2O4',
                                    'Li2O2', 'Li3Fe2(MoO4)3',
                                    'Li3Fe2P6(C5O27)2', 'Li1.5Si0.5']
        for i in xrange(len(self.comp)):
            self.assertEqual(self.comp[i]
                             .get_reduced_composition_and_factor()[0],
                             Composition
                             .from_formula(correct_reduced_formulas[i]))

    def test_reduced_formula(self):
        correct_reduced_formulas = ['Li3Fe2(PO4)3', 'Li3FePO5', 'LiMn2O4',
                                    'Li2O2', 'Li3Fe2(MoO4)3',
                                    'Li3Fe2P6(C5O27)2', 'Li1.5Si0.5']
        all_formulas = [c.reduced_formula for c in self.comp]
        self.assertEqual(all_formulas, correct_reduced_formulas)

    def test_num_atoms(self):
        correct_num_atoms = [20, 10, 7, 8, 20, 75, 2]
        all_natoms = [c.num_atoms for c in self.comp]
        self.assertEqual(all_natoms, correct_num_atoms)

    def test_weight(self):
        correct_weights = [417.427086, 187.63876199999999, 180.81469, 91.7616,
                           612.3258, 1302.430172, 24.454250000000002]
        all_weights = [c.weight for c in self.comp]
        self.assertAlmostEqual(all_weights, correct_weights, 5)

    def test_get_atomic_fraction(self):
        correct_at_frac = {"Li": 0.15, "Fe": 0.1, "P": 0.15, "O": 0.6}
        for el in ["Li", "Fe", "P", "O"]:
            self.assertEqual(self.comp[0].get_atomic_fraction(Element(el)),
                             correct_at_frac[el],
                             "Wrong computed atomic fractions")
        self.assertEqual(self.comp[0].get_atomic_fraction(Element("S")), 0,
                         "Wrong computed atomic fractions")

    def test_anonymized_formula(self):
        expected_formulas = ['A2B3C3D12', 'ABC3D5', 'AB2C4', 'A2B2',
                             'A2B3C3D12', 'A2B3C6D10E54', 'A0.5B1.5']
        for i in xrange(len(self.comp)):
            self.assertEqual(self.comp[i].anonymized_formula,
                             expected_formulas[i])

    def test_get_wt_fraction(self):
        correct_wt_frac = {"Li": 0.0498841610868, "Fe": 0.267567687258,
                           "P": 0.222604831158, "O": 0.459943320496}
        for el in ["Li", "Fe", "P", "O"]:
            self.assertAlmostEqual(correct_wt_frac[el],
                                   self.comp[0].get_wt_fraction(Element(el)),
                                   5, "Wrong computed weight fraction")
        self.assertEqual(self.comp[0].get_wt_fraction(Element("S")), 0,
                         "Wrong computed weight fractions")

    def test_from_dict(self):
        sym_dict = {"Fe": 6, "O": 8}
        self.assertEqual(Composition.from_dict(sym_dict).reduced_formula,
                         "Fe3O4",
                         "Creation form sym_amount dictionary failed!")

    def test_to_dict(self):
        c = Composition.from_dict({'Fe': 4, 'O': 6})
        d = c.to_dict
        correct_dict = {'Fe': 4.0, 'O': 6.0}
        self.assertEqual(d['Fe'], correct_dict['Fe'])
        self.assertEqual(d['O'], correct_dict['O'])
        correct_dict = {'Fe': 2.0, 'O': 3.0}
        d = c.to_reduced_dict
        self.assertEqual(d['Fe'], correct_dict['Fe'])
        self.assertEqual(d['O'], correct_dict['O'])

    def test_add(self):
        self.assertEqual((self.comp[0] + self.comp[2]).formula,
                         "Li4 Mn2 Fe2 P3 O16",
                         "Incorrect composition after addition!")
        self.assertEqual((self.comp[3] + {"Fe": 4, "O": 4}).formula,
                         "Li4 Fe4 O8", "Incorrect composition after addition!")

    def test_sub(self):
        self.assertEqual((self.comp[0]
                          - Composition.from_formula("Li2O")).formula,
                         "Li1 Fe2 P3 O11",
                         "Incorrect composition after addition!")
        self.assertEqual((self.comp[0] - {"Fe": 2, "O": 3}).formula,
                         "Li3 P3 O9")

    def test_mul(self):
        self.assertEqual((self.comp[0] * 4).formula, "Li12 Fe8 P12 O48",
                         "Incorrect composition after addition!")

    def test_equals(self):
        random_z = random.randint(1, 92)
        fixed_el = Element.from_Z(random_z)
        other_z = random.randint(1, 92)
        while other_z == random_z:
            other_z = random.randint(1, 92)
        comp1 = Composition({fixed_el: 1, Element.from_Z(other_z): 0})
        other_z = random.randint(1, 92)
        while other_z == random_z:
            other_z = random.randint(1, 92)
        comp2 = Composition({fixed_el: 1, Element.from_Z(other_z): 0})
        self.assertEqual(comp1, comp2,
                         "Composition equality test failed. " +
                         "%s should be equal to %s" % (comp1.formula,
                                                       comp2.formula))
        self.assertEqual(comp1.__hash__(), comp2.__hash__(),
                         "Hashcode equality test failed!")

    def test_equality(self):
        self.assertTrue(self.comp[0].__eq__(self.comp[0]))
        self.assertFalse(self.comp[0].__eq__(self.comp[1]))
        self.assertFalse(self.comp[0].__ne__(self.comp[0]))
        self.assertTrue(self.comp[0].__ne__(self.comp[1]))

    def test_get_fractional_composition(self):
        for c in self.comp:
            self.assertAlmostEqual(c.get_fractional_composition().num_atoms, 1)

if __name__ == '__main__':
    unittest.main()
