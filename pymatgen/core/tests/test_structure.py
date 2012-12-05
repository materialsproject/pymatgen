#!/usr/bin/python

import unittest

from pymatgen.core.periodic_table import Element, Specie
from pymatgen.core.structure import Structure, Molecule, StructureError
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
        self.assertEqual(str(s.composition), 'Mg2+0.8 O2-1')

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
Site: C (0.0000, 0.0000, 0.0000)
Site: H (0.0000, 0.0000, 1.0890)
Site: H (1.0267, 0.0000, -0.3630)
Site: H (-0.5134, -0.8892, -0.3630)
Site: H (-0.5134, 0.8892, -0.3630)"""
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


if __name__ == '__main__':
    unittest.main()
