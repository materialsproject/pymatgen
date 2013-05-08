#!/usr/bin/python

import unittest

from pymatgen.core.periodic_table import Element, Specie
from pymatgen.core.composition import Composition
from pymatgen.core.operations import SymmOp
from pymatgen.core.structure import IStructure, Structure, IMolecule, \
    StructureError, Molecule
from pymatgen.core.lattice import Lattice
import numpy as np
import random


class IStructureTest(unittest.TestCase):

    def setUp(self):
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        self.lattice = Lattice([[3.8401979337, 0.00, 0.00],
                                [1.9200989668, 3.3257101909, 0.00],
                                [0.00, -2.2171384943, 3.1355090603]])
        self.struct = IStructure(self.lattice, ["Si"] * 2, coords)
        self.assertEqual(len(self.struct), 2,
                         "Wrong number of sites in structure!")
        self.assertTrue(self.struct.is_ordered)
        self.assertTrue(self.struct.ntypat==1)
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0., 0, 0.0000001])
        self.assertRaises(StructureError, IStructure, self.lattice,
                          ["Si"] * 2, coords, True)
        self.propertied_structure = IStructure(
            self.lattice, ["Si"] * 2, coords,
            site_properties={'magmom': [5, -5]})

    def test_bad_structure(self):
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        coords.append([0.75, 0.5, 0.75])
        self.assertRaises(StructureError, IStructure, self.lattice,
                          ["Si"] * 3, coords, validate_proximity=True)

    def test_volume_and_density(self):
        self.assertAlmostEqual(self.struct.volume, 40.04, 2, "Volume wrong!")
        self.assertAlmostEqual(self.struct.density, 2.33, 2,
                               "Incorrect density")

    def test_specie_init(self):
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        s = IStructure(self.lattice, [{Specie('O', -2): 1.0},
                                      {Specie('Mg', 2): 0.8}], coords)
        self.assertEqual(str(s.composition), 'Mg2+0.8 O2-1')

    def test_get_sorted_structure(self):
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        s = IStructure(self.lattice, ["O", "Li"], coords,
                       site_properties={'charge': [-2, 1]})
        sorted_s = s.get_sorted_structure()
        self.assertEqual(sorted_s[0].species_and_occu, Composition("Li"))
        self.assertEqual(sorted_s[1].species_and_occu, Composition("O"))
        self.assertEqual(sorted_s[0].charge, 1)
        self.assertEqual(sorted_s[1].charge, -2)

    def test_fractional_occupations(self):
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        s = IStructure(self.lattice, [{'O': 1.0}, {'Mg': 0.8}],
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
        struct = IStructure(self.lattice, [{si: 0.5, mn: 0.5}, {si: 0.5}],
                            coords)
        self.assertIn("lattice", struct.to_dict)
        self.assertIn("sites", struct.to_dict)
        d = self.propertied_structure.to_dict
        self.assertEqual(d['sites'][0]['properties']['magmom'], 5)
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        s = IStructure(self.lattice, [{Specie('O', -2,
                                              properties={"spin": 3}): 1.0},
                                      {Specie('Mg', 2,
                                              properties={"spin": 2}): 0.8}],
                       coords, site_properties={'magmom': [5, -5]})
        d = s.to_dict
        self.assertEqual(d['sites'][0]['properties']['magmom'], 5)
        self.assertEqual(d['sites'][0]['species'][0]['properties']['spin'], 3)

    def test_from_dict(self):
        d = self.propertied_structure.to_dict
        s = IStructure.from_dict(d)
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
                                                  'oxidation_state': -2,
                                                  'properties': {'spin': 3},
                                                  'element': 'O'}],
                        'label': 'O2-', 'xyz': [0.0, 0.0, 0.0]},
                       {'properties': {'magmom': -5},
                        'abc': [0.75, 0.5, 0.75],
                        'occu': 0.8, 'species': [{'occu': 0.8,
                                                  'oxidation_state': 2,
                                                  'properties': {'spin': 2},
                                                  'element': 'Mg'}],
                        'label': 'Mg2+:0.800',
                        'xyz': [3.8401979336749994, 1.2247250003039056e-06,
                                2.351631795225]}]}
        s = IStructure.from_dict(d)
        self.assertEqual(s[0].magmom, 5)
        self.assertEqual(s[0].specie.spin, 3)
        self.assertEqual(type(s), IStructure)

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

        structure = IStructure(self.lattice, ["O", "Si"], coords,
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
        struct = IStructure(self.lattice, ["Si"] * 2, coords)
        coords2 = list()
        coords2.append([0, 0, 0])
        coords2.append([0.5, 0.5, 0.5])
        struct2 = IStructure(self.struct.lattice, ["Si"] * 2, coords2)
        int_s = struct.interpolate(struct2, 10)
        for s in int_s:
            self.assertIsNotNone(s, "Interpolation Failed!")
        self.assertTrue(np.array_equal(int_s[1][1].frac_coords,
                                       [0.725, 0.5, 0.725]))

        badlattice = [[1, 0.00, 0.00], [0, 1, 0.00], [0.00, 0, 1]]
        struct2 = IStructure(badlattice, ["Si"] * 2, coords2)
        self.assertRaises(ValueError, struct.interpolate, struct2)

        coords2 = list()
        coords2.append([0, 0, 0])
        coords2.append([0.5, 0.5, 0.5])
        struct2 = IStructure(self.struct.lattice, ["Si", "Fe"], coords2)
        self.assertRaises(ValueError, struct.interpolate, struct2)

    def test_get_primitive_structure(self):
        coords = [[0, 0, 0], [0.5, 0.5, 0], [0, 0.5, 0.5], [0.5, 0, 0.5]]
        fcc_ag = IStructure(Lattice.cubic(4.09), ["Ag"] * 4, coords)
        self.assertEqual(len(fcc_ag.get_primitive_structure()), 1)
        coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
        bcc_li = IStructure(Lattice.cubic(4.09), ["Li"] * 2, coords)
        self.assertEqual(len(bcc_li.get_primitive_structure()), 1)

    def test_primitive_structure_volume_check(self):
        l = Lattice.tetragonal(10, 30)
        coords = [[0.5, 0.8, 0], [0.5, 0.2, 0],
                  [0.5, 0.8, 0.333], [0.5, 0.5, 0.333],
                  [0.5, 0.5, 0.666], [0.5, 0.2, 0.666]]
        s = IStructure(l, ["Ag"] * 6, coords)
        sprim = s.get_primitive_structure(tolerance=0.1)
        self.assertEqual(len(sprim), 6)

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


class StructureTest(unittest.TestCase):

    def setUp(self):
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                           [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])
        self.structure = Structure(lattice, ["Si", "Si"], coords)

    def test_non_hash(self):
        self.assertRaises(TypeError, dict, [(self.structure, 1)])

    def test_append_insert_remove_replace(self):
        s = self.structure
        s.insert(1, "O", [0.5, 0.5, 0.5])
        self.assertEqual(s.formula, "Si2 O1")
        self.assertTrue(s.ntypat==2)
        s.remove(2)
        self.assertEqual(s.formula, "Si1 O1")
        s.append("N", [0.25, 0.25, 0.25])
        self.assertEqual(s.formula, "Si1 N1 O1")
        self.assertTrue(s.ntypat==3)
        s.replace(0, "Ge")
        self.assertEqual(s.formula, "Ge1 N1 O1")
        s.replace_species({"Ge": "Si"})
        self.assertEqual(s.formula, "Si1 N1 O1")
        self.assertTrue(s.ntypat==3)

        s.replace_species({"Si": {"Ge": 0.5, "Si": 0.5}})
        self.assertEqual(s.formula, "Si0.5 Ge0.5 N1 O1")
        #this should change the .5Si .5Ge sites to .75Si .25Ge
        s.replace_species({"Ge": {"Ge": 0.5, "Si": 0.5}})
        self.assertEqual(s.formula, "Si0.75 Ge0.25 N1 O1")

        # In this case, s.ntypat is ambiguous. 
        # for the time being, we raise AttributeError.
        with self.assertRaises(AttributeError):
            s.ntypat

        s.remove_species(["Si"])
        self.assertEqual(s.formula, "Ge0.25 N1 O1")

        s.remove_sites([1, 2])
        self.assertEqual(s.formula, "Ge0.25")

    def test_add_site_property(self):
        s = self.structure
        s.add_site_property("charge", [4.1, -5])
        self.assertEqual(s[0].charge, 4.1)
        self.assertEqual(s[1].charge, -5)
        s.add_site_property("magmom", [3, 2])
        self.assertEqual(s[0].charge, 4.1)
        self.assertEqual(s[0].magmom, 3)

    def test_perturb(self):
        d = 0.1
        pre_perturbation_sites = self.structure.sites[:]
        self.structure.perturb(distance=d)
        post_perturbation_sites = self.structure.sites

        for i, x in enumerate(pre_perturbation_sites):
            self.assertAlmostEqual(x.distance(post_perturbation_sites[i]), d,
                                   3, "Bad perturbation distance")

    def test_add_oxidation_states(self):
        oxidation_states = {"Si": -4}
        self.structure.add_oxidation_state_by_element(oxidation_states)
        for site in self.structure:
            for k in site.species_and_occu.keys():
                self.assertEqual(k.oxi_state, oxidation_states[k.symbol],
                                 "Wrong oxidation state assigned!")
        oxidation_states = {"Fe": 2}
        self.assertRaises(ValueError,
                          self.structure.add_oxidation_state_by_element,
                          oxidation_states)
        self.structure.add_oxidation_state_by_site([2, -4])
        self.assertEqual(self.structure[0].specie.oxi_state, 2)
        self.assertRaises(ValueError,
                          self.structure.add_oxidation_state_by_site,
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
        s_specie.remove_oxidation_states()
        self.assertEqual(s_elem, s_specie, "Oxidation state remover "
                                           "failed")

    def test_apply_operation(self):
        op = SymmOp.from_axis_angle_and_translation([0, 0, 1], 90)
        self.structure.apply_operation(op)
        self.assertTrue(np.allclose(
            [[0.000000, 3.840198, 0.000000],
             [-3.325710, 1.920099, 0.000000],
             [2.217138, -0.000000, 3.135509]], self.structure.lattice.matrix))

    def test_apply_strain(self):
        self.structure.apply_strain(0.01)
        self.assertAlmostEqual(
            self.structure.lattice.abc,
            (3.8785999130369997, 3.878600984287687, 3.8785999130549516))

    def test_translate_sites(self):
        self.structure.translate_sites([0, 1], [0.5, 0.5, 0.5],
                                       frac_coords=True)
        self.assertTrue(np.array_equal(self.structure.frac_coords[0],
                                       np.array([0.5, 0.5, 0.5])))

        self.structure.translate_sites([0], [0.5, 0.5, 0.5],
                                       frac_coords=False)
        self.assertTrue(np.allclose(self.structure.cart_coords[0],
                                    [3.38014845, 1.05428585, 2.06775453]))

    def test_make_supercell(self):
        self.structure.make_supercell([2, 1, 1])
        self.assertEqual(self.structure.formula, "Si4")
        self.structure.make_supercell([[1, 0, 0], [2, 1, 0], [0, 0, 1]])
        self.assertEqual(self.structure.formula, "Si4")
        self.structure.make_supercell(2)
        self.assertEqual(self.structure.formula, "Si32")
        self.assertTrue(np.allclose(self.structure.lattice.abc,
                                    [15.360792, 35.195996, 7.680396]))

    def test_to_from_dict(self):
        d = self.structure.to_dict
        s2 = Structure.from_dict(d)
        self.assertEqual(type(s2), Structure)


class IMoleculeTest(unittest.TestCase):

    def setUp(self):
        coords = [[0.000000, 0.000000, 0.000000],
                  [0.000000, 0.000000, 1.089000],
                  [1.026719, 0.000000, -0.363000],
                  [-0.513360, -0.889165, -0.363000],
                  [-0.513360, 0.889165, -0.363000]]
        self.coords = coords
        self.mol = Molecule(["C", "H", "H", "H", "H"], coords)

    def test_bad_molecule(self):
        coords = [[0.000000, 0.000000, 0.000000],
                  [0.000000, 0.000000, 1.089000],
                  [1.026719, 0.000000, -0.363000],
                  [-0.513360, -0.889165, -0.363000],
                  [-0.513360, 0.889165, -0.363000],
                  [-0.513360, 0.889165, -0.36301]]
        self.assertRaises(StructureError, Molecule,
                          ["C", "H", "H", "H", "H", "H"], coords,
                          validate_proximity=True)

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
Charge = 0, Spin Mult = 1
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
        self.assertEqual(propertied_mol, mol)
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

    def test_prop(self):
        self.assertEqual(self.mol.charge, 0)
        self.assertEqual(self.mol.spin_multiplicity, 1)
        self.assertEqual(self.mol.nelectrons, 10)
        self.assertTrue(np.allclose(self.mol.center_of_mass, np.zeros(3),
                                    atol=1e-7))
        self.assertRaises(ValueError, Molecule, ["C", "H", "H", "H", "H"],
                          self.coords, charge=1, spin_multiplicity=1)
        mol = Molecule(["C", "H", "H", "H", "H"], self.coords, charge=1)
        self.assertEqual(mol.spin_multiplicity, 2)
        self.assertEqual(mol.nelectrons, 9)

        #Triplet O2
        mol = IMolecule(["O"] * 2, [[0, 0, 0], [0, 0, 1.2]],
                        spin_multiplicity=3)
        self.assertEqual(mol.spin_multiplicity, 3)

    def test_equal(self):
        mol = IMolecule(["C", "H", "H", "H", "H"], self.coords, charge=1)
        self.assertNotEqual(mol, self.mol)

    def test_get_centered_molecule(self):
        mol = IMolecule(["O"] * 2, [[0, 0, 0], [0, 0, 1.2]],
                        spin_multiplicity=3)
        centered = mol.get_centered_molecule()
        self.assertFalse(np.allclose(mol.center_of_mass, np.zeros(3),
                                     atol=1e-7))
        self.assertTrue(np.allclose(centered.center_of_mass, np.zeros(3),
                                    atol=1e-7))

    def test_to_from_dict(self):
        d = self.mol.to_dict
        mol2 = IMolecule.from_dict(d)
        self.assertEqual(type(mol2), IMolecule)


class MoleculeTest(unittest.TestCase):

    def setUp(self):
        coords = [[0.000000, 0.000000, 0.000000],
                  [0.000000, 0.000000, 1.089000],
                  [1.026719, 0.000000, -0.363000],
                  [-0.513360, -0.889165, -0.363000],
                  [-0.513360, 0.889165, -0.363000]]
        self.mol = Molecule(["C", "H", "H", "H", "H"], coords)

    def test_insert_remove_append(self):
        mol = self.mol
        mol.insert(1, "O", [0.5, 0.5, 0.5])
        self.assertEqual(mol.formula, "H4 C1 O1")
        mol.remove(2)
        self.assertEqual(mol.formula, "H3 C1 O1")
        mol.set_charge_and_spin(0)
        self.assertEqual(mol.spin_multiplicity, 2)
        mol.append("N", [0.25, 0.25, 0.25])
        self.assertEqual(mol.formula, "H3 C1 N1 O1")
        self.assertRaises(TypeError, dict, [(mol, 1)])
        mol.remove_sites([0, 1])
        self.assertEqual(mol.formula, "H3 N1")

    def test_translate_sites(self):
        self.mol.translate_sites([0, 1], [0.5, 0.5, 0.5])
        self.assertTrue(np.array_equal(self.mol.cart_coords[0],
                                       np.array([0.5, 0.5, 0.5])))

    def test_replace(self):
        self.mol.replace(0, "Ge")
        self.assertEqual(self.mol.formula, "Ge1 H4")

        self.mol.replace_species({Element("Ge"): {Element("Ge"): 0.5,
                                                  Element("Si"): 0.5}})
        self.assertEqual(self.mol.formula, "Si0.5 Ge0.5 H4")

        #this should change the .5Si .5Ge sites to .75Si .25Ge
        self.mol.replace_species({Element("Ge"): {Element("Ge"): 0.5,
                                                  Element("Si"): 0.5}})
        self.assertEqual(self.mol.formula, "Si0.75 Ge0.25 H4")

        d = 0.1
        pre_perturbation_sites = self.mol.sites[:]
        self.mol.perturb(distance=d)
        post_perturbation_sites = self.mol.sites

        for i, x in enumerate(pre_perturbation_sites):
            self.assertAlmostEqual(x.distance(post_perturbation_sites[i]), d,
                                   3, "Bad perturbation distance")

    def test_add_site_property(self):
        self.mol.add_site_property("charge", [4.1, -2, -2, -2, -2])
        self.assertEqual(self.mol[0].charge, 4.1)
        self.assertEqual(self.mol[1].charge, -2)

        self.mol.add_site_property("magmom", [3, 2, 2, 2, 2])
        self.assertEqual(self.mol[0].charge, 4.1)
        self.assertEqual(self.mol[0].magmom, 3)

    def test_to_from_dict(self):
        d = self.mol.to_dict
        mol2 = Molecule.from_dict(d)
        self.assertEqual(type(mol2), Molecule)

    def test_apply_operation(self):
        op = SymmOp.from_axis_angle_and_translation([0, 0, 1], 90)
        self.mol.apply_operation(op)
        self.assertTrue(np.allclose(self.mol[2].coords,
                                    [0.000000, 1.026719, -0.363000]))


    def test_substitute(self):
        coords = [[0.000000, 0.000000, 1.08],
                  [0.000000, 0.000000, 0.000000],
                  [1.026719, 0.000000, -0.363000],
                  [-0.513360, -0.889165, -0.363000],
                  [-0.513360, 0.889165, -0.363000]]
        sub = Molecule(["X", "C", "H", "H", "H"], coords)
        self.mol.substitute(1, sub)
        self.assertAlmostEqual(self.mol.get_distance(0, 4), 1.54)
        sub = Molecule(["X", "F"], [[0, 0, 0], [0, 0, 1.11]])
        self.mol.substitute(2, sub)
        self.assertAlmostEqual(self.mol.get_distance(0, 7), 1.35)
        oh = Molecule(["X", "O", "H"],
                      [[0, 0.780362, -.456316], [0, 0, .114079],
                       [0, -.780362, -.456316]])
        self.mol.substitute(1, oh)
        self.assertAlmostEqual(self.mol.get_distance(0, 7), 1.43)

if __name__ == '__main__':
    unittest.main()
