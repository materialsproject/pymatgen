# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

from pymatgen.util.testing import PymatgenTest
from pymatgen.core.periodic_table import Element, Specie
from pymatgen.core.composition import Composition
from pymatgen.core.operations import SymmOp
from pymatgen.core.structure import IStructure, Structure, IMolecule, \
    StructureError, Molecule
from pymatgen.core.lattice import Lattice
import random
import warnings
import os


class IStructureTest(PymatgenTest):

    def setUp(self):
        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        self.lattice = Lattice([[3.8401979337, 0.00, 0.00],
                                [1.9200989668, 3.3257101909, 0.00],
                                [0.00, -2.2171384943, 3.1355090603]])
        self.struct = IStructure(self.lattice, ["Si"] * 2, coords)
        self.assertEqual(len(self.struct), 2,
                         "Wrong number of sites in structure!")
        self.assertTrue(self.struct.is_ordered)
        self.assertTrue(self.struct.ntypesp == 1)
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
        #these shouldn't raise an error
        IStructure(self.lattice, ["Si"] * 2, coords[:2], True)
        IStructure(self.lattice, ["Si"], coords[:1], True)


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
        self.assertEqual(s.composition.formula, 'Mg0.8 O1')

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
        s = IStructure(self.lattice, ["Se", "C", "Se", "C"],
                       [[0] * 3, [0.5] * 3, [0.25] * 3, [0.75] * 3])
        self.assertEqual([site.specie.symbol
                          for site in s.get_sorted_structure()],
                         ["C", "C", "Se", "Se"])


    def test_fractional_occupations(self):
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        s = IStructure(self.lattice, [{'O': 1.0}, {'Mg': 0.8}],
                       coords)
        self.assertEqual(s.composition.formula, 'Mg0.8 O1')
        self.assertFalse(s.is_ordered)

    def test_get_distance(self):
        self.assertAlmostEqual(self.struct.get_distance(0, 1), 2.35, 2,
                               "Distance calculated wrongly!")
        pt = [0.9, 0.9, 0.8]
        self.assertAlmostEqual(self.struct[0].distance_from_point(pt),
                               1.50332963784, 2,
                               "Distance calculated wrongly!")

    def test_as_dict(self):
        si = Specie("Si", 4)
        mn = Element("Mn")
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        struct = IStructure(self.lattice, [{si: 0.5, mn: 0.5}, {si: 0.5}],
                            coords)
        self.assertIn("lattice", struct.as_dict())
        self.assertIn("sites", struct.as_dict())
        d = self.propertied_structure.as_dict()
        self.assertEqual(d['sites'][0]['properties']['magmom'], 5)
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        s = IStructure(self.lattice, [{Specie('O', -2,
                                              properties={"spin": 3}): 1.0},
                                      {Specie('Mg', 2,
                                              properties={"spin": 2}): 0.8}],
                       coords, site_properties={'magmom': [5, -5]})
        d = s.as_dict()
        self.assertEqual(d['sites'][0]['properties']['magmom'], 5)
        self.assertEqual(d['sites'][0]['species'][0]['properties']['spin'], 3)

    def test_from_dict(self):
        d = self.propertied_structure.as_dict()
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
            self.assertEqual(int_s[0].lattice, s.lattice)
        self.assertArrayEqual(int_s[1][1].frac_coords, [0.725, 0.5, 0.725])

        badlattice = [[1, 0.00, 0.00], [0, 1, 0.00], [0.00, 0, 1]]
        struct2 = IStructure(badlattice, ["Si"] * 2, coords2)
        self.assertRaises(ValueError, struct.interpolate, struct2)

        coords2 = list()
        coords2.append([0, 0, 0])
        coords2.append([0.5, 0.5, 0.5])
        struct2 = IStructure(self.struct.lattice, ["Si", "Fe"], coords2)
        self.assertRaises(ValueError, struct.interpolate, struct2)

        # Test autosort feature.
        s1 = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3),
                                       ["Fe"], [[0, 0, 0]])
        s1.pop(0)
        s2 = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3),
                                       ["Fe"], [[0, 0, 0]])
        s2.pop(2)
        random.shuffle(s2)

        for s in s1.interpolate(s2, autosort_tol=0.5):
            self.assertArrayAlmostEqual(s1[0].frac_coords, s[0].frac_coords)
            self.assertArrayAlmostEqual(s1[2].frac_coords, s[2].frac_coords)

        # Make sure autosort has no effect on simpler interpolations,
        # and with shuffled sites.
        s1 = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3),
                                       ["Fe"], [[0, 0, 0]])
        s2 = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3),
                                       ["Fe"], [[0, 0, 0]])
        s2[0] = "Fe", [0.01, 0.01, 0.01]
        random.shuffle(s2)

        for s in s1.interpolate(s2, autosort_tol=0.5):
            self.assertArrayAlmostEqual(s1[1].frac_coords, s[1].frac_coords)
            self.assertArrayAlmostEqual(s1[2].frac_coords, s[2].frac_coords)
            self.assertArrayAlmostEqual(s1[3].frac_coords, s[3].frac_coords)

    def test_interpolate_lattice(self):
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        struct = IStructure(self.lattice, ["Si"] * 2, coords)
        coords2 = list()
        coords2.append([0, 0, 0])
        coords2.append([0.5, 0.5, 0.5])
        l2 = Lattice.from_lengths_and_angles([3,4,4], [100,100,70])
        struct2 = IStructure(l2, ["Si"] * 2, coords2)
        int_s = struct.interpolate(struct2, 2, interpolate_lattices=True)
        self.assertArrayAlmostEqual(struct.lattice.abc,
                                    int_s[0].lattice.abc)
        self.assertArrayAlmostEqual(struct.lattice.angles,
                                    int_s[0].lattice.angles)
        self.assertArrayAlmostEqual(struct2.lattice.abc,
                                    int_s[2].lattice.abc)
        self.assertArrayAlmostEqual(struct2.lattice.angles,
                                    int_s[2].lattice.angles)
        int_angles = [(a + struct2.lattice.angles[i]) / 2
                      for i, a in enumerate(struct.lattice.angles)]
        self.assertArrayAlmostEqual(int_angles,
                                    int_s[1].lattice.angles)

    def test_get_primitive_structure(self):
        coords = [[0, 0, 0], [0.5, 0.5, 0], [0, 0.5, 0.5], [0.5, 0, 0.5]]
        fcc_ag = IStructure(Lattice.cubic(4.09), ["Ag"] * 4, coords)
        self.assertEqual(len(fcc_ag.get_primitive_structure()), 1)
        coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
        bcc_li = IStructure(Lattice.cubic(4.09), ["Li"] * 2, coords)
        bcc_prim = bcc_li.get_primitive_structure()
        self.assertEqual(len(bcc_prim), 1)
        self.assertAlmostEqual(bcc_prim.lattice.alpha, 109.47122, 3)

        coords = [[0] * 3, [0.5] * 3, [0.25] * 3, [0.26] * 3]
        s = IStructure(Lattice.cubic(4.09), ["Ag"] * 4, coords)
        self.assertEqual(len(s.get_primitive_structure()), 4)

    def test_primitive_cell_site_merging(self):
        l = Lattice.cubic(10)
        coords = [[0, 0, 0], [0, 0, 0.5],
                  [0, 0, 0.26], [0, 0, 0.74]]
        sp = ['Ag', 'Ag', 'Be', 'Be']
        s = Structure(l, sp, coords)
        dm = s.get_primitive_structure().distance_matrix
        self.assertArrayAlmostEqual(dm, [[0, 2.5], [2.5, 0]])

    def test_primitive_on_large_supercell(self):
        coords = [[0, 0, 0], [0.5, 0.5, 0], [0, 0.5, 0.5], [0.5, 0, 0.5]]
        fcc_ag = Structure(Lattice.cubic(4.09), ["Ag"] * 4, coords)
        fcc_ag.make_supercell([2, 2, 2])
        fcc_ag_prim = fcc_ag.get_primitive_structure()
        self.assertEqual(len(fcc_ag_prim), 1)
        self.assertAlmostEqual(fcc_ag_prim.volume, 17.10448225)

    def test_primitive_positions(self):
        coords = [[0, 0, 0], [0.3, 0.35, 0.45]]
        s = Structure(Lattice.from_parameters(1,2,3,50,66,88), ["Ag"] * 2, coords)

        a = [[-1,2,-3], [3,2,-4], [1,0,-1]]
        b = [[4, 0, 0], [1, 1, 0], [3, 0, 1]]
        c = [[2, 0, 0], [1, 3, 0], [1, 1, 1]]

        for sc_matrix in [c]:
            sc = s.copy()
            sc.make_supercell(sc_matrix)
            prim = sc.get_primitive_structure(0.01)

            self.assertEqual(len(prim), 2)
            self.assertAlmostEqual(prim.distance_matrix[0,1], 1.0203432356739286)

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
        nn = s.get_neighbors_in_shell(s[0].frac_coords, 2, 4,
                                       include_index=True)
        self.assertEqual(len(nn), 47)
        self.assertEqual(nn[0][-1], 0)

        r = random.uniform(3, 6)
        all_nn = s.get_all_neighbors(r, True)
        for i in range(len(s)):
            self.assertEqual(len(all_nn[i]), len(s.get_neighbors(s[i], r)))

        for site, nns in zip(s, all_nn):
            for nn in nns:
                self.assertTrue(nn[0].is_periodic_image(s[nn[2]]))
                d = sum((site.coords - nn[0].coords) ** 2) ** 0.5
                self.assertAlmostEqual(d, nn[1])

        s = Structure(Lattice.cubic(1), ['Li'], [[0,0,0]])
        s.make_supercell([2,2,2])
        self.assertEqual(sum(map(len, s.get_all_neighbors(3))), 976)


    def test_get_all_neighbors_outside_cell(self):
        s = Structure(Lattice.cubic(2), ['Li', 'Li', 'Li', 'Si'],
                      [[3.1] * 3, [0.11] * 3, [-1.91] * 3, [0.5] * 3])
        all_nn = s.get_all_neighbors(0.2, True)
        for site, nns in zip(s, all_nn):
            for nn in nns:
                self.assertTrue(nn[0].is_periodic_image(s[nn[2]]))
                d = sum((site.coords - nn[0].coords) ** 2) ** 0.5
                self.assertAlmostEqual(d, nn[1])
        self.assertEqual(list(map(len, all_nn)), [2, 2, 2, 0])

    def test_get_dist_matrix(self):
        ans = [[0., 2.3516318],
               [2.3516318, 0.]]
        self.assertArrayAlmostEqual(self.struct.distance_matrix, ans)

    def test_to_from_file_string(self):
        for fmt in ["cif", "json", "poscar", "cssr"]:
            s = self.struct.to(fmt=fmt)
            self.assertIsNotNone(s)
            ss = IStructure.from_str(s, fmt=fmt)
            self.assertArrayAlmostEqual(
                ss.lattice.lengths_and_angles,
                self.struct.lattice.lengths_and_angles, decimal=5)
            self.assertArrayAlmostEqual(ss.frac_coords, self.struct.frac_coords)
            self.assertIsInstance(ss, IStructure)

        self.struct.to(filename="POSCAR.testing")
        self.assertTrue(os.path.exists("POSCAR.testing"))
        os.remove("POSCAR.testing")

        self.struct.to(filename="Si_testing.yaml")
        self.assertTrue(os.path.exists("Si_testing.yaml"))
        s = Structure.from_file("Si_testing.yaml")
        self.assertEqual(s, self.struct)
        os.remove("Si_testing.yaml")


class StructureTest(PymatgenTest):

    def setUp(self):
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                           [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])
        self.structure = Structure(lattice, ["Si", "Si"], coords)

    def test_mutable_sequence_methods(self):
        s = self.structure
        s[0] = "Fe"
        self.assertEqual(s.formula, "Fe1 Si1")
        s[0] = "Fe", [0.5, 0.5, 0.5]
        self.assertEqual(s.formula, "Fe1 Si1")
        self.assertArrayAlmostEqual(s[0].frac_coords, [0.5, 0.5, 0.5])
        s.reverse()
        self.assertEqual(s[0].specie, Element("Si"))
        self.assertArrayAlmostEqual(s[0].frac_coords, [0.75, 0.5, 0.75])
        s[0] = {"Mn": 0.5}
        self.assertEqual(s.formula, "Mn0.5 Fe1")
        del s[1]
        self.assertEqual(s.formula, "Mn0.5")
        s[0] = "Fe", [0.9, 0.9, 0.9], {"magmom": 5}
        self.assertEqual(s.formula, "Fe1")
        self.assertEqual(s[0].magmom, 5)

    def test_non_hash(self):
        self.assertRaises(TypeError, dict, [(self.structure, 1)])

    def test_sort(self):
        s = self.structure
        s[0] = "F"
        s.sort()
        self.assertEqual(s[0].species_string, "Si")
        self.assertEqual(s[1].species_string, "F")
        s.sort(key=lambda site: site.species_string)
        self.assertEqual(s[0].species_string, "F")
        self.assertEqual(s[1].species_string, "Si")
        s.sort(key=lambda site: site.species_string, reverse=True)
        self.assertEqual(s[0].species_string, "Si")
        self.assertEqual(s[1].species_string, "F")

    def test_append_insert_remove_replace(self):
        s = self.structure
        s.insert(1, "O", [0.5, 0.5, 0.5])
        self.assertEqual(s.formula, "Si2 O1")
        self.assertTrue(s.ntypesp == 2)
        self.assertTrue(s.symbol_set == ("Si", "O"))
        self.assertTrue(s.indices_from_symbol("Si") == (0,2))
        self.assertTrue(s.indices_from_symbol("O") == (1,))
        del s[2]
        self.assertEqual(s.formula, "Si1 O1")
        self.assertTrue(s.indices_from_symbol("Si") == (0,))
        self.assertTrue(s.indices_from_symbol("O") == (1,))
        s.append("N", [0.25, 0.25, 0.25])
        self.assertEqual(s.formula, "Si1 N1 O1")
        self.assertTrue(s.ntypesp == 3)
        self.assertTrue(s.symbol_set == ("Si", "O", "N"))
        self.assertTrue(s.indices_from_symbol("Si") == (0,))
        self.assertTrue(s.indices_from_symbol("O") == (1,))
        self.assertTrue(s.indices_from_symbol("N") == (2,))
        s[0] = "Ge"
        self.assertEqual(s.formula, "Ge1 N1 O1")
        self.assertTrue(s.symbol_set == ("Ge", "O", "N"))
        s.replace_species({"Ge": "Si"})
        self.assertEqual(s.formula, "Si1 N1 O1")
        self.assertTrue(s.ntypesp == 3)

        s.replace_species({"Si": {"Ge": 0.5, "Si": 0.5}})
        self.assertEqual(s.formula, "Si0.5 Ge0.5 N1 O1")
        #this should change the .5Si .5Ge sites to .75Si .25Ge
        s.replace_species({"Ge": {"Ge": 0.5, "Si": 0.5}})
        self.assertEqual(s.formula, "Si0.75 Ge0.25 N1 O1")

        # In this case, s.ntypesp is ambiguous.
        # for the time being, we raise AttributeError.
        with self.assertRaises(AttributeError):
            s.ntypesp

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

    def test_propertied_structure(self):
        #Make sure that site properties are set to None for missing values.
        s = self.structure
        s.add_site_property("charge", [4.1, -5])
        s.append("Li", [0.3, 0.3 ,0.3])
        self.assertEqual(len(s.site_properties["charge"]), 3)

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
        s = self.structure.copy()
        s.apply_operation(op)
        self.assertArrayAlmostEqual(
            s.lattice.matrix,
            [[0.000000, 3.840198, 0.000000],
             [-3.325710, 1.920099, 0.000000],
             [2.217138, -0.000000, 3.135509]], 5)

        op = SymmOp([[1, 1, 0, 0.5], [1, 0, 0, 0.5], [0, 0, 1, 0.5],
                     [0, 0, 0, 1]])
        s = self.structure.copy()
        s.apply_operation(op, fractional=True)
        self.assertArrayAlmostEqual(
            s.lattice.matrix,
            [[5.760297, 3.325710, 0.000000],
             [3.840198, 0.000000, 0.000000],
             [0.000000, -2.217138, 3.135509]], 5)

    def test_apply_strain(self):
        s = self.structure
        initial_coord = s[1].coords
        s.apply_strain(0.01)
        self.assertAlmostEqual(
            s.lattice.abc,
            (3.8785999130369997, 3.878600984287687, 3.8785999130549516))
        self.assertArrayAlmostEqual(s[1].coords, initial_coord * 1.01)
        a1, b1, c1 = s.lattice.abc
        s.apply_strain([0.1, 0.2, 0.3])
        a2, b2, c2 = s.lattice.abc
        self.assertAlmostEqual(a2 / a1, 1.1)
        self.assertAlmostEqual(b2 / b1, 1.2)
        self.assertAlmostEqual(c2 / c1, 1.3)

    def test_scale_lattice(self):
        initial_coord = self.structure[1].coords
        self.structure.scale_lattice(self.structure.volume * 1.01 ** 3)
        self.assertArrayAlmostEqual(
            self.structure.lattice.abc,
            (3.8785999130369997, 3.878600984287687, 3.8785999130549516))
        self.assertArrayAlmostEqual(self.structure[1].coords,
            initial_coord * 1.01)

    def test_translate_sites(self):
        self.structure.translate_sites([0, 1], [0.5, 0.5, 0.5],
                                       frac_coords=True)
        self.assertArrayEqual(self.structure.frac_coords[0],
                              [0.5, 0.5, 0.5])

        self.structure.translate_sites([0], [0.5, 0.5, 0.5],
                                       frac_coords=False)
        self.assertArrayAlmostEqual(self.structure.cart_coords[0],
                                    [3.38014845, 1.05428585, 2.06775453])

        self.structure.translate_sites([0], [0.5, 0.5, 0.5],
                                       frac_coords=True, to_unit_cell=False)
        self.assertArrayAlmostEqual(self.structure.frac_coords[0],
                                    [1.00187517, 1.25665291, 1.15946374])

    def test_make_supercell(self):
        self.structure.make_supercell([2, 1, 1])
        self.assertEqual(self.structure.formula, "Si4")
        self.structure.make_supercell([[1, 0, 0], [2, 1, 0], [0, 0, 1]])
        self.assertEqual(self.structure.formula, "Si4")
        self.structure.make_supercell(2)
        self.assertEqual(self.structure.formula, "Si32")
        self.assertArrayAlmostEqual(self.structure.lattice.abc,
                                    [15.360792, 35.195996, 7.680396], 5)

    def test_disordered_supercell_primitive_cell(self):
        l = Lattice.cubic(2)
        f = [[0.5, 0.5, 0.5]]
        sp = [{'Si': 0.54738}]
        s = Structure(l, sp, f)
        #this supercell often breaks things
        s.make_supercell([[0,-1,1],[-1,1,0],[1,1,1]])
        self.assertEqual(len(s.get_primitive_structure()), 1)

    def test_another_supercell(self):
        #this is included b/c for some reason the old algo was failing on it
        s = self.structure.copy()
        s.make_supercell([[0, 2, 2], [2, 0, 2], [2, 2, 0]])
        self.assertEqual(s.formula, "Si32")
        s = self.structure.copy()
        s.make_supercell([[0, 2, 0], [1, 0, 0], [0, 0, 1]])
        self.assertEqual(s.formula, "Si4")

    def test_to_from_dict(self):
        d = self.structure.as_dict()
        s2 = Structure.from_dict(d)
        self.assertEqual(type(s2), Structure)

    def test_propertied_structure_mod(self):
        prop_structure = Structure(
            self.structure.lattice, ["Si"] * 2, self.structure.frac_coords,
            site_properties={'magmom': [5, -5]})
        prop_structure.append("C", [0.25, 0.25, 0.25])
        d = prop_structure.as_dict()
        with warnings.catch_warnings(record=True) as w:
            # Cause all warnings to always be triggered.
            warnings.simplefilter("always")
            s2 = Structure.from_dict(d)
            self.assertEqual(len(w), 1)
            self.assertEqual(
                str(w[0].message),
                'Not all sites have property magmom. Missing values are set '
                'to None.')

    def test_to_from_file_string(self):
        for fmt in ["cif", "json", "poscar", "cssr", "yaml", "xsf"]:
            s = self.structure.to(fmt=fmt)
            self.assertIsNotNone(s)
            ss = Structure.from_str(s, fmt=fmt)
            self.assertArrayAlmostEqual(
                ss.lattice.lengths_and_angles,
                self.structure.lattice.lengths_and_angles, decimal=5)
            self.assertArrayAlmostEqual(ss.frac_coords,
                                        self.structure.frac_coords)
            self.assertIsInstance(ss, Structure)

        self.structure.to(filename="POSCAR.testing")
        self.assertTrue(os.path.exists("POSCAR.testing"))
        os.remove("POSCAR.testing")

        self.structure.to(filename="structure_testing.json")
        self.assertTrue(os.path.exists("structure_testing.json"))
        s = Structure.from_file("structure_testing.json")
        self.assertEqual(s, self.structure)
        os.remove("structure_testing.json")

    def test_from_spacegroup(self):
        s1 = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3), ["Li", "O"],
                                       [[0.25, 0.25, 0.25], [0, 0, 0]])
        self.assertEqual(s1.formula, "Li8 O4")
        s2 = Structure.from_spacegroup(225, Lattice.cubic(3), ["Li", "O"],
                                       [[0.25, 0.25, 0.25], [0, 0, 0]])
        self.assertEqual(s1, s2)

        s2 = Structure.from_spacegroup(225, Lattice.cubic(3), ["Li", "O"],
                                       [[0.25, 0.25, 0.25], [0, 0, 0]],
                                       site_properties={"charge": [1, -2]})
        self.assertEqual(sum(s2.site_properties["charge"]), 0)

        s = Structure.from_spacegroup("Pm-3m", Lattice.cubic(3), ["Cs", "Cl"],
                                      [[0, 0, 0], [0.5, 0.5, 0.5]])
        self.assertEqual(s.formula, "Cs1 Cl1")

        self.assertRaises(ValueError, Structure.from_spacegroup,
                          "Pm-3m", Lattice.tetragonal(1, 3), ["Cs", "Cl"],
                          [[0, 0, 0], [0.5, 0.5, 0.5]])

    def test_merge_sites(self):
        species = [{'Ag': 0.5}, {'Cl': 0.25}, {'Cl': 0.1},
                   {'Ag': 0.5}, {'F': 0.15}, {'F': 0.1}]
        coords = [[0, 0, 0], [0.5, 0.5, 0.5], [0.5, 0.5, 0.5],
                  [0, 0, 0], [0.5, 0.5, 1.501], [0.5, 0.5, 1.501]]
        s = Structure(Lattice.cubic(1), species, coords)
        s.merge_sites()
        self.assertEqual(s[0].specie.symbol, 'Ag')
        self.assertEqual(s[1].species_and_occu,
                         Composition({'Cl': 0.35, 'F': 0.25}))
        self.assertArrayAlmostEqual(s[1].frac_coords, [.5, .5, .5005])


class IMoleculeTest(PymatgenTest):

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
        ans = """Full Formula (H4 C1)
Reduced Formula: H4C
Charge = 0, Spin Mult = 1
Sites (5)
0 C     0.000000     0.000000     0.000000
1 H     0.000000     0.000000     1.089000
2 H     1.026719     0.000000    -0.363000
3 H    -0.513360    -0.889165    -0.363000
4 H    -0.513360     0.889165    -0.363000"""
        self.assertEqual(self.mol.__str__(), ans)
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

    def test_get_boxed_structure(self):
        s = self.mol.get_boxed_structure(9, 9, 9)
        # C atom should be in center of box.
        self.assertArrayAlmostEqual(s[4].frac_coords,
                                    [0.50000001,  0.5,  0.5])
        self.assertArrayAlmostEqual(s[1].frac_coords,
                                    [0.6140799, 0.5,  0.45966667])
        self.assertRaises(ValueError, self.mol.get_boxed_structure, 1, 1, 1)
        s2 = self.mol.get_boxed_structure(5, 5, 5, (2, 3, 4))
        self.assertEqual(len(s2), 24 * 5)
        self.assertEqual(s2.lattice.abc, (10, 15, 20))

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
        self.assertArrayAlmostEqual(self.mol.distance_matrix, ans)

    def test_break_bond(self):
        (mol1, mol2) = self.mol.break_bond(0, 1)
        self.assertEqual(mol1.formula, "H3 C1")
        self.assertEqual(mol2.formula, "H1")

    def test_prop(self):
        self.assertEqual(self.mol.charge, 0)
        self.assertEqual(self.mol.spin_multiplicity, 1)
        self.assertEqual(self.mol.nelectrons, 10)
        self.assertArrayAlmostEqual(self.mol.center_of_mass, [0, 0, 0])
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
        self.assertArrayAlmostEqual(centered.center_of_mass, [0, 0, 0])

    def test_to_from_dict(self):
        d = self.mol.as_dict()
        mol2 = IMolecule.from_dict(d)
        self.assertEqual(type(mol2), IMolecule)
        propertied_mol = Molecule(["C", "H", "H", "H", "H"], self.coords,
                                  charge=1,
                                  site_properties={'magmom':
                                                   [0.5, -0.5, 1, 2, 3]})
        d = propertied_mol.as_dict()
        self.assertEqual(d['sites'][0]['properties']['magmom'], 0.5)
        mol = Molecule.from_dict(d)
        self.assertEqual(propertied_mol, mol)
        self.assertEqual(mol[0].magmom, 0.5)
        self.assertEqual(mol.formula, "H4 C1")
        self.assertEqual(mol.charge, 1)

    def test_to_from_file_string(self):
        for fmt in ["xyz", "json", "g03", "yaml"]:
            s = self.mol.to(fmt=fmt)
            self.assertIsNotNone(s)
            m = IMolecule.from_str(s, fmt=fmt)
            self.assertEqual(m, self.mol)
            self.assertIsInstance(m, IMolecule)

        self.mol.to(filename="CH4_testing.xyz")
        self.assertTrue(os.path.exists("CH4_testing.xyz"))
        os.remove("CH4_testing.xyz")
        self.mol.to(filename="CH4_testing.yaml")
        self.assertTrue(os.path.exists("CH4_testing.yaml"))
        mol = Molecule.from_file("CH4_testing.yaml")
        self.assertEqual(self.mol, mol)
        os.remove("CH4_testing.yaml")


class MoleculeTest(PymatgenTest):

    def setUp(self):
        coords = [[0.000000, 0.000000, 0.000000],
                  [0.000000, 0.000000, 1.089000],
                  [1.026719, 0.000000, -0.363000],
                  [-0.513360, -0.889165, -0.363000],
                  [-0.513360, 0.889165, -0.363000]]
        self.mol = Molecule(["C", "H", "H", "H", "H"], coords)

    def test_mutable_sequence_methods(self):
        s = self.mol
        s[1] = ("F", [0.5, 0.5, 0.5])
        self.assertEqual(s.formula, "H3 C1 F1")
        self.assertArrayAlmostEqual(s[1].coords, [0.5, 0.5, 0.5])
        s.reverse()
        self.assertEqual(s[0].specie, Element("H"))
        self.assertArrayAlmostEqual(s[0].coords,
                                    [-0.513360, 0.889165, -0.363000])
        del s[1]
        self.assertEqual(s.formula, "H2 C1 F1")
        s[3] = "N", [0,0,0], {"charge": 4}
        self.assertEqual(s.formula, "H2 N1 F1")
        self.assertEqual(s[3].charge, 4)

    def test_insert_remove_append(self):
        mol = self.mol
        mol.insert(1, "O", [0.5, 0.5, 0.5])
        self.assertEqual(mol.formula, "H4 C1 O1")
        del mol[2]
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
        self.assertArrayEqual(self.mol.cart_coords[0],
                              [0.5, 0.5, 0.5])

    def test_replace(self):
        self.mol[0] = "Ge"
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
        d = self.mol.as_dict()
        mol2 = Molecule.from_dict(d)
        self.assertEqual(type(mol2), Molecule)

    def test_apply_operation(self):
        op = SymmOp.from_axis_angle_and_translation([0, 0, 1], 90)
        self.mol.apply_operation(op)
        self.assertArrayAlmostEqual(self.mol[2].coords,
                                    [0.000000, 1.026719, -0.363000])

    def test_substitute(self):
        coords = [[0.000000, 0.000000, 1.08],
                  [0.000000, 0.000000, 0.000000],
                  [1.026719, 0.000000, -0.363000],
                  [-0.513360, -0.889165, -0.363000],
                  [-0.513360, 0.889165, -0.363000]]
        sub = Molecule(["X", "C", "H", "H", "H"], coords)
        self.mol.substitute(1, sub)
        self.assertAlmostEqual(self.mol.get_distance(0, 4), 1.54)
        f = Molecule(["X", "F"], [[0, 0, 0], [0, 0, 1.11]])
        self.mol.substitute(2, f)
        self.assertAlmostEqual(self.mol.get_distance(0, 7), 1.35)
        oh = Molecule(["X", "O", "H"],
                      [[0, 0.780362, -.456316], [0, 0, .114079],
                       [0, -.780362, -.456316]])
        self.mol.substitute(1, oh)
        self.assertAlmostEqual(self.mol.get_distance(0, 7), 1.43)
        self.mol.substitute(3, "methyl")
        self.assertEqual(self.mol.formula, "H7 C3 O1 F1")
        coords = [[0.00000, 1.40272, 0.00000],
                  [0.00000, 2.49029, 0.00000],
                  [-1.21479, 0.70136, 0.00000],
                  [-2.15666, 1.24515, 0.00000],
                  [-1.21479, -0.70136, 0.00000],
                  [-2.15666, -1.24515, 0.00000],
                  [0.00000, -1.40272, 0.00000],
                  [0.00000, -2.49029, 0.00000],
                  [1.21479, -0.70136, 0.00000],
                  [2.15666, -1.24515, 0.00000],
                  [1.21479, 0.70136, 0.00000],
                  [2.15666, 1.24515, 0.00000]]
        benzene = Molecule(["C", "H", "C", "H", "C", "H", "C", "H", "C", "H",
                            "C", "H"], coords)
        benzene.substitute(1, sub)
        self.assertEqual(benzene.formula, "H8 C7")
        #Carbon attached should be in plane.
        self.assertAlmostEqual(benzene[11].coords[2], 0)

    def test_to_from_file_string(self):
        for fmt in ["xyz", "json", "g03"]:
            s = self.mol.to(fmt=fmt)
            self.assertIsNotNone(s)
            m = Molecule.from_str(s, fmt=fmt)
            self.assertEqual(m, self.mol)
            self.assertIsInstance(m, Molecule)

        self.mol.to(filename="CH4_testing.xyz")
        self.assertTrue(os.path.exists("CH4_testing.xyz"))
        os.remove("CH4_testing.xyz")


if __name__ == '__main__':
    import unittest
    unittest.main()
