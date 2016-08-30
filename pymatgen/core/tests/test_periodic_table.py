# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import unittest2 as unittest
import pickle

from pymatgen.util.testing import PymatgenTest
from pymatgen.core.periodic_table import Element, Specie, DummySpecie, get_el_sp
from pymatgen.core.composition import Composition
from copy import deepcopy


class ElementTestCase(PymatgenTest):

    def test_init(self):
        self.assertEqual("Fe", Element("Fe").symbol, "Fe test failed")

        fictional_symbols = ["D", "T", "Zebra"]

        for sym in fictional_symbols:
            self.assertRaises(ValueError, Element, sym)

        # Test caching
        self.assertEqual(id(Element("Fe")), id(Element("Fe")))

    def test_dict(self):
        fe = Element.Fe
        d = fe.as_dict()
        self.assertEqual(fe, Element.from_dict(d))

    def test_block(self):
        testsets = {"O": "p", "Fe": "d", "Li": "s", "U": "f", "Er": "f",
                    "Lu": "d", "Lr": "d"}
        for k, v in testsets.items():
            self.assertEqual(Element(k).block, v)

    def test_full_electronic_structure(self):
        testsets = {"O": [(1, "s", 2), (2, "s", 2), (2, "p", 4)],
                    "Fe": [(1, "s", 2), (2, "s", 2), (2, "p", 6), (3, "s", 2),
                           (3, "p", 6), (3, "d", 6), (4, "s", 2)],
                    "Li": [(1, "s", 2), (2, "s", 1)],
                    "U": [(1, "s", 2), (2, "s", 2), (2, "p", 6), (3, "s", 2),
                          (3, "p", 6), (3, "d", 10), (4, "s", 2), (4, "p", 6),
                          (4, "d", 10), (5, "s", 2), (5, "p", 6), (4, "f", 14),
                          (5, "d", 10), (6, "s", 2), (6, "p", 6), (5, "f", 3),
                          (6, "d", 1), (7, "s", 2)]}
        for k, v in testsets.items():
            self.assertEqual(Element(k).full_electronic_structure, v)

    def test_attributes(self):
        is_true = {("Xe", "Kr"): "is_noble_gas",
                   ("Fe", "Ni"): "is_transition_metal",
                   ("Li", "Cs"): "is_alkali",
                   ("Ca", "Mg"): "is_alkaline",
                   ("F", "Br", "I"): "is_halogen",
                   ("La",): "is_lanthanoid",
                   ("U", "Pu"): "is_actinoid",
                   ("Si", "Ge"): "is_metalloid",
                   ("O", "Te"): "is_chalcogen"}

        for k, v in is_true.items():
            for sym in k:
                self.assertTrue(getattr(Element(sym), v), sym + " is false")

        keys = ["name", "mendeleev_no", "atomic_mass",
                "electronic_structure", "X", "atomic_radius",
                "min_oxidation_state", "max_oxidation_state",
                "electrical_resistivity", "velocity_of_sound", "reflectivity",
                "refractive_index", "poissons_ratio", "molar_volume",
                "thermal_conductivity", "melting_point", "boiling_point",
                "liquid_range", "critical_temperature",
                "superconduction_temperature",
                "bulk_modulus", "youngs_modulus", "brinell_hardness",
                "rigidity_modulus", "mineral_hardness",
                "vickers_hardness", "density_of_solid",
                "coefficient_of_linear_thermal_expansion", "oxidation_states",
                "common_oxidation_states", "average_ionic_radius",
                "ionic_radii"]

        # Test all elements up to Uranium
        for i in range(1, 93):
            el = Element.from_Z(i)
            d = el.data
            for k in keys:
                k_str = k.capitalize().replace("_", " ")
                if k_str in d and (not str(d[k_str]).startswith("no data")):
                    self.assertIsNotNone(getattr(el, k))
            el = Element.from_Z(i)
            if len(el.oxidation_states) > 0:
                self.assertEqual(max(el.oxidation_states),
                                 el.max_oxidation_state)
                self.assertEqual(min(el.oxidation_states),
                                 el.min_oxidation_state)

        self.assertRaises(ValueError, Element.from_Z, 1000)


    def test_oxidation_states(self):
        el = Element.Fe
        self.assertEqual(el.oxidation_states, (-2, -1, 1, 2, 3, 4, 5, 6))
        self.assertEqual(el.common_oxidation_states, (2, 3))

    def test_deepcopy(self):
        el1 = Element.Fe
        el2 = Element.Na
        ellist = [el1, el2]
        self.assertEqual(ellist, deepcopy(ellist),
                         "Deepcopy operation doesn't produce exact copy")

    def test_radii(self):
        el = Element.Pd
        self.assertEqual(el.atomic_radius, 1.40)
        self.assertEqual(el.atomic_radius_calculated, 1.69)
        self.assertEqual(el.van_der_waals_radius, 1.63)

    def test_data(self):
        self.assertEqual(Element.Pd.data["Atomic radius"], 1.4)
        al = Element.Al
        val = al.thermal_conductivity
        self.assertEqual(val, 235)
        self.assertEqual(str(val.unit), "W K^-1 m^-1")
        val = al.electrical_resistivity
        self.assertEqual(val, 2.7e-08)
        self.assertEqual(str(val.unit), "m ohm")

    def test_sort(self):
        els = [Element.Se, Element.C]
        self.assertEqual(sorted(els), [Element.C, Element.Se])

    def test_pickle(self):
        el1 = Element.Fe
        o = pickle.dumps(el1)
        self.assertEqual(el1, pickle.loads(o))

        #Test all elements up to Uranium
        for i in range(1, 93):
            self.serialize_with_pickle(Element.from_Z(i), test_eq=True)

    def test_print_periodic_table(self):
        Element.print_periodic_table()


class SpecieTestCase(PymatgenTest):

    def setUp(self):
        self.specie1 = Specie.from_string("Fe2+")
        self.specie2 = Specie("Fe", 3)
        self.specie3 = Specie("Fe", 2)
        self.specie4 = Specie("Fe", 2, {"spin": 5})

    def test_init(self):
        self.assertRaises(ValueError, Specie, "Fe", 2, {"magmom": 5})

    def test_cached(self):
        specie5 = Specie("Fe", 2)
        self.assertEqual(id(specie5), id(self.specie3))

    def test_ionic_radius(self):
        self.assertEqual(self.specie2.ionic_radius, 78.5 / 100)
        self.assertEqual(self.specie3.ionic_radius, 92 / 100)
        self.assertAlmostEqual(Specie("Mn", 4).ionic_radius, 0.67)

    def test_eq(self):
        self.assertEqual(self.specie1, self.specie3,
                         "Static and actual constructor gives unequal result!")
        self.assertNotEqual(self.specie1, self.specie2,
                            "Fe2+ should not be equal to Fe3+")
        self.assertNotEqual(self.specie4, self.specie3)
        self.assertFalse(self.specie1 == Element("Fe"))
        self.assertFalse(Element("Fe") == self.specie1)

    def test_cmp(self):
        self.assertLess(self.specie1, self.specie2, "Fe2+ should be < Fe3+")
        self.assertLess(Specie("C", 1), Specie("Se", 1))

    def test_attr(self):
        self.assertEqual(self.specie1.Z, 26,
                         "Z attribute for Fe2+ should be = Element Fe.")
        self.assertEqual(self.specie4.spin, 5)

    def test_deepcopy(self):
        el1 = Specie("Fe", 4)
        el2 = Specie("Na", 1)
        ellist = [el1, el2]
        self.assertEqual(ellist, deepcopy(ellist),
                         "Deepcopy operation doesn't produce exact copy.")

    def test_pickle(self):
        self.assertEqual(self.specie1, pickle.loads(pickle.dumps(self.specie1)))
        for i in range(1, 5):
            self.serialize_with_pickle(getattr(self, "specie%d" % i) , test_eq=True)

    def test_get_crystal_field_spin(self):
        self.assertEqual(Specie("Fe", 2).get_crystal_field_spin(), 4)
        self.assertEqual(Specie("Fe", 3).get_crystal_field_spin(), 5)
        self.assertEqual(Specie("Fe", 4).get_crystal_field_spin(), 4)
        self.assertEqual(Specie("Co", 3).get_crystal_field_spin(
            spin_config="low"), 0)
        self.assertEqual(Specie("Co", 4).get_crystal_field_spin(
            spin_config="low"), 1)
        self.assertEqual(Specie("Ni", 3).get_crystal_field_spin(
            spin_config="low"), 1)
        self.assertEqual(Specie("Ni", 4).get_crystal_field_spin(
            spin_config="low"), 0)

        self.assertRaises(AttributeError,
                          Specie("Li", 1).get_crystal_field_spin)
        self.assertRaises(AttributeError,
                          Specie("Ge", 4).get_crystal_field_spin)
        self.assertRaises(AttributeError,
                          Specie("H", 1).get_crystal_field_spin)
        self.assertRaises(AttributeError,
                          Specie("Fe", 10).get_crystal_field_spin)
        self.assertRaises(ValueError, Specie("Fe", 2).get_crystal_field_spin,
                          "hex")

        s = Specie("Co", 3).get_crystal_field_spin("tet", spin_config="low")
        self.assertEqual(s, 2)

    def test_sort(self):
        els = map(get_el_sp, ["N3-", "Si4+", "Si3+"])
        self.assertEqual(sorted(els), [Specie("Si", 3), Specie("Si", 4),
                                       Specie("N", -3)])

    def test_to_from_string(self):
        fe3 = Specie("Fe", 3, {"spin": 5})
        self.assertEqual(str(fe3), "Fe3+spin=5")
        fe = Specie.from_string("Fe3+spin=5")
        self.assertEqual(fe.spin, 5)
        mo0 = Specie("Mo", 0, {"spin": 5})
        self.assertEqual(str(mo0), "Mo0+spin=5")
        mo = Specie.from_string("Mo0+spin=4")
        self.assertEqual(mo.spin, 4)


class DummySpecieTestCase(unittest.TestCase):

    def test_init(self):
        self.specie1 = DummySpecie("X")
        self.assertRaises(ValueError, DummySpecie, "Xe")
        self.assertRaises(ValueError, DummySpecie, "Xec")
        self.assertRaises(ValueError, DummySpecie, "Vac")
        self.specie2 = DummySpecie("X", 2, {"spin": 3})
        self.assertEqual(self.specie2.spin, 3)

    def test_cached(self):
        sp1 = DummySpecie("X", 2)
        sp2 = DummySpecie("X", 2)
        self.assertEqual(id(sp1), id(sp2))

    def test_eq(self):
        self.assertFalse(DummySpecie("Xg") == DummySpecie("Xh"))
        self.assertFalse(DummySpecie("Xg") == DummySpecie("Xg", 3))
        self.assertTrue(DummySpecie("Xg", 3) == DummySpecie("Xg", 3))

    def test_from_string(self):
        sp = DummySpecie.from_string("X")
        self.assertEqual(sp.oxi_state, 0)
        sp = DummySpecie.from_string("X2+")
        self.assertEqual(sp.oxi_state, 2)
        sp = DummySpecie.from_string("X2+spin=5")
        self.assertEqual(sp.oxi_state, 2)
        self.assertEqual(sp.spin, 5)

    def test_pickle(self):
        el1 = DummySpecie("X", 3)
        o = pickle.dumps(el1)
        self.assertEqual(el1, pickle.loads(o))

    def test_sort(self):
        r = sorted([Element.Fe, DummySpecie("X")])
        self.assertEqual(r, [DummySpecie("X"), Element.Fe])
        self.assertTrue(DummySpecie("X", 3) < DummySpecie("X", 4))

    def test_safe_from_composition(self):
        c = Composition({'Xa': 1, 'Fe': 1})
        self.assertEqual(DummySpecie.safe_from_composition(c).symbol, 'Xb')
        self.assertEqual(DummySpecie.safe_from_composition(c, 1).symbol, 'Xb')


class FuncTest(unittest.TestCase):

    def test_get_el_sp(self):
        self.assertEqual(get_el_sp("Fe2+"), Specie("Fe", 2))
        self.assertEqual(get_el_sp("3"), Element.Li)
        self.assertEqual(get_el_sp("3.0"), Element.Li)
        self.assertEqual(get_el_sp("U"), Element.U)
        self.assertEqual(get_el_sp("X2+"), DummySpecie("X", 2))
        self.assertEqual(get_el_sp("Mn3+"), Specie("Mn", 3))

if __name__ == "__main__":
    unittest.main()
