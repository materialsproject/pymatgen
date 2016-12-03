# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals
import unittest2 as unittest
import os
import json

import numpy as np

from pymatgen import Lattice, Structure, Specie
from pymatgen.transformations.standard_transformations import \
    OxidationStateDecorationTransformation, SubstitutionTransformation, \
    OrderDisorderedStructureTransformation
from pymatgen.transformations.advanced_transformations import \
    SuperTransformation, EnumerateStructureTransformation, \
    MultipleSubstitutionTransformation, ChargeBalanceTransformation, \
    SubstitutionPredictorTransformation, MagOrderingTransformation, \
    DopingTransformation, _find_codopant
from monty.os.path import which
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.energy_models import IsingModel
from pymatgen.util.testing import PymatgenTest

"""
Created on Jul 24, 2012
"""


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Jul 24, 2012"


test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


def get_table():
    """
    Loads a lightweight lambda table for use in unit tests to reduce
    initialization time, and make unit tests insensitive to changes in the
    default lambda table.
    """
    data_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                            'test_files', 'struct_predictor')
    json_file = os.path.join(data_dir, 'test_lambda.json')
    with open(json_file) as f:
        lambda_table = json.load(f)
    return lambda_table


enumlib_present = which('enum.x') and which('makestr.x')


class SuperTransformationTest(unittest.TestCase):

    def test_apply_transformation(self):
        tl = [SubstitutionTransformation({"Li+": "Na+"}),
              SubstitutionTransformation({"Li+": "K+"})]
        t = SuperTransformation(tl)
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.375, 0.375, 0.375])
        coords.append([.5, .5, .5])
        coords.append([0.875, 0.875, 0.875])
        coords.append([0.125, 0.125, 0.125])
        coords.append([0.25, 0.25, 0.25])
        coords.append([0.625, 0.625, 0.625])
        coords.append([0.75, 0.75, 0.75])

        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                           [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])
        struct = Structure(lattice, ["Li+", "Li+", "Li+", "Li+", "Li+", "Li+",
                                     "O2-", "O2-"], coords)
        s = t.apply_transformation(struct, return_ranked_list=True)

        for s_and_t in s:
            self.assertEqual(s_and_t['transformation']
                             .apply_transformation(struct),
                             s_and_t['structure'])

    @unittest.skipIf(not enumlib_present, "enum_lib not present.")
    def test_apply_transformation_mult(self):
        #Test returning multiple structures from each transformation.
        disord = Structure(np.eye(3) * 4.209, [{"Cs+": 0.5, "K+": 0.5}, "Cl-"],
                           [[0, 0, 0], [0.5, 0.5, 0.5]])
        disord.make_supercell([2, 2, 1])


        tl = [EnumerateStructureTransformation(),
              OrderDisorderedStructureTransformation()]
        t = SuperTransformation(tl, nstructures_per_trans=10)
        self.assertEqual(len(t.apply_transformation(disord,
                                                    return_ranked_list=20)), 8)
        t = SuperTransformation(tl)
        self.assertEqual(len(t.apply_transformation(disord,
                                                    return_ranked_list=20)), 2)


class MultipleSubstitutionTransformationTest(unittest.TestCase):

    def test_apply_transformation(self):
        sub_dict = {1: ["Na", "K"]}
        t = MultipleSubstitutionTransformation("Li+", 0.5, sub_dict, None)
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.75, 0.75])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.25, 0.25, 0.25])
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                           [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])
        struct = Structure(lattice, ["Li+", "Li+", "O2-", "O2-"], coords)
        self.assertEqual(len(t.apply_transformation(struct,
                                                    return_ranked_list=True)),
                         2)


class ChargeBalanceTransformationTest(unittest.TestCase):

    def test_apply_transformation(self):
        t = ChargeBalanceTransformation('Li+')
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.375, 0.375, 0.375])
        coords.append([.5, .5, .5])
        coords.append([0.875, 0.875, 0.875])
        coords.append([0.125, 0.125, 0.125])
        coords.append([0.25, 0.25, 0.25])
        coords.append([0.625, 0.625, 0.625])
        coords.append([0.75, 0.75, 0.75])

        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                           [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])
        struct = Structure(lattice, ["Li+", "Li+", "Li+", "Li+", "Li+", "Li+",
                                     "O2-", "O2-"], coords)
        s = t.apply_transformation(struct)

        self.assertAlmostEqual(s.charge, 0, 5)


@unittest.skipIf(not enumlib_present, "enum_lib not present.")
class EnumerateStructureTransformationTest(unittest.TestCase):

    def test_apply_transformation(self):
        enum_trans = EnumerateStructureTransformation(refine_structure=True)
        p = Poscar.from_file(os.path.join(test_dir, 'POSCAR.LiFePO4'),
                             check_for_POTCAR=False)
        struct = p.structure
        expected_ans = [1, 3, 1]
        for i, frac in enumerate([0.25, 0.5, 0.75]):
            trans = SubstitutionTransformation({'Fe': {'Fe': frac}})
            s = trans.apply_transformation(struct)
            oxitrans = OxidationStateDecorationTransformation(
                {'Li': 1, 'Fe': 2, 'P': 5, 'O': -2})
            s = oxitrans.apply_transformation(s)
            alls = enum_trans.apply_transformation(s, 100)
            self.assertEqual(len(alls), expected_ans[i])
            self.assertIsInstance(trans.apply_transformation(s), Structure)
            for s in alls:
                self.assertIn("energy", s)

        #make sure it works for non-oxidation state decorated structure
        trans = SubstitutionTransformation({'Fe': {'Fe': 0.5}})
        s = trans.apply_transformation(struct)
        alls = enum_trans.apply_transformation(s, 100)
        self.assertEqual(len(alls), 3)
        self.assertIsInstance(trans.apply_transformation(s), Structure)
        for s in alls:
            self.assertNotIn("energy", s)

    def test_to_from_dict(self):
        trans = EnumerateStructureTransformation()
        d = trans.as_dict()
        trans = EnumerateStructureTransformation.from_dict(d)
        self.assertEqual(trans.symm_prec, 0.1)


class SubstitutionPredictorTransformationTest(unittest.TestCase):
    def test_apply_transformation(self):
        t = SubstitutionPredictorTransformation(threshold=1e-3, alpha=-5,
                                                lambda_table=get_table())
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.75, 0.75])
        coords.append([0.5, 0.5, 0.5])
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                           [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])
        struct = Structure(lattice, ['O2-', 'Li1+', 'Li1+'], coords)

        outputs = t.apply_transformation(struct, return_ranked_list=True)
        self.assertEqual(len(outputs), 4, 'incorrect number of structures')

    def test_as_dict(self):
        t = SubstitutionPredictorTransformation(threshold=2, alpha=-2,
                                                lambda_table=get_table())
        d = t.as_dict()
        t = SubstitutionPredictorTransformation.from_dict(d)
        self.assertEqual(t.threshold, 2,
                         'incorrect threshold passed through dict')
        self.assertEqual(t._substitutor.p.alpha, -2,
                         'incorrect alpha passed through dict')


@unittest.skipIf(not enumlib_present, "enum_lib not present.")
class MagOrderingTransformationTest(PymatgenTest):

    def test_apply_transformation(self):
        trans = MagOrderingTransformation({"Fe": 5})
        p = Poscar.from_file(os.path.join(test_dir, 'POSCAR.LiFePO4'),
                             check_for_POTCAR=False)
        s = p.structure
        alls = trans.apply_transformation(s, 10)
        self.assertEqual(len(alls), 3)
        f = SpacegroupAnalyzer(alls[0]["structure"], 0.1)
        self.assertEqual(f.get_space_group_number(), 31)

        model = IsingModel(5, 5)
        trans = MagOrderingTransformation({"Fe": 5},
                                          energy_model=model)
        alls2 = trans.apply_transformation(s, 10)
        #Ising model with +J penalizes similar neighbor magmom.
        self.assertNotEqual(alls[0]["structure"], alls2[0]["structure"])
        self.assertEqual(alls[0]["structure"], alls2[2]["structure"])

        s = self.get_structure('Li2O')
        #Li2O doesn't have magnetism of course, but this is to test the
        # enumeration.
        trans = MagOrderingTransformation({"Li+": 1}, max_cell_size=3)
        alls = trans.apply_transformation(s, 100)
        self.assertEqual(len(alls), 10)

    def test_ferrimagnetic(self):
        trans = MagOrderingTransformation({"Fe": 5}, 0.75, max_cell_size=1)
        p = Poscar.from_file(os.path.join(test_dir, 'POSCAR.LiFePO4'),
                             check_for_POTCAR=False)
        s = p.structure
        alls = trans.apply_transformation(s, 10)
        self.assertEqual(len(alls), 2)

    def test_as_from_dict(self):
        trans = MagOrderingTransformation({"Fe": 5}, 0.75)
        d = trans.as_dict()
        #Check json encodability
        s = json.dumps(d)
        trans = MagOrderingTransformation.from_dict(d)
        self.assertEqual(trans.mag_species_spin, {"Fe": 5})
        from pymatgen.analysis.energy_models import SymmetryModel
        self.assertIsInstance(trans.energy_model, SymmetryModel)

    def test_zero_spin_case(self):
        #ensure that zero spin case maintains sites and formula
        s = self.get_structure('Li2O')
        trans = MagOrderingTransformation({"Li+": 0.0}, 0.5)
        alls = trans.apply_transformation(s)
        #Ensure s does not have a spin property
        self.assertFalse('spin' in s.sites[0].specie._properties)
        #ensure sites are assigned a spin property in alls
        self.assertTrue('spin' in alls.sites[0].specie._properties)


@unittest.skipIf(not enumlib_present, "enum_lib not present.")
class DopingTransformationTest(PymatgenTest):

    def test_apply_transformation(self):
        structure = PymatgenTest.get_structure("LiFePO4")
        t = DopingTransformation("Ca2+", min_length=10)
        ss = t.apply_transformation(structure, 100)
        self.assertEqual(len(ss), 1)

        t = DopingTransformation("Al3+", min_length=15, ionic_radius_tol=0.1)
        ss = t.apply_transformation(structure, 100)
        self.assertEqual(len(ss), 0)

        # Aliovalent doping with vacancies
        for dopant, nstructures in [("Al3+", 4), ("N3-", 420), ("Cl-", 16)]:
            t = DopingTransformation(dopant, min_length=4, alio_tol=1,
                                     max_structures_per_enum=1000)
            ss = t.apply_transformation(structure, 1000)
            self.assertEqual(len(ss), nstructures)
            for d in ss:
                self.assertEqual(d["structure"].charge, 0)

        # Aliovalent doping with codopant
        for dopant, nstructures in [("Al3+", 3), ("N3-", 60), ("Cl-", 60)]:
            t = DopingTransformation(dopant, min_length=4, alio_tol=1,
                                     codopant=True,
                                     max_structures_per_enum=1000)
            ss = t.apply_transformation(structure, 1000)
            self.assertEqual(len(ss), nstructures)
            if __name__ == '__main__':
                for d in ss:
                    self.assertEqual(d["structure"].charge, 0)

        # Make sure compensation is done with lowest oxi state
        structure = PymatgenTest.get_structure("SrTiO3")
        t = DopingTransformation("Nb5+", min_length=5, alio_tol=1,
                                 max_structures_per_enum=1000,
                                 allowed_doping_species=["Ti4+"])
        ss = t.apply_transformation(structure, 1000)
        self.assertEqual(len(ss), 3)
        for d in ss:
            self.assertEqual(d["structure"].formula, "Sr7 Ti6 Nb2 O24")

    def test_as_from_dict(self):
        trans = DopingTransformation("Al3+", min_length=5, alio_tol=1,
                                     codopant=False, max_structures_per_enum=1)
        d = trans.as_dict()
        # Check json encodability
        s = json.dumps(d)
        trans = DopingTransformation.from_dict(d)
        self.assertEqual(str(trans.dopant), "Al3+")
        self.assertEqual(trans.max_structures_per_enum, 1)

    def test_find_codopant(self):
        self.assertEqual(_find_codopant(Specie("Fe", 2), 1), Specie("Cu", 1))
        self.assertEqual(_find_codopant(Specie("Fe", 2), 3), Specie("In", 3))

if __name__ == "__main__":
    import logging
    logging.basicConfig(level=logging.INFO)
    unittest.main()
