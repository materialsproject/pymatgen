# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals
import unittest
import os
import warnings
from pymatgen.command_line.enumlib_caller import EnumlibAdaptor
from pymatgen import Element, Structure
from pymatgen.transformations.standard_transformations import \
    SubstitutionTransformation
from monty.os.path import which
from pymatgen.transformations.site_transformations import \
    RemoveSitesTransformation
from pymatgen.util.testing import PymatgenTest
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

"""
Created on Jul 22, 2012
"""

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Jul 22, 2012"

enum_cmd = which('enum.x') or which('multienum.x')
makestr_cmd = which('makestr.x') or which('makeStr.x') or which('makeStr.py')
enumlib_present = enum_cmd and makestr_cmd
test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


@unittest.skipIf(not enumlib_present, "enum_lib not present.")
class EnumlibAdaptorTest(PymatgenTest):
    def test_init(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                    'test_files')
            struct = self.get_structure("LiFePO4")
            subtrans = SubstitutionTransformation({'Li': {'Li': 0.5}})
            adaptor = EnumlibAdaptor(subtrans.apply_transformation(struct), 1, 2)
            adaptor.run()
            structures = adaptor.structures
            self.assertEqual(len(structures), 86)
            for s in structures:
                self.assertAlmostEqual(
                    s.composition.get_atomic_fraction(Element("Li")), 0.5 / 6.5)
            adaptor = EnumlibAdaptor(subtrans.apply_transformation(struct), 1, 2,
                                     refine_structure=True)
            adaptor.run()
            structures = adaptor.structures
            self.assertEqual(len(structures), 52)

            subtrans = SubstitutionTransformation({'Li': {'Li': 0.25}})
            adaptor = EnumlibAdaptor(subtrans.apply_transformation(struct), 1, 1,
                                     refine_structure=True)
            adaptor.run()
            structures = adaptor.structures
            self.assertEqual(len(structures), 1)
            for s in structures:
                self.assertAlmostEqual(s.composition
                                       .get_atomic_fraction(Element("Li")),
                                       0.25 / 6.25)

            # Make sure it works for completely disordered structures.
            struct = Structure([[10, 0, 0], [0, 10, 0], [0, 0, 10]], [{'Fe': 0.5}],
                               [[0, 0, 0]])
            adaptor = EnumlibAdaptor(struct, 1, 2)
            adaptor.run()
            self.assertEqual(len(adaptor.structures), 3)

            # Make sure it works properly when symmetry is broken by ordered sites.
            struct = self.get_structure("LiFePO4")
            subtrans = SubstitutionTransformation({'Li': {'Li': 0.25}})
            s = subtrans.apply_transformation(struct)
            # REmove some ordered sites to break symmetry.
            removetrans = RemoveSitesTransformation([4, 7])
            s = removetrans.apply_transformation(s)
            adaptor = EnumlibAdaptor(s, 1, 1, enum_precision_parameter=0.01)
            adaptor.run()
            structures = adaptor.structures
            self.assertEqual(len(structures), 4)

            struct = Structure([[3, 0, 0], [0, 3, 0], [0, 0, 3]],
                               [{"Si": 0.5}] * 2, [[0, 0, 0], [0.5, 0.5, 0.5]])
            adaptor = EnumlibAdaptor(struct, 1, 3, enum_precision_parameter=0.01)
            adaptor.run()
            structures = adaptor.structures
            self.assertEqual(len(structures), 10)

            struct = Structure.from_file(
                os.path.join(test_dir, "EnumerateTest.json"))
            adaptor = EnumlibAdaptor(struct, 1, 1)
            adaptor.run()
            structures = adaptor.structures
            self.assertEqual(len(structures), 2)

    def test_partial_disorder(self):
        s = Structure.from_file(filename=os.path.join(test_dir, "garnet.cif"))
        a = SpacegroupAnalyzer(s, 0.1)
        prim = a.find_primitive()
        s = prim.copy()
        s["Al3+"] = {"Al3+": 0.5, "Ga3+": 0.5}
        adaptor = EnumlibAdaptor(s, 1, 1, enum_precision_parameter=0.01)
        adaptor.run()
        structures = adaptor.structures
        self.assertEqual(len(structures), 7)
        for s in structures:
            self.assertEqual(s.formula, 'Ca12 Al4 Ga4 Si12 O48')
        s = prim.copy()
        s["Ca2+"] = {"Ca2+": 1/3, "Mg2+": 2/3}
        adaptor = EnumlibAdaptor(s, 1, 1, enum_precision_parameter=0.01)
        adaptor.run()
        structures = adaptor.structures
        self.assertEqual(len(structures), 20)
        for s in structures:
            self.assertEqual(s.formula, 'Ca4 Mg8 Al8 Si12 O48')

        s = prim.copy()
        s["Si4+"] = {"Si4+": 1/3, "Ge4+": 2/3}
        adaptor = EnumlibAdaptor(s, 1, 1, enum_precision_parameter=0.01)
        adaptor.run()
        structures = adaptor.structures
        self.assertEqual(len(structures), 18)
        for s in structures:
            self.assertEqual(s.formula, 'Ca12 Al8 Si4 Ge8 O48')


if __name__ == '__main__':
    unittest.main()
