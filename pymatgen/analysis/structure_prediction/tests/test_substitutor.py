# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest
import os
import json

from pymatgen.core.periodic_table import Specie
from pymatgen.core.composition import Composition
from pymatgen.analysis.structure_prediction.substitutor import Substitutor
from pymatgen.util.testing import PymatgenTest


def get_table():
    """
    Loads a lightweight lambda table for use in unit tests to reduce
    initialization time, and make unit tests insensitive to changes in the
    default lambda table.
    """
    data_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                            'test_files', "struct_predictor")
    json_file = os.path.join(data_dir, 'test_lambda.json')
    with open(json_file) as f:
        lambda_table = json.load(f)
    return lambda_table


class SubstitutorTest(PymatgenTest):

    def setUp(self):
        self.s = Substitutor(threshold=1e-3, lambda_table=get_table(),
                             alpha= -5.)

    def test_substitutor(self):
        s_list = [Specie('O', -2), Specie('Li', 1)]
        subs = self.s.pred_from_list(s_list)
        self.assertEqual(len(subs), 4
                         , 'incorrect number of substitutions')
        c = Composition({'O2-': 1, 'Li1+': 2})
        subs = self.s.pred_from_comp(c)
        self.assertEqual(len(subs), 4
                         , 'incorrect number of substitutions')

        structures = [{"structure": PymatgenTest.get_structure("Li2O"),
                       "id": "pmgtest"}]
        subs = self.s.pred_from_structures(["Na+", "O2-"], structures)
        self.assertEqual(subs[0].formula, "Na2 O1")

    def test_as_dict(self):
        Substitutor.from_dict(self.s.as_dict())


if __name__ == "__main__":
    unittest.main()
