from __future__ import annotations

import json

from pymatgen.analysis.structure_prediction.substitutor import Substitutor
from pymatgen.core import Composition, Species
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

TEST_DIR = f"{TEST_FILES_DIR}/analysis/struct_predictor"


def get_table():
    """
    Loads a lightweight lambda table for use in unit tests to reduce
    initialization time, and make unit tests insensitive to changes in the
    default lambda table.
    """
    json_path = f"{TEST_DIR}/test_lambda.json"
    with open(json_path) as file:
        return json.load(file)


class TestSubstitutor(PymatgenTest):
    def setUp(self):
        self.substitutor = Substitutor(threshold=1e-3, lambda_table=get_table(), alpha=-5.0)

    def test_substitutor(self):
        s_list = [Species("O", -2), Species("Li", 1)]
        subs = self.substitutor.pred_from_list(s_list)
        assert len(subs) == 4, "incorrect number of substitutions"
        comp = Composition({"O2-": 1, "Li1+": 2})
        subs = self.substitutor.pred_from_comp(comp)
        assert len(subs) == 4, "incorrect number of substitutions"

        structures = [{"structure": PymatgenTest.get_structure("Li2O"), "id": "pmgtest"}]
        subs = self.substitutor.pred_from_structures(["Na+", "O2-"], structures)
        assert subs[0].formula == "Na2 O1"

    def test_as_dict(self):
        Substitutor.from_dict(self.substitutor.as_dict())
