# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


from __future__ import annotations

import os
import warnings

from pymatgen.alchemy.filters import ContainsSpecieFilter
from pymatgen.alchemy.transmuters import CifTransmuter, PoscarTransmuter
from pymatgen.transformations.advanced_transformations import SuperTransformation
from pymatgen.transformations.standard_transformations import (
    OrderDisorderedStructureTransformation,
    RemoveSpeciesTransformation,
    SubstitutionTransformation,
)
from pymatgen.util.testing import PymatgenTest


class CifTransmuterTest(PymatgenTest):
    def setUp(self):
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_init(self):
        trans = []
        trans.append(SubstitutionTransformation({"Fe": "Mn", "Fe2+": "Mn2+"}))
        tsc = CifTransmuter.from_filenames([os.path.join(self.TEST_FILES_DIR, "MultiStructure.cif")], trans)
        assert len(tsc) == 2
        expected_ans = {"Mn", "O", "Li", "P"}
        for s in tsc:
            els = {el.symbol for el in s.final_structure.composition.elements}
            assert expected_ans == els


class PoscarTransmuterTest(PymatgenTest):
    def test_init(self):
        trans = []
        trans.append(SubstitutionTransformation({"Fe": "Mn"}))
        tsc = PoscarTransmuter.from_filenames(
            [os.path.join(self.TEST_FILES_DIR, "POSCAR"), os.path.join(self.TEST_FILES_DIR, "POSCAR")], trans
        )
        assert len(tsc) == 2
        expected_ans = {"Mn", "O", "P"}
        for s in tsc:
            els = {el.symbol for el in s.final_structure.composition.elements}
            assert expected_ans == els

    def test_transmuter(self):
        tsc = PoscarTransmuter.from_filenames([os.path.join(self.TEST_FILES_DIR, "POSCAR")])
        tsc.append_transformation(RemoveSpeciesTransformation("O"))
        assert len(tsc[0].final_structure) == 8

        tsc.append_transformation(SubstitutionTransformation({"Fe": {"Fe2+": 0.25, "Mn3+": 0.75}, "P": "P5+"}))
        tsc.append_transformation(OrderDisorderedStructureTransformation(), extend_collection=50)
        assert len(tsc) == 4

        t = SuperTransformation(
            [
                SubstitutionTransformation({"Fe2+": "Mg2+"}),
                SubstitutionTransformation({"Fe2+": "Zn2+"}),
                SubstitutionTransformation({"Fe2+": "Be2+"}),
            ]
        )
        tsc.append_transformation(t, extend_collection=True)
        assert len(tsc) == 12
        for x in tsc:
            # should be 4 trans + starting structure
            assert len(x) == 5, "something might be wrong with the number of transformations in the history"

        # test the filter
        tsc.apply_filter(ContainsSpecieFilter(["Zn2+", "Be2+", "Mn4+"], strict_compare=True, AND=False))
        assert len(tsc) == 8
        assert tsc.transformed_structures[0].as_dict()["history"][-1]["@class"] == "ContainsSpecieFilter"

        tsc.apply_filter(ContainsSpecieFilter(["Be2+"]))
        assert len(tsc) == 4

        # Test set_parameter and add_tag.
        tsc.set_parameter("para1", "hello")
        assert tsc.transformed_structures[0].as_dict()["other_parameters"]["para1"] == "hello"
        tsc.add_tags(["world", "universe"])
        assert tsc.transformed_structures[0].as_dict()["other_parameters"]["tags"] == ["world", "universe"]


if __name__ == "__main__":
    import unittest

    unittest.main()
