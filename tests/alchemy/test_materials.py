from __future__ import annotations

import json

import pytest

from pymatgen.alchemy.filters import ContainsSpecieFilter
from pymatgen.alchemy.materials import TransformedStructure
from pymatgen.core import SETTINGS
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.transformations.standard_transformations import (
    PartialRemoveSpecieTransformation,
    SubstitutionTransformation,
    SupercellTransformation,
)
from pymatgen.util.provenance import StructureNL
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest


class TestTransformedStructure(PymatgenTest):
    def setUp(self):
        structure = PymatgenTest.get_structure("LiFePO4")
        self.structure = structure
        trans = [SubstitutionTransformation({"Li": "Na"})]
        self.trans = TransformedStructure(structure, trans)

    def test_append_transformation(self):
        trafo = SubstitutionTransformation({"Fe": "Mn"})
        self.trans.append_transformation(trafo)
        assert self.trans.final_structure.composition.reduced_formula == "NaMnPO4"
        assert len(self.trans.structures) == 3
        coords = []
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        lattice = [
            [3.8401979337, 0, 0],
            [1.9200989668, 3.3257101909, 0],
            [0, -2.2171384943, 3.1355090603],
        ]
        struct = Structure(lattice, ["Si4+", "Si4+"], coords)
        ts = TransformedStructure(struct, [])
        ts.append_transformation(SupercellTransformation.from_scaling_factors(2, 1, 1))
        alt = ts.append_transformation(
            PartialRemoveSpecieTransformation("Si4+", 0.5, algo=PartialRemoveSpecieTransformation.ALGO_COMPLETE), 5
        )
        assert len(alt) == 2

    def test_append_filter(self):
        f3 = ContainsSpecieFilter(["O2-"], strict_compare=True, AND=False)
        self.trans.append_filter(f3)

    def test_get_vasp_input(self):
        SETTINGS["PMG_VASP_PSP_DIR"] = TEST_FILES_DIR
        potcar = self.trans.get_vasp_input(MPRelaxSet)["POTCAR"]
        assert "\n".join(p.symbol for p in potcar) == "Na_pv\nFe_pv\nP\nO"
        assert len(self.trans.structures) == 2

    def test_final_structure(self):
        assert self.trans.final_structure.composition.reduced_formula == "NaFePO4"

    def test_from_dict(self):
        with open(f"{TEST_FILES_DIR}/transformations.json") as file:
            dct = json.load(file)
        dct["other_parameters"] = {"tags": ["test"]}
        ts = TransformedStructure.from_dict(dct)
        ts.other_parameters["author"] = "Will"
        ts.append_transformation(SubstitutionTransformation({"Fe": "Mn"}))
        assert ts.final_structure.composition.reduced_formula == "MnPO4"
        assert ts.other_parameters == {"author": "Will", "tags": ["test"]}

    def test_undo_and_redo_last_change(self):
        trans = [
            SubstitutionTransformation({"Li": "Na"}),
            SubstitutionTransformation({"Fe": "Mn"}),
        ]
        ts = TransformedStructure(self.structure, trans)
        assert ts.final_structure.composition.reduced_formula == "NaMnPO4"
        ts.undo_last_change()
        assert ts.final_structure.composition.reduced_formula == "NaFePO4"
        ts.undo_last_change()
        assert ts.final_structure.composition.reduced_formula == "LiFePO4"
        with pytest.raises(IndexError, match="No more changes to undo"):
            ts.undo_last_change()
        ts.redo_next_change()
        assert ts.final_structure.composition.reduced_formula == "NaFePO4"
        ts.redo_next_change()
        assert ts.final_structure.composition.reduced_formula == "NaMnPO4"
        with pytest.raises(IndexError, match="No more changes to redo"):
            ts.redo_next_change()
        # Make sure that this works with filters.
        f3 = ContainsSpecieFilter(["O2-"], strict_compare=True, AND=False)
        ts.append_filter(f3)
        ts.undo_last_change()
        ts.redo_next_change()

    def test_as_dict(self):
        self.trans.set_parameter("author", "will")
        d = self.trans.as_dict()
        assert "last_modified" in d
        assert "history" in d
        assert "author" in d["other_parameters"]
        assert Structure.from_dict(d).formula == "Na4 Fe4 P4 O16"

    def test_snl(self):
        self.trans.set_parameter("author", "will")
        with pytest.warns(UserWarning) as warns:
            snl = self.trans.to_snl([("will", "will@test.com")])
        assert len(warns) == 1, "Warning not raised on type conversion with other_parameters"

        ts = TransformedStructure.from_snl(snl)
        assert ts.history[-1]["@class"] == "SubstitutionTransformation"

        hist = ("testname", "testURL", {"test": "testing"})
        snl = StructureNL(ts.final_structure, [("will", "will@test.com")], history=[hist])
        snl = TransformedStructure.from_snl(snl).to_snl([("notwill", "notwill@test.com")])
        assert snl.history == [hist]
        assert snl.authors == [("notwill", "notwill@test.com")]
