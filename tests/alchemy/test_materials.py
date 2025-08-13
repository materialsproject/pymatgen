from __future__ import annotations

from copy import deepcopy

import orjson
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
from pymatgen.util.testing import FAKE_POTCAR_DIR, TEST_FILES_DIR, MatSciTest

TEST_DIR = f"{TEST_FILES_DIR}/alchemy"


class TestTransformedStructure(MatSciTest):
    def setup_method(self):
        structure = MatSciTest.get_structure("LiFePO4")
        self.structure = structure
        trafos = [SubstitutionTransformation({"Li": "Na"})]
        self.trans = TransformedStructure(structure, trafos)

    def test_append_transformation(self):
        trafo = SubstitutionTransformation({"Fe": "Mn"})
        self.trans.append_transformation(trafo)
        assert self.trans.final_structure.reduced_formula == "NaMnPO4"
        assert len(self.trans.structures) == 3
        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        lattice = [
            [3.8401979337, 0, 0],
            [1.9200989668, 3.3257101909, 0],
            [0, -2.2171384943, 3.1355090603],
        ]
        struct = Structure(lattice, ["Si4+", "Si4+"], coords)
        t_struct = TransformedStructure(struct, [])
        t_struct.append_transformation(SupercellTransformation.from_scaling_factors(2, 1, 1))
        alt = t_struct.append_transformation(
            PartialRemoveSpecieTransformation("Si4+", 0.5, algo=PartialRemoveSpecieTransformation.ALGO_COMPLETE),
            5,
        )
        assert len(alt) == 2

    def test_append_filter(self):
        f3 = ContainsSpecieFilter(["O2-"], strict_compare=True, AND=False)
        self.trans.append_filter(f3)

    def test_get_vasp_input(self):
        SETTINGS["PMG_VASP_PSP_DIR"] = FAKE_POTCAR_DIR
        potcar = self.trans.get_vasp_input(MPRelaxSet)["POTCAR"]
        assert "\n".join(p.symbol for p in potcar) == "Na_pv\nFe_pv\nP\nO"
        assert len(self.trans.structures) == 2

    def test_final_structure(self):
        assert self.trans.final_structure.reduced_formula == "NaFePO4"
        # https://github.com/materialsproject/pymatgen/pull/3617
        assert isinstance(deepcopy(self.trans), TransformedStructure)

    def test_from_dict(self):
        with open(f"{TEST_DIR}/transformations.json", "rb") as file:
            dct = orjson.loads(file.read())
        dct["other_parameters"] = {"tags": ["test"]}
        t_struct = TransformedStructure.from_dict(dct)
        t_struct.other_parameters["author"] = "Will"
        t_struct.append_transformation(SubstitutionTransformation({"Fe": "Mn"}))
        assert t_struct.final_structure.reduced_formula == "MnPO4"
        assert t_struct.other_parameters == {"author": "Will", "tags": ["test"]}

    def test_undo_and_redo_last_change(self):
        trafos = [
            SubstitutionTransformation({"Li": "Na"}),
            SubstitutionTransformation({"Fe": "Mn"}),
        ]
        t_struct = TransformedStructure(self.structure, trafos)
        assert t_struct.final_structure.reduced_formula == "NaMnPO4"
        t_struct.undo_last_change()
        assert t_struct.final_structure.reduced_formula == "NaFePO4"
        t_struct.undo_last_change()
        assert t_struct.final_structure.reduced_formula == "LiFePO4"
        with pytest.raises(IndexError, match="No more changes to undo"):
            t_struct.undo_last_change()
        t_struct.redo_next_change()
        assert t_struct.final_structure.reduced_formula == "NaFePO4"
        t_struct.redo_next_change()
        assert t_struct.final_structure.reduced_formula == "NaMnPO4"
        with pytest.raises(IndexError, match="No more changes to redo"):
            t_struct.redo_next_change()
        # Make sure that this works with filters.
        f3 = ContainsSpecieFilter(["O2-"], strict_compare=True, AND=False)
        t_struct.append_filter(f3)
        t_struct.undo_last_change()
        t_struct.redo_next_change()

    def test_set_parameter(self):
        trans = self.trans.set_parameter("author", "will")
        assert trans.other_parameters["author"] == "will"
        assert trans is self.trans

    def test_as_dict(self):
        self.trans.set_parameter("author", "will")
        dct = self.trans.as_dict()
        assert "last_modified" in dct
        assert "history" in dct
        assert "author" in dct["other_parameters"]
        assert Structure.from_dict(dct).formula == "Na4 Fe4 P4 O16"

    def test_snl(self):
        self.trans.set_parameter("author", "will")
        with pytest.warns(UserWarning, match="discarded during type conversion to SNL") as warns:
            struct_nl = self.trans.to_snl([("will", "will@test.com")])

        assert len(warns) >= 1, f"Warning not raised on type conversion with other_parameters {len(warns)=}"
        assert (
            str(warns[0].message)
            == "Data in TransformedStructure.other_parameters discarded during type conversion to SNL"
        )

        t_struct = TransformedStructure.from_snl(struct_nl)
        assert t_struct.history[-1]["@class"] == "SubstitutionTransformation"

        hist = ("testname", "testURL", {"test": "testing"})
        struct_nl = StructureNL(t_struct.final_structure, [("will", "will@test.com")], history=[hist])
        t_struct = TransformedStructure.from_snl(struct_nl).to_snl([("notwill", "notwill@test.com")])
        assert t_struct.history == [hist]
        assert t_struct.authors == [("notwill", "notwill@test.com")]
