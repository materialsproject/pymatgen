"""Unit tests for StructureNL (SNL) format."""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
from unittest import TestCase

import numpy as np
import pytest

from pymatgen.core.structure import Molecule, Structure
from pymatgen.util.provenance import Author, HistoryNode, StructureNL

__author__ = "Anubhav Jain"
__credits__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Anubhav Jain"
__email__ = "ajain@lbl.gov"
__date__ = "2/14/13"


class TestStructureNL(TestCase):
    def setUp(self):
        # set up a Structure
        self.struct = Structure(np.eye(3, 3) * 3, ["Fe"], [[0, 0, 0]])
        self.s2 = Structure(np.eye(3, 3) * 3, ["Al"], [[0, 0, 0]])
        self.mol = Molecule(["He"], [[0, 0, 0]])
        # set up BibTeX strings
        self.matproj = "@misc{MaterialsProject,\ntitle = {{Materials Project}},\nurl = {http://materialsproject.org}\n}"
        self.pmg = (
            "@article{Ong2013,\n author = {Ong, "
            "Shyue Ping and Richards, William Davidson and Jain, "
            "Anubhav and Hautier, Geoffroy and Kocher, "
            "Michael and Cholia, Shreyas and Gunter, Dan and Chevrier,"
            " Vincent L. and Persson, Kristin A. and Ceder, Gerbrand},"
            "\n doi = {10.1016/j.commatsci.2012.10.028},"
            "\n issn = {09270256},\n journal = {Computational "
            "Materials Science},\n month = feb,\n pages = {314--319},"
            "\n publisher = {Elsevier B.V.},"
            "\n title = {{Python Materials Genomics (pymatgen): A "
            "robust, open-source python library for materials "
            "analysis}},\n url = {http://linkinghub.elsevier"
            ".com/retrieve/pii/S0927025612006295},\n volume = {68},"
            "\n year = {2013}\n}"
        )
        repeat = "REPEAT" * 10000
        self.superlong = f"@misc{{SuperLong,\ntitle = {{{repeat}}}}}"
        self.unicode_title = "@misc{Unicode_Title,\ntitle = {{A \u73ab is a rose}}}"
        self.junk = "This is junk text, not a BibTeX reference"

        # set up remarks
        self.remark_fail = [
            "This is a really long remark that is clearly invalid and must fail, don't you agree? It would be silly "
            "to allow remarks that went on forever and ever."
        ]

        # set up some authors
        self.hulk = [{"name": "Hulk", "email": "hulk@avengers.com"}]
        self.america = "Captain America <captainamerica@avengers.com>"
        self.thor = [("Thor", "thor@avengers.com")]
        self.duo = "Iron Man <ironman@avengers.com>, Black Widow <blackwidow@avengers.com>"

        # set up HistoryNodes
        self.valid_node = HistoryNode("DB 1", "www.db1URLgoeshere.com", {"db1_id": 12424})
        self.valid_node2 = {
            "name": "DB 2",
            "url": "www.db2URLgoeshere.com",
            "description": {"db2_id": 12424},
        }
        self.invalid_node = {"name": "DB 3", "url": "http://www.db3isnotavalidnode.com"}

    def test_authors(self):
        struct_nl = StructureNL(self.struct, self.hulk, references=self.pmg)
        assert struct_nl.authors[0].name == "Hulk"
        assert struct_nl.authors[0].email == "hulk@avengers.com"

        struct_nl = StructureNL(self.struct, self.america, references=self.pmg)
        assert struct_nl.authors[0].name == "Captain America"
        assert struct_nl.authors[0].email == "captainamerica@avengers.com"

        struct_nl = StructureNL(self.struct, self.thor, references=self.pmg)
        assert struct_nl.authors[0].name == "Thor"
        assert struct_nl.authors[0].email == "thor@avengers.com"

        struct_nl = StructureNL(self.struct, self.duo, references=self.pmg)
        assert struct_nl.authors[0].name == "Iron Man"
        assert struct_nl.authors[0].email == "ironman@avengers.com"
        assert struct_nl.authors[1].name == "Black Widow"
        assert struct_nl.authors[1].email == "blackwidow@avengers.com"
        StructureNL(self.struct, self.hulk, references=self.pmg)

    def test_references(self):
        # An empty string should be OK
        StructureNL(self.struct, self.hulk, references="")

        # An empty list should not work
        with pytest.raises(
            TypeError,
            match="Invalid format for SNL reference! Should be empty string or BibTeX string.",
        ):
            StructureNL(self.struct, self.hulk, references=[])

        # junk reference should not work
        with pytest.raises(ValueError, match="Invalid format for SNL reference! Should be BibTeX string."):
            StructureNL(self.struct, self.hulk, references=self.junk)

        # good references should be ok
        StructureNL(self.struct, self.hulk, references=self.pmg)

        # unicode references should work
        StructureNL(self.struct, self.hulk, references=self.unicode_title)

        # multi-line references should be OK
        StructureNL(self.struct, self.hulk, references=f"{self.matproj}\n{self.pmg}")

        # super long references are bad
        with pytest.raises(ValueError, match="The BibTeX string must be fewer than 20000 chars, you have 60028"):
            StructureNL(self.struct, self.hulk, references=self.superlong)

    def test_history_nodes(self):
        struct_nl = StructureNL(self.struct, self.hulk, history=[self.valid_node])
        assert struct_nl.history[0].name == "DB 1"
        assert struct_nl.history[0].url == "www.db1URLgoeshere.com"
        assert struct_nl.history[0].description == {"db1_id": 12424}

        struct_nl = StructureNL(self.struct, self.hulk, history=[self.valid_node, self.valid_node2])
        assert struct_nl.history[1].name == "DB 2"
        assert struct_nl.history[1].url == "www.db2URLgoeshere.com"
        assert struct_nl.history[1].description == {"db2_id": 12424}

        # invalid nodes should not work
        with pytest.raises(KeyError, match="description"):
            StructureNL(self.struct, self.hulk, history=[self.invalid_node])

        # too many nodes should not work
        n_nodes = 1000
        with pytest.raises(ValueError, match=f"A maximum of 100 History nodes are supported, you have {n_nodes}!"):
            StructureNL(self.struct, self.hulk, history=[self.valid_node] * n_nodes)

    def test_data(self):
        # Structure data is OK due to PMGEncoder/Decoder
        struct_nl = StructureNL(self.struct, self.hulk, data={"_structure": self.s2})
        assert struct_nl.data["_structure"] == self.s2, "Data storage is broken"
        with pytest.raises(ValueError, match="data must contain properly namespaced data with keys starting "):
            StructureNL(self.struct, self.hulk, data={"bad_key": 1})

    def test_remarks(self):
        struct_nl = StructureNL(self.struct, self.hulk, remarks="string format")
        assert struct_nl.remarks[0] == "string format"
        with pytest.raises(ValueError, match="The remark exceeds the maximum size of 140 characters: 150"):
            StructureNL(self.struct, self.hulk, remarks=self.remark_fail)

    def test_eq(self):
        # test basic Equal()
        created_at = datetime.now(tz=timezone.utc)
        struct_nl = StructureNL(
            self.struct,
            self.hulk,
            ["test_project"],
            self.pmg,
            ["remark1"],
            {"_my_data": self.s2},
            [self.valid_node, self.valid_node2],
            created_at,
        )
        struct_nl2 = StructureNL(
            self.struct,
            self.hulk,
            ["test_project"],
            self.pmg,
            ["remark1"],
            {"_my_data": self.s2},
            [self.valid_node, self.valid_node2],
            created_at,
        )
        assert struct_nl == struct_nl2

        # change the created at date, now they are no longer equal
        created_at = datetime.now(tz=timezone.utc) + timedelta(days=-1)
        snl_new_date = StructureNL(
            self.struct,
            self.hulk,
            ["test_project"],
            self.pmg,
            ["remark1"],
            {"_my_data": self.s2},
            [self.valid_node, self.valid_node2],
            created_at,
        )
        assert struct_nl != snl_new_date, "__eq__() method is broken! false positive"

        # or try a different structure, those should not be equal
        snl_diff_struct = StructureNL(
            self.s2,
            self.hulk,
            ["test_project"],
            self.pmg,
            ["remark1"],
            {"_my_data": self.s2},
            [self.valid_node, self.valid_node2],
            created_at,
        )
        assert struct_nl != snl_diff_struct, "__eq__() method is broken! false positive"

    def test_as_from_dict(self):
        # no complicated objects in the 'data' or 'nodes' field
        struct_nl = StructureNL(
            self.struct,
            self.hulk,
            ["test_project"],
            self.pmg,
            ["remark1"],
            {"_my_data": "string"},
            [self.valid_node, self.valid_node2],
        )
        round_trip_from_dict = StructureNL.from_dict(struct_nl.as_dict())
        assert struct_nl == round_trip_from_dict
        # complicated objects in the 'data' and 'nodes' field
        complicated_node = {
            "name": "complicated node",
            "url": "www.complicatednodegoeshere.com",
            "description": {"structure": self.s2},
        }
        struct_nl = StructureNL(
            self.struct,
            self.hulk,
            ["test_project"],
            self.pmg,
            ["remark1"],
            {"_my_data": {"structure": self.s2}},
            [complicated_node, self.valid_node],
        )
        round_trip_from_dict = StructureNL.from_dict(struct_nl.as_dict())
        assert (
            struct_nl == round_trip_from_dict
        ), "to/from dict is broken when object embedding is used! Apparently MontyEncoding is broken..."

        # Test molecule
        mol_nl = StructureNL(self.mol, self.hulk, references=self.pmg)
        round_trip_from_dict = StructureNL.from_dict(mol_nl.as_dict())
        assert mol_nl == round_trip_from_dict

    def test_from_structures(self):
        s1 = Structure(np.eye(3) * 5, ["Fe"], [[0, 0, 0]])
        s2 = Structure(np.eye(3) * 5, ["Mn"], [[0, 0, 0]])
        remarks = ["unittest"]
        authors = "Test User <test@materialsproject.com>"
        snl_list = StructureNL.from_structures([s1, s2], authors, remarks=remarks)

        assert len(snl_list) == 2
        snl1 = snl_list[0]
        snl2 = snl_list[1]
        assert snl1.remarks == remarks
        assert snl2.remarks == remarks
        assert snl1.authors == [Author.parse_author(authors)]
        assert snl2.authors == [Author.parse_author(authors)]
