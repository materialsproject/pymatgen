# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals
import datetime
import unittest
import numpy as np

from pymatgen import Structure, Molecule
from pymatgen.util.provenance import StructureNL, HistoryNode, Author

"""
Unit tests for StructureNL (SNL) format
"""

__author__ = "Anubhav Jain"
__credits__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Anubhav Jain"
__email__ = "ajain@lbl.gov"
__date__ = "2/14/13"


class StructureNLCase(unittest.TestCase):

    def setUp(self):
        # set up a Structure
        self.s = Structure(np.eye(3, 3) * 3, ["Fe"], [[0, 0, 0]])
        self.s2 = Structure(np.eye(3, 3) * 3, ["Al"], [[0, 0, 0]])
        self.mol = Molecule(["He"], [[0, 0, 0]])
        # set up BibTeX strings
        self.matproj = "@misc{MaterialsProject,\ntitle = {{Materials " \
                       "Project}},\nurl = {http://www.materialsproject.org}\n}"
        self.pmg = "@article{Ong2013,\n author = {Ong, " \
                   "Shyue Ping and Richards, William Davidson and Jain, " \
                   "Anubhav and Hautier, Geoffroy and Kocher, " \
                   "Michael and Cholia, Shreyas and Gunter, Dan and Chevrier," \
                   " Vincent L. and Persson, Kristin A. and Ceder, Gerbrand}," \
                   "\n doi = {10.1016/j.commatsci.2012.10.028}," \
                   "\n issn = {09270256},\n journal = {Computational " \
                   "Materials Science},\n month = feb,\n pages = {314--319}," \
                   "\n publisher = {Elsevier B.V.}," \
                   "\n title = {{Python Materials Genomics (pymatgen): A " \
                   "robust, open-source python library for materials " \
                   "analysis}},\n url = {http://linkinghub.elsevier" \
                   ".com/retrieve/pii/S0927025612006295},\n volume = {68}," \
                   "\n year = {2013}\n}"
        repeat = "REPEAT" * 10000
        self.superlong = "@misc{SuperLong,\ntitle = {{" + repeat + "}}}"
        self.unicode_title = "@misc{Unicode_Title,\ntitle = {{A \u73ab is a rose}}}"
        self.junk = "This is junk text, not a BibTeX reference"

        # set up remarks
        self.remark_fail = ["This is a really long remark that is clearly invalid and must fail, don't you agree? It would be silly to allow remarks that went on forever and ever."]

        # set up some authors
        self.hulk = [{"name": "Hulk", "email": "hulk@avengers.com"}]
        self.america = "Captain America <captainamerica@avengers.com>"
        self.thor = [("Thor", "thor@avengers.com")]
        self.duo = "Iron Man <ironman@avengers.com>, " \
                   "Black Widow <blackwidow@avengers.com>"

        # set up HistoryNodes
        self.valid_node = HistoryNode("DB 1", "www.db1URLgoeshere.com",
                                      {"db1_id": 12424})
        self.valid_node2 = {"name": "DB 2", "url": "www.db2URLgoeshere.com",
                            "description": {"db2_id": 12424}}
        self.invalid_node = {"name": "DB 3",
                             "url": "http://www.db3isnotavalidnode.com"}

    def test_authors(self):
        a = StructureNL(self.s, self.hulk, references=self.pmg)
        self.assertEqual(a.authors[0].name, "Hulk")
        self.assertEqual(a.authors[0].email, "hulk@avengers.com")

        a = StructureNL(self.s, self.america, references=self.pmg)
        self.assertEqual(a.authors[0].name, "Captain America")
        self.assertEqual(a.authors[0].email, "captainamerica@avengers.com")

        a = StructureNL(self.s, self.thor, references=self.pmg)
        self.assertEqual(a.authors[0].name, "Thor")
        self.assertEqual(a.authors[0].email, "thor@avengers.com")

        a = StructureNL(self.s, self.duo, references=self.pmg)
        self.assertEqual(a.authors[0].name, "Iron Man")
        self.assertEqual(a.authors[0].email, "ironman@avengers.com")
        self.assertEqual(a.authors[1].name, "Black Widow")
        self.assertEqual(a.authors[1].email, "blackwidow@avengers.com")
        StructureNL(self.s, self.hulk, references=self.pmg)

    def test_references(self):
        # An empty string should be OK
        StructureNL(self.s, self.hulk, references="")

        # An empty list should not work
        self.assertRaises(ValueError, StructureNL, self.s, self.hulk,
                          references=[])

        # junk reference should not work
        self.assertRaises(ValueError, StructureNL, self.s, self.hulk,
                          references=self.junk)

        # good references should be ok
        StructureNL(self.s, self.hulk, references=self.pmg)

        # unicode references should work
        StructureNL(self.s, self.hulk, references=self.unicode_title)

        # multi-line references should be OK
        StructureNL(self.s, self.hulk,
                    references='\n'.join([self.matproj, self.pmg]))

        # super long references are bad
        self.assertRaises(ValueError, StructureNL, self.s, self.hulk,
                          references=self.superlong)

    def test_historynodes(self):
        a = StructureNL(self.s, self.hulk, history=[self.valid_node])
        self.assertEqual(a.history[0].name, "DB 1")
        self.assertEqual(a.history[0].url, "www.db1URLgoeshere.com")
        self.assertEqual(a.history[0].description, {"db1_id": 12424})

        a = StructureNL(self.s, self.hulk,
                        history=[self.valid_node, self.valid_node2])
        self.assertEqual(a.history[1].name, "DB 2")
        self.assertEqual(a.history[1].url, "www.db2URLgoeshere.com")
        self.assertEqual(a.history[1].description, {"db2_id": 12424})

        # invalid nodes should not work
        self.assertRaises(Exception, StructureNL, self.s, self.hulk,
                          history=[self.invalid_node])

        # too many nodes should not work
        self.assertRaises(ValueError, StructureNL, self.s, self.hulk,
                          history=[self.valid_node] * 1000)

    def test_data(self):
        # Structure data is OK due to PMGEncoder/Decoder
        a = StructureNL(self.s, self.hulk, data={"_structure": self.s2})
        self.assertEqual(a.data["_structure"], self.s2,
                         'Data storage is broken')
        self.assertRaises(ValueError, StructureNL, self.s, self.hulk,
                          data={"bad_key": 1})

    def test_remarks(self):
        a = StructureNL(self.s, self.hulk, remarks="string format")
        self.assertEqual(a.remarks[0], "string format")
        self.assertRaises(ValueError, StructureNL, self.s, self.hulk,
                          remarks=self.remark_fail)

    def test_eq(self):
        # test basic Equal()
        created_at = datetime.datetime.now()
        a = StructureNL(self.s, self.hulk, ['test_project'], self.pmg,
                        ['remark1'], {"_my_data": self.s2},
                        [self.valid_node, self.valid_node2], created_at)
        b = StructureNL(self.s, self.hulk, ['test_project'], self.pmg,
                        ['remark1'], {"_my_data": self.s2},
                        [self.valid_node, self.valid_node2], created_at)
        self.assertEqual(a, b, "__eq__() method is broken! false negative")

        # change the created at date, now they are no longer equal
        created_at = datetime.datetime.now() + datetime.timedelta(days=-1)
        c = StructureNL(self.s, self.hulk, ['test_project'], self.pmg,
                        ['remark1'], {"_my_data": self.s2},
                        [self.valid_node, self.valid_node2], created_at)
        self.assertNotEqual(a, c, "__eq__() method is broken! false positive")

        # or try a different structure, those should not be equal
        d = StructureNL(self.s2, self.hulk, ['test_project'], self.pmg,
                        ['remark1'], {"_my_data": self.s2},
                        [self.valid_node, self.valid_node2], created_at)
        self.assertNotEqual(a, d, "__eq__() method is broken! false positive")

    def test_to_from_dict(self):
        # no complicated objects in the 'data' or 'nodes' field
        a = StructureNL(self.s, self.hulk, ['test_project'], self.pmg,
                        ['remark1'], {"_my_data": "string"},
                        [self.valid_node, self.valid_node2])
        b = StructureNL.from_dict(a.as_dict())
        self.assertEqual(a, b)
        # complicated objects in the 'data' and 'nodes' field
        complicated_node = {"name": "complicated node",
                            "url": "www.complicatednodegoeshere.com",
                            "description": {"structure": self.s2}}
        a = StructureNL(self.s, self.hulk, ['test_project'], self.pmg,
                        ['remark1'], {"_my_data": {"structure": self.s2}},
                        [complicated_node, self.valid_node])
        b = StructureNL.from_dict(a.as_dict())
        self.assertEqual(a, b,
                         'to/from dict is broken when object embedding is '
                         'used! Apparently MontyEncoding is broken...')

        #Test molecule
        molnl = StructureNL(self.mol, self.hulk, references=self.pmg)
        b = StructureNL.from_dict(molnl.as_dict())
        self.assertEqual(molnl, b)

    def test_from_structures(self):
        s1 = Structure([[5, 0, 0], [0, 5, 0], [0, 0, 5]], ["Fe"], [[0, 0, 0]])
        s2 = Structure([[5, 0, 0], [0, 5, 0], [0, 0, 5]], ["Mn"], [[0, 0, 0]])
        remarks = ["unittest"]
        authors="Test User <test@materialsproject.com>"
        snl_list = StructureNL.from_structures([s1, s2], authors, remarks=remarks)

        self.assertEqual(len(snl_list), 2)
        snl1 = snl_list[0]
        snl2 = snl_list[1]
        self.assertEqual(snl1.remarks, remarks)
        self.assertEqual(snl2.remarks, remarks)
        self.assertEqual(snl1.authors, [Author.parse_author(authors)])
        self.assertEqual(snl2.authors, [Author.parse_author(authors)])


if __name__ == '__main__':
    unittest.main()
