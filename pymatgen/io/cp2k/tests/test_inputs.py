# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import unittest
from pathlib import Path
import numpy as np

from pymatgen.core.structure import Molecule, Structure
from pymatgen.io.cp2k.inputs import Coord, Cp2kInput, Keyword, KeywordList, Kind
from pymatgen.util.testing import PymatgenTest

Si_structure = Structure(
    lattice=[[0, 2.734364, 2.734364], [2.734364, 0, 2.734364], [2.734364, 2.734364, 0]],
    species=["Si", "Si"],
    coords=[[0, 0, 0], [0.25, 0.25, 0.25]],
)

nonsense_Structure = Structure(
    lattice=[[-1.0, -10.0, -100.0], [0.1, 0.01, 0.001], [7.0, 11.0, 21.0]],
    species=["X"],
    coords=[[-1, -1, -1]],
)

molecule = Molecule(species=["C", "H"], coords=[[0, 0, 0], [1, 1, 1]])


class InputTest(PymatgenTest):
    def setUp(self):
        self.TEST_FILES_DIR = Path.joinpath(self.TEST_FILES_DIR, "cp2k")
        self.ci = Cp2kInput.from_file(Path.joinpath(self.TEST_FILES_DIR, "cp2k.inp"))

    def test_basic_sections(self):
        s = """
        &GLOBAL
            RUN_TYPE ENERGY
            PROJECT_NAME CP2K ! default name
        &END
        """
        ci = Cp2kInput.from_string(s)
        self.assertEqual(ci["GLOBAL"]["RUN_TYPE"], Keyword("RUN_TYPE", "energy"))
        self.assertEqual(ci["GLOBAL"]["PROJECT_NAME"].description, "default name")
        self.assertMSONable(ci)

    def test_basic_keywords(self):
        kwd = Keyword("TEST1", 1, 2)
        self.assertEqual(kwd.values, (1, 2))
        kwd = Keyword("TEST2", [1, 2, 3])
        self.assertEqual(kwd.values, ([1, 2, 3],))
        kwd = Keyword("TEST3", "xyz", description="testing", units="Ha")
        self.assertEqual(kwd.description, "testing")
        self.assertIn("[Ha]", kwd.get_string())

    def test_coords(self):
        for strucs in [nonsense_Structure, Si_structure, molecule]:
            coords = Coord(strucs)
            self.assertEqual(len(strucs.symbol_set), len(coords.keywords))
            for c in coords.keywords.values():
                self.assertIsInstance(c, KeywordList)

    def test_kind(self):
        for s in [nonsense_Structure, Si_structure, molecule]:
            for spec in s.species:
                self.assertEqual(spec, Kind(spec).specie)

    def test_ci_file(self):
        # proper keyword retrieval
        self.assertEqual(
            self.ci["FORCE_EVAL"]["DFT"]["SCF"]["OT"]["MINIMIZER"],
            Keyword("MINIMIZER", "DIIS"),
        )

        # proper type retrieval
        self.assertIsInstance(self.ci["FORCE_EVAL"]["DFT"]["MGRID"]["NGRIDS"].values[0], int)
        # self.assertIsInstance(self.ci["FORCE_EVAL"]["SUBSYS"]["COORD"]["Si"], Sequence)
        self.assertIsInstance(self.ci["FORCE_EVAL"]["DFT"]["UKS"].values[0], bool)
        self.assertIsInstance(self.ci["FORCE_EVAL"]["DFT"]["QS"]["EPS_DEFAULT"].values[0], float)

        # description retrieval
        self.assertEqual(
            self.ci["FORCE_EVAL"]["SUBSYS"]["CELL"].description,
            "Input parameters needed to set up the CELL.",
        )
        self.assertEqual(
            self.ci["FORCE_EVAL"]["DFT"]["MGRID"]["CUTOFF"].description,
            "Cutoff in [Ry] for finest level of the MG.",
        )

    def test_odd_file(self):
        scramble = ""
        for s in self.ci.get_string():
            if np.random.rand(1) > 0.5:
                if s == "\t":
                    scramble += " "
                elif s == " ":
                    scramble += "  "
                elif s == "&":
                    scramble += s
                elif s == "\n":
                    scramble += s
                elif s.isalpha():
                    scramble += s.lower()
                else:
                    scramble += s
            else:
                scramble += s
        # Can you initialize from jumbled input
        # should be case insensitive and ignore
        # excessive white space or tabs
        ci = Cp2kInput.from_string(scramble)
        self.assertEqual(ci["FORCE_EVAL"]["DFT"]["UKS"], Keyword("UKS", True))
        self.assertEqual(
            [k.name.upper() for k in ci["FORCE_EVAL"]["DFT"]["BASIS_SET_FILE_NAME"]],
            ["BASIS_SET_FILE_NAME", "BASIS_SET_FILE_NAME"],
        )

    def test_preprocessor(self):
        self.assertTrue(self.ci.check("INCLUDE"))
        self.assertEqual(self.ci["INCLUDE"]["KEYWORD"], Keyword("KEYWORD", "VALUE"))
        self.assertEqual(self.ci["FORCE_EVAL"]["METHOD"], Keyword("METHOD", "QS"))
        self.assertEqual(self.ci["FORCE_EVAL"]["DFT"]["SCF"]["MAX_SCF"], Keyword("MAX_SCF", 20))

    def test_mongo(self):
        s = """
        &GLOBAL
            RUN_TYPE ENERGY
            PROJECT_NAME CP2K ! default name
        &END
        """
        s = Cp2kInput.from_string(s)
        s.inc({"GLOBAL": {"TEST": 1}})
        self.assertEqual(s["global"]["test"], Keyword("TEST", 1))

        s.unset({"GLOBAL": "RUN_TYPE"})
        self.assertFalse("RUN_TYPE" in s["global"].keywords)

        s.set({"GLOBAL": {"SUBSEC": {"TEST2": 2}, "SUBSEC2": {"Test2": 1}}})
        self.assertTrue(s.check("global/SUBSEC"))
        self.assertTrue(s.check("global/subsec2"))


if __name__ == "__main__":
    unittest.main()
