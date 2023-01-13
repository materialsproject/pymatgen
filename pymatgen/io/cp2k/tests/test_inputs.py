# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import annotations

import unittest
from pathlib import Path

import numpy as np

from pymatgen.core.structure import Molecule, Structure
from pymatgen.io.cp2k.inputs import (
    BasisFile,
    BasisInfo,
    Coord,
    Cp2kInput,
    GaussianTypeOrbitalBasisSet,
    GthPotential,
    Keyword,
    KeywordList,
    Kind,
    PotentialFile,
    PotentialInfo,
    Section,
    SectionList,
)
from pymatgen.util.testing import PymatgenTest

Si_structure = Structure(
    lattice=[[0, 2.734364, 2.734364], [2.734364, 0, 2.734364], [2.734364, 2.734364, 0]],
    species=["Si", "Si"],
    coords=[[0, 0, 0], [0.25, 0.25, 0.25]],
)

nonsense_Structure = Structure(
    lattice=[[-1.0, -10.0, -100.0], [0.1, 0.01, 0.001], [7.0, 11.0, 21.0]],
    species=["H"],
    coords=[[-1, -1, -1]],
)

molecule = Molecule(species=["C", "H"], coords=[[0, 0, 0], [1, 1, 1]])

basis = """
 H  SZV-MOLOPT-GTH SZV-MOLOPT-GTH-q1
 1
 2 0 0 7 1
     11.478000339908  0.024916243200
      3.700758562763  0.079825490000
      1.446884268432  0.128862675300
      0.716814589696  0.379448894600
      0.247918564176  0.324552432600
      0.066918004004  0.037148121400
      0.021708243634 -0.001125195500
"""
all_H = """
H ALLELECTRON ALL
    1    0    0
     0.20000000    0
"""
pot_H = """
H GTH-PBE-q1 GTH-PBE
    1
     0.20000000    2    -4.17890044     0.72446331
    0
"""


class BasisAndPotentialTest(PymatgenTest):
    def test_basis_info(self):
        # Ensure basis metadata can be read from string
        b = BasisInfo.from_string("cc-pc-DZVP-MOLOPT-q1-SCAN")
        self.assertEqual(b.valence, 2)
        self.assertEqual(b.molopt, True)
        self.assertEqual(b.electrons, 1)
        self.assertEqual(b.polarization, 1)
        self.assertEqual(b.cc, True)
        self.assertEqual(b.pc, True)
        self.assertEqual(b.xc, "SCAN")

        # Ensure one-way softmatching works
        b2 = BasisInfo.from_string("cc-pc-DZVP-MOLOPT-q1")
        self.assertTrue(b2.softmatch(b))
        self.assertFalse(b.softmatch(b2))

        b3 = BasisInfo.from_string("cpFIT3")
        self.assertEqual(b3.valence, 3)
        self.assertEqual(b3.polarization, 1)
        self.assertTrue(b3.contracted, True)

    def test_potential_info(self):
        # Ensure potential metadata can be read from string
        p = PotentialInfo.from_string("GTH-PBE-q1-NLCC")
        self.assertEqual(p.potential_type, "GTH")
        self.assertEqual(p.xc, "PBE")
        self.assertEqual(p.nlcc, True)

        # Ensure one-way softmatching works
        p2 = PotentialInfo.from_string("GTH-q1-NLCC")
        self.assertTrue(p2.softmatch(p))
        self.assertFalse(p.softmatch(p2))

    def test_basis(self):
        # Ensure cp2k formatted string can be read for data correctly
        molopt = GaussianTypeOrbitalBasisSet.from_string(basis)
        self.assertEqual(molopt.nexp, [7])
        # Basis file can read from strings
        bf = BasisFile.from_string(basis)
        for obj in [molopt, bf.objects[0]]:
            self.assertArrayAlmostEqual(
                obj.exponents[0],
                [
                    11.478000339908,
                    3.700758562763,
                    1.446884268432,
                    0.716814589696,
                    0.247918564176,
                    0.066918004004,
                    0.021708243634,
                ],
            )

        # Ensure keyword can be properly generated
        kw = molopt.get_keyword()
        self.assertEqual(kw.values[0], "SZV-MOLOPT-GTH")
        molopt.info.admm = True
        kw = molopt.get_keyword()
        self.assertArrayEqual(kw.values, ["AUX_FIT", "SZV-MOLOPT-GTH"])
        molopt.info.admm = False

    def test_potentials(self):
        # Ensure cp2k formatted string can be read for data correctly
        alle = GthPotential.from_string(all_H)
        self.assertEqual(alle.potential, "All Electron")
        pot = GthPotential.from_string(pot_H)
        self.assertEqual(pot.potential, "Pseudopotential")
        self.assertAlmostEqual(pot.r_loc, 0.2)
        self.assertAlmostEqual(pot.nexp_ppl, 2)
        self.assertArrayAlmostEqual(pot.c_exp_ppl, [-4.17890044, 0.72446331])

        # Basis file can read from strings
        pf = PotentialFile.from_string(pot_H)
        self.assertEqual(pf.objects[0], pot)

        # Ensure keyword can be properly generated
        kw = pot.get_keyword()
        self.assertEqual(kw.values[0], "GTH-PBE-q1")
        kw = alle.get_keyword()
        self.assertEqual(kw.values[0], "ALL")


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

    def test_sectionlist(self):
        s1 = Section("TEST")
        sl = SectionList(sections=[s1, s1])
        for s in sl:
            assert isinstance(s, Section)
        self.assertEqual(sl[0].name, "TEST")
        self.assertEqual(sl[1].name, "TEST")
        self.assertEqual(len(sl), 2)
        sl += s1
        self.assertEqual(len(sl), 3)

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
            for c in coords.keywords.values():
                self.assertTrue(isinstance(c, (Keyword, KeywordList)))

    def test_kind(self):
        for s in [nonsense_Structure, Si_structure, molecule]:
            for spec in s.species:
                self.assertEqual(spec, Kind(spec).specie)

    def test_ci_file(self):
        # proper type retrieval
        self.assertIsInstance(self.ci["FORCE_EVAL"]["DFT"]["MGRID"]["NGRIDS"].values[0], int)
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
        self.assertEqual(self.ci["FORCE_EVAL"]["DFT"]["SCF"]["MAX_SCF"], Keyword("MAX_SCF", 1))

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
