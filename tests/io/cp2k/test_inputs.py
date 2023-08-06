from __future__ import annotations

import numpy as np
from numpy.testing import assert_array_equal
from pytest import approx

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
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

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


class TestBasisAndPotential(PymatgenTest):
    def test_basis_info(self):
        # Ensure basis metadata can be read from string
        b = BasisInfo.from_str("cc-pc-DZVP-MOLOPT-q1-SCAN")
        assert b.valence == 2
        assert b.molopt
        assert b.electrons == 1
        assert b.polarization == 1
        assert b.cc
        assert b.pc
        assert b.xc == "SCAN"

        # Ensure one-way softmatching works
        b2 = BasisInfo.from_str("cc-pc-DZVP-MOLOPT-q1")
        assert b2.softmatch(b)
        assert not b.softmatch(b2)

        b3 = BasisInfo.from_str("cpFIT3")
        assert b3.valence == 3
        assert b3.polarization == 1
        assert b3.contracted, True

    def test_potential_info(self):
        # Ensure potential metadata can be read from string
        p = PotentialInfo.from_str("GTH-PBE-q1-NLCC")
        assert p.potential_type == "GTH"
        assert p.xc == "PBE"
        assert p.nlcc

        # Ensure one-way softmatching works
        p2 = PotentialInfo.from_str("GTH-q1-NLCC")
        assert p2.softmatch(p)
        assert not p.softmatch(p2)

    def test_basis(self):
        # Ensure cp2k formatted string can be read for data correctly
        mol_opt = GaussianTypeOrbitalBasisSet.from_str(basis)
        assert mol_opt.nexp == [7]
        # Basis file can read from strings
        bf = BasisFile.from_str(basis)
        for obj in [mol_opt, bf.objects[0]]:
            assert np.allclose(
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
        kw = mol_opt.get_keyword()
        assert kw.values[0] == "SZV-MOLOPT-GTH"  # noqa: PD011
        mol_opt.info.admm = True
        kw = mol_opt.get_keyword()
        assert_array_equal(kw.values, ["AUX_FIT", "SZV-MOLOPT-GTH"])
        mol_opt.info.admm = False

    def test_potentials(self):
        # Ensure cp2k formatted string can be read for data correctly
        all = GthPotential.from_str(all_H)
        assert all.potential == "All Electron"
        pot = GthPotential.from_str(pot_H)
        assert pot.potential == "Pseudopotential"
        assert pot.r_loc == approx(0.2)
        assert pot.nexp_ppl == approx(2)
        assert np.allclose(pot.c_exp_ppl, [-4.17890044, 0.72446331])

        # Basis file can read from strings
        pf = PotentialFile.from_str(pot_H)
        assert pf.objects[0] == pot

        # Ensure keyword can be properly generated
        kw = pot.get_keyword()
        assert kw.values[0] == "GTH-PBE-q1"  # noqa: PD011
        kw = all.get_keyword()
        assert kw.values[0] == "ALL"  # noqa: PD011


class TestInput(PymatgenTest):
    def setUp(self):
        self.ci = Cp2kInput.from_file(f"{TEST_FILES_DIR}/cp2k/cp2k.inp")

    def test_basic_sections(self):
        s = """
        &GLOBAL
            RUN_TYPE ENERGY
            PROJECT_NAME CP2K ! default name
        &END
        """
        ci = Cp2kInput.from_str(s)
        assert ci["GLOBAL"]["RUN_TYPE"] == Keyword("RUN_TYPE", "energy")
        assert ci["GLOBAL"]["PROJECT_NAME"].description == "default name"
        self.assert_msonable(ci)

    def test_sectionlist(self):
        s1 = Section("TEST")
        sl = SectionList(sections=[s1, s1])
        for s in sl:
            assert isinstance(s, Section)
        assert sl[0].name == "TEST"
        assert sl[1].name == "TEST"
        assert len(sl) == 2
        sl += s1
        assert len(sl) == 3

    def test_basic_keywords(self):
        kwd = Keyword("TEST1", 1, 2)
        assert kwd.values == (1, 2)  # noqa: PD011
        kwd = Keyword("TEST2", [1, 2, 3])
        assert kwd.values == ([1, 2, 3],)  # noqa: PD011
        kwd = Keyword("TEST3", "xyz", description="testing", units="Ha")
        assert kwd.description == "testing"
        assert "[Ha]" in kwd.get_string()

    def test_coords(self):
        for strucs in [nonsense_Structure, Si_structure, molecule]:
            coords = Coord(strucs)
            for c in coords.keywords.values():
                assert isinstance(c, (Keyword, KeywordList))

    def test_kind(self):
        for s in [nonsense_Structure, Si_structure, molecule]:
            for spec in s.species:
                assert spec == Kind(spec).specie

    def test_ci_file(self):
        # proper type retrieval
        assert isinstance(self.ci["FORCE_EVAL"]["DFT"]["MGRID"]["NGRIDS"].values[0], int)  # noqa: PD011
        assert isinstance(self.ci["FORCE_EVAL"]["DFT"]["UKS"].values[0], bool)  # noqa: PD011
        assert isinstance(self.ci["FORCE_EVAL"]["DFT"]["QS"]["EPS_DEFAULT"].values[0], float)  # noqa: PD011

        # description retrieval
        assert self.ci["FORCE_EVAL"]["SUBSYS"]["CELL"].description == "Input parameters needed to set up the CELL."
        assert (
            self.ci["FORCE_EVAL"]["DFT"]["MGRID"]["CUTOFF"].description == "Cutoff in [Ry] for finest level of the MG."
        )

    def test_odd_file(self):
        scramble = ""
        for s in self.ci.get_string():
            if np.random.rand(1) > 0.5:
                if s == "\t":
                    scramble += " "
                elif s == " ":
                    scramble += "  "
                elif s in ("&", "\n"):
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
        ci = Cp2kInput.from_str(scramble)
        assert ci["FORCE_EVAL"]["DFT"]["UKS"] == Keyword("UKS", True)
        assert [k.name.upper() for k in ci["FORCE_EVAL"]["DFT"]["BASIS_SET_FILE_NAME"]] == [
            "BASIS_SET_FILE_NAME",
            "BASIS_SET_FILE_NAME",
        ]

    def test_preprocessor(self):
        assert self.ci.check("INCLUDE")
        assert self.ci["INCLUDE"]["KEYWORD"] == Keyword("KEYWORD", "VALUE")
        assert self.ci["FORCE_EVAL"]["METHOD"] == Keyword("METHOD", "QS")
        assert self.ci["FORCE_EVAL"]["DFT"]["SCF"]["MAX_SCF"] == Keyword("MAX_SCF", 1)

    def test_mongo(self):
        s = """
        &GLOBAL
            RUN_TYPE ENERGY
            PROJECT_NAME CP2K ! default name
        &END
        """
        s = Cp2kInput.from_str(s)
        s.inc({"GLOBAL": {"TEST": 1}})
        assert s["global"]["test"] == Keyword("TEST", 1)

        s.unset({"GLOBAL": "RUN_TYPE"})
        assert "RUN_TYPE" not in s["global"].keywords

        s.set({"GLOBAL": {"SUBSEC": {"TEST2": 2}, "SUBSEC2": {"Test2": 1}}})
        assert s.check("global/SUBSEC")
        assert s.check("global/subsec2")
