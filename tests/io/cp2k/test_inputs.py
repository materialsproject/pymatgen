from __future__ import annotations

import numpy as np
import pytest
from numpy.testing import assert_allclose
from pytest import approx

from pymatgen.core.structure import Molecule, Structure
from pymatgen.io.cp2k.inputs import (
    BasisFile,
    BasisInfo,
    Coord,
    Cp2kInput,
    DataFile,
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
from pymatgen.util.testing import TEST_FILES_DIR, MatSciTest

TEST_DIR = f"{TEST_FILES_DIR}/io/cp2k"


BASIS_FILE_STR = """
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


class TestBasis(MatSciTest):
    def test_basis_info(self):
        # Ensure basis metadata can be read from string
        basis_info = BasisInfo.from_str("cc-pc-DZVP-MOLOPT-q1-SCAN")
        assert basis_info.valence == 2
        assert basis_info.molopt
        assert basis_info.electrons == 1
        assert basis_info.polarization == 1
        assert basis_info.cc
        assert basis_info.pc
        assert basis_info.xc == "SCAN"

        # Ensure one-way soft-matching works
        basis_info2 = BasisInfo.from_str("cc-pc-DZVP-MOLOPT-q1")
        assert basis_info2.softmatch(basis_info)
        assert not basis_info.softmatch(basis_info2)

        basis_info3 = BasisInfo.from_str("cpFIT3")
        assert basis_info3.valence == 3
        assert basis_info3.polarization == 1
        assert basis_info3.contracted

    def test_basis(self):
        # Ensure cp2k formatted string can be read for data correctly
        mol_opt = GaussianTypeOrbitalBasisSet.from_str(BASIS_FILE_STR)
        assert mol_opt.nexp == [7]
        # Basis file can read from strings
        bf = BasisFile.from_str(BASIS_FILE_STR)
        for obj in (mol_opt, bf.objects[0]):
            assert_allclose(
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
        kw: Keyword = mol_opt.get_keyword()
        assert kw.values[0] == "SZV-MOLOPT-GTH"
        mol_opt.info.admm = True
        kw = mol_opt.get_keyword()
        assert list(kw.values) == ["AUX_FIT", "SZV-MOLOPT-GTH"]
        mol_opt.info.admm = False


class TestPotential(MatSciTest):
    all_hydrogen_str = """
H ALLELECTRON ALL
    1    0    0
    0.20000000    0
"""

    pot_hydrogen_str = """
H GTH-PBE-q1 GTH-PBE
    1
    0.20000000    2    -4.17890044     0.72446331
    0
"""

    def test_potential_info(self):
        # Ensure potential metadata can be read from string
        pot_info = PotentialInfo.from_str("GTH-PBE-q1-NLCC")
        assert pot_info.potential_type == "GTH"
        assert pot_info.xc == "PBE"
        assert pot_info.nlcc

        # Ensure one-way soft-matching works
        pot_info2 = PotentialInfo.from_str("GTH-q1-NLCC")
        assert pot_info2.softmatch(pot_info)
        assert not pot_info.softmatch(pot_info2)

    def test_potentials(self):
        # Ensure cp2k formatted string can be read for data correctly
        h_all_elec = GthPotential.from_str(self.all_hydrogen_str)
        assert h_all_elec.potential == "All Electron"
        pot = GthPotential.from_str(self.pot_hydrogen_str)
        assert pot.potential == "Pseudopotential"
        assert pot.r_loc == approx(0.2)
        assert pot.nexp_ppl == approx(2)
        assert_allclose(pot.c_exp_ppl, [-4.17890044, 0.72446331])

        # Basis file can read from strings
        pot_file = PotentialFile.from_str(self.pot_hydrogen_str)
        assert pot_file.objects[0] == pot

        pot_file_path = self.tmp_path / "potential-file"
        pot_file_path.write_text(self.pot_hydrogen_str)
        pot_from_file = PotentialFile.from_file(pot_file_path)
        assert pot_file != pot_from_file  # unequal because pot_from_file has filename != None

        # Ensure keyword can be properly generated
        kw = pot.get_keyword()
        assert kw.values[0] == "GTH-PBE-q1"
        kw = h_all_elec.get_keyword()
        assert kw.values[0] == "ALL"


class TestCp2kInput(MatSciTest):
    si_struct = Structure(
        lattice=[
            [0, 2.734364, 2.734364],
            [2.734364, 0, 2.734364],
            [2.734364, 2.734364, 0],
        ],
        species=["Si", "Si"],
        coords=[[0, 0, 0], [0.25, 0.25, 0.25]],
    )

    invalid_struct = Structure(
        lattice=[[-1.0, -10.0, -100.0], [0.1, 0.01, 0.001], [7.0, 11.0, 21.0]],
        species=["H"],
        coords=[[-1, -1, -1]],
    )

    ch_mol = Molecule(species=["C", "H"], coords=[[0, 0, 0], [1, 1, 1]])

    cp2k_input_str = """
&GLOBAL
    RUN_TYPE ENERGY
    PROJECT_NAME CP2K ! default name
&END
"""

    def setup_method(self):
        self.ci = Cp2kInput.from_file(f"{TEST_DIR}/cp2k.inp")

    def test_basic_sections(self):
        cp2k_input = Cp2kInput.from_str(self.cp2k_input_str)
        assert cp2k_input["GLOBAL"]["RUN_TYPE"] == Keyword("RUN_TYPE", "energy")
        assert cp2k_input["GLOBAL"]["PROJECT_NAME"].description == "default name"
        self.assert_msonable(cp2k_input)

    def test_section_list(self):
        sec1 = Section("TEST")
        sec_list = SectionList(sections=[sec1, sec1])
        for s in sec_list:
            assert isinstance(s, Section)
        assert sec_list[0].name == "TEST"
        assert sec_list[1].name == "TEST"
        assert len(sec_list) == 2
        sec_list += sec1
        assert len(sec_list) == 3

    def test_basic_keywords(self):
        kwd = Keyword("TEST1", 1, 2)
        assert kwd.values == (1, 2)
        kwd = Keyword("TEST2", [1, 2, 3])
        assert kwd.values == ([1, 2, 3],)
        kwd = Keyword("TEST3", "xyz", description="testing", units="Ha")
        assert kwd.description == "testing"
        assert "[Ha]" in kwd.get_str()

    def test_coords(self):
        for struct in (self.invalid_struct, self.si_struct, self.ch_mol):
            coords = Coord(struct)
            for val in coords.keywords.values():
                assert isinstance(val, Keyword | KeywordList)

    def test_kind(self):
        for struct in (self.invalid_struct, self.si_struct, self.ch_mol):
            for spec in struct.species:
                assert spec == Kind(spec).specie

    def test_ci_file(self):
        # proper type retrieval
        assert isinstance(self.ci["FORCE_EVAL"]["DFT"]["MGRID"]["NGRIDS"].values[0], int)
        assert isinstance(self.ci["FORCE_EVAL"]["DFT"]["UKS"].values[0], bool)
        assert isinstance(self.ci["FORCE_EVAL"]["DFT"]["QS"]["EPS_DEFAULT"].values[0], float)

        # description retrieval
        assert (
            self.ci["FORCE_EVAL"]["SUBSYS"].description == '"/" used to break input file parser in Section start line'
        )
        assert self.ci["FORCE_EVAL"]["SUBSYS"]["CELL"].description == "Input parameters needed to set up the CELL."
        assert (
            self.ci["FORCE_EVAL"]["DFT"]["MGRID"]["CUTOFF"].description == "Cutoff in [Ry] for finest level of the MG."
        )
        assert self.ci["FORCE_EVAL"]["METHOD"].description is None

    def test_odd_file(self):
        scramble = ""
        rng = np.random.default_rng()
        for string in self.ci.get_str():
            if rng.choice((True, False)):
                if string == "\t":
                    scramble += " "
                elif string == " ":
                    scramble += "  "
                elif string in ("&", "\n"):
                    scramble += string
                elif string.isalpha():
                    scramble += string.lower()
                else:
                    scramble += string
            else:
                scramble += string
        # Can you initialize from jumbled input
        # should be case insensitive and ignore
        # excessive white space or tabs
        ci = Cp2kInput.from_str(scramble)
        assert ci["FORCE_EVAL"]["DFT"]["UKS"] == Keyword("UKS", True)  # noqa: FBT003
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
        cp2k_input = Cp2kInput.from_str(self.cp2k_input_str)
        cp2k_input.inc({"GLOBAL": {"TEST": 1}})
        assert cp2k_input["global"]["test"] == Keyword("TEST", 1)

        cp2k_input.unset({"GLOBAL": "RUN_TYPE"})
        assert "RUN_TYPE" not in cp2k_input["global"].keywords

        cp2k_input.set({"GLOBAL": {"SUBSEC": {"TEST2": 2}, "SUBSEC2": {"Test2": 1}}})
        assert cp2k_input.check("global/SUBSEC")
        assert cp2k_input.check("global/subsec2")


class TestDataFile(MatSciTest):
    def test_data_file(self):
        # make temp file with BASIS_FILE_STR
        data_file = self.tmp_path / "data-file"
        data_file.write_text(BASIS_FILE_STR)
        with pytest.raises(NotImplementedError):
            DataFile.from_file(data_file)
