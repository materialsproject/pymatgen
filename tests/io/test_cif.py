from __future__ import annotations

import unittest

import numpy as np
import pytest
from pytest import approx

from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Composition, DummySpecies, Element, Lattice, Species, Structure
from pymatgen.electronic_structure.core import Magmom
from pymatgen.io.cif import CifBlock, CifParser, CifWriter
from pymatgen.symmetry.structure import SymmetrizedStructure
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

try:
    import pybtex
except ImportError:
    pybtex = None


class TestCifBlock(PymatgenTest):
    def test_to_string(self):
        with open(f"{TEST_FILES_DIR}/Graphite.cif") as file:
            cif_str = file.read()
        cif_block = CifBlock.from_str(cif_str)
        cif_str_2 = str(CifBlock.from_str(str(cif_block)))
        cif_str = """data_53781-ICSD
_database_code_ICSD   53781
_audit_creation_date   2003-04-01
_audit_update_record   2013-02-01
_chemical_name_systematic   Carbon
_chemical_formula_structural   C
_chemical_formula_sum   C1
_chemical_name_structure_type   Graphite(2H)
_chemical_name_mineral   'Graphite 2H'
_exptl_crystal_density_diffrn   2.22
_publ_section_title   'Structure of graphite'
loop_
 _citation_id
 _citation_journal_full
 _citation_year
 _citation_journal_volume
 _citation_page_first
 _citation_page_last
 _citation_journal_id_ASTM
  primary  'Physical Review (1,1893-132,1963/141,1966-188,1969)'
  1917  10  661  696  PHRVAO
loop_
 _publ_author_name
  'Hull, A.W.'
_cell_length_a   2.47
_cell_length_b   2.47
_cell_length_c   6.8
_cell_angle_alpha   90.
_cell_angle_beta   90.
_cell_angle_gamma   120.
_cell_volume   35.93
_cell_formula_units_Z   4
_symmetry_space_group_name_H-M   'P 63/m m c'
_symmetry_Int_Tables_number   194
loop_
 _symmetry_equiv_pos_site_id
 _symmetry_equiv_pos_as_xyz
  1  'x, x-y, -z+1/2'
  2  '-x+y, y, -z+1/2'
  3  '-y, -x, -z+1/2'
  4  '-x+y, -x, -z+1/2'
  5  '-y, x-y, -z+1/2'
  6  'x, y, -z+1/2'
  7  '-x, -x+y, z+1/2'
  8  'x-y, -y, z+1/2'
  9  'y, x, z+1/2'
  10  'x-y, x, z+1/2'
  11  'y, -x+y, z+1/2'
  12  '-x, -y, z+1/2'
  13  '-x, -x+y, -z'
  14  'x-y, -y, -z'
  15  'y, x, -z'
  16  'x-y, x, -z'
  17  'y, -x+y, -z'
  18  '-x, -y, -z'
  19  'x, x-y, z'
  20  '-x+y, y, z'
  21  '-y, -x, z'
  22  '-x+y, -x, z'
  23  '-y, x-y, z'
  24  'x, y, z'
loop_
 _atom_type_symbol
 _atom_type_oxidation_number
  C0+  0
loop_
 _atom_site_label
 _atom_site_type_symbol
 _atom_site_symmetry_multiplicity
 _atom_site_Wyckoff_symbol
 _atom_site_fract_x
 _atom_site_fract_y
 _atom_site_fract_z
 _atom_site_B_iso_or_equiv
 _atom_site_occupancy
 _atom_site_attached_hydrogens
  C1  C0+  2  b  0  0  0.25  .  1.  0
  C2  C0+  2  c  0.3333  0.6667  0.25  .  1.  0"""
        for l1, l2, l3 in zip(str(cif_block).split("\n"), cif_str.split("\n"), cif_str_2.split("\n")):
            assert l1.strip() == l2.strip()
            assert l2.strip() == l3.strip()

    def test_double_quotes_and_underscore_data(self):
        cif_str = """data_test
_symmetry_space_group_name_H-M   "P -3 m 1"
_thing   '_annoying_data'"""
        cb = CifBlock.from_str(cif_str)
        assert cb["_symmetry_space_group_name_H-M"] == "P -3 m 1"
        assert cb["_thing"] == "_annoying_data"
        assert str(cb) == cif_str.replace('"', "'")

    def test_double_quoted_data(self):
        cif_str = """data_test
_thing   ' '_annoying_data''
_other   " "_more_annoying_data""
_more   ' "even more" ' """
        cb = CifBlock.from_str(cif_str)
        assert cb["_thing"] == " '_annoying_data'"
        assert cb["_other"] == ' "_more_annoying_data"'
        assert cb["_more"] == ' "even more" '

    def test_nested_fake_multiline_quotes(self):
        cif_str = """data_test
_thing
;
long quotes
 ;
 still in the quote
 ;
actually going to end now
;"""
        cb = CifBlock.from_str(cif_str)
        assert cb["_thing"] == " long quotes  ;  still in the quote  ; actually going to end now"

    def test_long_loop(self):
        data = {
            "_stuff1": ["A" * 30] * 2,
            "_stuff2": ["B" * 30] * 2,
            "_stuff3": ["C" * 30] * 2,
        }
        loops = [["_stuff1", "_stuff2", "_stuff3"]]
        cif_str = """data_test
loop_
 _stuff1
 _stuff2
 _stuff3
  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA  BBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA  BBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"""
        assert str(CifBlock(data, loops, "test")) == cif_str


class TestCifIO(PymatgenTest):
    def test_cif_parser(self):
        parser = CifParser(f"{TEST_FILES_DIR}/LiFePO4.cif")
        for struct in parser.parse_structures():
            assert struct.formula == "Li4 Fe4 P4 O16", "Incorrectly parsed cif."

        parser = CifParser(f"{TEST_FILES_DIR}/V2O3.cif")
        for struct in parser.parse_structures():
            assert struct.formula == "V4 O6"

        bibtex_str = """
@article{cifref0,
    author = "Andersson, G.",
    title = "Studies on vanadium oxides. I. Phase analysis",
    journal = "Acta Chemica Scandinavica (1-27,1973-42,1988)",
    volume = "8",
    year = "1954",
    pages = "1599--1606"
}
        """
        assert parser.get_bibtex_string().strip() == bibtex_str.strip()

        parser = CifParser(f"{TEST_FILES_DIR}/Li2O.cif")
        prim = parser.parse_structures()[0]
        assert prim.formula == "Li8 O4"
        conv = parser.parse_structures(primitive=False)[0]
        assert conv.formula == "Li8 O4"

        # test for disordered structures
        parser = CifParser(f"{TEST_FILES_DIR}/Li10GeP2S12.cif")
        for struct in parser.parse_structures():
            assert struct.formula == "Li20.2 Ge2.06 P3.94 S24", "Incorrectly parsed cif."

        with open(f"{TEST_FILES_DIR}/FePO4.cif") as cif_file:
            cif_str = cif_file.read()
        parser = CifParser.from_str(cif_str)
        struct = parser.parse_structures(primitive=False)[0]
        assert struct.formula == "Fe4 P4 O16"
        assert struct.lattice.a == approx(10.4117668699)
        assert struct.lattice.b == approx(6.06717187997)
        assert struct.lattice.c == approx(4.75948953998)
        assert struct.lattice.alpha == approx(91)
        assert struct.lattice.beta == approx(92)
        assert struct.lattice.gamma == approx(93)

        parser = CifParser(f"{TEST_FILES_DIR}/srycoo.cif")
        assert parser.parse_structures()[0].formula == "Sr11.2 Y4.8 Co16 O42"

        # Test with a decimal Xyz. This should parse as two atoms in
        # conventional cell if it is correct, one if not.
        parser = CifParser(f"{TEST_FILES_DIR}/Fe.cif")
        assert len(parser.parse_structures(primitive=False)[0]) == 2
        assert not parser.has_errors

    def test_get_symmetrized_structure(self):
        parser = CifParser(f"{TEST_FILES_DIR}/Li2O.cif")
        sym_structure = parser.parse_structures(primitive=False, symmetrized=True)[0]
        structure = parser.parse_structures(primitive=False, symmetrized=False)[0]
        assert isinstance(sym_structure, SymmetrizedStructure)
        assert structure == sym_structure
        assert sym_structure.equivalent_indices == [[0, 1, 2, 3], [4, 5, 6, 7, 8, 9, 10, 11]]
        assert set(sym_structure.labels) == {"O1", "Li1"}

    def test_site_symbol_preference(self):
        parser = CifParser(f"{TEST_FILES_DIR}/site_type_symbol_test.cif")
        assert parser.parse_structures()[0].formula == "Ge1.6 Sb1.6 Te4"

    def test_implicit_hydrogen(self):
        parser = CifParser(f"{TEST_FILES_DIR}/Senegalite_implicit_hydrogen.cif")
        for struct in parser.parse_structures():
            assert struct.formula == "Al8 P4 O32"
            assert sum(struct.site_properties["implicit_hydrogens"]) == 20
        assert (
            "Structure has implicit hydrogens defined, "
            "parsed structure unlikely to be suitable for use "
            "in calculations unless hydrogens added." in parser.warnings
        )
        parser = CifParser(f"{TEST_FILES_DIR}/cif_implicit_hydrogens_cod_1011130.cif")
        struct = parser.parse_structures()[0]
        assert (
            "Structure has implicit hydrogens defined, "
            "parsed structure unlikely to be suitable for use "
            "in calculations unless hydrogens added." in parser.warnings
        )

    def test_site_labels(self):
        parser = CifParser(f"{TEST_FILES_DIR}/garnet.cif")
        struct = parser.parse_structures(primitive=True)[0]

        # ensure structure has correct number of labels
        assert len(struct.labels) == len(struct)

        # ensure the labels are unique and match the expected site names
        expected_site_names = {"Al1", "Ca1", "O1", "Si1"}
        assert {*struct.labels} == expected_site_names

        # check label of each site
        for site, label in zip(struct, struct.labels):
            assert site.label == label
            # Ensure the site label starts with the site species name
            assert site.label.startswith(site.specie.name)

        # ensure multiple species with different names have correct labels
        parser2 = CifParser(f"{TEST_FILES_DIR}/Fe3O4.cif")
        struct2 = parser2.parse_structures(primitive=False)[0]

        expected_site_names2 = {*"O1 O2 O3 O4 O5 O6 O7 O8 Fe9 Fe10 Fe11 Fe12 Fe13 Fe14".split()}
        assert set(struct2.labels) == expected_site_names2

    def test_cif_writer_labeled(self):
        parser = CifParser(f"{TEST_FILES_DIR}/garnet.cif")
        struct = parser.parse_structures()[0]
        for idx, site in enumerate(struct):
            site.label = f"my_{site.specie.name}{idx}"
        writer = CifWriter(struct)

        parser2 = CifParser.from_str(str(writer))
        struct2 = parser2.parse_structures()[0]

        assert set(struct.labels) == set(struct2.labels)

    def test_cif_parser_springer_pauling(self):
        # Below are 10 tests for CIFs from the Springer Materials/Pauling file DBs.

        # Partial occupancy on sites, incorrect label, previously unparsable
        parser = CifParser(f"{TEST_FILES_DIR}/PF_sd_1928405.cif")
        for struct in parser.parse_structures():
            assert struct.formula == "Er1 Mn3.888 Fe2.112 Sn6"
        assert parser.has_errors

        # Partial occupancy on sites, previously parsed as an ordered structure
        parser = CifParser(f"{TEST_FILES_DIR}/PF_sd_1011081.cif")
        for struct in parser.parse_structures():
            assert struct.formula == "Zr0.4 Nb1.6"
        assert parser.has_errors

        # Partial occupancy on sites, incorrect label, previously unparsable
        parser = CifParser(f"{TEST_FILES_DIR}/PF_sd_1615854.cif")
        for idx, struct in enumerate(parser.parse_structures()):
            if idx == 0:
                assert struct.formula == "Na2 Al2 Si6 O16"
            else:
                assert struct.formula == "Na4 Al4 Si12 O32"
        assert parser.has_errors

        # Partial occupancy on sites, incorrect label, previously unparsable
        parser = CifParser(f"{TEST_FILES_DIR}/PF_sd_1622133.cif")
        for struct in parser.parse_structures():
            assert struct.formula == "Ca0.184 Mg13.016 Fe2.8 Si16 O48"
        assert parser.has_errors

        # Partial occupancy on sites, previously parsed as an ordered structure
        parser = CifParser(f"{TEST_FILES_DIR}/PF_sd_1908491.cif")
        for struct in parser.parse_structures():
            assert struct.formula == "Mn0.96 Zn1.04 Ga4 Se8"
        assert parser.has_errors

        # Partial occupancy on sites, incorrect label, previously unparsable
        parser = CifParser(f"{TEST_FILES_DIR}/PF_sd_1811457.cif")
        for struct in parser.parse_structures():
            assert struct.formula == "Ba8 Mg2.4 Zr0.8 Ta4.8 O24"
        assert parser.has_errors

        # Incomplete powder diffraction data, previously unparsable
        # This CIF file contains the molecular species "NH3" which is
        # parsed as "N" because the label is "N{x}" (x = 1,2,..) and the
        # corresponding symbol is "NH3". Since, the label and symbol are switched
        # in CIFs from Springer Materials/Pauling file DBs, CifParser parses the
        # element as "Nh" (Nihonium).
        parser = CifParser(f"{TEST_FILES_DIR}/PF_sd_1002871.cif")
        assert parser.parse_structures()[0].formula == "Cu2 Br4 Nh12"
        assert parser.parse_structures()[1].formula == "Cu2 Br8 Nh12"
        assert parser.has_errors

        # Incomplete powder diffraction data, previously unparsable
        parser = CifParser(f"{TEST_FILES_DIR}/PF_sd_1704003.cif")
        for struct in parser.parse_structures():
            assert struct.formula == "Rb4 Mn2 F12"
        assert parser.has_errors

        # Unparsable species 'OH/OH2', previously parsed as "O"
        parser = CifParser(f"{TEST_FILES_DIR}/PF_sd_1500382.cif")
        for struct in parser.parse_structures():
            assert struct.formula == "Mg6 B2 O6 F1.764"
        assert parser.has_errors

        # Unparsable species 'OH/OH2', previously parsed as "O"
        parser = CifParser(f"{TEST_FILES_DIR}/PF_sd_1601634.cif")
        for struct in parser.parse_structures():
            assert struct.formula == "Zn2.58 Fe1.38 As4 Pb2.04 O16"

    def test_cif_parser_cod(self):
        """Parsing problematic CIF files from the COD database."""
        # Symbol in capital letters
        parser = CifParser(f"{TEST_FILES_DIR}/Cod_2100513.cif")
        for struct in parser.parse_structures():
            assert struct.formula == "Ca4 Nb2.0 Al2 O12"

        # Label in capital letters
        parser = CifParser(f"{TEST_FILES_DIR}/Cod_4115344.cif")
        for struct in parser.parse_structures():
            assert struct.formula == "Mo8 P4 H120 C120 I8 O8"

    def test_parse_symbol(self):
        """
        Test the _parse_symbol function with several potentially
        problematic examples of symbols and labels.
        """
        test_cases = {
            "MgT": "Mg",
            "MgT1": "Mg",
            "H(46A)": "H",
            "O(M)": "O",
            "N(Am)": "N",
            "H1N2a": "H",
            "CO(1)": "Co",
            "Wat1": "O",
            "MgM2A": "Mg",
            "CaX": "Ca",
            "X1": "X",
            "X": "X",
            "OA1": "O",
            "NaA2": "Na",
            "O-H2": "O",
            "OD2": "O",
            "OW": "O",
            "SiT": "Si",
            "SiTet": "Si",
            "Na-Int": "Na",
            "CaD1": "Ca",
            "KAm": "K",
            "D+1": "D",
            "D": "D",
            "D1-": "D",
            "D4": "D",
            "D0": "D",
            "NH": "Nh",
            "NH2": "Nh",
            "NH3": "Nh",
            "SH": "S",
        }

        for e in Element:
            name = e.name
            test_cases[name] = name
            if len(name) == 2:
                test_cases[name.upper()] = name
                test_cases[f"{name.upper()}1"] = name
                test_cases[name.upper() + "A"] = name
            test_cases[f"{name}1"] = name
            test_cases[f"{name}2"] = name
            test_cases[f"{name}3"] = name
            test_cases[f"{name}1A"] = name

        special = {"Hw": "H", "Ow": "O", "Wat": "O", "wat": "O", "OH": "", "OH2": ""}
        test_cases.update(special)

        parser = CifParser(f"{TEST_FILES_DIR}/LiFePO4.cif")
        for sym, expected_symbol in test_cases.items():
            assert parser._parse_symbol(sym) == expected_symbol

    def test_cif_writer(self):
        filepath = f"{TEST_FILES_DIR}/POSCAR"
        struct = Structure.from_file(filepath)
        writer = CifWriter(struct, symprec=0.01)
        answer = """# generated using pymatgen
data_FePO4
_symmetry_space_group_name_H-M   Pnma
_cell_length_a   10.41176687
_cell_length_b   6.06717188
_cell_length_c   4.75948954
_cell_angle_alpha   90.00000000
_cell_angle_beta   90.00000000
_cell_angle_gamma   90.00000000
_symmetry_Int_Tables_number   62
_chemical_formula_structural   FePO4
_chemical_formula_sum   'Fe4 P4 O16'
_cell_volume   300.65685512
_cell_formula_units_Z   4
loop_
 _symmetry_equiv_pos_site_id
 _symmetry_equiv_pos_as_xyz
  1  'x, y, z'
  2  '-x, -y, -z'
  3  '-x+1/2, -y, z+1/2'
  4  'x+1/2, y, -z+1/2'
  5  'x+1/2, -y+1/2, -z+1/2'
  6  '-x+1/2, y+1/2, z+1/2'
  7  '-x, y+1/2, -z'
  8  'x, -y+1/2, z'
loop_
 _atom_site_type_symbol
 _atom_site_label
 _atom_site_symmetry_multiplicity
 _atom_site_fract_x
 _atom_site_fract_y
 _atom_site_fract_z
 _atom_site_occupancy
  Fe  Fe0  4  0.21872822  0.75000000  0.47486711  1
  P  P1  4  0.09461309  0.25000000  0.41824327  1
  O  O2  8  0.16570974  0.04607233  0.28538394  1
  O  O3  4  0.04337231  0.75000000  0.70713767  1
  O  O4  4  0.09664244  0.25000000  0.74132035  1"""
        for l1, l2 in zip(str(writer).split("\n"), answer.split("\n")):
            assert l1.strip() == l2.strip()

    def test_symmetrized(self):
        filepath = f"{TEST_FILES_DIR}/POSCAR"
        struct = Structure.from_file(filepath)
        writer = CifWriter(struct, symprec=0.1)

        cif = CifParser.from_str(str(writer))
        m = StructureMatcher()

        assert m.fit(cif.parse_structures()[0], struct)

        # for l1, l2 in zip(str(writer).split("\n"), answer.split("\n")):
        #     assert l1.strip() == l2.strip()

        struct = Structure.from_file(f"{TEST_FILES_DIR}/LiFePO4.cif")
        writer = CifWriter(struct, symprec=0.1)
        s2 = CifParser.from_str(str(writer)).parse_structures()[0]

        assert m.fit(struct, s2)

        struct = self.get_structure("Li2O")
        writer = CifWriter(struct, symprec=0.1)
        s2 = CifParser.from_str(str(writer)).parse_structures()[0]
        assert m.fit(struct, s2)

        # test angle tolerance.
        struct = Structure.from_file(f"{TEST_FILES_DIR}/LiFePO4.cif")
        writer = CifWriter(struct, symprec=0.1, angle_tolerance=0)
        d = next(iter(writer.cif_file.data.values()))
        assert d["_symmetry_Int_Tables_number"] == 14
        struct = Structure.from_file(f"{TEST_FILES_DIR}/LiFePO4.cif")
        writer = CifWriter(struct, symprec=0.1, angle_tolerance=2)
        d = next(iter(writer.cif_file.data.values()))
        assert d["_symmetry_Int_Tables_number"] == 62

    def test_disordered(self):
        si = Element("Si")
        n = Element("N")
        coords = []
        coords.extend((np.array([0, 0, 0]), np.array([0.75, 0.5, 0.75])))
        lattice = Lattice(
            [
                [3.8401979337, 0.00, 0.00],
                [1.9200989668, 3.3257101909, 0.00],
                [0.00, -2.2171384943, 3.1355090603],
            ]
        )
        struct = Structure(lattice, [si, {si: 0.5, n: 0.5}], coords)
        writer = CifWriter(struct)
        answer = """# generated using pymatgen
data_Si1.5N0.5
_symmetry_space_group_name_H-M   'P 1'
_cell_length_a   3.84019793
_cell_length_b   3.84019899
_cell_length_c   3.84019793
_cell_angle_alpha   119.99999086
_cell_angle_beta   90.00000000
_cell_angle_gamma   60.00000914
_symmetry_Int_Tables_number   1
_chemical_formula_structural   Si1.5N0.5
_chemical_formula_sum   'Si1.5 N0.5'
_cell_volume   40.04479464
_cell_formula_units_Z   1
loop_
 _symmetry_equiv_pos_site_id
 _symmetry_equiv_pos_as_xyz
  1  'x, y, z'
loop_
 _atom_site_type_symbol
 _atom_site_label
 _atom_site_symmetry_multiplicity
 _atom_site_fract_x
 _atom_site_fract_y
 _atom_site_fract_z
 _atom_site_occupancy
  Si  Si0  1  0.00000000  0.00000000  0.00000000  1
  Si  Si1  1  0.75000000  0.50000000  0.75000000  0.5
  N  N2  1  0.75000000  0.50000000  0.75000000  0.5"""

        for l1, l2 in zip(str(writer).split("\n"), answer.split("\n")):
            assert l1.strip() == l2.strip()

    def test_cif_writer_without_refinement(self):
        si2 = Structure.from_file(f"{TEST_FILES_DIR}/abinit/si.cif")

        writer = CifWriter(si2, symprec=1e-3, significant_figures=10, refine_struct=False)
        cif_str = str(writer)
        assert "Fd-3m" in cif_str
        same_si2 = CifParser.from_str(cif_str).parse_structures()[0]
        assert len(si2) == len(same_si2)

    def test_specie_cif_writer(self):
        si4 = Species("Si", 4)
        si3 = Species("Si", 3)
        n = DummySpecies("X", -3)
        coords = []
        coords.extend((np.array([0.5, 0.5, 0.5]), np.array([0.75, 0.5, 0.75]), np.array([0, 0, 0])))
        lattice = Lattice(
            [
                [3.8401979337, 0.00, 0.00],
                [1.9200989668, 3.3257101909, 0.00],
                [0.00, -2.2171384943, 3.1355090603],
            ]
        )
        struct = Structure(lattice, [n, {si3: 0.5, n: 0.5}, si4], coords)
        writer = CifWriter(struct)
        answer = """# generated using pymatgen
data_X1.5Si1.5
_symmetry_space_group_name_H-M   'P 1'
_cell_length_a   3.84019793
_cell_length_b   3.84019899
_cell_length_c   3.84019793
_cell_angle_alpha   119.99999086
_cell_angle_beta   90.00000000
_cell_angle_gamma   60.00000914
_symmetry_Int_Tables_number   1
_chemical_formula_structural   X1.5Si1.5
_chemical_formula_sum   'X1.5 Si1.5'
_cell_volume   40.04479464
_cell_formula_units_Z   1
loop_
 _symmetry_equiv_pos_site_id
 _symmetry_equiv_pos_as_xyz
  1  'x, y, z'
loop_
 _atom_type_symbol
 _atom_type_oxidation_number
  X3-  -3.0
  Si3+  3.0
  Si4+  4.0
loop_
 _atom_site_type_symbol
 _atom_site_label
 _atom_site_symmetry_multiplicity
 _atom_site_fract_x
 _atom_site_fract_y
 _atom_site_fract_z
 _atom_site_occupancy
  X3-  X0  1  0.50000000  0.50000000  0.50000000  1
  X3-  X1  1  0.75000000  0.50000000  0.75000000  0.5
  Si3+  Si2  1  0.75000000  0.50000000  0.75000000  0.5
  Si4+  Si3  1  0.00000000  0.00000000  0.00000000  1
"""
        for l1, l2 in zip(str(writer).split("\n"), answer.split("\n")):
            assert l1.strip() == l2.strip()

        # test that mixed valence works properly
        s2 = Structure.from_str(answer, "cif")
        assert struct.composition == s2.composition

    def test_primes(self):
        parser = CifParser(f"{TEST_FILES_DIR}/C26H16BeN2O2S2.cif")
        for struct in parser.parse_structures(primitive=False):
            assert struct.composition == 8 * Composition("C26H16BeN2O2S2")

    def test_missing_atom_site_type_with_oxi_states(self):
        parser = CifParser(f"{TEST_FILES_DIR}/P24Ru4H252C296S24N16.cif")
        comp = Composition({"S0+": 24, "Ru0+": 4, "H0+": 252, "C0+": 296, "N0+": 16, "P0+": 24})
        for struct in parser.parse_structures(primitive=False):
            assert struct.composition == comp

    def test_no_coords_or_species(self):
        string = """#generated using pymatgen
    data_Si1.5N1.5
    _symmetry_space_group_name_H-M   'P 1'
    _cell_length_a   3.84019793
    _cell_length_b   3.84019899
    _cell_length_c   3.84019793
    _cell_angle_alpha   119.99999086
    _cell_angle_beta   90.00000000
    _cell_angle_gamma   60.00000914
    _symmetry_Int_Tables_number   1
    _chemical_formula_structural   Si1.5N1.5
    _chemical_formula_sum   'Si1.5 N1.5'
    _cell_volume   40.0447946443
    _cell_formula_units_Z   0
    loop_
      _symmetry_equiv_pos_site_id
      _symmetry_equiv_pos_as_xyz
      1  'x, y, z'
    loop_
      _atom_type_symbol
      _atom_type_oxidation_number
       Si3+  3.0
       Si4+  4.0
       N3-  -3.0
    loop_
      _atom_site_type_symbol
      _atom_site_label
      _atom_site_symmetry_multiplicity
      _atom_site_fract_x
      _atom_site_fract_y
      _atom_site_fract_z
      _atom_site_occupancy
      ? ? ? ? ? ? ?
    """
        parser = CifParser.from_str(string)
        with pytest.raises(ValueError, match="Invalid CIF file with no structures"):
            parser.parse_structures()

    def test_get_lattice_from_lattice_type(self):
        cif_structure = """#generated using pymatgen
data_FePO4
_symmetry_space_group_name_H-M   Pnma
_cell_length_a   10.41176687
_cell_length_b   6.06717188
_cell_length_c   4.75948954
_chemical_formula_structural   FePO4
_chemical_formula_sum   'Fe4 P4 O16'
_cell_volume   300.65685512
_cell_formula_units_Z   4
_symmetry_cell_setting Orthorhombic
loop_
  _symmetry_equiv_pos_site_id
  _symmetry_equiv_pos_as_xyz
   1  'x, y, z'
loop_
  _atom_site_type_symbol
  _atom_site_label
  _atom_site_symmetry_multiplicity
  _atom_site_fract_x
  _atom_site_fract_y
  _atom_site_fract_z
  _atom_site_occupancy
  Fe  Fe1  1  0.218728  0.750000  0.474867  1
  Fe  Fe2  1  0.281272  0.250000  0.974867  1
  Fe  Fe3  1  0.718728  0.750000  0.025133  1
  Fe  Fe4  1  0.781272  0.250000  0.525133  1
  P  P5  1  0.094613  0.250000  0.418243  1
  P  P6  1  0.405387  0.750000  0.918243  1
  P  P7  1  0.594613  0.250000  0.081757  1
  P  P8  1  0.905387  0.750000  0.581757  1
  O  O9  1  0.043372  0.750000  0.707138  1
  O  O10  1  0.096642  0.250000  0.741320  1
  O  O11  1  0.165710  0.046072  0.285384  1
  O  O12  1  0.165710  0.453928  0.285384  1
  O  O13  1  0.334290  0.546072  0.785384  1
  O  O14  1  0.334290  0.953928  0.785384  1
  O  O15  1  0.403358  0.750000  0.241320  1
  O  O16  1  0.456628  0.250000  0.207138  1
  O  O17  1  0.543372  0.750000  0.792862  1
  O  O18  1  0.596642  0.250000  0.758680  1
  O  O19  1  0.665710  0.046072  0.214616  1
  O  O20  1  0.665710  0.453928  0.214616  1
  O  O21  1  0.834290  0.546072  0.714616  1
  O  O22  1  0.834290  0.953928  0.714616  1
  O  O23  1  0.903358  0.750000  0.258680  1
  O  O24  1  0.956628  0.250000  0.292862  1

"""
        parser = CifParser.from_str(cif_structure)
        s_test = parser.parse_structures(primitive=False)[0]
        filepath = f"{TEST_FILES_DIR}/POSCAR"
        struct = Structure.from_file(filepath)

        sm = StructureMatcher(stol=0.05, ltol=0.01, angle_tol=0.1)
        assert sm.fit(struct, s_test)

    def test_empty(self):
        # single line
        cb = CifBlock.from_str("data_mwe\nloop_\n_tag\n ''")
        assert cb.data["_tag"][0] == ""

        # multi line
        cb = CifBlock.from_str("data_mwe\nloop_\n_tag\n;\n;")
        assert cb.data["_tag"][0] == ""

        cb2 = CifBlock.from_str(str(cb))
        assert cb == cb2

    def test_bad_cif(self):
        filepath = f"{TEST_FILES_DIR}/bad_occu.cif"
        parser = CifParser(filepath)
        with pytest.raises(
            ValueError, match="No structure parsed for section 1 in CIF.\nSpecies occupancies sum to more than 1!"
        ):
            parser.parse_structures(on_error="raise")
        parser = CifParser(filepath, occupancy_tolerance=2)
        struct = parser.parse_structures()[0]
        assert struct[0].species["Al3+"] == approx(0.778)

    def test_one_line_symm(self):
        f = f"{TEST_FILES_DIR}/OneLineSymmP1.cif"
        parser = CifParser(f)
        struct = parser.parse_structures()[0]
        assert struct.formula == "Ga4 Pb2 O8"

    def test_no_symmops(self):
        f = f"{TEST_FILES_DIR}/nosymm.cif"
        parser = CifParser(f)
        struct = parser.parse_structures()[0]
        assert struct.formula == "H96 C60 O8"

    def test_dot_positions(self):
        f = f"{TEST_FILES_DIR}/ICSD59959.cif"
        parser = CifParser(f)
        struct = parser.parse_structures()[0]
        assert struct.formula == "K1 Mn1 F3"

    def test_replacing_finite_precision_frac_coords(self):
        cif = f"{TEST_FILES_DIR}/cif_finite_precision_frac_coord_error.cif"
        parser = CifParser(cif)
        warn_msg = "4 fractional coordinates rounded to ideal values to avoid issues with finite precision."
        with pytest.warns(UserWarning, match=warn_msg):
            struct = parser.parse_structures()[0]
        assert str(struct.composition) == "N5+72"
        assert warn_msg in parser.warnings

    def test_empty_deque(self):
        cif_str = """data_1526655
_journal_name_full
_space_group_IT_number           227
_symmetry_space_group_name_Hall  'F 4d 2 3 -1d'
_symmetry_space_group_name_H-M   'F d -3 m :1'
_cell_angle_alpha                90
_cell_angle_beta                 90
_cell_angle_gamma                90
_cell_formula_units_Z            8
_cell_length_a                   5.381
_cell_length_b                   5.381
_cell_length_c                   5.381
_cell_volume                     155.808
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
_atom_site_U_iso_or_equiv
Si1 Si 0 0 0 1 0.0
_iucr_refine_fcf_details
;
data_symmetries
loop_
  _space_group_symop_id
  _space_group_symop_operation_xyz
  1  x,y,z
  2  -x+1/2,y+1/2,-z+1/2
  3  -x,-y,-z
  4  x-1/2,-y-1/2,z-1/2
;"""
        parser = CifParser.from_str(cif_str)
        assert parser.parse_structures()[0].formula == "Si1"
        cif = """
data_1526655
_journal_name_full
_space_group_IT_number           227
_symmetry_space_group_name_Hall  'F 4d 2 3 -1d'
_symmetry_space_group_name_H-M   'F d -3 m :1'
_cell_angle_alpha                90
_cell_angle_beta                 90
_cell_angle_gamma                90
_cell_formula_units_Z            8
_cell_length_a                   5.381
_cell_length_b                   5.381
_cell_length_c                   5.381
_cell_volume                     155.808
_iucr_refine_fcf_details
;
data_symmetries
Some arbitrary multiline string
;
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
_atom_site_U_iso_or_equiv
Si1 Si 0 0 0 1 0.0
"""
        parser = CifParser.from_str(cif)
        with pytest.raises(ValueError, match="Invalid CIF file with no structures"):
            parser.parse_structures()

    def test_no_check_occu(self):
        with open(f"{TEST_FILES_DIR}/site_type_symbol_test.cif") as cif_file:
            cif_str = cif_file.read()
        cif_str = cif_str.replace("Te    Te 1.0000", "Te    Te 1.5000", 1)

        with pytest.raises(ValueError, match="Invalid CIF file with no structures"):
            # should fail without setting custom occupancy tolerance
            CifParser.from_str(cif_str).parse_structures()

        for tol in (1.5, 10):
            parser = CifParser.from_str(cif_str, occupancy_tolerance=tol)
            structs = parser.parse_structures(primitive=False, check_occu=False)[0]
            assert structs[0].species.as_dict()["Te"] == 1.5

    def test_cif_writer_write_file(self):
        struct1 = Structure.from_file(f"{TEST_FILES_DIR}/POSCAR")
        out_path = f"{self.tmp_path}/test.cif"
        CifWriter(struct1).write_file(out_path)
        read_structs = CifParser(out_path).parse_structures()
        assert len(read_structs) == 1
        assert struct1.matches(read_structs[0])

        # test write_file append mode='a'
        struct2 = Structure.from_file(f"{TEST_FILES_DIR}/Graphite.cif")
        CifWriter(struct2).write_file(out_path, mode="a")

        read_structs = CifParser(out_path).parse_structures()
        assert len(read_structs) == 2
        assert [x.formula for x in read_structs] == ["Fe4 P4 O16", "C4"]


class TestMagCif(PymatgenTest):
    def setUp(self):
        self.mcif = CifParser(f"{TEST_FILES_DIR}/magnetic.example.NiO.mcif")
        self.mcif_ncl = CifParser(f"{TEST_FILES_DIR}/magnetic.ncl.example.GdB4.mcif")
        self.mcif_incommensurate = CifParser(f"{TEST_FILES_DIR}/magnetic.incommensurate.example.Cr.mcif")
        self.mcif_disordered = CifParser(f"{TEST_FILES_DIR}/magnetic.disordered.example.CuMnO2.mcif")
        self.mcif_ncl2 = CifParser(f"{TEST_FILES_DIR}/Mn3Ge_IR2.mcif")

    def test_mcif_detection(self):
        assert self.mcif.feature_flags["magcif"]
        assert self.mcif_ncl.feature_flags["magcif"]
        assert self.mcif_incommensurate.feature_flags["magcif"]
        assert self.mcif_disordered.feature_flags["magcif"]
        assert not self.mcif.feature_flags["magcif_incommensurate"]
        assert not self.mcif_ncl.feature_flags["magcif_incommensurate"]
        assert self.mcif_incommensurate.feature_flags["magcif_incommensurate"]
        assert not self.mcif_disordered.feature_flags["magcif_incommensurate"]

    def test_parse_structures(self):
        # incommensurate structures not currently supported
        with pytest.raises(NotImplementedError, match="Incommensurate structures not currently supported"):
            self.mcif_incommensurate.parse_structures()

        # disordered magnetic structures not currently supported
        with pytest.raises(NotImplementedError, match="Disordered magnetic structures not currently supported"):
            self.mcif_disordered.parse_structures()

        # taken from self.mcif_ncl, removing explicit magnetic symmops
        # so that MagneticSymmetryGroup() has to be invoked
        mag_cif_str = """
data_5yOhtAoR

_space_group.magn_name_BNS     "P 4/m' b' m' "
_cell_length_a                 7.1316
_cell_length_b                 7.1316
_cell_length_c                 4.0505
_cell_angle_alpha              90.00
_cell_angle_beta               90.00
_cell_angle_gamma              90.00

loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Gd1 Gd 0.31746 0.81746 0.00000 1
B1 B 0.00000 0.00000 0.20290 1
B2 B 0.17590 0.03800 0.50000 1
B3 B 0.08670 0.58670 0.50000 1

loop_
_atom_site_moment_label
_atom_site_moment_crystalaxis_x
_atom_site_moment_crystalaxis_y
_atom_site_moment_crystalaxis_z
Gd1 5.05 5.05 0.0"""

        struct = self.mcif.parse_structures(primitive=False)[0]
        assert struct.formula == "Ni32 O32"
        assert Magmom.are_collinear(struct.site_properties["magmom"])

        # example with non-collinear spin
        s_ncl = self.mcif_ncl.parse_structures(primitive=False)[0]
        s_ncl_from_msg = CifParser.from_str(mag_cif_str).parse_structures(primitive=False)[0]
        assert s_ncl.formula == "Gd4 B16"
        assert not Magmom.are_collinear(s_ncl.site_properties["magmom"])

        assert s_ncl.matches(s_ncl_from_msg)

    def test_write(self):
        with open(f"{TEST_FILES_DIR}/GdB4-writer-ref.mcif") as file:
            cw_ref_string = file.read()
        s_ncl = self.mcif_ncl.parse_structures(primitive=False)[0]

        cw = CifWriter(s_ncl, write_magmoms=True)
        assert str(cw) == cw_ref_string

        # from list-type magmoms
        list_magmoms = [list(m) for m in s_ncl.site_properties["magmom"]]

        # float magmoms (magnitude only)
        float_magmoms = [float(m) for m in s_ncl.site_properties["magmom"]]

        s_ncl.add_site_property("magmom", list_magmoms)
        cw = CifWriter(s_ncl, write_magmoms=True)
        assert str(cw) == cw_ref_string

        s_ncl.add_site_property("magmom", float_magmoms)
        cw = CifWriter(s_ncl, write_magmoms=True)

        with open(f"{TEST_FILES_DIR}/GdB4-str-magnitudes-ref.mcif") as file:
            cw_ref_string_magnitudes = file.read()

        assert str(cw).strip() == cw_ref_string_magnitudes.strip()
        # test we're getting correct magmoms in ncl case
        s_ncl2 = self.mcif_ncl2.parse_structures()[0]
        list_magmoms = [list(m) for m in s_ncl2.site_properties["magmom"]]
        assert list_magmoms[0][0] == 0.0
        assert list_magmoms[0][1] == approx(5.9160793408726366)
        assert list_magmoms[1][0] == approx(-5.1234749999999991)
        assert list_magmoms[1][1] == approx(2.9580396704363183)

        # test creating a structure without oxidation state doesn't raise errors
        s_manual = Structure(Lattice.cubic(4.2), ["Cs", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        s_manual.add_spin_by_site([1, -1])
        cw = CifWriter(s_manual, write_magmoms=True)

        # check oxidation state
        with open(f"{TEST_FILES_DIR}/CsCl-manual-oxi-ref.mcif") as file:
            cw_manual_oxi_string = file.read()
        s_manual.add_oxidation_state_by_site([1, 1])
        cw = CifWriter(s_manual, write_magmoms=True)
        assert str(cw) == cw_manual_oxi_string

    @unittest.skipIf(pybtex is None, "pybtex not present")
    def test_bibtex(self):
        ref_bibtex_string = """@article{cifref0,
    author = "Blanco, J.A.",
    journal = "PHYSICAL REVIEW B",
    volume = "73",
    year = "2006",
    pages = "?--?"
}
"""
        assert self.mcif_ncl.get_bibtex_string() == ref_bibtex_string
