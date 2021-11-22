# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import unittest
import warnings

import numpy as np

from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core.composition import Composition
from pymatgen.core.lattice import Lattice
from pymatgen.core.periodic_table import DummySpecies, Element, Species
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.core import Magmom
from pymatgen.io.cif import CifBlock, CifParser, CifWriter
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.symmetry.structure import SymmetrizedStructure
from pymatgen.util.testing import PymatgenTest

try:
    import pybtex
except ImportError:
    pybtex = None


class CifBlockTest(PymatgenTest):
    def test_to_string(self):
        with open(self.TEST_FILES_DIR / "Graphite.cif") as f:
            s = f.read()
        c = CifBlock.from_string(s)
        cif_str_2 = str(CifBlock.from_string(str(c)))
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
        for l1, l2, l3 in zip(str(c).split("\n"), cif_str.split("\n"), cif_str_2.split("\n")):
            self.assertEqual(l1.strip(), l2.strip())
            self.assertEqual(l2.strip(), l3.strip())

    def test_double_quotes_and_underscore_data(self):
        cif_str = """data_test
_symmetry_space_group_name_H-M   "P -3 m 1"
_thing   '_annoying_data'"""
        cb = CifBlock.from_string(cif_str)
        self.assertEqual(cb["_symmetry_space_group_name_H-M"], "P -3 m 1")
        self.assertEqual(cb["_thing"], "_annoying_data")
        self.assertEqual(str(cb), cif_str.replace('"', "'"))

    def test_double_quoted_data(self):
        cif_str = """data_test
_thing   ' '_annoying_data''
_other   " "_more_annoying_data""
_more   ' "even more" ' """
        cb = CifBlock.from_string(cif_str)
        self.assertEqual(cb["_thing"], " '_annoying_data'")
        self.assertEqual(cb["_other"], ' "_more_annoying_data"')
        self.assertEqual(cb["_more"], ' "even more" ')

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
        cb = CifBlock.from_string(cif_str)
        self.assertEqual(
            cb["_thing"],
            " long quotes  ;  still in the quote  ; actually going to end now",
        )

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
        self.assertEqual(str(CifBlock(data, loops, "test")), cif_str)


class CifIOTest(PymatgenTest):
    def test_CifParser(self):
        parser = CifParser(self.TEST_FILES_DIR / "LiFePO4.cif")
        for s in parser.get_structures(True):
            self.assertEqual(s.formula, "Li4 Fe4 P4 O16", "Incorrectly parsed cif.")

        parser = CifParser(self.TEST_FILES_DIR / "V2O3.cif")
        for s in parser.get_structures(True):
            self.assertEqual(s.formula, "V4 O6")

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
        self.assertEqual(parser.get_bibtex_string().strip(), bibtex_str.strip())

        parser = CifParser(self.TEST_FILES_DIR / "Li2O.cif")
        prim = parser.get_structures(True)[0]
        self.assertEqual(prim.formula, "Li2 O1")
        conv = parser.get_structures(False)[0]
        self.assertEqual(conv.formula, "Li8 O4")

        # test for disordered structures
        parser = CifParser(self.TEST_FILES_DIR / "Li10GeP2S12.cif")
        for s in parser.get_structures(True):
            self.assertEqual(s.formula, "Li20.2 Ge2.06 P3.94 S24", "Incorrectly parsed cif.")
        cif_str = r"""#\#CIF1.1
##########################################################################
#               Crystallographic Information Format file
#               Produced by PyCifRW module
#
#  This is a CIF file.  CIF has been adopted by the International
#  Union of Crystallography as the standard for data archiving and
#  transmission.
#
#  For information on this file format, follow the CIF links at
#  http://www.iucr.org
##########################################################################

data_FePO4
_symmetry_space_group_name_H-M          'P 1'
_cell_length_a                          10.4117668699
_cell_length_b                          6.06717187997
_cell_length_c                          4.75948953998
loop_ # sometimes this is in a loop (incorrectly)
_cell_angle_alpha
91.0
_cell_angle_beta                        92.0
_cell_angle_gamma                       93.0
_chemical_name_systematic               'Generated by pymatgen'
_symmetry_Int_Tables_number             1
_chemical_formula_structural            FePO4
_chemical_formula_sum                   'Fe4 P4 O16'
_cell_volume                            300.65685512
_cell_formula_units_Z                   4
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
  _atom_site_attached_hydrogens
  _atom_site_B_iso_or_equiv
  _atom_site_occupancy
    Fe  Fe1  1  0.218728  0.750000  0.474867  0  .  1
    Fe  JJ2  1  0.281272  0.250000  0.974867  0  .  1
    # there's a typo here, parser should read the symbol from the
    # _atom_site_type_symbol
    Fe  Fe3  1  0.718728  0.750000  0.025133  0  .  1
    Fe  Fe4  1  0.781272  0.250000  0.525133  0  .  1
    P  P5  1  0.094613  0.250000  0.418243  0  .  1
    P  P6  1  0.405387  0.750000  0.918243  0  .  1
    P  P7  1  0.594613  0.250000  0.081757  0  .  1
    P  P8  1  0.905387  0.750000  0.581757  0  .  1
    O  O9  1  0.043372  0.750000  0.707138  0  .  1
    O  O10  1  0.096642  0.250000  0.741320  0  .  1
    O  O11  1  0.165710  0.046072  0.285384  0  .  1
    O  O12  1  0.165710  0.453928  0.285384  0  .  1
    O  O13  1  0.334290  0.546072  0.785384  0  .  1
    O  O14  1  0.334290  0.953928  0.785384  0  .  1
    O  O15  1  0.403358  0.750000  0.241320  0  .  1
    O  O16  1  0.456628  0.250000  0.207138  0  .  1
    O  O17  1  0.543372  0.750000  0.792862  0  .  1
    O  O18  1  0.596642  0.250000  0.758680  0  .  1
    O  O19  1  0.665710  0.046072  0.214616  0  .  1
    O  O20  1  0.665710  0.453928  0.214616  0  .  1
    O  O21  1  0.834290  0.546072  0.714616  0  .  1
    O  O22  1  0.834290  0.953928  0.714616  0  .  1
    O  O23  1  0.903358  0.750000  0.258680  0  .  1
    O  O24  1  0.956628  0.250000  0.292862  0  .  1

"""
        parser = CifParser.from_string(cif_str)
        struct = parser.get_structures(primitive=False)[0]
        self.assertEqual(struct.formula, "Fe4 P4 O16")
        self.assertAlmostEqual(struct.lattice.a, 10.4117668699)
        self.assertAlmostEqual(struct.lattice.b, 6.06717187997)
        self.assertAlmostEqual(struct.lattice.c, 4.75948953998)
        self.assertAlmostEqual(struct.lattice.alpha, 91)
        self.assertAlmostEqual(struct.lattice.beta, 92)
        self.assertAlmostEqual(struct.lattice.gamma, 93)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            parser = CifParser(self.TEST_FILES_DIR / "srycoo.cif")
        self.assertEqual(parser.get_structures()[0].formula, "Sr5.6 Y2.4 Co8 O21")

        # Test with a decimal Xyz. This should parse as two atoms in
        # conventional cell if it is correct, one if not.
        parser = CifParser(self.TEST_FILES_DIR / "Fe.cif")
        self.assertEqual(len(parser.get_structures(primitive=False)[0]), 2)
        self.assertFalse(parser.has_errors)

    def test_get_symmetrized_structure(self):
        parser = CifParser(self.TEST_FILES_DIR / "Li2O.cif")
        sym_structure = parser.get_structures(primitive=False, symmetrized=True)[0]
        structure = parser.get_structures(primitive=False, symmetrized=False)[0]
        self.assertIsInstance(sym_structure, SymmetrizedStructure)
        self.assertEqual(structure, sym_structure)
        self.assertEqual(sym_structure.equivalent_indices, [[0, 1, 2, 3], [4, 5, 6, 7, 8, 9, 10, 11]])

    def test_site_symbol_preference(self):
        parser = CifParser(self.TEST_FILES_DIR / "site_type_symbol_test.cif")
        self.assertEqual(parser.get_structures()[0].formula, "Ge0.4 Sb0.4 Te1")

    def test_implicit_hydrogen(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            parser = CifParser(self.TEST_FILES_DIR / "Senegalite_implicit_hydrogen.cif")
            for s in parser.get_structures():
                self.assertEqual(s.formula, "Al8 P4 O32")
                self.assertEqual(sum(s.site_properties["implicit_hydrogens"]), 20)
            self.assertIn(
                "Structure has implicit hydrogens defined, "
                "parsed structure unlikely to be suitable for use "
                "in calculations unless hydrogens added.",
                parser.warnings,
            )
            parser = CifParser(self.TEST_FILES_DIR / "cif_implicit_hydrogens_cod_1011130.cif")
            s = parser.get_structures()[0]
            self.assertIn(
                "Structure has implicit hydrogens defined, "
                "parsed structure unlikely to be suitable for use "
                "in calculations unless hydrogens added.",
                parser.warnings,
            )

    def test_CifParserSpringerPauling(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # Below are 10 tests for CIFs from the Springer Materials/Pauling file DBs.

            # Partial occupancy on sites, incorrect label, previously unparsable
            parser = CifParser(self.TEST_FILES_DIR / "PF_sd_1928405.cif")
            for s in parser.get_structures(True):
                self.assertEqual(s.formula, "Er1 Mn3.888 Fe2.112 Sn6")
            self.assertTrue(parser.has_errors)

            # Partial occupancy on sites, previously parsed as an ordered structure
            parser = CifParser(self.TEST_FILES_DIR / "PF_sd_1011081.cif")
            for s in parser.get_structures(True):
                self.assertEqual(s.formula, "Zr0.2 Nb0.8")
            self.assertTrue(parser.has_errors)

            # Partial occupancy on sites, incorrect label, previously unparsable
            parser = CifParser(self.TEST_FILES_DIR / "PF_sd_1615854.cif")
            for s in parser.get_structures(True):
                self.assertEqual(s.formula, "Na2 Al2 Si6 O16")
            self.assertTrue(parser.has_errors)

            # Partial occupancy on sites, incorrect label, previously unparsable
            parser = CifParser(self.TEST_FILES_DIR / "PF_sd_1622133.cif")
            for s in parser.get_structures(True):
                self.assertEqual(s.formula, "Ca0.184 Mg13.016 Fe2.8 Si16 O48")
            self.assertTrue(parser.has_errors)

            # Partial occupancy on sites, previously parsed as an ordered structure
            parser = CifParser(self.TEST_FILES_DIR / "PF_sd_1908491.cif")
            for s in parser.get_structures(True):
                self.assertEqual(s.formula, "Mn0.48 Zn0.52 Ga2 Se4")
            self.assertTrue(parser.has_errors)

            # Partial occupancy on sites, incorrect label, previously unparsable
            parser = CifParser(self.TEST_FILES_DIR / "PF_sd_1811457.cif")
            for s in parser.get_structures(True):
                self.assertEqual(s.formula, "Ba2 Mg0.6 Zr0.2 Ta1.2 O6")
            self.assertTrue(parser.has_errors)

            # Incomplete powder diffraction data, previously unparsable
            # This CIF file contains the molecular species "NH3" which is
            # parsed as "N" because the label is "N{x}" (x = 1,2,..) and the
            # corresponding symbol is "NH3". Since, the label and symbol are switched
            # in CIFs from Springer Materials/Pauling file DBs, CifParser parses the
            # element as "Nh" (Nihonium).
            parser = CifParser(self.TEST_FILES_DIR / "PF_sd_1002871.cif")
            self.assertEqual(parser.get_structures(True)[0].formula, "Cu1 Br2 Nh6")
            self.assertEqual(parser.get_structures(True)[1].formula, "Cu1 Br4 Nh6")
            self.assertTrue(parser.has_errors)

            # Incomplete powder diffraction data, previously unparsable
            parser = CifParser(self.TEST_FILES_DIR / "PF_sd_1704003.cif")
            for s in parser.get_structures():
                self.assertEqual(s.formula, "Rb4 Mn2 F12")
            self.assertTrue(parser.has_errors)

            # Unparsable species 'OH/OH2', previously parsed as "O"
            parser = CifParser(self.TEST_FILES_DIR / "PF_sd_1500382.cif")
            for s in parser.get_structures():
                self.assertEqual(s.formula, "Mg6 B2 O6 F1.764")
            self.assertTrue(parser.has_errors)

            # Unparsable species 'OH/OH2', previously parsed as "O"
            parser = CifParser(self.TEST_FILES_DIR / "PF_sd_1601634.cif")
            for s in parser.get_structures():
                self.assertEqual(s.formula, "Zn1.29 Fe0.69 As2 Pb1.02 O8")

    def test_CifParserCod(self):
        """
        Parsing problematic cif files from the COD database
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Symbol in capital letters
            parser = CifParser(self.TEST_FILES_DIR / "Cod_2100513.cif")
            for s in parser.get_structures(True):
                self.assertEqual(s.formula, "Ca4 Nb2.0 Al2 O12")

            # Label in capital letters
            parser = CifParser(self.TEST_FILES_DIR / "Cod_4115344.cif")
            for s in parser.get_structures(True):
                self.assertEqual(s.formula, "Mo4 P2 H60 C60 I4 O4")

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
                test_cases[name.upper() + str(1)] = name
                test_cases[name.upper() + "A"] = name
            test_cases[name + str(1)] = name
            test_cases[name + str(2)] = name
            test_cases[name + str(3)] = name
            test_cases[name + str(1) + "A"] = name

        special = {"Hw": "H", "Ow": "O", "Wat": "O", "wat": "O", "OH": "", "OH2": ""}
        test_cases.update(special)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            parser = CifParser(self.TEST_FILES_DIR / "LiFePO4.cif")
            for sym, expected_symbol in test_cases.items():
                self.assertEqual(parser._parse_symbol(sym), expected_symbol)

    def test_CifWriter(self):
        filepath = self.TEST_FILES_DIR / "POSCAR"
        poscar = Poscar.from_file(filepath)
        writer = CifWriter(poscar.structure, symprec=0.01)
        ans = """# generated using pymatgen
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
        for l1, l2 in zip(str(writer).split("\n"), ans.split("\n")):
            self.assertEqual(l1.strip(), l2.strip())

    def test_symmetrized(self):
        filepath = self.TEST_FILES_DIR / "POSCAR"
        poscar = Poscar.from_file(filepath, check_for_POTCAR=False)
        writer = CifWriter(poscar.structure, symprec=0.1)

        cif = CifParser.from_string(str(writer))
        m = StructureMatcher()

        self.assertTrue(m.fit(cif.get_structures()[0], poscar.structure))

        # for l1, l2 in zip(str(writer).split("\n"), ans.split("\n")):
        #     self.assertEqual(l1.strip(), l2.strip())

        s = Structure.from_file(self.TEST_FILES_DIR / "LiFePO4.cif")
        writer = CifWriter(s, symprec=0.1)
        s2 = CifParser.from_string(str(writer)).get_structures()[0]

        self.assertTrue(m.fit(s, s2))

        s = self.get_structure("Li2O")
        writer = CifWriter(s, symprec=0.1)
        s2 = CifParser.from_string(str(writer)).get_structures()[0]
        self.assertTrue(m.fit(s, s2))

        # test angle tolerance.
        s = Structure.from_file(self.TEST_FILES_DIR / "LiFePO4.cif")
        writer = CifWriter(s, symprec=0.1, angle_tolerance=0)
        d = list(writer.ciffile.data.values())[0]
        self.assertEqual(d["_symmetry_Int_Tables_number"], 14)
        s = Structure.from_file(self.TEST_FILES_DIR / "LiFePO4.cif")
        writer = CifWriter(s, symprec=0.1, angle_tolerance=2)
        d = list(writer.ciffile.data.values())[0]
        self.assertEqual(d["_symmetry_Int_Tables_number"], 62)

    def test_disordered(self):
        si = Element("Si")
        n = Element("N")
        coords = []
        coords.append(np.array([0, 0, 0]))
        coords.append(np.array([0.75, 0.5, 0.75]))
        lattice = Lattice(
            np.array(
                [
                    [3.8401979337, 0.00, 0.00],
                    [1.9200989668, 3.3257101909, 0.00],
                    [0.00, -2.2171384943, 3.1355090603],
                ]
            )
        )
        struct = Structure(lattice, [si, {si: 0.5, n: 0.5}], coords)
        writer = CifWriter(struct)
        ans = """# generated using pymatgen
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

        for l1, l2 in zip(str(writer).split("\n"), ans.split("\n")):
            self.assertEqual(l1.strip(), l2.strip())

    def test_cifwrite_without_refinement(self):
        si2 = Structure.from_file(self.TEST_FILES_DIR / "abinit" / "si.cif")
        writer = CifWriter(si2, symprec=1e-3, significant_figures=10, refine_struct=False)
        s = str(writer)
        assert "Fd-3m" in s
        same_si2 = CifParser.from_string(s).get_structures()[0]
        assert len(si2) == len(same_si2)

    def test_specie_cifwriter(self):
        si4 = Species("Si", 4)
        si3 = Species("Si", 3)
        n = DummySpecies("X", -3)
        coords = []
        coords.append(np.array([0.5, 0.5, 0.5]))
        coords.append(np.array([0.75, 0.5, 0.75]))
        coords.append(np.array([0, 0, 0]))
        lattice = Lattice(
            np.array(
                [
                    [3.8401979337, 0.00, 0.00],
                    [1.9200989668, 3.3257101909, 0.00],
                    [0.00, -2.2171384943, 3.1355090603],
                ]
            )
        )
        struct = Structure(lattice, [n, {si3: 0.5, n: 0.5}, si4], coords)
        writer = CifWriter(struct)
        ans = """# generated using pymatgen
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
        for l1, l2 in zip(str(writer).split("\n"), ans.split("\n")):
            self.assertEqual(l1.strip(), l2.strip())

        # test that mixed valence works properly
        s2 = Structure.from_str(ans, "cif")
        self.assertEqual(struct.composition, s2.composition)

    def test_primes(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            parser = CifParser(self.TEST_FILES_DIR / "C26H16BeN2O2S2.cif")
            for s in parser.get_structures(False):
                self.assertEqual(s.composition, 8 * Composition("C26H16BeN2O2S2"))

    def test_missing_atom_site_type_with_oxistates(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            parser = CifParser(self.TEST_FILES_DIR / "P24Ru4H252C296S24N16.cif")
            c = Composition({"S0+": 24, "Ru0+": 4, "H0+": 252, "C0+": 296, "N0+": 16, "P0+": 24})
            for s in parser.get_structures(False):
                self.assertEqual(s.composition, c)

    def test_no_coords_or_species(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
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
            parser = CifParser.from_string(string)
            self.assertRaises(ValueError, parser.get_structures)

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
        cp = CifParser.from_string(cif_structure)
        s_test = cp.get_structures(False)[0]
        filepath = self.TEST_FILES_DIR / "POSCAR"
        poscar = Poscar.from_file(filepath)
        s_ref = poscar.structure

        sm = StructureMatcher(stol=0.05, ltol=0.01, angle_tol=0.1)
        self.assertTrue(sm.fit(s_ref, s_test))

    def test_empty(self):
        # single line
        cb = CifBlock.from_string("data_mwe\nloop_\n_tag\n ''")
        self.assertEqual(cb.data["_tag"][0], "")

        # multi line
        cb = CifBlock.from_string("data_mwe\nloop_\n_tag\n;\n;")
        self.assertEqual(cb.data["_tag"][0], "")

        cb2 = CifBlock.from_string(str(cb))
        self.assertEqual(cb, cb2)

    def test_bad_cif(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            f = self.TEST_FILES_DIR / "bad_occu.cif"
            p = CifParser(f)
            self.assertRaises(ValueError, p.get_structures)
            p = CifParser(f, occupancy_tolerance=2)
            s = p.get_structures()[0]
            self.assertAlmostEqual(s[0].species["Al3+"], 0.5)

    def test_one_line_symm(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            f = self.TEST_FILES_DIR / "OneLineSymmP1.cif"
            p = CifParser(f)
            s = p.get_structures()[0]
            self.assertEqual(s.formula, "Ga4 Pb2 O8")

    def test_no_symmops(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            f = self.TEST_FILES_DIR / "nosymm.cif"
            p = CifParser(f)
            s = p.get_structures()[0]
            self.assertEqual(s.formula, "H96 C60 O8")

    def test_dot_positions(self):
        f = self.TEST_FILES_DIR / "ICSD59959.cif"
        p = CifParser(f)
        s = p.get_structures()[0]
        self.assertEqual(s.formula, "K1 Mn1 F3")

    def test_replacing_finite_precision_frac_coords(self):
        f = self.TEST_FILES_DIR / "cif_finite_precision_frac_coord_error.cif"
        with warnings.catch_warnings():
            p = CifParser(f)
            s = p.get_structures()[0]
            self.assertEqual(str(s.composition), "N5+24")
            self.assertIn(
                "Some fractional co-ordinates rounded to ideal values to avoid issues with finite precision.",
                p.warnings,
            )

    def test_empty_deque(self):
        s = """data_1526655
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
        p = CifParser.from_string(s)
        self.assertEqual(p.get_structures()[0].formula, "Si1")
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
        p = CifParser.from_string(cif)
        self.assertRaises(ValueError, p.get_structures)


class MagCifTest(PymatgenTest):
    def setUp(self):
        warnings.filterwarnings("ignore")
        self.mcif = CifParser(self.TEST_FILES_DIR / "magnetic.example.NiO.mcif")
        self.mcif_ncl = CifParser(self.TEST_FILES_DIR / "magnetic.ncl.example.GdB4.mcif")
        self.mcif_incom = CifParser(self.TEST_FILES_DIR / "magnetic.incommensurate.example.Cr.mcif")
        self.mcif_disord = CifParser(self.TEST_FILES_DIR / "magnetic.disordered.example.CuMnO2.mcif")
        self.mcif_ncl2 = CifParser(self.TEST_FILES_DIR / "Mn3Ge_IR2.mcif")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_mcif_detection(self):
        self.assertTrue(self.mcif.feature_flags["magcif"])
        self.assertTrue(self.mcif_ncl.feature_flags["magcif"])
        self.assertTrue(self.mcif_incom.feature_flags["magcif"])
        self.assertTrue(self.mcif_disord.feature_flags["magcif"])
        self.assertFalse(self.mcif.feature_flags["magcif_incommensurate"])
        self.assertFalse(self.mcif_ncl.feature_flags["magcif_incommensurate"])
        self.assertTrue(self.mcif_incom.feature_flags["magcif_incommensurate"])
        self.assertFalse(self.mcif_disord.feature_flags["magcif_incommensurate"])

    def test_get_structures(self):
        # incommensurate structures not currently supported
        self.assertRaises(NotImplementedError, self.mcif_incom.get_structures)

        # disordered magnetic structures not currently supported
        self.assertRaises(NotImplementedError, self.mcif_disord.get_structures)

        # taken from self.mcif_ncl, removing explicit magnetic symmops
        # so that MagneticSymmetryGroup() has to be invoked
        magcifstr = """
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

        s = self.mcif.get_structures(primitive=False)[0]
        self.assertEqual(s.formula, "Ni32 O32")
        self.assertTrue(Magmom.are_collinear(s.site_properties["magmom"]))

        # example with non-collinear spin
        s_ncl = self.mcif_ncl.get_structures(primitive=False)[0]
        s_ncl_from_msg = CifParser.from_string(magcifstr).get_structures(primitive=False)[0]
        self.assertEqual(s_ncl.formula, "Gd4 B16")
        self.assertFalse(Magmom.are_collinear(s_ncl.site_properties["magmom"]))

        self.assertTrue(s_ncl.matches(s_ncl_from_msg))

    def test_write(self):
        cw_ref_string = """# generated using pymatgen
data_GdB4
_symmetry_space_group_name_H-M   'P 1'
_cell_length_a   7.13160000
_cell_length_b   7.13160000
_cell_length_c   4.05050000
_cell_angle_alpha   90.00000000
_cell_angle_beta   90.00000000
_cell_angle_gamma   90.00000000
_symmetry_Int_Tables_number   1
_chemical_formula_structural   GdB4
_chemical_formula_sum   'Gd4 B16'
_cell_volume   206.00729003
_cell_formula_units_Z   4
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
  Gd  Gd0  1  0.31746000  0.81746000  0.00000000  1.0
  Gd  Gd1  1  0.18254000  0.31746000  0.00000000  1.0
  Gd  Gd2  1  0.81746000  0.68254000  0.00000000  1.0
  Gd  Gd3  1  0.68254000  0.18254000  0.00000000  1.0
  B  B4  1  0.00000000  0.00000000  0.20290000  1.0
  B  B5  1  0.50000000  0.50000000  0.79710000  1.0
  B  B6  1  0.00000000  0.00000000  0.79710000  1.0
  B  B7  1  0.50000000  0.50000000  0.20290000  1.0
  B  B8  1  0.17590000  0.03800000  0.50000000  1.0
  B  B9  1  0.96200000  0.17590000  0.50000000  1.0
  B  B10  1  0.03800000  0.82410000  0.50000000  1.0
  B  B11  1  0.67590000  0.46200000  0.50000000  1.0
  B  B12  1  0.32410000  0.53800000  0.50000000  1.0
  B  B13  1  0.82410000  0.96200000  0.50000000  1.0
  B  B14  1  0.53800000  0.67590000  0.50000000  1.0
  B  B15  1  0.46200000  0.32410000  0.50000000  1.0
  B  B16  1  0.08670000  0.58670000  0.50000000  1.0
  B  B17  1  0.41330000  0.08670000  0.50000000  1.0
  B  B18  1  0.58670000  0.91330000  0.50000000  1.0
  B  B19  1  0.91330000  0.41330000  0.50000000  1.0
loop_
 _atom_site_moment_label
 _atom_site_moment_crystalaxis_x
 _atom_site_moment_crystalaxis_y
 _atom_site_moment_crystalaxis_z
  Gd0  5.05000000  5.05000000  0.00000000
  Gd1  -5.05000000  5.05000000  0.00000000
  Gd2  5.05000000  -5.05000000  0.00000000
  Gd3  -5.05000000  -5.05000000  0.00000000
"""
        s_ncl = self.mcif_ncl.get_structures(primitive=False)[0]

        cw = CifWriter(s_ncl, write_magmoms=True)
        self.assertEqual(cw.__str__(), cw_ref_string)

        # from list-type magmoms
        list_magmoms = [list(m) for m in s_ncl.site_properties["magmom"]]

        # float magmoms (magnitude only)
        float_magmoms = [float(m) for m in s_ncl.site_properties["magmom"]]

        s_ncl.add_site_property("magmom", list_magmoms)
        cw = CifWriter(s_ncl, write_magmoms=True)
        self.assertEqual(cw.__str__(), cw_ref_string)

        s_ncl.add_site_property("magmom", float_magmoms)
        cw = CifWriter(s_ncl, write_magmoms=True)

        cw_ref_string_magnitudes = """# generated using pymatgen
data_GdB4
_symmetry_space_group_name_H-M   'P 1'
_cell_length_a   7.13160000
_cell_length_b   7.13160000
_cell_length_c   4.05050000
_cell_angle_alpha   90.00000000
_cell_angle_beta   90.00000000
_cell_angle_gamma   90.00000000
_symmetry_Int_Tables_number   1
_chemical_formula_structural   GdB4
_chemical_formula_sum   'Gd4 B16'
_cell_volume   206.00729003
_cell_formula_units_Z   4
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
  Gd  Gd0  1  0.31746000  0.81746000  0.00000000  1.0
  Gd  Gd1  1  0.18254000  0.31746000  0.00000000  1.0
  Gd  Gd2  1  0.81746000  0.68254000  0.00000000  1.0
  Gd  Gd3  1  0.68254000  0.18254000  0.00000000  1.0
  B  B4  1  0.00000000  0.00000000  0.20290000  1.0
  B  B5  1  0.50000000  0.50000000  0.79710000  1.0
  B  B6  1  0.00000000  0.00000000  0.79710000  1.0
  B  B7  1  0.50000000  0.50000000  0.20290000  1.0
  B  B8  1  0.17590000  0.03800000  0.50000000  1.0
  B  B9  1  0.96200000  0.17590000  0.50000000  1.0
  B  B10  1  0.03800000  0.82410000  0.50000000  1.0
  B  B11  1  0.67590000  0.46200000  0.50000000  1.0
  B  B12  1  0.32410000  0.53800000  0.50000000  1.0
  B  B13  1  0.82410000  0.96200000  0.50000000  1.0
  B  B14  1  0.53800000  0.67590000  0.50000000  1.0
  B  B15  1  0.46200000  0.32410000  0.50000000  1.0
  B  B16  1  0.08670000  0.58670000  0.50000000  1.0
  B  B17  1  0.41330000  0.08670000  0.50000000  1.0
  B  B18  1  0.58670000  0.91330000  0.50000000  1.0
  B  B19  1  0.91330000  0.41330000  0.50000000  1.0
loop_
 _atom_site_moment_label
 _atom_site_moment_crystalaxis_x
 _atom_site_moment_crystalaxis_y
 _atom_site_moment_crystalaxis_z
  Gd0  0.00000000  0.00000000  7.14177849
  Gd1  0.00000000  0.00000000  7.14177849
  Gd2  0.00000000  0.00000000  -7.14177849
  Gd3  0.00000000  0.00000000  -7.14177849
"""
        self.assertEqual(cw.__str__().strip(), cw_ref_string_magnitudes.strip())
        # test we're getting correct magmoms in ncl case
        s_ncl2 = self.mcif_ncl2.get_structures()[0]
        list_magmoms = [list(m) for m in s_ncl2.site_properties["magmom"]]
        self.assertEqual(list_magmoms[0][0], 0.0)
        self.assertAlmostEqual(list_magmoms[0][1], 5.9160793408726366)
        self.assertAlmostEqual(list_magmoms[1][0], -5.1234749999999991)
        self.assertAlmostEqual(list_magmoms[1][1], 2.9580396704363183)

        # test creating an structure without oxidation state doesn't raise errors
        s_manual = Structure(Lattice.cubic(4.2), ["Cs", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        s_manual.add_spin_by_site([1, -1])
        cw = CifWriter(s_manual, write_magmoms=True)

        # check oxidation state
        cw_manual_oxi_string = """# generated using pymatgen
data_CsCl
_symmetry_space_group_name_H-M   'P 1'
_cell_length_a   4.20000000
_cell_length_b   4.20000000
_cell_length_c   4.20000000
_cell_angle_alpha   90.00000000
_cell_angle_beta   90.00000000
_cell_angle_gamma   90.00000000
_symmetry_Int_Tables_number   1
_chemical_formula_structural   CsCl
_chemical_formula_sum   'Cs1 Cl1'
_cell_volume   74.08800000
_cell_formula_units_Z   1
loop_
 _symmetry_equiv_pos_site_id
 _symmetry_equiv_pos_as_xyz
  1  'x, y, z'
loop_
 _atom_type_symbol
 _atom_type_oxidation_number
  Cs+  1.0
  Cl+  1.0
loop_
 _atom_site_type_symbol
 _atom_site_label
 _atom_site_symmetry_multiplicity
 _atom_site_fract_x
 _atom_site_fract_y
 _atom_site_fract_z
 _atom_site_occupancy
  Cs+  Cs0  1  0.00000000  0.00000000  0.00000000  1
  Cl+  Cl1  1  0.50000000  0.50000000  0.50000000  1
loop_
 _atom_site_moment_label
 _atom_site_moment_crystalaxis_x
 _atom_site_moment_crystalaxis_y
 _atom_site_moment_crystalaxis_z
"""
        s_manual.add_oxidation_state_by_site([1, 1])
        cw = CifWriter(s_manual, write_magmoms=True)
        self.assertEqual(cw.__str__(), cw_manual_oxi_string)

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
        self.assertEqual(self.mcif_ncl.get_bibtex_string(), ref_bibtex_string)


if __name__ == "__main__":
    unittest.main()
