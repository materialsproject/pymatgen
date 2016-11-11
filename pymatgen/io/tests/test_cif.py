# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals, division, print_function

import unittest2 as unittest
import os
import warnings

import numpy as np

from pymatgen.io.cif import CifParser, CifWriter, CifBlock
from pymatgen.io.vasp.inputs import Poscar
from pymatgen import Element, Specie, Lattice, Structure, Composition, DummySpecie
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.util.testing import PymatgenTest

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class CifBlockTest(unittest.TestCase):

    def test_to_string(self):
        with open(os.path.join(test_dir, 'Graphite.cif')) as f:
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
        for l1, l2, l3 in zip(str(c).split("\n"), cif_str.split("\n"),
                          cif_str_2.split("\n")):
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
        self.assertEqual(cb["_thing"], " long quotes  ;  still in the quote"
                                       "  ; actually going to end now")

    def test_long_loop(self):
        data = {'_stuff1': ['A' * 30] * 2,
                '_stuff2': ['B' * 30] * 2,
                '_stuff3': ['C' * 30] * 2}
        loops = [['_stuff1', '_stuff2', '_stuff3']]
        cif_str = """data_test
loop_
 _stuff1
 _stuff2
 _stuff3
  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA  BBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA  BBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"""
        self.assertEqual(str(CifBlock(data, loops, 'test')), cif_str)


class CifIOTest(PymatgenTest):

    def test_CifParser(self):
        parser = CifParser(os.path.join(test_dir, 'LiFePO4.cif'))
        for s in parser.get_structures(True):
            self.assertEqual(s.formula, "Li4 Fe4 P4 O16",
                             "Incorrectly parsed cif.")

        parser = CifParser(os.path.join(test_dir, 'V2O3.cif'))
        for s in parser.get_structures(True):
            self.assertEqual(s.formula, "V4 O6")

        parser = CifParser(os.path.join(test_dir, 'Li2O.cif'))
        prim = parser.get_structures(True)[0]
        self.assertEqual(prim.formula, "Li2 O1")
        conv = parser.get_structures(False)[0]
        self.assertEqual(conv.formula, "Li8 O4")

        #test for disordered structures
        parser = CifParser(os.path.join(test_dir, 'Li10GeP2S12.cif'))
        for s in parser.get_structures(True):
            self.assertEqual(s.formula, "Li20.2 Ge2.06 P3.94 S24",
                             "Incorrectly parsed cif.")
        cif_str = """#\#CIF1.1
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
        self.assertEqual(struct.formula, "Fe4 P4 O16",
                         "Incorrectly parsed cif.")
        self.assertAlmostEqual(struct.lattice.a, 10.4117668699)
        self.assertAlmostEqual(struct.lattice.b, 6.06717187997)
        self.assertAlmostEqual(struct.lattice.c, 4.75948953998)
        self.assertAlmostEqual(struct.lattice.alpha, 91)
        self.assertAlmostEqual(struct.lattice.beta, 92)
        self.assertAlmostEqual(struct.lattice.gamma, 93)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            parser = CifParser(os.path.join(test_dir, 'srycoo.cif'))
        self.assertEqual(parser.get_structures()[0].formula,
                         "Sr5.6 Y2.4 Co8 O21")

        # Test with a decimal Xyz. This should parse as two atoms in
        # conventional cell if it is correct, one if not.
        parser = CifParser(os.path.join(test_dir, "Fe.cif"))
        self.assertEqual(len(parser.get_structures(primitive=False)[0]), 2)

        # Below are 10 tests for CIFs from the Springer Materials/Pauling file DBs.

        # Partial occupancy on sites, incorrect label, previously unparsable
        parser = CifParser(os.path.join(test_dir, 'PF_sd_1928405.cif'))
        for s in parser.get_structures(True):
            self.assertEqual(s.formula, "Er1 Mn3.888 Fe2.112 Sn6")

        # Partial occupancy on sites, previously parsed as an ordered structure
        parser = CifParser(os.path.join(test_dir, 'PF_sd_1011081.cif'))
        for s in parser.get_structures(True):
            self.assertEqual(s.formula, "Zr0.2 Nb0.8")

        # Partial occupancy on sites, incorrect label, previously unparsable
        parser = CifParser(os.path.join(test_dir, 'PF_sd_1615854.cif'))
        for s in parser.get_structures(True):
            self.assertEqual(s.formula, "Na2 Al2 Si6 O16")

        # Partial occupancy on sites, incorrect label, previously unparsable
        parser = CifParser(os.path.join(test_dir, 'PF_sd_1622133.cif'))
        for s in parser.get_structures(True):
            self.assertEqual(s.formula, "Ca0.184 Mg13.016 Fe2.8 Si16 O48")

        # Partial occupancy on sites, previously parsed as an ordered structure
        parser = CifParser(os.path.join(test_dir, 'PF_sd_1908491.cif'))
        for s in parser.get_structures(True):
            self.assertEqual(s.formula, "Mn0.48 Zn0.52 Ga2 Se4")

        # Partial occupancy on sites, incorrect label, previously unparsable
        parser = CifParser(os.path.join(test_dir, 'PF_sd_1811457.cif'))
        for s in parser.get_structures(True):
            self.assertEqual(s.formula, "Ba2 Mg0.6 Zr0.2 Ta1.2 O6")

        # Incomplete powder diffraction data, previously unparsable
        # This CIF file contains the molecular species "NH3" which is
        # parsed as "N" because the label is "N{x}" (x = 1,2,..) and the
        # corresponding symbol is "NH3". Since, the label and symbol are switched
        # in CIFs from Springer Materials/Pauling file DBs, CifParser parses the
        # element as "N".
        parser = CifParser(os.path.join(test_dir, 'PF_sd_1002871.cif'))
        self.assertEqual(parser.get_structures(True)[0].formula, "Cu1 Br2 N6")
        self.assertEqual(parser.get_structures(True)[1].formula, "Cu1 Br4 N6")

        # Incomplete powder diffraction data, previously unparsable
        parser = CifParser(os.path.join(test_dir, 'PF_sd_1704003.cif'))
        for s in parser.get_structures():
            self.assertEqual(s.formula, "Rb4 Mn2 F12")

        # Unparsable species 'OH/OH2', previously parsed as "O"
        parser = CifParser(os.path.join(test_dir, 'PF_sd_1500382.cif'))
        for s in parser.get_structures():
            self.assertEqual(s.formula, "Mg6 B2 O6 F1.764")

        # Unparsable species 'OH/OH2', previously parsed as "O"
        parser = CifParser(os.path.join(test_dir, 'PF_sd_1601634.cif'))
        for s in parser.get_structures():
            self.assertEqual(s.formula, "Zn1.29 Fe0.69 As2 Pb1.02 O8")

    def test_CifWriter(self):
        filepath = os.path.join(test_dir, 'POSCAR')
        poscar = Poscar.from_file(filepath)
        writer = CifWriter(poscar.structure, symprec=0.001)
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
  Fe  Fe1  4  0.218728  0.250000  0.525133  1
  P  P2  4  0.094613  0.750000  0.581757  1
  O  O3  8  0.165710  0.546072  0.714616  1
  O  O4  4  0.043372  0.250000  0.292862  1
  O  O5  4  0.096642  0.750000  0.258680  1"""
        for l1, l2 in zip(str(writer).split("\n"), ans.split("\n")):
            self.assertEqual(l1.strip(), l2.strip())

    def test_symmetrized(self):
        filepath = os.path.join(test_dir, 'POSCAR')
        poscar = Poscar.from_file(filepath)
        writer = CifWriter(poscar.structure, symprec=0.1)
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
  Fe  Fe1  4  0.218728  0.250000  0.525133  1
  P  P2  4  0.094613  0.750000  0.581757  1
  O  O3  8  0.165710  0.546072  0.714616  1
  O  O4  4  0.043372  0.250000  0.292862  1
  O  O5  4  0.096642  0.750000  0.258680  1"""
        for l1, l2 in zip(str(writer).split("\n"), ans.split("\n")):
            self.assertEqual(l1.strip(), l2.strip())

        ans = """# generated using pymatgen
data_LiFePO4
_symmetry_space_group_name_H-M   Pnma
_cell_length_a   10.41037000
_cell_length_b   6.06577000
_cell_length_c   4.74480000
_cell_angle_alpha   90.00000000
_cell_angle_beta   90.00000000
_cell_angle_gamma   90.00000000
_symmetry_Int_Tables_number   62
_chemical_formula_structural   LiFePO4
_chemical_formula_sum   'Li4 Fe4 P4 O16'
_cell_volume   299.619458734
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
  Li  Li1  4  0.000000  0.000000  0.000000  1.0
  Fe  Fe2  4  0.218845  0.750000  0.474910  1.0
  P  P3  4  0.094445  0.250000  0.417920  1.0
  O  O4  8  0.165815  0.044060  0.286540  1.0
  O  O5  4  0.043155  0.750000  0.708460  1.0
  O  O6  4  0.096215  0.250000  0.741480  1.0
"""
        s = Structure.from_file(os.path.join(test_dir, 'LiFePO4.cif'))
        writer = CifWriter(s, symprec=0.1)
        s2 = CifParser.from_string(str(writer)).get_structures()[0]
        m = StructureMatcher()
        self.assertTrue(m.fit(s, s2))

        s = self.get_structure("Li2O")
        writer = CifWriter(s, symprec=0.1)
        ans = """# generated using pymatgen
data_Li2O
_symmetry_space_group_name_H-M   Fm-3m
_cell_length_a   4.61000000
_cell_length_b   4.61000000
_cell_length_c   4.61000000
_cell_angle_alpha   90.00000000
_cell_angle_beta   90.00000000
_cell_angle_gamma   90.00000000
_symmetry_Int_Tables_number   225
_chemical_formula_structural   Li2O
_chemical_formula_sum   'Li8 O4'
_cell_volume   97.972181
_cell_formula_units_Z   4
loop_
 _symmetry_equiv_pos_site_id
 _symmetry_equiv_pos_as_xyz
  1  'x, y, z'
  2  '-x, -y, -z'
  3  'z, y, -x'
  4  '-z, -y, x'
  5  '-x, y, -z'
  6  'x, -y, z'
  7  '-z, y, x'
  8  'z, -y, -x'
  9  'x, -y, -z'
  10  '-x, y, z'
  11  'z, -y, x'
  12  '-z, y, -x'
  13  '-x, -y, z'
  14  'x, y, -z'
  15  '-z, -y, -x'
  16  'z, y, x'
  17  'y, -z, -x'
  18  '-y, z, x'
  19  'y, x, -z'
  20  '-y, -x, z'
  21  'y, z, x'
  22  '-y, -z, -x'
  23  'y, -x, z'
  24  '-y, x, -z'
  25  '-y, z, -x'
  26  'y, -z, x'
  27  '-y, -x, -z'
  28  'y, x, z'
  29  '-y, -z, x'
  30  'y, z, -x'
  31  '-y, x, z'
  32  'y, -x, -z'
  33  '-z, x, -y'
  34  'z, -x, y'
  35  'x, z, -y'
  36  '-x, -z, y'
  37  'z, -x, -y'
  38  '-z, x, y'
  39  '-x, -z, -y'
  40  'x, z, y'
  41  'z, x, y'
  42  '-z, -x, -y'
  43  '-x, z, y'
  44  'x, -z, -y'
  45  '-z, -x, y'
  46  'z, x, -y'
  47  'x, -z, y'
  48  '-x, z, -y'
  49  'x+1/2, y+1/2, z'
  50  '-x+1/2, -y+1/2, -z'
  51  'z+1/2, y+1/2, -x'
  52  '-z+1/2, -y+1/2, x'
  53  '-x+1/2, y+1/2, -z'
  54  'x+1/2, -y+1/2, z'
  55  '-z+1/2, y+1/2, x'
  56  'z+1/2, -y+1/2, -x'
  57  'x+1/2, -y+1/2, -z'
  58  '-x+1/2, y+1/2, z'
  59  'z+1/2, -y+1/2, x'
  60  '-z+1/2, y+1/2, -x'
  61  '-x+1/2, -y+1/2, z'
  62  'x+1/2, y+1/2, -z'
  63  '-z+1/2, -y+1/2, -x'
  64  'z+1/2, y+1/2, x'
  65  'y+1/2, -z+1/2, -x'
  66  '-y+1/2, z+1/2, x'
  67  'y+1/2, x+1/2, -z'
  68  '-y+1/2, -x+1/2, z'
  69  'y+1/2, z+1/2, x'
  70  '-y+1/2, -z+1/2, -x'
  71  'y+1/2, -x+1/2, z'
  72  '-y+1/2, x+1/2, -z'
  73  '-y+1/2, z+1/2, -x'
  74  'y+1/2, -z+1/2, x'
  75  '-y+1/2, -x+1/2, -z'
  76  'y+1/2, x+1/2, z'
  77  '-y+1/2, -z+1/2, x'
  78  'y+1/2, z+1/2, -x'
  79  '-y+1/2, x+1/2, z'
  80  'y+1/2, -x+1/2, -z'
  81  '-z+1/2, x+1/2, -y'
  82  'z+1/2, -x+1/2, y'
  83  'x+1/2, z+1/2, -y'
  84  '-x+1/2, -z+1/2, y'
  85  'z+1/2, -x+1/2, -y'
  86  '-z+1/2, x+1/2, y'
  87  '-x+1/2, -z+1/2, -y'
  88  'x+1/2, z+1/2, y'
  89  'z+1/2, x+1/2, y'
  90  '-z+1/2, -x+1/2, -y'
  91  '-x+1/2, z+1/2, y'
  92  'x+1/2, -z+1/2, -y'
  93  '-z+1/2, -x+1/2, y'
  94  'z+1/2, x+1/2, -y'
  95  'x+1/2, -z+1/2, y'
  96  '-x+1/2, z+1/2, -y'
  97  'x+1/2, y, z+1/2'
  98  '-x+1/2, -y, -z+1/2'
  99  'z+1/2, y, -x+1/2'
  100  '-z+1/2, -y, x+1/2'
  101  '-x+1/2, y, -z+1/2'
  102  'x+1/2, -y, z+1/2'
  103  '-z+1/2, y, x+1/2'
  104  'z+1/2, -y, -x+1/2'
  105  'x+1/2, -y, -z+1/2'
  106  '-x+1/2, y, z+1/2'
  107  'z+1/2, -y, x+1/2'
  108  '-z+1/2, y, -x+1/2'
  109  '-x+1/2, -y, z+1/2'
  110  'x+1/2, y, -z+1/2'
  111  '-z+1/2, -y, -x+1/2'
  112  'z+1/2, y, x+1/2'
  113  'y+1/2, -z, -x+1/2'
  114  '-y+1/2, z, x+1/2'
  115  'y+1/2, x, -z+1/2'
  116  '-y+1/2, -x, z+1/2'
  117  'y+1/2, z, x+1/2'
  118  '-y+1/2, -z, -x+1/2'
  119  'y+1/2, -x, z+1/2'
  120  '-y+1/2, x, -z+1/2'
  121  '-y+1/2, z, -x+1/2'
  122  'y+1/2, -z, x+1/2'
  123  '-y+1/2, -x, -z+1/2'
  124  'y+1/2, x, z+1/2'
  125  '-y+1/2, -z, x+1/2'
  126  'y+1/2, z, -x+1/2'
  127  '-y+1/2, x, z+1/2'
  128  'y+1/2, -x, -z+1/2'
  129  '-z+1/2, x, -y+1/2'
  130  'z+1/2, -x, y+1/2'
  131  'x+1/2, z, -y+1/2'
  132  '-x+1/2, -z, y+1/2'
  133  'z+1/2, -x, -y+1/2'
  134  '-z+1/2, x, y+1/2'
  135  '-x+1/2, -z, -y+1/2'
  136  'x+1/2, z, y+1/2'
  137  'z+1/2, x, y+1/2'
  138  '-z+1/2, -x, -y+1/2'
  139  '-x+1/2, z, y+1/2'
  140  'x+1/2, -z, -y+1/2'
  141  '-z+1/2, -x, y+1/2'
  142  'z+1/2, x, -y+1/2'
  143  'x+1/2, -z, y+1/2'
  144  '-x+1/2, z, -y+1/2'
  145  'x, y+1/2, z+1/2'
  146  '-x, -y+1/2, -z+1/2'
  147  'z, y+1/2, -x+1/2'
  148  '-z, -y+1/2, x+1/2'
  149  '-x, y+1/2, -z+1/2'
  150  'x, -y+1/2, z+1/2'
  151  '-z, y+1/2, x+1/2'
  152  'z, -y+1/2, -x+1/2'
  153  'x, -y+1/2, -z+1/2'
  154  '-x, y+1/2, z+1/2'
  155  'z, -y+1/2, x+1/2'
  156  '-z, y+1/2, -x+1/2'
  157  '-x, -y+1/2, z+1/2'
  158  'x, y+1/2, -z+1/2'
  159  '-z, -y+1/2, -x+1/2'
  160  'z, y+1/2, x+1/2'
  161  'y, -z+1/2, -x+1/2'
  162  '-y, z+1/2, x+1/2'
  163  'y, x+1/2, -z+1/2'
  164  '-y, -x+1/2, z+1/2'
  165  'y, z+1/2, x+1/2'
  166  '-y, -z+1/2, -x+1/2'
  167  'y, -x+1/2, z+1/2'
  168  '-y, x+1/2, -z+1/2'
  169  '-y, z+1/2, -x+1/2'
  170  'y, -z+1/2, x+1/2'
  171  '-y, -x+1/2, -z+1/2'
  172  'y, x+1/2, z+1/2'
  173  '-y, -z+1/2, x+1/2'
  174  'y, z+1/2, -x+1/2'
  175  '-y, x+1/2, z+1/2'
  176  'y, -x+1/2, -z+1/2'
  177  '-z, x+1/2, -y+1/2'
  178  'z, -x+1/2, y+1/2'
  179  'x, z+1/2, -y+1/2'
  180  '-x, -z+1/2, y+1/2'
  181  'z, -x+1/2, -y+1/2'
  182  '-z, x+1/2, y+1/2'
  183  '-x, -z+1/2, -y+1/2'
  184  'x, z+1/2, y+1/2'
  185  'z, x+1/2, y+1/2'
  186  '-z, -x+1/2, -y+1/2'
  187  '-x, z+1/2, y+1/2'
  188  'x, -z+1/2, -y+1/2'
  189  '-z, -x+1/2, y+1/2'
  190  'z, x+1/2, -y+1/2'
  191  'x, -z+1/2, y+1/2'
  192  '-x, z+1/2, -y+1/2'
loop_
 _atom_type_symbol
 _atom_type_oxidation_number
  Li+  1.0
  O2-  -2.0
loop_
 _atom_site_type_symbol
 _atom_site_label
 _atom_site_symmetry_multiplicity
 _atom_site_fract_x
 _atom_site_fract_y
 _atom_site_fract_z
 _atom_site_occupancy
  Li+  Li1  8  0.250000  0.250000  0.250000  1.0
  O2-  O2  4  0.000000  0.000000  0.000000  1.0"""

        for l1, l2 in zip(str(writer).split("\n"), ans.split("\n")):
            self.assertEqual(l1.strip(), l2.strip())

    def test_disordered(self):
        si = Element("Si")
        n = Element("N")
        coords = list()
        coords.append(np.array([0, 0, 0]))
        coords.append(np.array([0.75, 0.5, 0.75]))
        lattice = Lattice(np.array([[3.8401979337, 0.00, 0.00],
                                    [1.9200989668, 3.3257101909, 0.00],
                                    [0.00, -2.2171384943, 3.1355090603]]))
        struct = Structure(lattice, [si, {si:0.5, n:0.5}], coords)
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
_cell_volume   40.0447946443
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
Si  Si1  1  0.000000  0.000000  0.000000  1
Si  Si2  1  0.750000  0.500000  0.750000  0.5
N  N3  1  0.750000  0.500000  0.750000  0.5

    """
        for l1, l2 in zip(str(writer).split("\n"), ans.split("\n")):
            self.assertEqual(l1.strip(), l2.strip())

    def test_specie_cifwriter(self):
        si4 = Specie("Si", 4)
        si3 = Specie("Si", 3)
        n = DummySpecie("X", -3)
        coords = list()
        coords.append(np.array([0.5, 0.5, 0.5]))
        coords.append(np.array([0.75, 0.5, 0.75]))
        coords.append(np.array([0, 0, 0]))
        lattice = Lattice(np.array([[3.8401979337, 0.00, 0.00],
                                    [1.9200989668, 3.3257101909, 0.00],
                                    [0.00, -2.2171384943, 3.1355090603]]))
        struct = Structure(lattice, [n, {si3:0.5, n:0.5}, si4], coords)
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
_cell_volume   40.0447946443
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
  X3-  X1  1  0.500000  0.500000  0.500000  1
  X3-  X2  1  0.750000  0.500000  0.750000  0.5
  Si3+  Si3  1  0.750000  0.500000  0.750000  0.5
  Si4+  Si4  1  0.000000  0.000000  0.000000  1

"""
        for l1, l2 in zip(str(writer).split("\n"), ans.split("\n")):
            self.assertEqual(l1.strip(), l2.strip())

    def test_primes(self):
        parser = CifParser(os.path.join(test_dir, 'C26H16BeN2O2S2.cif'))
        for s in parser.get_structures(False):
            self.assertEqual(s.composition, 8 * Composition('C26H16BeN2O2S2'))

    def test_missing_atom_site_type_with_oxistates(self):
        parser = CifParser(os.path.join(test_dir, 'P24Ru4H252C296S24N16.cif'))
        c = Composition({'S0+': 24, 'Ru0+': 4, 'H0+': 252, 'C0+': 296,
                         'N0+': 16, 'P0+': 24})
        for s in parser.get_structures(False):
            self.assertEqual(s.composition, c)

    def test_no_coords_or_species(self):
        string=  """#generated using pymatgen
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
        filepath = os.path.join(test_dir, 'POSCAR')
        poscar = Poscar.from_file(filepath)
        s_ref = poscar.structure

        sm = StructureMatcher(stol=0.05, ltol=0.01, angle_tol=0.1)
        self.assertTrue(sm.fit(s_ref, s_test))

    def test_empty(self):
        # single line
        cb = CifBlock.from_string("data_mwe\nloop_\n_tag\n ''")
        self.assertEqual(cb.data['_tag'][0], '')

        # multi line
        cb = CifBlock.from_string("data_mwe\nloop_\n_tag\n;\n;")
        self.assertEqual(cb.data['_tag'][0], '')

        cb2 = CifBlock.from_string(str(cb))
        self.assertEqual(cb, cb2)

    def test_bad_cif(self):
        f = os.path.join(test_dir, "bad_occu.cif")
        p = CifParser(f)
        self.assertRaises(ValueError, p.get_structures)
        p = CifParser(f, occupancy_tolerance=2)
        s = p.get_structures()[0]
        self.assertAlmostEqual(s[0].species_and_occu["Al3+"], 0.5)

    def test_one_line_symm(self):
        f = os.path.join(test_dir, "OneLineSymmP1.cif")
        p = CifParser(f)
        s = p.get_structures()[0]
        self.assertEqual(s.formula, "Ga4 Pb2 O8")

    def test_no_symmops(self):
        f = os.path.join(test_dir, "nosymm.cif")
        p = CifParser(f)
        s = p.get_structures()[0]
        self.assertEqual(s.formula, "H96 C60 O8")

    def test_dot_positions(self):
        f = os.path.join(test_dir, "ICSD59959.cif")
        p = CifParser(f)
        s = p.get_structures()[0]
        self.assertEqual(s.formula, "K1 Mn1 F3")


if __name__ == '__main__':
    unittest.main()
