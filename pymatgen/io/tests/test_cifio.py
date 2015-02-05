# coding: utf-8

from __future__ import unicode_literals

import unittest
import os
import warnings

import numpy as np

from pymatgen.io.cifio import CifParser, CifWriter, CifBlock
from pymatgen.io.vaspio.vasp_input import Poscar
from pymatgen import Element, Specie, Lattice, Structure, Composition

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')

class CifBlockTest(unittest.TestCase):

    def test_to_string(self):
        with open(os.path.join(test_dir, 'Graphite.cif')) as f:
            s = f.read()
        c = CifBlock.from_string(s)
        cif_str_2 = str(CifBlock.from_string(str(c)))
        cif_str = """data_28417-ICSD
_database_code_ICSD   28417
_audit_creation_date   1980-01-01
_audit_update_record   2003-04-01
_chemical_name_systematic   Carbon
_chemical_formula_structural   C
_chemical_formula_sum   C1
_chemical_name_mineral   'Graphite, nitrated'
_exptl_crystal_density_diffrn   1.36
_publ_section_title
'Order-disorder transformations in graphite nitrates'
loop_
  _citation_id
  _citation_journal_full
  _citation_year
  _citation_journal_volume
  _citation_page_first
  _citation_page_last
  _citation_journal_id_ASTM
   primary
;
Proceedings of the Royal Society of London, Series A: Mathematical and
Physical Sciences (76,1906-)
;
   1966  291  324  339  PRLAAZ
loop_
  _publ_author_name
   'Nixon, D.E.'
   'Parry, G.S.'
   'Ubbelohde, A.R.'
_cell_length_a   2.46
_cell_length_b   2.46
_cell_length_c   33.45
_cell_angle_alpha   90.
_cell_angle_beta   90.
_cell_angle_gamma   120.
_cell_volume   175.31
_cell_formula_units_Z   12
_symmetry_space_group_name_H-M   'R -3 m H'
_symmetry_Int_Tables_number   166
loop_
  _symmetry_equiv_pos_site_id
  _symmetry_equiv_pos_as_xyz
   1  'x-y, -y, -z'
   2  '-x, -x+y, -z'
   3  'y, x, -z'
   4  'x-y, x, -z'
   5  'y, -x+y, -z'
   6  '-x, -y, -z'
   7  '-x+y, y, z'
   8  'x, x-y, z'
   9  '-y, -x, z'
   10  '-x+y, -x, z'
   11  '-y, x-y, z'
   12  'x, y, z'
   13  'x-y+2/3, -y+1/3, -z+1/3'
   14  '-x+2/3, -x+y+1/3, -z+1/3'
   15  'y+2/3, x+1/3, -z+1/3'
   16  'x-y+2/3, x+1/3, -z+1/3'
   17  'y+2/3, -x+y+1/3, -z+1/3'
   18  '-x+2/3, -y+1/3, -z+1/3'
   19  '-x+y+2/3, y+1/3, z+1/3'
   20  'x+2/3, x-y+1/3, z+1/3'
   21  '-y+2/3, -x+1/3, z+1/3'
   22  '-x+y+2/3, -x+1/3, z+1/3'
   23  '-y+2/3, x-y+1/3, z+1/3'
   24  'x+2/3, y+1/3, z+1/3'
   25  'x-y+1/3, -y+2/3, -z+2/3'
   26  '-x+1/3, -x+y+2/3, -z+2/3'
   27  'y+1/3, x+2/3, -z+2/3'
   28  'x-y+1/3, x+2/3, -z+2/3'
   29  'y+1/3, -x+y+2/3, -z+2/3'
   30  '-x+1/3, -y+2/3, -z+2/3'
   31  '-x+y+1/3, y+2/3, z+2/3'
   32  'x+1/3, x-y+2/3, z+2/3'
   33  '-y+1/3, -x+2/3, z+2/3'
   34  '-x+y+1/3, -x+2/3, z+2/3'
   35  '-y+1/3, x-y+2/3, z+2/3'
   36  'x+1/3, y+2/3, z+2/3'
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
   C1  C0+  6  c  0  0  0.05  0.  1.  0
   C2  C0+  6  c  0  0  0.283  0.  1.  0"""
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


class CifIOTest(unittest.TestCase):

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
90.0
_cell_angle_beta                        90.0
_cell_angle_gamma                       90.0
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
        self.assertEqual(parser.get_structures()[0].formula, "Fe4 P4 O16",
                         "Incorrectly parsed cif.")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            parser = CifParser(os.path.join(test_dir, 'srycoo.cif'))
        self.assertEqual(parser.get_structures()[0].formula,
                         "Sr5.6 Y2.4 Co8 O21")


    def test_CifWriter(self):
        filepath = os.path.join(test_dir, 'POSCAR')
        poscar = Poscar.from_file(filepath)
        writer = CifWriter(poscar.structure, find_spacegroup=True)
        ans = """#generated using pymatgen
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
        for l1, l2 in zip(str(writer).split("\n"), ans.split("\n")):
            self.assertEqual(l1.strip(), l2.strip())

    def test_symmetrized(self):
        filepath = os.path.join(test_dir, 'POSCAR')
        poscar = Poscar.from_file(filepath)
        writer = CifWriter(poscar.structure, find_spacegroup=True,
                           symprec=0.1)
        ans = """#generated using pymatgen
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
  Fe  Fe1  4  0.218728  0.750000  0.474867  1
  P  P2  4  0.094613  0.250000  0.418243  1
  O  O3  4  0.043372  0.750000  0.707138  1
  O  O4  4  0.096642  0.250000  0.741320  1
  O  O5  8  0.165710  0.046072  0.285384  1"""
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
        ans = """#generated using pymatgen
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
        n = Specie("N", -3)
        coords = list()
        coords.append(np.array([0, 0, 0]))
        coords.append(np.array([0.75, 0.5, 0.75]))
        coords.append(np.array([0.5, 0.5, 0.5]))
        lattice = Lattice(np.array([[3.8401979337, 0.00, 0.00],
                                    [1.9200989668, 3.3257101909, 0.00],
                                    [0.00, -2.2171384943, 3.1355090603]]))
        struct = Structure(lattice, [si4, {si3:0.5, n:0.5}, n], coords)
        writer = CifWriter(struct)
        ans = """#generated using pymatgen
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
  Si4+  Si1  1  0.000000  0.000000  0.000000  1
  Si3+  Si2  1  0.750000  0.500000  0.750000  0.5
  N3-  N3  1  0.750000  0.500000  0.750000  0.5
  N3-  N4  1  0.500000  0.500000  0.500000  1
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
        self.assertEqual(len(parser.get_structures()), 0)

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
        from pymatgen.analysis.structure_matcher import StructureMatcher
        sm = StructureMatcher(stol=0.05, ltol=0.01, angle_tol=0.1)
        self.assertTrue(sm.fit(s_ref, s_test))

if __name__ == '__main__':
    unittest.main()
