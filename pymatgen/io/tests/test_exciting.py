# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import unittest
import os
import xml.etree.cElementTree as ET
from pymatgen.io.exciting import ExcitingInput
from pymatgen.core import Structure, Lattice

"""
Tests the module to parse exciting input
Created 1st of December, 2016
"""

__author__ = "Christian Vorwerk"
__copyright__ = "Copyright 2016"
__version__ = "0.1"
__maintainer__ = "Christian Vorwerk"
__email__ = "vorwerk@physik.hu-berlin.de"
__date__ = "Dec 01, 2016"

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class ExcitingInputTest(unittest.TestCase):
    def test_fromfile(self):
        # Test for the import of a structure directly from an exciting
        # input file
        filepath = os.path.join(test_dir, 'input1.xml')
        excin = ExcitingInput.from_file(filepath)
        lattice = [[0.0, 2.81, 2.81], [2.81, 0.0, 2.81], [2.81, 2.81, 0.0]]
        atoms = ['Na', 'Cl']
        fraccoords = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
        self.assertEqual(lattice, excin.structure.lattice.matrix.tolist())
        self.assertEqual(atoms,
                         [site.specie.symbol for site in excin.structure])
        self.assertEqual(fraccoords, [site.frac_coords.tolist() for site in
                                      excin.structure])

    def test_writestring(self):
        # Test for the string export of s atructure into the exciting input xml schema
        input_string = (
            '<input xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
            'xsi:noNamespaceSchemaLocation="http://xml.exciting-code.org/excitinginput'
            '.xsd">\n  <title>Na4 Cl4</title>\n  <structure speciespath="./">\n    '
            '<crystal scale="1.8897543768634038">\n      <basevect>      5.62000000'
            '       0.00000000       0.00000000</basevect>\n      <basevect>      '
            '0.00000000       5.62000000       0.00000000</basevect>\n      '
            '<basevect>      0.00000000       0.00000000       5.62000000</basevect>'
            '\n    </crystal>\n    <species speciesfile="Na.xml">\n      <atom coord='
            '"      0.00000000       0.00000000       0.00000000" />\n      <atom coor'
            'd="      0.50000000       0.50000000       0.00000000" />\n      <atom co'
            'ord="      0.50000000       0.00000000       0.50000000" />\n      <atom '
            'coord="      0.00000000       0.50000000       0.50000000" />\n    </spec'
            'ies>\n    <species speciesfile="Cl.xml">\n      <atom coord="      0.5000'
            '0000       0.00000000       0.00000000" />\n      <atom coord="      0.00'
            '000000       0.50000000       0.00000000" />\n      <atom coord="      0.'
            '00000000       0.00000000       0.50000000" />\n      <atom coord="      '
            '0.50000000       0.50000000       0.50000000" />\n    </species>\n  </str'
            'ucture>\n</input>\n')
        lattice = Lattice.cubic('5.62')
        structure = Structure(lattice, ['Na', 'Na', 'Na', 'Na',
                                        'Cl', 'Cl', 'Cl', 'Cl'],
                              [[0, 0, 0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5],
                               [0.0, 0.5, 0.5],
                               [0.5, 0.0, 0.0], [0.0, 0.5, 0.0],
                               [0.0, 0.0, 0.5], [0.5, 0.5, 0.5]])
        excin = ExcitingInput(structure)
        for l1, l2 in zip(input_string.split("\n"),
                          excin.write_string('unchanged').split("\n")):
            if not l1.strip().startswith("<crystal scale"):
                self.assertEqual(l1, l2)

    def test_writebandstr(self):
        filepath = os.path.join(test_dir, 'CsI3Pb.cif')
        structure = Structure.from_file(filepath)
        excin = ExcitingInput(structure)
        string = excin.write_string('primitive', bandstr=True)
        bandstr = string.split('<properties>')[1].split('</properties>')[0]

        coord = []
        label = []
        coord_ref = [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [0.5, 0.5, 0.0],
                     [0.0, 0.5, 0.0],
                     [0.0, 0.0, 0.0], [0.0, 0.0, 0.5], [0.5, 0.0, 0.5],
                     [0.5, 0.5, 0.5],
                     [0.0, 0.5, 0.5], [0.0, 0.0, 0.5], [0.0, 0.5, 0.0],
                     [0.0, 0.5, 0.5],
                     [0.5, 0.0, 0.5], [0.5, 0.0, 0.0], [0.5, 0.5, 0.0],
                     [0.5, 0.5, 0.5]]
        label_ref = ['GAMMA', 'X', 'S', 'Y', 'GAMMA', 'Z', 'U', 'R', 'T', 'Z',
                     'Y', 'T',
                     'U', 'X', 'S', 'R']
        root = ET.fromstring(bandstr)
        for plot1d in root.iter('plot1d'):
            for point in plot1d.iter('point'):
                coord.append([float(i) for i in point.get('coord').split()])
                label.append(point.get('label'))
        self.assertEqual(label, label_ref)
        self.assertEqual(coord, coord_ref)


if __name__ == "__main__":
    unittest.main()
