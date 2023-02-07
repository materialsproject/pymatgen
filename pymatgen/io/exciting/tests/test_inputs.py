# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import annotations

import os
import unittest
import xml.etree.ElementTree as ET

from pymatgen.core import Lattice, Structure
from pymatgen.io.exciting import ExcitingInput
from pymatgen.util.testing import PymatgenTest

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


class ExcitingInputTest(PymatgenTest):
    def test_fromfile(self):
        # Test for the import of a structure directly from an exciting
        # input file
        filepath = os.path.join(PymatgenTest.TEST_FILES_DIR, "input_exciting1.xml")
        excin = ExcitingInput.from_file(filepath)
        lattice = [[0.0, 2.81, 2.81], [2.81, 0.0, 2.81], [2.81, 2.81, 0.0]]
        atoms = ["Na", "Cl"]
        fraccoords = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
        self.assertArrayAlmostEqual(lattice, excin.structure.lattice.matrix.tolist())
        assert atoms == [site.specie.symbol for site in excin.structure]
        assert fraccoords == [site.frac_coords.tolist() for site in excin.structure]

    def test_writestring(self):
        # Test for the string export of s atructure into the exciting input xml schema
        input_string = (
            '<input xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
            'xsi:noNamespaceSchemaLocation="http://xml.exciting-code.org/excitinginput'
            '.xsd">\n  <title>Na4 Cl4</title>\n  <structure speciespath="./">\n    '
            '<crystal scale="1.8897543768634038">\n      <basevect>      5.62000000'
            "       0.00000000       0.00000000</basevect>\n      <basevect>      "
            "0.00000000       5.62000000       0.00000000</basevect>\n      "
            "<basevect>      0.00000000       0.00000000       5.62000000</basevect>"
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
            "ucture>\n</input>\n"
        )
        lattice = Lattice.cubic("5.62")
        structure = Structure(
            lattice,
            ["Na", "Na", "Na", "Na", "Cl", "Cl", "Cl", "Cl"],
            [
                [0, 0, 0],
                [0.5, 0.5, 0.0],
                [0.5, 0.0, 0.5],
                [0.0, 0.5, 0.5],
                [0.5, 0.0, 0.0],
                [0.0, 0.5, 0.0],
                [0.0, 0.0, 0.5],
                [0.5, 0.5, 0.5],
            ],
        )
        excin = ExcitingInput(structure)
        for l1, l2 in zip(input_string.split("\n"), excin.write_string("unchanged").split("\n")):
            if not l1.strip().startswith("<crystal scale"):
                assert l1.strip() == l2.strip()

    def test_writebandstr(self):
        filepath = os.path.join(PymatgenTest.TEST_FILES_DIR, "CsI3Pb.cif")
        structure = Structure.from_file(filepath)
        excin = ExcitingInput(structure)
        string = excin.write_string("primitive", bandstr=True)
        bandstr = string.split("<properties>")[1].split("</properties>")[0]

        coord = []
        label = []
        coord_ref = [
            [0.0, 0.0, 0.0],
            [0.5, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.0, 0.5, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.5],
            [0.5, 0.0, 0.5],
            [0.5, 0.5, 0.5],
            [0.0, 0.5, 0.5],
            [0.0, 0.0, 0.5],
            [0.0, 0.5, 0.0],
            [0.0, 0.5, 0.5],
            [0.5, 0.0, 0.5],
            [0.5, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.5, 0.5, 0.5],
        ]
        label_ref = [
            "GAMMA",
            "X",
            "S",
            "Y",
            "GAMMA",
            "Z",
            "U",
            "R",
            "T",
            "Z",
            "Y",
            "T",
            "U",
            "X",
            "S",
            "R",
        ]
        root = ET.fromstring(bandstr)
        for plot1d in root.iter("plot1d"):
            for point in plot1d.iter("point"):
                coord.append([float(i) for i in point.get("coord").split()])
                label.append(point.get("label"))
        assert label == label_ref
        assert coord == coord_ref

    def test_paramdict(self):
        coords = [[0.0, 0.0, 0.0], [0.75, 0.5, 0.75]]
        lattice = Lattice.from_parameters(a=3.84, b=3.84, c=3.84, alpha=120, beta=90, gamma=60)
        struct = Structure(lattice, ["Si", "Si"], coords)
        paradir = {
            "grst": {
                "do": "fromscratch",
                "ngridk": "8 8 8",
                "xctype": "GGA_PBE_SOL",
                "gmaxvr": "14.0",
            },
            "xs": {
                "xstype": "BSE",
                "ngridk": "4 4 4",
                "ngridq": "4 4 4",
                "nempty": "30",
                "gqmax": "3.0",
                "broad": "0.07",
                "tevout": "true",
                "energywindow": {"intv": "0.0 1.0", "points": "1200"},
                "screening": {"screentype": "full", "nempty": "100"},
                "BSE": {"bsetype": "singlet", "nstlbse": "1 5 1 4"},
            },
        }

        test_input = ExcitingInput(struct)
        test_string = test_input.write_string("unchanged", **paradir)

        # read reference file
        filepath = os.path.join(PymatgenTest.TEST_FILES_DIR, "input_exciting2.xml")
        tree = ET.parse(filepath)
        root = tree.getroot()
        ref_string = ET.tostring(root, encoding="unicode")

        assert ref_string.strip() == test_string.strip()


if __name__ == "__main__":
    unittest.main()
