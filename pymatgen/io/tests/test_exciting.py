# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

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

import unittest
import os
from pymatgen.io.exciting import excitingInput
from pymatgen.core import Structure, Lattice

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')
class ExcitingInputTest(unittest.TestCase):
    def test_fromfile(self):
        # Test for the import of a structure directly from an exciting
        # input file
        filepath=os.path.join(test_dir,'input1.xml')
        excin=excitingInput.from_file(filepath)
        lattice=[[5.62, 0.0, 0.0],[0.0, 5.62, 0.0],[0.0, 0.0, 5.62]]
        atoms=['Na', 'Na', 'Na', 'Na', 'Cl', 'Cl', 'Cl', 'Cl']
        fraccoords=[[0, 0, 0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5], 
                    [0.0, 0.5, 0.5], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0], 
                    [0.0, 0.0, 0.5], [0.5, 0.5, 0.5]]
        self.assertEqual(lattice, excin.lattice.matrix.tolist())
        self.assertEqual(atoms, [site.specie.symbol for site in excin])
        self.assertEqual(fraccoords,[site.frac_coords.tolist() for site in excin.sites])

    def test_writestring(self):
        # Test for the string export of s atructure into the exciting input xml schema
        input_string="""<input xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="http://xml.exciting-code.org/excitinginput.xsd">
  <title>Na4 Cl4</title>
  <structure speciespath="./">
    <crystal scale="1.8897543768634038">
      <basevect>      0.00000000       2.81000000       2.81000000</basevect>
      <basevect>      2.81000000       0.00000000       2.81000000</basevect>
      <basevect>      2.81000000       2.81000000       0.00000000</basevect>
    </crystal>
    <species speciesfile="Na.xml">
      <atom coord="      0.00000000       0.00000000       0.00000000" />
    </species>
    <species speciesfile="Cl.xml">
      <atom coord="      0.50000000       0.50000000       0.50000000" />
    </species>
  </structure>
  <properties>
    <bandstructure>
      <plot1d>
        <path steps="100">
          <point coord="      0.00000000       0.00000000       0.00000000" label="GAMMA" />
          <point coord="      0.50000000       0.00000000       0.50000000" label="X" />
          <point coord="      0.50000000       0.25000000       0.75000000" label="W" />
          <point coord="      0.37500000       0.37500000       0.75000000" label="K" />
          <point coord="      0.00000000       0.00000000       0.00000000" label="GAMMA" />
          <point coord="      0.50000000       0.50000000       0.50000000" label="L" />
          <point coord="      0.62500000       0.25000000       0.62500000" label="U" />
          <point coord="      0.50000000       0.25000000       0.75000000" label="W" />
          <point coord="      0.50000000       0.50000000       0.50000000" label="L" />
          <point coord="      0.37500000       0.37500000       0.75000000" label="K" />
        </path>
      </plot1d>
      <plot1d>
        <path steps="100">
          <point coord="      0.62500000       0.25000000       0.62500000" label="U" />
          <point coord="      0.50000000       0.00000000       0.50000000" label="X" />
        </path>
      </plot1d>
    </bandstructure>
  </properties>
</input>"""
        lattice=Lattice.cubic('5.62')
        structure=Structure(lattice,['Na','Na','Na','Na',
            'Cl','Cl','Cl','Cl'],
            [[0,0,0],[0.5,0.5,0.0],[0.5,0.0,0.5],[0.0,0.5,0.5],
            [0.5,0.0,0.0],[0.0,0.5,0.0],[0.0,0.0,0.5],[0.5,0.5,0.5]])
         excin=excitingInput(structure)
         assertEqual(input_string, excin.write_string('primitive', False, True))
        
