# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import unittest
import os
import json
import warnings

from pymatgen import SETTINGS

from pymatgen.core.structure import Structure
from pymatgen.transformations.standard_transformations import \
    SubstitutionTransformation, PartialRemoveSpecieTransformation, \
    SupercellTransformation
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.alchemy.filters import ContainsSpecieFilter
from pymatgen.alchemy.materials import TransformedStructure
from pymatgen.util.provenance import StructureNL
from pymatgen.util.testing import PymatgenTest

'''
Created on Mar 5, 2012
'''


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Mar 5, 2012"



test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class TransformedStructureTest(PymatgenTest):

    def setUp(self):
        structure = PymatgenTest.get_structure("LiFePO4")
        self.structure = structure
        trans = [SubstitutionTransformation({"Li": "Na"})]
        self.trans = TransformedStructure(structure, trans)

    def test_append_transformation(self):
        t = SubstitutionTransformation({"Fe": "Mn"})
        self.trans.append_transformation(t)
        self.assertEqual("NaMnPO4",
                         self.trans.final_structure.composition
                         .reduced_formula)
        self.assertEqual(len(self.trans.structures), 3)
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        lattice = [[3.8401979337, 0.00, 0.00],
                   [1.9200989668, 3.3257101909, 0.00],
                   [0.00, -2.2171384943, 3.1355090603]]
        struct = Structure(lattice, ["Si4+", "Si4+"], coords)
        ts = TransformedStructure(struct, [])
        ts.append_transformation(SupercellTransformation
                                 .from_scaling_factors(2, 1, 1))
        alt = ts.append_transformation(
            PartialRemoveSpecieTransformation(
                'Si4+', 0.5,
                algo=PartialRemoveSpecieTransformation.ALGO_COMPLETE), 5)
        self.assertEqual(len(alt), 2)

    def test_append_filter(self):
        f3 = ContainsSpecieFilter(['O2-'], strict_compare=True, AND=False)
        self.trans.append_filter(f3)

    def test_get_vasp_input(self):
        SETTINGS["PMG_VASP_PSP_DIR"] = os.path.abspath(
                os.path.join(os.path.dirname(__file__), "..", "..", "..",
                             "test_files"))
        potcar = self.trans.get_vasp_input(MPRelaxSet)['POTCAR']
        self.assertEqual("Na_pv\nFe_pv\nP\nO",
                         "\n".join([p.symbol for p in potcar]))
        self.assertEqual(len(self.trans.structures), 2)

    def test_final_structure(self):
        self.assertEqual("NaFePO4", self.trans.final_structure.composition
                         .reduced_formula)

    def test_from_dict(self):
        d = json.load(open(os.path.join(test_dir, 'transformations.json'),
                           'r'))
        d['other_parameters'] = {'tags': ['test']}
        ts = TransformedStructure.from_dict(d)
        ts.other_parameters['author'] = 'Will'
        ts.append_transformation(SubstitutionTransformation({"Fe": "Mn"}))
        self.assertEqual("MnPO4",
                         ts.final_structure.composition.reduced_formula)
        self.assertEqual(ts.other_parameters, {'author': 'Will',
                                               'tags': ['test']})

    def test_undo_and_redo_last_change(self):
        trans = [SubstitutionTransformation({"Li": "Na"}),
                 SubstitutionTransformation({"Fe": "Mn"})]
        ts = TransformedStructure(self.structure, trans)
        self.assertEqual("NaMnPO4",
                         ts.final_structure.composition.reduced_formula)
        ts.undo_last_change()
        self.assertEqual("NaFePO4",
                         ts.final_structure.composition.reduced_formula)
        ts.undo_last_change()
        self.assertEqual("LiFePO4",
                         ts.final_structure.composition.reduced_formula)
        self.assertRaises(IndexError, ts.undo_last_change)
        ts.redo_next_change()
        self.assertEqual("NaFePO4",
                         ts.final_structure.composition.reduced_formula)
        ts.redo_next_change()
        self.assertEqual("NaMnPO4",
                         ts.final_structure.composition.reduced_formula)
        self.assertRaises(IndexError, ts.redo_next_change)
        #Make sure that this works with filters.
        f3 = ContainsSpecieFilter(['O2-'], strict_compare=True, AND=False)
        ts.append_filter(f3)
        ts.undo_last_change()
        ts.redo_next_change()

    def test_as_dict(self):
        self.trans.set_parameter('author', 'will')
        d = self.trans.as_dict()
        self.assertIn('last_modified', d)
        self.assertIn('history', d)
        self.assertIn('version', d)
        self.assertIn('author', d['other_parameters'])
        self.assertEqual(Structure.from_dict(d).formula, 'Na4 Fe4 P4 O16')

    def test_snl(self):
        self.trans.set_parameter('author', 'will')
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            snl = self.trans.to_snl([('will', 'will@test.com')])
            self.assertEqual(len(w), 1, 'Warning not raised on type conversion '
                             'with other_parameters')
        ts = TransformedStructure.from_snl(snl)
        self.assertEqual(ts.history[-1]['@class'], 'SubstitutionTransformation')
        
        h = ('testname', 'testURL', {'test' : 'testing'})
        snl = StructureNL(ts.final_structure,[('will', 'will@test.com')], 
                          history = [h])
        snl = TransformedStructure.from_snl(snl).to_snl([('notwill', 
                                                          'notwill@test.com')])
        self.assertEqual(snl.history, [h])
        self.assertEqual(snl.authors, [('notwill', 'notwill@test.com')])
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
