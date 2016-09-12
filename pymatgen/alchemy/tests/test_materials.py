# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

'''
Created on Mar 5, 2012
'''


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Mar 5, 2012"

import unittest2 as unittest
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
from pymatgen.matproj.snl import StructureNL

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class TransformedStructureTest(unittest.TestCase):

    def setUp(self):
        structure_dict = {
            "lattice": {"a": 4.754150115,
                        "volume": 302.935463898643,
                        "c": 10.462573348,
                        "b": 6.090300362,
                        "matrix": [[4.754150115, 0.0, 0.0],
                                   [0.0, 6.090300362, 0.0],
                                   [0.0, 0.0, 10.462573348]],
                        "alpha": 90.0,
                        "beta": 90.0,
                        "gamma": 90.0},
            "sites": [{"occu": 1, "abc": [0.0, 0.0, 0.0],
                       "xyz": [0.0, 0.0, 0.0],
                       "species": [{"occu": 1, "element": "Li"}],
                       "label": "Li"}, {"occu": 1,
                                        "abc": [0.5000010396179928, 0.0,
                                                0.5000003178950235],
                                        "xyz": [2.37708, 0.0, 5.23129],
                                        "species": [{"occu": 1,
                                                     "element": "Li"}],
                                        "label": "Li"},
                                    {"occu": 1, "abc": [0.0,
                                                        0.49999997028061194,
                                                        0.0],
                                     "xyz": [0.0, 3.04515, 0.0],
                                     "species": [{"occu": 1, "element": "Li"}], "label": "Li"}, {"occu": 1, "abc": [0.5000010396179928, 0.49999997028061194, 0.5000003178950235], "xyz": [2.37708, 3.04515, 5.23129], "species": [{"occu": 1, "element": "Li"}], "label": "Li"}, {"occu": 1, "abc": [0.7885825876997996, 0.5473161916279229, 0.3339168944194627], "xyz": [3.74904, 3.33332, 3.4936300000000005], "species": [{"occu": 1, "element": "O"}], "label": "O"}, {"occu": 1, "abc": [0.2114173881108085, 0.452683748933301, 0.6660827855827808], "xyz": [1.00511, 2.75698, 6.968940000000001], "species": [{"occu": 1, "element": "O"}], "label": "O"}, {"occu": 1, "abc": [0.7114184277288014, 0.5473161916279229, 0.8339172123144861], "xyz": [3.38219, 3.33332, 8.72492], "species": [{"occu": 1, "element": "O"}], "label": "O"}, {"occu": 1, "abc": [0.7885825876997996, 0.9526820772587701, 0.3339168944194627], "xyz": [3.74904, 5.8021199999999995, 3.4936300000000005], "species": [{"occu": 1, "element": "O"}], "label": "O"}, {"occu": 1, "abc": [0.28858365150718424, 0.047317863302453654, 0.16608342347556082], "xyz": [1.37197, 0.28818, 1.73766], "species": [{"occu": 1, "element": "O"}], "label": "O"}, {"occu": 1, "abc": [0.7440972443925447, 0.25000080611787734, 0.09613791622232937], "xyz": [3.537549999999999, 1.52258, 1.00585], "species": [{"occu": 1, "element": "O"}], "label": "O"}, {"occu": 1, "abc": [0.28858365150718424, 0.452683748933301, 0.16608342347556082], "xyz": [1.37197, 2.75698, 1.73766], "species": [{"occu": 1, "element": "O"}], "label": "O"}, {"occu": 1, "abc": [0.2114173881108085, 0.047317863302453654, 0.6660827855827808], "xyz": [1.00511, 0.28818, 6.968940000000001], "species": [{"occu": 1, "element": "O"}], "label": "O"}, {"occu": 1, "abc": [0.2559006279926859, 0.7499991344433464, 0.9038627195677177], "xyz": [1.21659, 4.56772, 9.45673], "species": [{"occu": 1, "element": "O"}], "label": "O"}, {"occu": 1, "abc": [0.7559016676106785, 0.25000080611787734, 0.5961372783295493], "xyz": [3.5936699999999986, 1.52258, 6.2371300000000005], "species": [{"occu": 1, "element": "O"}], "label": "O"}, {"occu": 1, "abc": [0.7939989080466804, 0.7499991344433464, 0.5421304884886912], "xyz": [3.77479, 4.56772, 5.67208], "species": [{"occu": 1, "element": "O"}], "label": "O"}, {"occu": 1, "abc": [0.24409830819992942, 0.7499991344433464, 0.40386240167269416], "xyz": [1.16048, 4.56772, 4.22544], "species": [{"occu": 1, "element": "O"}], "label": "O"}, {"occu": 1, "abc": [0.7060021073819206, 0.7499991344433464, 0.04213017059366761], "xyz": [3.35644, 4.56772, 0.44079000000000007], "species": [{"occu": 1, "element": "O"}], "label": "O"}, {"occu": 1, "abc": [0.2939978684286875, 0.25000080611787734, 0.9578695094085758], "xyz": [1.3977099999999996, 1.52258, 10.02178], "species": [{"occu": 1, "element": "O"}], "label": "O"}, {"occu": 1, "abc": [0.20600106776392774, 0.25000080611787734, 0.4578701473013559], "xyz": [0.9793599999999998, 1.52258, 4.7905], "species": [{"occu": 1, "element": "O"}], "label": "O"}, {"occu": 1, "abc": [0.7114184277288014, 0.9526820772587701, 0.8339172123144861], "xyz": [3.38219, 5.8021199999999995, 8.72492], "species": [{"occu": 1, "element": "O"}], "label": "O"}, {"occu": 1, "abc": [0.5793611756830275, 0.7499991344433464, 0.9051119342269868], "xyz": [2.75437, 4.56772, 9.4698], "species": [{"occu": 1, "element": "P"}], "label": "P"}, {"occu": 1, "abc": [0.9206377363201961, 0.7499991344433464, 0.40511161633196324], "xyz": [4.37685, 4.56772, 4.23851], "species": [{"occu": 1, "element": "P"}], "label": "P"}, {"occu": 1, "abc": [0.42063880012758065, 0.25000080611787734, 0.09488774577525667], "xyz": [1.9997799999999994, 1.52258, 0.99277], "species": [{"occu": 1, "element": "P"}], "label": "P"}, {"occu": 1, "abc": [0.07936223949041206, 0.25000080611787734, 0.5948880636702801], "xyz": [0.3773, 1.52258, 6.22406], "species": [{"occu": 1, "element": "P"}], "label": "P"}, {"occu": 1, "abc": [0.021860899947623972, 0.7499991344433464, 0.7185507570598875], "xyz": [0.10393, 4.56772, 7.517890000000001], "species": [{"occu": 1, "element": "Fe"}], "label": "Fe"}, {"occu": 1, "abc": [0.478135932819614, 0.7499991344433464, 0.21855043916486389], "xyz": [2.27313, 4.56772, 2.2866], "species": [{"occu": 1, "element": "Fe"}], "label": "Fe"}, {"occu": 1, "abc": [0.9781369724376069, 0.25000080611787734, 0.2814489229423561], "xyz": [4.65021, 1.52258, 2.9446800000000004], "species": [{"occu": 1, "element": "Fe"}], "label": "Fe"}, {"occu": 1, "abc": [0.5218619395656168, 0.25000080611787734, 0.7814492408373795], "xyz": [2.48101, 1.52258, 8.17597], "species": [{"occu": 1, "element": "Fe"}], "label": "Fe"}]}
        structure = Structure.from_dict(structure_dict)
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
        SETTINGS["VASP_PSP_DIR"] = os.path.abspath(
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
