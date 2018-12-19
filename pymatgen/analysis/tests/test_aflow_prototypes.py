# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


from pymatgen.util.testing import PymatgenTest
from pymatgen.analysis.aflow_prototypes import AflowPrototypeMatcher

class AflowPrototypeMatcherTest(PymatgenTest):

    def test_prototype_matching(self):

        af = AflowPrototypeMatcher()

        struct = self.get_structure('Sn')
        prototype = af.get_prototypes(struct)[0]

        self.assertDictEqual(prototype['tags'], {
            'aflow': 'A_cF8_227_a',
            'mineral': 'diamond',
            'pearson': 'cF8',
            'strukturbericht': 'A4'
        })

        struct = self.get_structure('CsCl')
        prototype = af.get_prototypes(struct)[0]

        self.assertDictEqual(prototype['tags'], {
            'aflow': 'AB_cP2_221_b_a',
            'mineral': '',
            'pearson': 'cP2',
            'strukturbericht': 'B2'
        })

        struct = self.get_structure('Li2O')
        prototype = af.get_prototypes(struct)[0]

        self.assertDictEqual(prototype['tags'], {
            'aflow': 'AB2_cF12_225_a_c',
            'mineral': 'Fluorite',
            'pearson': 'cF12',
            'strukturbericht': 'C1'
        })
