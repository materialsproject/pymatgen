#!/usr/bin/env python

'''
Created on Apr 30, 2012
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Apr 30, 2012"

import unittest

from pymatgen.core.structure import Structure, Molecule
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.transformations.standard_transformations import IdentityTransformation
import json
from pymatgen.serializers.json_coders import PMGJSONEncoder, PMGJSONDecoder, MSONable


class MSONableTest(unittest.TestCase):

    def setUp(self):
        class GoodMSONClass(MSONable):

            def __init__(self, a, b):
                self.a = a
                self.b = b

            @property
            def to_dict(self):
                d = {'a':self.a, 'b':self.b}
                return d

            @staticmethod
            def from_dict(d):
                return GoodMSONClass(d['a'], d['b'])

        self.good_cls = GoodMSONClass

        class BadMSONClass(MSONable):

            def __init__(self, a, b):
                self.a = a
                self.b = b

            @property
            def to_dict(self):
                d = {'a':self.a, 'b':self.b}
                return d

        self.bad_cls = BadMSONClass

    def test_to_from_dict(self):
        obj = self.good_cls("Hello", "World")
        d = obj.to_dict
        self.assertIsNotNone(d)
        self.good_cls.from_dict(d)
        obj = self.bad_cls("Hello", "World")
        d = obj.to_dict
        self.assertIsNotNone(d)
        self.assertRaises(NotImplementedError, self.bad_cls.from_dict, d)

    def test_to_json(self):
        obj = self.good_cls("Hello", "World")
        self.assertIsNotNone(obj.to_json)


class PMGJSONTest(unittest.TestCase):

    def test_core(self):
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        lattice = [[ 3.8401979337, 0.00, 0.00], [1.9200989668, 3.3257101909, 0.00], [0.00, -2.2171384943, 3.1355090603]]
        struct = Structure(lattice, ["Si4+", "Si4+"], coords)
        objs = [struct, struct[0], struct.lattice, struct[0].specie, struct.composition]
        for o in  objs:
            jsonstr = json.dumps(o, cls=PMGJSONEncoder)
            d = json.loads(jsonstr, cls=PMGJSONDecoder)
            self.assertEqual(type(d), type(o))

        mol = Molecule(["O", "O"], coords)
        objs = [mol, mol[0]]
        for o in  objs:
            jsonstr = json.dumps(o, cls=PMGJSONEncoder)
            d = json.loads(jsonstr, cls=PMGJSONDecoder)
            self.assertEqual(type(d), type(o))

        #Check dict of things
        o = {'structure':struct, "molecule":mol}
        jsonstr = json.dumps(o, cls=PMGJSONEncoder)
        d = json.loads(jsonstr, cls=PMGJSONDecoder)
        self.assertEqual(type(d['structure']), Structure)
        self.assertEqual(type(d['molecule']), Molecule)

    def test_entry(self):
        enc = PMGJSONEncoder()
        dec = PMGJSONDecoder()

        entry = ComputedEntry("Fe2O3", 2.3)
        jsonstr = enc.encode(entry)
        d = dec.decode(jsonstr)
        self.assertEqual(type(d), ComputedEntry)

        #Check list of entries
        entries = [entry, entry, entry]
        jsonstr = enc.encode(entries)
        d = dec.decode(jsonstr)
        for i in d:
            self.assertEqual(type(i), ComputedEntry)
        self.assertEqual(len(d), 3)

    def test_transformations(self):
        trans = IdentityTransformation()
        jsonstr = json.dumps(trans, cls=PMGJSONEncoder)
        d = json.loads(jsonstr, cls=PMGJSONDecoder)
        self.assertEqual(type(d), IdentityTransformation)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
