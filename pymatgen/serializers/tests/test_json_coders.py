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
from pymatgen.serializers.json_coders import PMGJSONEncoder, PMGJSONDecoder

class PMGJSONTest(unittest.TestCase):

    def test_core(self):
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        lattice = [[ 3.8401979337, 0.00, 0.00], [1.9200989668, 3.3257101909, 0.00], [0.00, -2.2171384943, 3.1355090603]]
        struct = Structure(lattice, ["Si4+", "Si4+"], coords)
        objs = [struct, struct[0], struct.lattice, struct[0].specie]
        for o in  objs:
            enc = PMGJSONEncoder()
            jsonstr = enc.encode(o)
            dec = PMGJSONDecoder()
            d = dec.decode(jsonstr)
            self.assertEqual(type(d), type(o))

        mol = Molecule(["O", "O"], coords)
        objs = [mol, mol[0]]
        for o in  objs:
            enc = PMGJSONEncoder()
            jsonstr = enc.encode(o)
            dec = PMGJSONDecoder()
            d = dec.decode(jsonstr)
            self.assertEqual(type(d), type(o))

        #Check dict of things
        o = {'structure':struct, "molecule":mol}
        enc = PMGJSONEncoder()
        jsonstr = enc.encode(o)
        dec = PMGJSONDecoder()
        d = dec.decode(jsonstr)
        self.assertEqual(type(d['structure']), Structure)
        self.assertEqual(type(d['molecule']), Molecule)

    def test_entry(self):
        entry = ComputedEntry("Fe2O3", 2.3)
        enc = PMGJSONEncoder()
        jsonstr = enc.encode(entry)
        dec = PMGJSONDecoder()
        d = dec.decode(jsonstr)
        self.assertEqual(type(d), ComputedEntry)

        #Check list of entries
        entries = [entry, entry, entry]
        enc = PMGJSONEncoder()
        jsonstr = enc.encode(entries)
        dec = PMGJSONDecoder()
        d = dec.decode(jsonstr)
        for i in d:
            self.assertEqual(type(i), ComputedEntry)
        self.assertEqual(len(d), 3)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
