# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
Created on Apr 30, 2012
"""


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Apr 30, 2012"

import unittest2 as unittest

from pymatgen.core.structure import Structure, Molecule
from pymatgen.entries.computed_entries import ComputedEntry
import json

from monty.json import MontyEncoder, MontyDecoder
import datetime


class MontyTest(unittest.TestCase):

    def test_core(self):
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        lattice = [[3.8401979337, 0.00, 0.00],
                   [1.9200989668, 3.3257101909, 0.00],
                   [0.00, -2.2171384943, 3.1355090603]]
        struct = Structure(lattice, ["Si4+", "Si4+"], coords)
        objs = [struct, struct[0], struct.lattice, struct[0].species_and_occu,
                struct.composition]
        for o in objs:
            jsonstr = json.dumps(o, cls=MontyEncoder)
            d = json.loads(jsonstr, cls=MontyDecoder)
            self.assertEqual(type(d), type(o))

        mol = Molecule(["O", "O"], coords)
        objs = [mol, mol[0]]
        for o in objs:
            jsonstr = json.dumps(o, cls=MontyEncoder)
            d = json.loads(jsonstr, cls=MontyDecoder)
            self.assertEqual(type(d), type(o))

        #Check dict of things
        o = {'structure': struct, "molecule": mol}
        jsonstr = json.dumps(o, cls=MontyEncoder)
        d = json.loads(jsonstr, cls=MontyDecoder)
        self.assertEqual(type(d['structure']), Structure)
        self.assertEqual(type(d['molecule']), Molecule)

    def test_entry(self):
        enc = MontyEncoder()
        dec = MontyDecoder()

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

    def test_datetime(self):
        dt = datetime.datetime.now()
        jsonstr = json.dumps(dt, cls=MontyEncoder)
        d = json.loads(jsonstr, cls=MontyDecoder)
        self.assertEqual(type(d), datetime.datetime)
        self.assertEqual(dt, d)
        #Test a nested datetime.
        a = {'dt': dt, "a": 1}
        jsonstr = json.dumps(a, cls=MontyEncoder)
        d = json.loads(jsonstr, cls=MontyDecoder)
        self.assertEqual(type(d["dt"]), datetime.datetime)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
