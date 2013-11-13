#!/usr/bin/env python
from __future__ import division, print_function

from pymatgen.util.testing import PymatgenTest
from pymatgen.io.abinitio.utils import *

class RpnTest(PymatgenTest):

    def test_mongodb_like_conditions(self):
        class Foo(object):
            one = 1.0
            two = 2.0
            three = 3.0
            four = 4.0

        map_res = [
            ( {"one": 1.0}, True),
            ( {"one": {"$eq": 1.0}}, True),
            ( {"one": {"$eq": "one"}}, True),
            ( {"one": {"$ne": "two"}}, True),
            ( {"one": {"$ne": 1.0}}, False),
            ( {"four": {"$divisible": 2.0}}, True),
            ( {"four": {"$divisible": 3.0}}, False),
            ( {"two": {"$gt": "one"}}, True ),
            ( {"$and": [ {"one": 1.0}, {"two": {"$lt": 3}}]}, True),
            ( {"$or": [ {"$not": {"one": 1.0}}, {"two": {"$lt": 20}}]}, True),
            ( {"$not": {"$and": [ {"$not": {"one": 1.0}}, {"two": {"$lt": 3}}] }}, True),
        ]

        for map, res in map_res:
            print("map", map)
            rpn = map2rpn(map, obj=Foo)
            print("rpn", rpn)
            self.assertTrue(res == evaluate_rpn(rpn))


if __name__ == '__main__':
    import unittest
    unittest.main()
