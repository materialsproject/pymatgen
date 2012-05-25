#!/usr/bin/env python
# encoding: utf-8
"""
test_etl.py
"""
__author__ = "Dan Gunter"
__copyright__ = "Copyright 2012, LBNL"
__maintainer__ = "Dan Gunter"
__email__ = "dkgunter@lbl.gov"
__date__ = "<DATE>"
__rcsid__ = "$Id$"

import pymongo
import time
import unittest

from pymatgen.etl import etl

class test_etl(unittest.TestCase):
    def setUp(self):
        self.uniq = str(int(time.time() * 1e6))
        self.dropdb = [ ]
    
    def tearDown(self):
        conn = pymongo.Connection()
        for dbname in self.dropdb:
            conn.drop_database(dbname)
                
    def testNoop(self):
        conf = { 
            "sources" : [{
                "collection" : "src1",
                "module" : ".etl.tests.test_etl",
                "class" : "NoAuth",
                }],
            "target": {
                "db" : "test_" + self.uniq,
                "collection" : "tgt"
            }
        }
        self.dropdb.append(conf["target"]["db"])
        #
        runner = etl.ETLRunner(conf)
        count = runner.run()
        self.assert_(count == 1)

class NoAuth(etl.ETLBase):
    def extract_transform_load(self):
        #print("tgt={}".format(self.tgt))
        self.tgt.insert({"ok":1})

if __name__ == '__main__':
    unittest.main()