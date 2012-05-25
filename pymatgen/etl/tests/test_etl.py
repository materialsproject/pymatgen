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

import os
import pymongo
import tempfile
import time
import unittest

from pymatgen.etl import etl

class test_etl(unittest.TestCase):
    def setUp(self):
        self.uniq = str(int(time.time() * 1e6))
        self.dropdb = [ ]
        self.rmtemp = [ ]
        
    def tearDown(self):
        # drop temporary databases
        conn = pymongo.Connection()
        for dbname in self.dropdb:
            conn.drop_database(dbname)
        # remove temporary files
        for fname in self.rmtemp:
            try:
                os.unlink(fname)
            except OSError:
                pass

    # No-op tests
    # ------------
    
    def test_noop_noyaml(self):
        "no-op with dictionary config"
        conf = { 
            "sources" : [{
                "collection" : "src1",
                "module" : ".etl.tests.test_etl",
                "class" : "Noop",
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
        self.assert_(count == len(runner)*Noop.NUM_INSERTS)

    _yaml_noop_conf = """
sources:
    - collection: src1
      module: .etl.tests.test_etl
      class: Noop
    - collection: src2
      module: .etl.tests.test_etl
      class: Noop
target:
    db: {db}
    collection: tgt
"""
        
    def test_noop_yaml_str(self):
        "no-op with YAML string config"
        db = "test_" + self.uniq
        self.dropdb.append(db)
        conf = self._yaml_noop_conf.format(db=db)
        #
        runner = etl.ETLRunner(conf)
        count = runner.run()
        self.assert_(count == len(runner)*Noop.NUM_INSERTS)

    def test_noop_yaml_file(self):
        "no-op with YAML file config"
        db = "test_" + self.uniq
        self.dropdb.append(db)
        conf_str = self._yaml_noop_conf.format(db=db)
        fd, conf_name = tempfile.mkstemp(suffix=".conf")
        self.rmtemp.append(conf_name)
        os.write(fd, conf_str)
        conf = os.fdopen(fd, "r")
        conf.seek(0)
        #
        runner = etl.ETLRunner(conf)
        count = runner.run()
        self.assert_(count == len(runner)*Noop.NUM_INSERTS)
    
class Noop(etl.ETLBase):
    """No-op ETL class, except to insert a record to make sure
    the database is really accessible.
    """
    NUM_INSERTS = 1
    def extract_transform_load(self):
        #print("tgt={}".format(self.tgt))
        self.tgt.insert({"ok":1})

if __name__ == '__main__':
    unittest.main()