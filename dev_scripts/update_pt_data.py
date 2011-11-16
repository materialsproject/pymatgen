#!/usr/bin/env python

'''
Developer script to convert yaml periodic table to json format.
Created on Nov 15, 2011
'''

from __future__ import division

__author__="Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Nov 15, 2011"

import json
import yaml

def generate_json_from_yaml():
    with open('periodic_table.yaml', 'r') as f:
        data = yaml.load(f)
    
    with open('periodic_table.json', 'w') as f:
        json.dump(data, f)
    
def test_yaml():
    with open('periodic_table.yaml', 'r') as f:
        data = yaml.load(f)

def test_json():
    with open('periodic_table.json', 'r') as f:
        data = json.load(f)

def test_yaml_json():
    from timeit import Timer
    t = Timer("test_yaml()", "from __main__ import test_yaml")
    print t.timeit(1)
    t = Timer("test_json()", "from __main__ import test_json")
    print t.timeit(1)
    