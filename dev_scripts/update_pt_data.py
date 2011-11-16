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
import re

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
    
def parse_oxi_state():
    with open('periodic_table.yaml', 'r') as f:
        data = yaml.load(f)
    f = open('oxidation_states.txt', 'r')
    oxidata = f.read()
    f.close()
    oxidata = re.sub('[\n\r]', '', oxidata)
    patt = re.compile('<tr>(.*?)</tr>', re.MULTILINE)
    
    for m in patt.finditer(oxidata):
        line = m.group(1)
        line = re.sub('</td>', '', line)
        line = re.sub('(<td>)+','<td>', line)
        line = re.sub('</*a[^>]*>','', line)
        el = None
        oxistates = []
        common_oxi = []
        for tok in re.split('<td>',line.strip()):
            m2 = re.match("<b>([A-Z][a-z]*)</b>", tok)
            if m2:
                el = m2.group(1)
            else:
                m3 = re.match("(<b>)*([\+\-]\d)(</b>)*", tok)
                if m3:
                    oxistates.append(int(m3.group(2)))
                    if m3.group(1):
                        common_oxi.append(int(m3.group(2)))
        if el in data:
            data[el]['oxidation_states'] = oxistates
            data[el]['common_oxidation_states'] = common_oxi
        else:
            print el
    with open('periodic_table2.yaml', 'w') as f:
        yaml.dump(data, f)
        
generate_json_from_yaml()