#!/usr/bin/env python

'''
Developer script to convert yaml periodic table to json format.
Created on Nov 15, 2011
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Nov 15, 2011"

import json
import yaml
import re


def test_yaml():
    with open('periodic_table.yaml', 'r') as f:
        data = yaml.load(f)
        print data


def test_json():
    with open('periodic_table.json', 'r') as f:
        data = json.load(f)
        print data


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
        line = re.sub('(<td>)+', '<td>', line)
        line = re.sub('</*a[^>]*>', '', line)
        el = None
        oxistates = []
        common_oxi = []
        for tok in re.split('<td>', line.strip()):
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
            del data[el]['Max oxidation state']
            del data[el]['Min oxidation state']
            del data[el]['Oxidation_states']
            del data[el]['Common_oxidation_states']
            data[el]['Oxidation states'] = oxistates
            data[el]['Common oxidation states'] = common_oxi
        else:
            print el
    with open('periodic_table2.yaml', 'w') as f:
        yaml.dump(data, f)


def parse_ionic_radii():
    with open('periodic_table.yaml', 'r') as f:
        data = yaml.load(f)
    f = open('ionic_radii.csv', 'r')
    radiidata = f.read()
    f.close()
    radiidata = radiidata.split("\r")
    header = radiidata[0].split(",")
    for i in xrange(1, len(radiidata)):
        line = radiidata[i]
        toks = line.strip().split(",")
        suffix = ""
        name = toks[1]
        if len(name.split(" ")) > 1:
            suffix = "_" + name.split(" ")[1]
        el = toks[2]

        ionic_radii = {}
        for j in xrange(3, len(toks)):
            m = re.match("^\s*([0-9\.]+)", toks[j])
            if m:
                ionic_radii[int(header[j])] = float(m.group(1))

        if el in data:
            data[el]['Ionic_radii' + suffix] = ionic_radii
            if suffix == '_hs':
                data[el]['Ionic_radii'] = ionic_radii
        else:
            print el
    with open('periodic_table2.yaml', 'w') as f:
        yaml.dump(data, f)


def parse_radii():
    with open('periodic_table.yaml', 'r') as f:
        data = yaml.load(f)
    f = open('radii.csv', 'r')
    radiidata = f.read()
    f.close()
    radiidata = radiidata.split("\r")
    header = radiidata[0].split(",")
    for i in xrange(1, len(radiidata)):
        line = radiidata[i]
        toks = line.strip().split(",")
        el = toks[1]
        try:
            atomic_radii = float(toks[3]) / 100
        except:
            atomic_radii = toks[3]

        try:
            atomic_radii_calc = float(toks[4]) / 100
        except:
            atomic_radii_calc = toks[4]

        try:
            vdw_radii = float(toks[5]) / 100
        except:
            vdw_radii = toks[5]

        if el in data:
            data[el]['Atomic radius'] = atomic_radii
            data[el]['Atomic radius calculated'] = atomic_radii_calc
            data[el]['Van der waals radius'] = vdw_radii
        else:
            print el
    with open('periodic_table2.yaml', 'w') as f:
        yaml.dump(data, f)
    with open('periodic_table.json', 'w') as f:
        json.dump(data, f)


def update_ionic_radii():
    with open('periodic_table.yaml', 'r') as f:
        data = yaml.load(f)

    for el, d in data.items():
        if "Ionic_radii" in d:
            d["Ionic radii"] = {k: v / 100
                                for k, v in d["Ionic_radii"].items()}
            del d["Ionic_radii"]
        if "Ionic_radii_hs" in d:
            d["Ionic radii hs"] = {k: v / 100
                                   for k, v in d["Ionic_radii_hs"].items()}
            del d["Ionic_radii_hs"]
        if "Ionic_radii_ls" in d:
            d["Ionic radii ls"] = {k: v / 100
                                   for k, v in d["Ionic_radii_ls"].items()}
            del d["Ionic_radii_ls"]
    with open('periodic_table2.yaml', 'w') as f:
        yaml.dump(data, f)
    with open('periodic_table.json', 'w') as f:
        json.dump(data, f)


def gen_periodic_table():
    with open('periodic_table.yaml', 'r') as f:
        data = yaml.load(f)

    with open('periodic_table.json', 'w') as f:
        json.dump(data, f)

if __name__ == "__main__":
    gen_periodic_table()
