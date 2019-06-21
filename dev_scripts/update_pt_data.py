#!/usr/bin/env python

'''
Developer script to convert yaml periodic table to json format.
Created on Nov 15, 2011
'''

from __future__ import division
import json
from itertools import product

import ruamel.yaml as yaml
import re

from monty.serialization import loadfn

from pymatgen import Element
from pymatgen.core.periodic_table import get_el_sp

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Nov 15, 2011"


def test_yaml():
    with open('periodic_table.yaml', 'r') as f:
        data = yaml.load(f)
        print(data)


def test_json():
    with open('periodic_table.json', 'r') as f:
        data = json.load(f)
        print(data)


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
            print(el)
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
    for i in range(1, len(radiidata)):
        line = radiidata[i]
        toks = line.strip().split(",")
        suffix = ""
        name = toks[1]
        if len(name.split(" ")) > 1:
            suffix = "_" + name.split(" ")[1]
        el = toks[2]

        ionic_radii = {}
        for j in range(3, len(toks)):
            m = re.match("^\s*([0-9\.]+)", toks[j])
            if m:
                ionic_radii[int(header[j])] = float(m.group(1))

        if el in data:
            data[el]['Ionic_radii' + suffix] = ionic_radii
            if suffix == '_hs':
                data[el]['Ionic_radii'] = ionic_radii
        else:
            print(el)
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
    for i in range(1, len(radiidata)):
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
            print(el)
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


def parse_shannon_radii():
    with open('periodic_table.yaml', 'r') as f:
        data = yaml.load(f)
    from openpyxl import load_workbook
    import collections
    wb = load_workbook('Shannon Radii.xlsx')
    print(wb.get_sheet_names())
    sheet = wb["Sheet1"]
    i = 2
    radii = collections.defaultdict(dict)
    while sheet["E%d" % i].value:
        if sheet["A%d" % i].value:
            el = sheet["A%d" % i].value
        if sheet["B%d" % i].value:
            charge = int(sheet["B%d" % i].value)
            radii[el][charge] = dict()
        if sheet["C%d" % i].value:
            cn = sheet["C%d" % i].value
            if cn not in radii[el][charge]:
                radii[el][charge][cn] = dict()

        if sheet["D%d" % i].value is not None:
            spin = sheet["D%d" % i].value
        else:
            spin = ""
        # print("%s - %d - %s" % (el, charge, cn))

        radii[el][charge][cn][spin] = {
            "crystal_radius": float(sheet["E%d" % i].value),
            "ionic_radius": float(sheet["F%d" % i].value),
        }
        i += 1

    for el in radii.keys():
        if el in data:
            data[el]["Shannon radii"] = dict(radii[el])

    with open('periodic_table.yaml', 'w') as f:
        yaml.safe_dump(data, f)
    with open('periodic_table.json', 'w') as f:
        json.dump(data, f)


def gen_periodic_table():
    with open('periodic_table.yaml', 'r') as f:
        data = yaml.load(f)

    with open('periodic_table.json', 'w') as f:
        json.dump(data, f)


def gen_iupac_ordering():
    periodic_table = loadfn("periodic_table.json")
    order = [([18], range(6, 0, -1)),  # noble gasses
             ([1], range(7, 1, -1)),  # alkali metals
             ([2], range(7, 1, -1)),  # alkali earth metals
             (range(17, 2, -1), [9]),  # actinides
             (range(17, 2, -1), [8]),  # lanthanides
             ([3], (5, 4)),  # Y, Sc
             ([4], (6, 5, 4)),  # Hf -> Ti
             ([5], (6, 5, 4)),  # Ta -> V
             ([6], (6, 5, 4)),  # W -> Cr
             ([7], (6, 5, 4)),  # Re -> Mn
             ([8], (6, 5, 4)),  # Os -> Fe
             ([9], (6, 5, 4)),  # Ir -> Co
             ([10], (6, 5, 4)),  # Pt -> Ni
             ([11], (6, 5, 4)),  # Au -> Cu
             ([12], (6, 5, 4)),  # Hg -> Zn
             ([13], range(6, 1, -1)),  # Tl -> B
             ([14], range(6, 1, -1)),  # Pb -> C
             ([15], range(6, 1, -1)),  # Bi -> N
             ([1], [1]),  # Hydrogen
             ([16], range(6, 1, -1)),  # Po -> O
             ([17], range(6, 1, -1))]  # At -> F

    order = sum([list(product(x, y)) for x, y in order], [])
    iupac_ordering_dict = dict(zip(
        [Element.from_row_and_group(row, group) for group, row in order],
        range(len(order))))

    # first clean periodic table of any IUPAC ordering
    for el in periodic_table:
        periodic_table[el].pop('IUPAC ordering', None)

    # now add iupac ordering
    for el in periodic_table:
        if 'IUPAC ordering' in periodic_table[el]:
            # sanity check that we don't cover the same element twice
            raise KeyError("IUPAC ordering already exists for {}".format(el))

        periodic_table[el]['IUPAC ordering'] = iupac_ordering_dict[get_el_sp(el)]


if __name__ == "__main__":
    parse_shannon_radii()
    #gen_periodic_table()
