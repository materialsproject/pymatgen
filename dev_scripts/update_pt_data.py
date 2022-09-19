#!/usr/bin/env python

"""
Developer script to convert yaml periodic table to json format.
Created on Nov 15, 2011
"""

import json
import re
from itertools import product

import ruamel.yaml as yaml
from monty.serialization import dumpfn, loadfn

from pymatgen.core import Element
from pymatgen.core.periodic_table import get_el_sp


def test_yaml():
    with open("periodic_table.yaml") as f:
        data = yaml.load(f)
        print(data)


def test_json():
    with open("periodic_table.json") as f:
        data = json.load(f)
        print(data)


def parse_oxi_state():
    with open("periodic_table.yaml") as f:
        data = yaml.load(f)
    with open("oxidation_states.txt") as f:
        oxidata = f.read()
    oxidata = re.sub("[\n\r]", "", oxidata)
    patt = re.compile("<tr>(.*?)</tr>", re.MULTILINE)

    for m in patt.finditer(oxidata):
        line = m.group(1)
        line = re.sub("</td>", "", line)
        line = re.sub("(<td>)+", "<td>", line)
        line = re.sub("</*a[^>]*>", "", line)
        el = None
        oxistates = []
        common_oxi = []
        for tok in re.split("<td>", line.strip()):
            m2 = re.match(r"<b>([A-Z][a-z]*)</b>", tok)
            if m2:
                el = m2.group(1)
            else:
                m3 = re.match(r"(<b>)*([\+\-]\d)(</b>)*", tok)
                if m3:
                    oxistates.append(int(m3.group(2)))
                    if m3.group(1):
                        common_oxi.append(int(m3.group(2)))
        if el in data:
            del data[el]["Max oxidation state"]
            del data[el]["Min oxidation state"]
            del data[el]["Oxidation_states"]
            del data[el]["Common_oxidation_states"]
            data[el]["Oxidation states"] = oxistates
            data[el]["Common oxidation states"] = common_oxi
        else:
            print(el)
    with open("periodic_table2.yaml", "w") as f:
        yaml.dump(data, f)


def parse_ionic_radii():
    with open("periodic_table.yaml") as f:
        data = yaml.load(f)
    with open("ionic_radii.csv") as f:
        radiidata = f.read()
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
            m = re.match(r"^\s*([0-9\.]+)", toks[j])
            if m:
                ionic_radii[int(header[j])] = float(m.group(1))

        if el in data:
            data[el]["Ionic_radii" + suffix] = ionic_radii
            if suffix == "_hs":
                data[el]["Ionic_radii"] = ionic_radii
        else:
            print(el)
    with open("periodic_table2.yaml", "w") as f:
        yaml.dump(data, f)


def parse_radii():
    with open("periodic_table.yaml") as f:
        data = yaml.load(f)
    with open("radii.csv") as f:
        radiidata = f.read()
    radiidata = radiidata.split("\r")

    for i in range(1, len(radiidata)):
        line = radiidata[i]
        toks = line.strip().split(",")
        el = toks[1]
        try:
            atomic_radii = float(toks[3]) / 100
        except Exception:
            atomic_radii = toks[3]

        try:
            atomic_radii_calc = float(toks[4]) / 100
        except Exception:
            atomic_radii_calc = toks[4]

        try:
            vdw_radii = float(toks[5]) / 100
        except Exception:
            vdw_radii = toks[5]

        if el in data:
            data[el]["Atomic radius"] = atomic_radii
            data[el]["Atomic radius calculated"] = atomic_radii_calc
            data[el]["Van der waals radius"] = vdw_radii
        else:
            print(el)
    with open("periodic_table2.yaml", "w") as f:
        yaml.dump(data, f)
    with open("periodic_table.json", "w") as f:
        json.dump(data, f)


def update_ionic_radii():
    with open("periodic_table.yaml") as f:
        data = yaml.load(f)

    for el, d in data.items():
        if "Ionic_radii" in d:
            d["Ionic radii"] = {k: v / 100 for k, v in d["Ionic_radii"].items()}
            del d["Ionic_radii"]
        if "Ionic_radii_hs" in d:
            d["Ionic radii hs"] = {k: v / 100 for k, v in d["Ionic_radii_hs"].items()}
            del d["Ionic_radii_hs"]
        if "Ionic_radii_ls" in d:
            d["Ionic radii ls"] = {k: v / 100 for k, v in d["Ionic_radii_ls"].items()}
            del d["Ionic_radii_ls"]
    with open("periodic_table2.yaml", "w") as f:
        yaml.dump(data, f)
    with open("periodic_table.json", "w") as f:
        json.dump(data, f)


def parse_shannon_radii():
    with open("periodic_table.yaml") as f:
        data = yaml.load(f)
    import collections

    from openpyxl import load_workbook

    wb = load_workbook("Shannon Radii.xlsx")
    print(wb.get_sheet_names())
    sheet = wb["Sheet1"]
    i = 2
    radii = collections.defaultdict(dict)
    while sheet[f"E{i}"].value:
        if sheet[f"A{i}"].value:
            el = sheet[f"A{i}"].value
        if sheet[f"B{i}"].value:
            charge = int(sheet[f"B{i}"].value)
            radii[el][charge] = dict()
        if sheet[f"C{i}"].value:
            cn = sheet[f"C{i}"].value
            if cn not in radii[el][charge]:
                radii[el][charge][cn] = dict()

        if sheet[f"D{i}"].value is not None:
            spin = sheet[f"D{i}"].value
        else:
            spin = ""

        radii[el][charge][cn][spin] = {
            "crystal_radius": float(sheet[f"E{i}"].value),
            "ionic_radius": float(sheet[f"F{i}"].value),
        }
        i += 1

    for el in radii:
        if el in data:
            data[el]["Shannon radii"] = dict(radii[el])

    with open("periodic_table.yaml", "w") as f:
        yaml.safe_dump(data, f)
    with open("periodic_table.json", "w") as f:
        json.dump(data, f)


def gen_periodic_table():
    with open("periodic_table.yaml") as f:
        data = yaml.load(f)

    with open("periodic_table.json", "w") as f:
        json.dump(data, f)


def gen_iupac_ordering():
    periodic_table = loadfn("periodic_table.json")
    order = [
        ([18], range(6, 0, -1)),  # noble gasses
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
        ([17], range(6, 1, -1)),
    ]  # At -> F

    order = sum((list(product(x, y)) for x, y in order), [])
    iupac_ordering_dict = dict(zip([Element.from_row_and_group(row, group) for group, row in order], range(len(order))))

    # first clean periodic table of any IUPAC ordering
    for el in periodic_table:
        periodic_table[el].pop("IUPAC ordering", None)

    # now add iupac ordering
    for el in periodic_table:
        if "IUPAC ordering" in periodic_table[el]:
            # sanity check that we don't cover the same element twice
            raise KeyError(f"IUPAC ordering already exists for {el}")

        periodic_table[el]["IUPAC ordering"] = iupac_ordering_dict[get_el_sp(el)]


def add_electron_affinities():
    """
    Update the periodic table data file with electron affinities.
    """
    import requests
    from bs4 import BeautifulSoup

    req = requests.get("https://en.wikipedia.org/wiki/Electron_affinity_(data_page)")
    soup = BeautifulSoup(req.text, "html.parser")
    for t in soup.find_all("table"):
        if "Hydrogen" in t.text:
            break
    data = []
    for tr in t.find_all("tr"):
        row = []
        for td in tr.find_all("td"):
            row.append(td.get_text().strip())
        data.append(row)
    data.pop(0)
    ea = {int(r[0]): float(re.sub(r"[\s\(\)]", "", r[3].strip("()[]"))) for r in data}
    assert set(ea).issuperset(range(1, 93))  # Ensure that we have data for up to U.
    print(ea)
    pt = loadfn("../pymatgen/core/periodic_table.json")
    for k, v in pt.items():
        v["Electron affinity"] = ea.get(Element(k).Z, None)
    dumpfn(pt, "../pymatgen/core/periodic_table.json")


def add_ionization_energies():
    """
    Update the periodic table data file with ground level and ionization energies from NIST.
    """
    import collections

    from bs4 import BeautifulSoup

    with open("NIST Atomic Ionization Energies Output.html") as f:
        soup = BeautifulSoup(f.read(), "html.parser")
    for t in soup.find_all("table"):
        if "Hydrogen" in t.text:
            break
    data = collections.defaultdict(list)
    for tr in t.find_all("tr"):
        row = [td.get_text().strip() for td in tr.find_all("td")]
        if row:
            Z = int(row[0])
            val = re.sub(r"\s", "", row[8].strip("()[]"))
            if val == "":
                val = None
            else:
                val = float(val)
            data[Z].append(val)
    print(data)
    print(data[51])
    assert set(data).issuperset(range(1, 93))  # Ensure that we have data for up to U.
    pt = loadfn("../pymatgen/core/periodic_table.json")
    for k, v in pt.items():
        del v["Ionization energy"]
        v["Ionization energies"] = data.get(Element(k).Z, [])
    dumpfn(pt, "../pymatgen/core/periodic_table.json")


if __name__ == "__main__":
    # parse_shannon_radii()
    # add_ionization_energies()
    add_electron_affinities()
    # gen_periodic_table()
