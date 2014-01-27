# coding: utf-8
import xml.etree.ElementTree as et
from pymatgen import Structure
read = False
read_data = False
structures = []
with open("vasprun.xml") as f:
    for ev, el in et.iterparse(f, events=["start", "end"]):
        if el.tag == "structure" and ev == "start":
            read = True
            structure = {}
        elif el.tag == "structure" and ev == "end":
            structures.append(structure)
            read = False
        elif read:
            if ev == "start" and el.attrib.get("name", "") in ["basis", "positions"]:
                data = []
                read_data = True
            elif ev == "end" and el.attrib.get("name", "") in ["basis", "positions"]:
                data = "".join(filter(lambda x: x is not None, data))
                structure[el.attrib.get("name", "")] = map(float, data.strip().split())
                read_data = False
            elif read_data:
                data.append(el.text)


#            data[el.tag + el.attrib.get("name", "")].append(el.text)
 
print structures[-1]       
#final = ''.join(structures[-1])
#data = map(float, final.strip().split())
#print data