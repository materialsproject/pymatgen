
from pymatgen.io.vasp import Vasprun

def parse_xml():
    v = Vasprun("../test_files/vasprun.xml")

if __name__ == "__main__":
    import timeit
    print(timeit.timeit("parse_xml()", setup="from __main__ import parse_xml",
        number=1))
