import pymatgen
from pymatgen.io.cifio import CifParser



struct = CifParser('/Users/dooheeyou/Desktop/pymatgen/pymatgen/phonons/aluminum.cif')
print struct.get_structures(primitive=False)



