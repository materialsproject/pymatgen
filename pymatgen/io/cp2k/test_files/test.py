from pymatgen.io.cp2k.inputs import *

ci = Cp2kInput.from_file('cp2k.inp')
print(ci)