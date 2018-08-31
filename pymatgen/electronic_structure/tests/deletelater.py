import os
from pymatgen.electronic_structure.cohp import CompleteCohp, Cohp, IcohpValue, IcohpCollection

structure = os.path.join("../../../test_files/cohp", "POSCAR.orbitalwise")
filepath = os.path.join("../../../test_files/cohp", "COHPCAR.lobster.notot.orbitalwise")
cohp_notot = CompleteCohp.from_file("lobster",filename=filepath,structure_file=structure)


