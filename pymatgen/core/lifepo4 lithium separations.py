

from fractions import gcd
import math
import numpy as np
import itertools

from pymatgen.core.structure import Structure
from pymatgen import Lattice, Structure
from
import numpy as np
import copy
import scipy.cluster.hierarchy

from pymatgen.io.smartio import CifParser
from pymatgen import write_structure

Li_Fe_PO4  = CifParser("/media/sharedfolder/LiFePO4/LiFePO4.cif")
lifepo4 = (Li_Fe_PO4.get_structures(primitive = False)[0])

l001 = Slab(lifepo4, [0, 0, 1], 4.787, 4.787, 0.001)
fileName = "/media/sharedfolder/LiFePO4/LiFePO4.cif"

Name = fileName[:-4] + "_%s_%s_%s_%s_" \
                          %(str(l001.miller_index),
                            l001.min_slab_size,
                            l001.min_vac_size,
                            l001.temp_thresh)
fileType = ".cif"

# For visual debugging
for iii in range(0, len(l001.enum)):
    name_num = str(iii)
    newFile = Name + name_num + fileType
    write_structure(l001.enum[iii], newFile)
print(len(l001.enum))
