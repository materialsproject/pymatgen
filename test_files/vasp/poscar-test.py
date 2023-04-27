import numpy as np
from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar

poscar = Poscar.from_file('POSCAR.conventional')
Ni_fcc_conventional = poscar.structure

poscar = Poscar.from_file('POSCAR.bulk')
Ni_fcc_bulk = poscar.structure

# Fixed the bottom layer of atoms along c direction
# False: fixed, True: mobile
fixed_indices = Ni_fcc_bulk.cart_coords[:,2] > Ni_fcc_conventional.lattice.c*0.5

poscar = Poscar(
    Ni_fcc_bulk,
    selective_dynamics=np.tile(fixed_indices.reshape(-1, 1), [1,3])
)

poscar.write_file('POSCAR')