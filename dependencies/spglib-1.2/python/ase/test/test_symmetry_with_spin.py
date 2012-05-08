#!/usr/bin/evn python

import sys
from atoms import Atoms
from vasp import read_vasp
from pyspglib import spglib
import numpy as np

def get_magmom(text):
    magmom = []
    for numxmag in text.split():
        num, mag = numxmag.split('*')
        magmom += [float(mag)] * int(num)
    return magmom

def parse_incar(filename):
    for line in open(filename):
        for conf in line.split(';'):
            if 'MAGMOM' in conf:
                return get_magmom(conf.split('=')[1])

cell = read_vasp(sys.argv[1]) # POSCAR
magmoms = parse_incar(sys.argv[2]) # INCAR
cell.set_magnetic_moments(magmoms)
symmetry = spglib.get_symmetry(cell, symprec=1e-3)
print len(symmetry['rotations'])
