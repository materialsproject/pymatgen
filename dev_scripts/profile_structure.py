#!/usr/bin/env python

from pymatgen.io.vaspio import Poscar
import cProfile
import pstats
import os

p = Poscar.from_file("../test_files/POSCAR.LiFePO4", check_for_POTCAR=False)
s = p.structure

def test():
    nn = s.get_sites_in_sphere([0, 0, 0], 20)
    print len(nn)

cProfile.run('test()', 'testprof')
p = pstats.Stats('testprof')
p.sort_stats('cumulative').print_stats(20)
os.remove("testprof")
