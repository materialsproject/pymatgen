#!/usr/bin/env python

from pymatgen.io.vaspio import Poscar
import cProfile
import pstats
import os
import logging

logging.basicConfig(level=logging.DEBUG)

p = Poscar.from_file("../test_files/POSCAR.LiFePO4", check_for_POTCAR=False)
s = p.structure

def test():
    nn = s.get_all_neighbors(20)
    print len(nn)

def chgcar_test():
    from pymatgen.io.vaspio import Chgcar
    c = Chgcar.from_file("../test_files/CHGCAR.noncubic")
    print c.get_integrated_diff(1, 2.5, 3)

cProfile.run('chgcar_test()', 'testprof')
p = pstats.Stats('testprof')
p.sort_stats('cumulative').print_stats(20)
os.remove("testprof")
