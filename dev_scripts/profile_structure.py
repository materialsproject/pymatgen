#!/usr/bin/env python

from pymatgen.io.vasp import Poscar
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
    from pymatgen.io.vasp import Chgcar
    c = Chgcar.from_file("../test_files/CHGCAR.noncubic")
    print c.get_integrated_diff(1, 2.5, 3)

def vasprun_test():
    from pymatgen.io.vasp import Vasprun
    v = Vasprun("../test_files/vasprun.xml")
    print v.final_energy

def matcher_test():
    p = Poscar.from_file("../test_files/POSCAR.Li2O")
    s = p.structure
    from pymatgen.analysis.structure_matcher import StructureMatcher
    print StructureMatcher().fit(s, s)


cProfile.run('matcher_test()', 'testprof')
p = pstats.Stats('testprof')
p.sort_stats('cumulative').print_stats(20)
os.remove("testprof")
