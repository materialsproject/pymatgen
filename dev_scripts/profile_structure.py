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
    nn = s.get_sites_in_sphere([0, 0, 0], 20)
    print len(nn)

def test2():
    from pymatgen.analysis.structure_fitter import StructureFitter
    fitter = StructureFitter(s, s)
    print fitter.fit_found

cProfile.run('test2()', 'testprof')
p = pstats.Stats('testprof')
p.sort_stats('cumulative').print_stats(20)
os.remove("testprof")
